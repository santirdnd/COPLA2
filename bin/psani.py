#!/usr/bin/env python3

# psani is a Python re-implementation, with additional features, of the two-way ANI
# of ani.rb from the Enveomics collection (Luis M. Rodriguez-R, Artistic License 2.0,
# https://github.com/lmrodriguezr/enveomics). Distributed under the GPL-3.0 license.

import os
import re
import sys
import gzip
import json
import math
import shutil
import pathlib
import argparse
import tempfile
import subprocess
import multiprocessing


prog = os.path.basename(sys.argv[0])

# Version of the database JSON format. Increase it when its keys change. Databases
# created with psani 1.0 do not store it and are taken as format 1
DB_FORMAT_VERSION = 2
# Keys of the database JSON
DB_REQUIRED = ('format_version', 'db_name', 'input_file', 'window', 'step',
               'ji_sketch_db', 'ji_sketchsize64', 'ji_kmerlen', 'ji_minhashtype', 'ji_dens',
               'ji_bbits', 'ji_randseed', 'seqs')
# Decimals of the ANI, its SD and the coverage in the dist output. A single SNP changes
# the ANI of a plasmid of length L by about 100/L, whatever the window and step, so 5
# decimals tell it apart from 100 up to about 20 Mb
OUTPUT_DECIMALS = 5

regex_notAlpha = re.compile(r'[^A-Za-z]')
regex_allN = re.compile(r'^[Nn]+$')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='BLAST-based ANI computation with optional Jaccard pre-screening using BinDash')
    parser.add_argument('--version', action='version', version='%(prog)s 2.0.0')
    subparsers = parser.add_subparsers(
        title='subcommands', dest='subcommand', required=True)

    parser_makedb = subparsers.add_parser('makedb', aliases=['mk'],
        help='create database of sequence fragments for ANI computation')
    parser_makedb.add_argument('input', type=str,
        help='sequence file (or a file formatted as a file of filenames). Gzip-compressed sequences are supported')
    parser_makedb.add_argument('dbdir', type=str,
        help='database directory to create (database JSON, fwa/ and kmer/)')
    parser_makedb.add_argument('--absolute-paths', action='store_true',
        help='write absolute paths in the database JSON. Default: paths relative to the database directory, so the database can be moved')
    parser_makedb.add_argument('-d', '--db-name', type=str,
        help='database name for input sequences. Default: the basename of the database directory')
    parser_makedb.add_argument('--seq-id-source', type=str,
        choices=['header', 'filename'], default='header',
        help='source for the ID of each sequence: filename or the first FASTA header. Default: header')
    parser_makedb.add_argument('-w', '--window', type=int,
        default=1000,
        help='sliding window fragment size. Default: 1000')
    parser_makedb.add_argument('-s', '--step', type=int,
        default=200,
        help='sliding window step. Default: 200')
    parser_makedb.add_argument('--ji-sketchsize64', type=int,
        default=32,
        help='number of sketches (in multiples of 64) for JI calculation. Use 0 to deactivate the JI shortcut. Default: 32 (i.e. 2048 sketches)')
    parser_makedb.add_argument('--ji-kmerlen', type=int,
        default=21,
        help='kmer length to avoid ANI comparisons by precomputing their JI. Default: 21')
    parser_makedb.add_argument('--ji-minhashtype', type=int,
        choices=[-1, 0, 1, 2],
        default=2,
        help='Type of minhash. Default: 2 (for more info see Bindash options)')
    parser_makedb.add_argument('--ji-dens', type=int,
        choices=[1, 2],
        default=1,
        help='densification strategy for minhashtype 2. Default: 1 (for more info see Bindash options)')
    parser_makedb.add_argument('--ji-bbits', type=int,
        default=14,
        help='number of bits kept for JI calculation as in b-bits minhash. Default: 14 (for more info see Bindash options)')
    parser_makedb.add_argument('--ji-randseed', type=int,
        default=41,
        help='random seed for the minhash sketches. Query and reference databases must use the same seed. Default: 41 (for more info see Bindash options)')
    parser_makedb.add_argument('-n', '--n-jobs', type=int,
        default=0,
        help='number of processes to spawn or threads to use. Default: 0 (use all available CPUs; for more info see os.sched_getaffinity())')
    parser_makedb.add_argument('--chunksize', type=int,
        default=10,
        help='chunksize of parallel tasks. Default: 10 (for more info see multiprocessing.Pool().imap_unordered())')
    parser_makedb.add_argument('-q', '--quiet', action='store_true',
        help='do not report the progress on stderr (warnings and errors are still reported)')
    parser_makedb.add_argument('--force', action='store_true',
        help='overwrite an existing database in the database directory (database JSON, fwa/ and kmer/). Other files in the directory are kept')
    parser_makedb.set_defaults(func=makedb)

    parser_dist = subparsers.add_parser('dist', aliases=['di'],
        help='compute ANI comparisons between sets of genomes (databases must be generated using the %(prog)s makedb subcommand)')
    parser_dist.add_argument('query_db', type=str,
        help='query database (JSON file defining the query database)')
    parser_dist.add_argument('reference_db', type=str,
        default = None, nargs='?',
        help='reference database (JSON file defining the reference database). If not present, %(prog)s will compare the query against itself and not perform redundant comparisons (upper triangular matrix)')
    parser_dist.add_argument('-o', '--output', type=str,
        default='-',
        help='ANI output file. Default: - (/dev/stdout)')
    parser_dist.add_argument('--header', dest='header', action='store_true',
        help='write output file with a header. Default: True')
    parser_dist.add_argument('--no-header', dest='header', action='store_false',
        help='write output file without header')
    parser_dist.set_defaults(header=True)
    parser_dist.add_argument('-i', '--fragment-id', type=float,
        default=0.7,
        help='minimum identity for a fragment to be considered homologous. Default: 0.7')
    parser_dist.add_argument('-l', '--fragment-len', type=float,
        default=0.7,
        help='minimum percentage of fragment length to be considered homologous. If value > 1, it will be considered as the minimum length in bp. Default: 0.7')
    parser_dist.add_argument('--correct-ani', dest='correct_ani', action='store_true',
        help='correct the identity of each homologous fragment as identity_corr = 100 - (100 - identity) / 0.8621 before averaging, which lowers the reported ANI. The fragments are selected with the uncorrected identity (see -i). Carried over from ani.rb of the Enveomics collection, where it is applied by default. Default: False')
    parser_dist.add_argument('--no-correct-ani', dest='correct_ani', action='store_false',
        help='report ANI values without post-hoc correction')
    parser_dist.set_defaults(correct_ani=False)
    parser_dist.add_argument('-c', '--coverage', type=float,
        default=0.5,
        help='minimum fraction of the smallest genome of each comparison. If value >= 1, it will be ceiled and considered as the minimum number of fragments. Default: 0.5')
    parser_dist.add_argument('-r', '--size-ratio', type=float,
        default=0,
        help='ratio of genome sizes allowed for an ANI comparison to be considered valid. Set to 0 to deactivate this constraint. Default: 0')
    parser_dist.add_argument('-m', '--mut-threshold', type=float,
        default=0.3,
        help='mutational distance threshold to avoid computing all ANI comparisons. Set to 1 to deactivate this check. Default: 0.3 (i.e., only compute ANI comparisons for genomes with %%id >= 70 [JI >= 0.001])')
    parser_dist.add_argument('-n', '--n-jobs', type=int,
        default=0,
        help='number of processes to spawn or threads to use. Default: 0 (use all available CPUs; for more info see os.sched_getaffinity())')
    parser_dist.add_argument('--chunksize', type=int,
        default=3,
        help='chunksize of parallel tasks. Default: 3 (for more info see multiprocessing.Pool().imap_unordered())')
    parser_dist.add_argument('-q', '--quiet', action='store_true',
        help='do not report the progress on stderr (warnings and errors are still reported)')
    parser_dist.add_argument('--force', action='store_true',
        help='overwrite the output file if it already exists')
    parser_dist.set_defaults(func=dist)

    return parser.parse_args()


def get_available_cpus() -> int:
    """
    Return the number of CPUs available to this process.

    On Linux, respect CPU affinity (cgroups, sched_setaffinity, etc.).
    On other systems, fallback to os.cpu_count().
    Always return at least 1.
    """
    # Linux with sched_getaffinity
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except (OSError, AttributeError):
            pass

    # Portable fallback
    count = os.cpu_count()
    return count if (count is not None and count > 0) else 1


def progress_reporter(label, total, quiet):
    # Report the progress of the parallel tasks on stderr every 10% of them, starting with
    # 0% so the start of each stage and its number of tasks are shown
    step = max(1, math.ceil(total / 10))

    def report(done):
        if (not quiet) and ((done % step == 0) or (done == total)):
            print(f"{prog}: {label}: {done}/{total} ({100 * done // total}%)", file=sys.stderr, flush=True)

    if total > 0:
        report(0)
    return report


def validate_arguments(args):
    available_cpus = get_available_cpus()

    # Common arguments for makedb and dist subcommands
    if args.subcommand in ('makedb', 'mk', 'dist', 'di'):
        if args.n_jobs < 0:
            print(f"{prog}: error: number of spawned processes or threads ({args.n_jobs}) must be >= 0", file=sys.stderr)
            sys.exit(2)

        if args.n_jobs == 0:
            args.n_jobs = available_cpus

        if args.n_jobs > available_cpus:
            print(f"{prog}: warning: required number of spawned processes or threads ({args.n_jobs}) is greater than the available CPUs ({available_cpus})", file=sys.stderr)
            print(f"{prog}: warning cont.: program execution will continue, but performance may degrade", file=sys.stderr)

        if args.chunksize <= 0:
            print(f"{prog}: error: chunksize ({args.chunksize}) must be > 0", file=sys.stderr)
            sys.exit(2)

    # makedb specific arguments
    if args.subcommand in ('makedb', 'mk'):
        if args.window <= 0:
            print(f"{prog}: error: window size ({args.window}) must be >= 1bp", file=sys.stderr)
            sys.exit(2)

        if (args.step <= 0) or (args.step > args.window):
            print(f"{prog}: error: step size ({args.step}) must be in the [1, {args.window}) range, being {args.window} the window size", file=sys.stderr)
            sys.exit(2)

        if not os.path.isfile(args.input):
            print(f"{prog}: error: {args.input} is not a valid file", file=sys.stderr)
            sys.exit(2)

        if os.path.basename(os.path.abspath(args.dbdir)) == '':
            print(f"{prog}: error: {args.dbdir} is not a valid database directory", file=sys.stderr)
            sys.exit(2)

        if (args.db_name is not None) and (os.sep in args.db_name):
            print(f"{prog}: error: {args.db_name} should not contain folder separator characters ({os.sep})", file=sys.stderr)
            sys.exit(2)

        if args.ji_sketchsize64 < 0:
            print(f"{prog}: error: number of sketches ({args.ji_sketchsize64}) must be >= 0", file=sys.stderr)
            sys.exit(2)

        # The other JI options are only used when sketches are computed
        if args.ji_sketchsize64 != 0:
            if args.ji_kmerlen <= 0:
                print(f"{prog}: error: kmer length ({args.ji_kmerlen}) must be >= 1bp", file=sys.stderr)
                sys.exit(2)

            if args.ji_minhashtype != 2:
                args.ji_dens = None

            # TODO: Check the range of posible values of bbits
            if args.ji_bbits <= 0:
                print(f"{prog}: error: number of bbits ({args.ji_bbits}) must be > 0", file=sys.stderr)
                sys.exit(2)

    # dist specific arguments
    if args.subcommand in ('dist', 'di'):
        if not os.path.isfile(args.query_db):
            print(f"{prog}: error: {args.query_db} is not a valid file", file=sys.stderr)
            sys.exit(2)

        if (args.reference_db is not None) and (not os.path.isfile(args.reference_db)):
            print(f"{prog}: error: {args.reference_db} is not a valid file", file=sys.stderr)
            sys.exit(2)

        if (args.output != '-') and os.path.exists(args.output) and (not args.force):
            print(f"{prog}: error: output file {args.output} already exists. Use --force to overwrite it", file=sys.stderr)
            sys.exit(2)

        if (args.fragment_id <= 0) or (args.fragment_id > 1):
            print(f"{prog}: error: fragment id percentage ({args.fragment_id}) must be in the (0, 1] range", file=sys.stderr)
            sys.exit(2)

        if args.fragment_len <= 0:
            print(f"{prog}: error: fragment length ({args.fragment_len}) must be > 0", file=sys.stderr)
            sys.exit(2)

        if args.coverage <= 0:
            print(f"{prog}: error: small genome coverage ({args.coverage}) must be > 0", file=sys.stderr)
            sys.exit(2)

        if (args.size_ratio != 0) and (args.size_ratio < 1):
            print(f"{prog}: error: ratio of genome sizes ({args.size_ratio}) must be 0 or >= 1", file=sys.stderr)
            sys.exit(2)

        if (args.mut_threshold < 0) or (args.mut_threshold > 1):
            print(f"{prog}: error: mutational distance threshold ({args.mut_threshold}) must be in the [0, 1] range", file=sys.stderr)
            sys.exit(2)


def is_gzip_file(fname):
    # Look for the magic number of Gzip files
    with open(fname, 'rb') as fh:
        # return fh.read(2) == b'\x1f\x8b'
        return fh.read(3) == b'\x1f\x8b\x08'


def is_fasta_file(fname):
    if is_gzip_file(fname):
        fh = gzip.open(fname, 'rt')
    else:
        fh = open(fname, 'rt')
    # TODO: Check if it is a nucleotide fasta file, not only a fasta file. It seems overkill though
    char = fh.read(1)
    fh.close()

    return char == '>'


def get_sequence_id_from_filename(fname):
    seq_id = os.path.basename(fname)

    exts = ['.gz']
    for ext in exts:
        if seq_id.endswith(ext):
            seq_id = seq_id[0:len(seq_id)-len(ext)]
    exts = ['.fna', '.fa', '.fasta']
    for ext in exts:
        if seq_id.endswith(ext):
            seq_id = seq_id[0:len(seq_id)-len(ext)]
            break

    return seq_id


def get_sequence_id_from_fasta_header(fname):
    if is_gzip_file(fname):
        fh = gzip.open(fname, 'rt')
    else:
        fh = open(fname, 'rt')
    line = fh.readline()
    fh.close()

    seq_id = line.strip().split(' ')[0][1:]
    return seq_id


def get_sequence_id(fname, seq_id_source):
    if seq_id_source == 'header':
        return get_sequence_id_from_fasta_header(fname)
    else:
        return get_sequence_id_from_filename(fname)


def get_number_of_contigs(fname):
    contigs = 0

    if is_gzip_file(fname):
        fh = gzip.open(fname, 'rt')
    else:
        fh = open(fname, 'rt')
    for line in fh:
        if line.startswith('>'):
            contigs += 1
    fh.close()

    return contigs


def get_databases_from_json(fname):
    with open(fname, 'rt') as fh:
        db_desc = json.load(fh)

    if db_desc.get('format_version', 1) != DB_FORMAT_VERSION:
        print(f"{prog}: error: {fname} was created with a different psani version (database format {db_desc.get('format_version', 1)}, expected {DB_FORMAT_VERSION}). Create the database again with psani makedb", file=sys.stderr)
        sys.exit(1)
    missing = [key for key in DB_REQUIRED if key not in db_desc]
    if missing:
        print(f"{prog}: error: {fname} is not a valid psani database (missing: {', '.join(missing)})", file=sys.stderr)
        sys.exit(1)

    # Relative paths are relative to the database directory, where the JSON is. Absolute
    # paths (makedb --absolute-paths) are kept as they are by os.path.join
    db_dir = os.path.dirname(os.path.abspath(fname))
    db_desc['input_file'] = os.path.join(db_dir, db_desc['input_file'])
    if db_desc['ji_sketch_db'] != '':
        db_desc['ji_sketch_db'] = os.path.join(db_dir, db_desc['ji_sketch_db'])
    for desc in db_desc['seqs'].values():
        for key in ('fna', 'fwa'):
            if desc[key] != '':
                desc[key] = os.path.join(db_dir, desc[key])

    return db_desc


def run_bindash_sketch(fofn, sketch_db, sketches_64, kmer_len, minhashtype, dens, bbits, randseed, n_jobs):
    # TODO: Don't wait for BinDash to finish. Write its output log in real time
    cmd = ['bindash', 'sketch', '--listfname='+fofn, '--outfname='+sketch_db, '--kmerlen='+str(kmer_len), '--sketchsize64='+str(sketches_64),'--minhashtype='+str(minhashtype), '--bbits='+str(bbits), '--randseed='+str(randseed), '--nthreads='+str(n_jobs)]
    # The densification strategy only applies to minhashtype 2 (dens is None otherwise)
    if dens is not None:
        cmd.append('--dens='+str(dens))
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{prog}: error: {sketch_db} sketch database could not be created", file=sys.stderr)
        print(exc, file=sys.stderr)
        if (exc.stdout != ''):
            print(exc.stdout, file=sys.stdout)
        if (exc.stderr != ''):
            print(exc.stderr, file=sys.stderr)
        sys.exit(1)
    else:
        if (cp.stdout != ''):
            print(cp.stdout, file=sys.stdout)


def run_bindash_dist(query_sketches, target_sketches, mut_thres, n_jobs):
    if target_sketches == '':
        cmd = ['bindash', 'dist', '--mthres='+str(mut_thres), '--nthreads='+str(n_jobs), query_sketches]
    else:
        cmd = ['bindash', 'dist', '--mthres='+str(mut_thres), '--nthreads='+str(n_jobs), query_sketches, target_sketches]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{prog}: error: mutational distances could not be computed using the JI", file=sys.stderr)
        print(exc, file=sys.stderr)
        if (exc.stdout != ''):
            print(exc.stdout, file=sys.stdout)
        if (exc.stderr != ''):
            print(exc.stderr, file=sys.stderr)
        sys.exit(1)
    else:
        output = cp.stdout.split("\n")
        if output[-1] == '':
            del output[-1]

        # Add output for self-comparisons when target_sketches is not provided
        # Format: id {TAB} id {TAB} 0.0000e+00 {TAB} 0.0000e+00 {TAB} 2048/2048
        if target_sketches == '':
            with open(query_sketches+'.txt', 'rt') as fh:
                next(fh)
                for line in fh:
                    (id, _) = line.split("\t", 1)
                    # TODO: Get number of scketches as a parameter and use it to replace the JI field
                    output.append(f"{id}\t{id}\t0.0000e+00\t0.0000e+00\t1.0")

        return output


def run_makeblastdb(fname):
    cmd = ['makeblastdb', '-in', fname, '-dbtype', 'nucl']
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{prog}: error: {fname} BLAST database could not be created", file=sys.stderr)
        print(exc, file=sys.stderr)
        if (exc.stdout != ''):
            print(exc.stdout, file=sys.stdout)
        if (exc.stderr != ''):
            print(exc.stderr, file=sys.stderr)
        # This runs in a worker process of the pool, where sys.exit() would leave the
        # main process waiting forever. The exception is raised again in the main process
        raise RuntimeError(f"{fname} BLAST database could not be created")


def run_blast(qry, ref_db):
    # Return the BLASTn output lines, or None and the error message if BLASTn fails
    cmd = ['blastn', '-query', qry, '-db', ref_db, '-dust', 'no', '-max_target_seqs', '1', '-outfmt', '6']
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        # Keep the error lines, without the BLASTn warnings (e.g. the one about
        # -max_target_seqs 1, always printed), in a single line
        lines = [l.strip() for l in (getattr(exc, 'stderr', None) or str(exc)).splitlines() if l.strip()]
        errors = [l for l in lines if not l.startswith('Warning:')] or lines
        return (None, '; '.join(errors))
    else:
        return (cp.stdout.split("\n"), '')


def makedb_task(params):
    fname = params[0]
    id = params[1]
    outdir = params[2]
    window = params[3]
    step = params[4]

    if is_gzip_file(fname):
        fh = gzip.open(fname, 'rt')
    else:
        fh = open(fname, 'rt')

    # Based on ani.rb from the Enveomics collection suite
    # Create fasta file of fragments using a sliding window algorithm. Discard last
    # fragment if shorter than window size. Discard all fragments composed only of Ns
    size = 0
    contigs = -1
    #disc = 0
    buffer = ''
    pos = 0
    idx = -1
    fwa = os.path.join(outdir, id) + '.fwa'
    fh_fwa = None
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            #contig = line.split()[0][1:]
            contigs += 1
            #disc += len(buffer)
            buffer = ''
            pos = 0
        else:
            line = regex_notAlpha.sub('', line)
            size += len(line)
            buffer += line

            while(len(buffer) >= window):
                frg = buffer[0:window]
                if not regex_allN.search(frg):
                    idx += 1
                    if not fh_fwa:
                        # TODO: Check how to make makeblastdb work with compressed files
                        fh_fwa = open(fwa, 'wt')
                    fh_fwa.write(f">{id}~{idx} {contigs} {pos}\n{frg}\n")
                #else:
                #    disc += len(frg)
                buffer = buffer[step:]
                pos += step
    if fh_fwa:
        fh_fwa.close()
    fh.close()

    # If no fragments were produced (e.g., all contigs shorter than window),
    # do not create an empty/missing BLAST database.
    if idx < 0:
        return (id, '', size, contigs+1, 0)

    # Create the BLAST database of sliding window fragments
    run_makeblastdb(fwa)

    return (id, fwa, size, contigs+1, idx+1)


def makedb(args):
    dbdir = os.path.normpath(args.dbdir)

    if args.db_name:
        db_name = args.db_name
    else:
        db_name = os.path.basename(os.path.abspath(dbdir))

    # Do not overwrite a previous database unless requested. Only the files created by
    # makedb are removed, as the database directory may contain other files
    outputs = [os.path.join(dbdir, db_name + '.json'), os.path.join(dbdir, 'fwa'), os.path.join(dbdir, 'kmer')]
    existing = [path for path in outputs if os.path.lexists(path)]
    if existing and not args.force:
        print(f"{prog}: error: {', '.join(existing)} already exist(s). Use --force to overwrite", file=sys.stderr)
        sys.exit(2)
    for path in existing:
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.remove(path)

    db_desc = {'format_version': DB_FORMAT_VERSION, 'db_name': db_name, 'input_file': args.input, 'window': args.window, 'step': args.step, 'ji_sketch_db': '', 'ji_sketchsize64': 0, 'ji_kmerlen': 0, 'ji_minhashtype': 0, 'ji_dens': 0, 'ji_bbits': 0, 'ji_randseed': 0, 'seqs':{}}

    # Verify that the query sequences are FASTA files and collect their metadata
    if is_fasta_file(args.input):
        fname = args.input
        id = get_sequence_id(fname, args.seq_id_source)
        db_desc['seqs'][id] = {'fna': fname, 'fwa': '', 'size': 0, 'contigs': 0, 'fragments': 0}
    else:
        with open(args.input, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                tmp = line.split("\t")
                if len(tmp) >= 2:
                    # Assume the FOFN format is: fname - ID - ...
                    fname = tmp[0].strip()
                    id = tmp[1].strip()
                    if os.path.isfile(fname) and is_fasta_file(fname):
                        db_desc['seqs'][id] = {'fna': fname, 'fwa': '', 'size': 0, 'contigs': 0, 'fragments': 0}
                elif len(tmp) == 1:
                    fname = tmp[0].strip()
                    id = get_sequence_id(fname, args.seq_id_source)
                    if os.path.isfile(fname) and is_fasta_file(fname):
                        db_desc['seqs'][id] = {'fna': fname, 'fwa': '', 'size': 0, 'contigs': 0, 'fragments': 0}
                else:
                    continue
    if len(db_desc['seqs']) == 0:
        print(f"{prog}: error: no valid FASTA file(s) found as input", file=sys.stderr)
        sys.exit(1)

    # Create output folders
    pathlib.Path(dbdir).mkdir(parents=True, exist_ok=True)
    dbdir_fwa = os.path.join(dbdir, 'fwa')
    pathlib.Path(dbdir_fwa).mkdir(parents=True, exist_ok=True)

    # If required create the sketch database to compute the JI shortcut
    if args.ji_sketchsize64 != 0:
        dbdir_kmer = os.path.join(dbdir, 'kmer')
        pathlib.Path(dbdir_kmer).mkdir(parents=True, exist_ok=True)

        # Create temporary FOFN with the correct ID of the sequences for BinDash
        fh_tmp = tempfile.NamedTemporaryFile(mode='wt', delete=False)
        for id in db_desc['seqs']:
            # BinDash requires a FOFN with three fields: fname - ID - contigs
            db_desc['seqs'][id]['contigs'] = get_number_of_contigs(db_desc['seqs'][id]['fna'])
            fh_tmp.write(f"{db_desc['seqs'][id]['fna']}\t{id}\t{db_desc['seqs'][id]['contigs']}\n")
        fh_tmp.close()

        sketch_db = os.path.join(dbdir_kmer, db_name + '.bdsh')

        if not args.quiet:
            print(f"{prog}: computing sketches with BinDash", file=sys.stderr, flush=True)
        run_bindash_sketch(fh_tmp.name, sketch_db, args.ji_sketchsize64, args.ji_kmerlen, args.ji_minhashtype, args.ji_dens, args.ji_bbits, args.ji_randseed, args.n_jobs)
        os.remove(fh_tmp.name)
        db_desc['ji_sketch_db'] = sketch_db
        db_desc['ji_sketchsize64'] = args.ji_sketchsize64
        db_desc['ji_kmerlen'] = args.ji_kmerlen
        db_desc['ji_minhashtype'] = args.ji_minhashtype
        db_desc['ji_dens'] = args.ji_dens
        db_desc['ji_bbits'] = args.ji_bbits
        db_desc['ji_randseed'] = args.ji_randseed

    # Parallel database creation
    with multiprocessing.Pool(args.n_jobs) as pool:
        task_params = []
        for id, desc in db_desc['seqs'].items():
            task_params.append((desc['fna'], id, dbdir_fwa, args.window, args.step))
        tasks = pool.imap_unordered(makedb_task, task_params, chunksize=args.chunksize)
        report = progress_reporter('fragment databases', len(task_params), args.quiet)
        done = 0
        while True:
            # Abort without writing the database JSON if a fragment database fails
            try:
                task_output = next(tasks)
            except StopIteration:
                break
            except RuntimeError as exc:
                print(f"{prog}: error: {exc}. Database not created", file=sys.stderr)
                sys.exit(1)
            done += 1
            report(done)
            if task_output != None:
                seq_id = task_output[0]
                db_desc['seqs'][seq_id]['fwa'] = task_output[1]
                db_desc['seqs'][seq_id]['size'] = task_output[2]
                db_desc['seqs'][seq_id]['contigs'] = task_output[3]
                db_desc['seqs'][seq_id]['fragments'] = task_output[4]

                if db_desc['seqs'][seq_id]['fragments'] == 0:
                    print(
                        f"{prog}: warning: sequence '{seq_id}' produced 0 fragments "
                        f"(all contigs shorter than window={args.window} or only-N fragments). "
                        "It will be skipped in ANI comparisons.",
                        file=sys.stderr
                    )

    # Store the paths relative to the database directory, so the database can be moved,
    # or absolute if requested
    store_path = os.path.abspath if args.absolute_paths else (lambda path: os.path.relpath(path, dbdir))
    db_desc['input_file'] = store_path(db_desc['input_file'])
    if db_desc['ji_sketch_db'] != '':
        db_desc['ji_sketch_db'] = store_path(db_desc['ji_sketch_db'])
    for desc in db_desc['seqs'].values():
        for key in ('fna', 'fwa'):
            if desc[key] != '':
                desc[key] = store_path(desc[key])

    # TODO: Try to embed this in the last for loop, writing each desc when the task_ouput is ready. This way it could serve as a log of the work already done
    db_json = os.path.join(dbdir, db_name + '.json')
    with open(db_json, 'wt') as fh_json:
        fh_json.write(json.dumps(db_desc, ensure_ascii=False, indent=2))


def two_way_ani(blast_src_trg, blast_trg_src, fragment_pident, fragment_len, correct, min_fragments):
    # Based on ani.rb from the Enveomics collection suite
    # For every pair of windowed genomes perform two BLASTn comparisons, one of the
    # first genome against the second and the other for the reverse comparison
    # Not sure about the biologic reason behind the identity correction, just leave it as is

    rbh = {}
    identities = []
    for (i, blast) in ((1, blast_src_trg), (2, blast_trg_src)):
        src_seen = set()
        for line in blast:
            if line != '':
                # BLASTn output format 6: qseqid - sseqid - pident - length - ...
                row = line.split("\t", 4)
                if ((row[0] not in src_seen) and (float(row[2]) >= fragment_pident) and (int(row[3]) >= fragment_len)):
                    src_seen.add(row[0])
                    identity_corr = 100 - (100 - float(row[2]))/(0.8621 if correct else 1.0)
                    if (i == 1):
                        # Save the best hits on the first iterarion
                        rbh[row[0]] = row[1]
                    else:
                        # Keep the identities on the second iteration only for reciprocal best hits
                        if ((row[1] in rbh) and (rbh[row[1]] == row[0])):
                            identities.append(identity_corr)
    n2 = len(identities)
    if n2 >= min_fragments:
        ani = sum(identities) / n2
        # Two-pass variance: the mean of the squares minus the square of the mean loses
        # precision with identities close to 100, and gave a tiny non-zero (or negative)
        # variance when all identities are equal
        ani_var = sum((x - ani) * (x - ani) for x in identities) / n2
        return (ani, ani_var, n2)
    else:
        return (None, None, None)


def fragments_in_length(length, window, step):
    # Number of sliding-window fragments of a sequence of this length (at least 1)
    if length <= window:
        return 1
    return 1 + int((length - window) // step)


def dist_task(params):
    qry_path = params[0]
    ref_path = params[1]
    qry_id = params[2]
    ref_id = params[3]
    qry_size = params[4]
    ref_size = params[5]
    qry_fragments = params[6]
    ref_fragments = params[7]
    window = params[8]
    step = params[9]
    fragment_pident = params[10]
    fragment_len = params[11]
    correct = params[12]
    coverage = params[13]
    extra_output = params[14]

    # Compute the minimum number of fragments required for this particular ANI comparison,
    # from the length of the smaller sequence (including its N)
    if coverage >= 1:
        min_fragments = int(coverage)
    else:
        min_fragments = fragments_in_length(coverage * min(qry_size, ref_size), window, step)
    if min_fragments < 1:
        print(f"{prog}: warning: minimum number of shared fragments computed for {qry_id}:{ref_id} comparison is < 1. Value changed to 1", file=sys.stderr)
        min_fragments = 1

    # Compute both BLASTn needed for this particular ANI comparison. Both are needed for the
    # reciprocal best hits, so if one fails the pair is not reported (None tells dist() to
    # count it as failed)
    (blast_qry_ref, error_qry_ref) = run_blast(qry_path, ref_path)
    (blast_ref_qry, error_ref_qry) = run_blast(ref_path, qry_path)
    if (blast_qry_ref is None) or (blast_ref_qry is None):
        error = error_qry_ref or error_ref_qry
        print(f"{prog}: warning: BLASTn failed for the {qry_id}:{ref_id} comparison, which is not "
              f"reported: {error}", file=sys.stderr)
        return None
    (ani_fwd, ani_var_fwd, n2_fwd) = two_way_ani(blast_qry_ref, blast_ref_qry, fragment_pident, fragment_len, correct, min_fragments)
    (ani_rev, ani_var_rev, n2_rev) = two_way_ani(blast_ref_qry, blast_qry_ref, fragment_pident, fragment_len, correct, min_fragments)

    # Both directions count the same reciprocal best hits (fragment pairs in which each one
    # is the best hit of the other), so they reach the minimum number of fragments together
    # and give the same number; only the identities, taken from one BLASTn or the other,
    # differ
    if ani_fwd is None:
        # No homology detected (empty output, not enough shared fragments)
        return ''
    n2 = n2_fwd
    # Sample mean and variance of the combined samples (Ref.: https://math.stackexchange.com/questions/2971315)
    ani = ((n2_fwd * ani_fwd) + (n2_rev * ani_rev)) / (n2_fwd + n2_rev)
    ani_var = ((((n2_fwd - 1) * ani_var_fwd) + ((n2_rev - 1) * ani_var_rev)) / (n2_fwd + n2_rev - 1)) + \
              ((n2_fwd * n2_rev * (ani_fwd - ani_rev) * (ani_fwd - ani_rev)) / ((n2_fwd + n2_rev) * (n2_fwd + n2_rev - 1)))
    ani_sd = math.sqrt(ani_var)
    # Coverage of the smaller sequence, computed as the threshold above, so that it can be
    # compared with -c/--coverage
    cov = n2 / fragments_in_length(min(qry_size, ref_size), window, step)
    # Values rounded to OUTPUT_DECIMALS decimals
    ani = f"{ani:.{OUTPUT_DECIMALS}f}"
    ani_sd = f"{ani_sd:.{OUTPUT_DECIMALS}f}"
    cov = f"{cov:.{OUTPUT_DECIMALS}f}"
    if extra_output == '':
        return f"{qry_id}\t{ref_id}\t{qry_size}\t{ref_size}\t{qry_fragments}\t{ref_fragments}\t{min_fragments}\t{n2}\t{ani}\t{ani_sd}\t{cov}\n"
    else:
        return f"{qry_id}\t{ref_id}\t{qry_size}\t{ref_size}\t{qry_fragments}\t{ref_fragments}\t{min_fragments}\t{n2}\t{ani}\t{ani_sd}\t{cov}\t{extra_output}\n"


def dist(args):
    query_db = get_databases_from_json(args.query_db)
    if len(query_db['seqs']) == 0:
        print(f"{prog}: error: no valid query databases found as input", file=sys.stderr)
        sys.exit(1)

    if args.reference_db is not None:
        reference_db = get_databases_from_json(args.reference_db)
        if len(reference_db['seqs']) == 0:
            print(f"{prog}: error: no valid reference databases found as input", file=sys.stderr)
            sys.exit(1)
        if query_db['window'] != reference_db['window']:
            print(f"{prog}: error: query and reference databases were built using different window parameters", file=sys.stderr)
            sys.exit(1)
        if query_db['step'] != reference_db['step']:
            print(f"{prog}: error: query and reference databases were built using different step parameters", file=sys.stderr)
            sys.exit(1)
    else:
        reference_db = query_db

    # The fragments are those of the databases (the same window and step in both)
    window = query_db['window']
    step = query_db['step']
    if args.fragment_len > 1:
        if args.fragment_len > window:
            print(f"{prog}: error: fragment length ({args.fragment_len}) must be <= {window}, being {window} the window size of the databases", file=sys.stderr)
            sys.exit(2)
        fragment_len = math.ceil(args.fragment_len)
    else:
        fragment_len = math.ceil(window * args.fragment_len)
    # Minimum identity in percent, the scale of the BLAST pident column
    fragment_pident = 100 * args.fragment_id

    # Compute JI similarities for JI shortcut
    if (query_db['ji_sketch_db'] != '') and (reference_db['ji_sketch_db'] != '') and (args.mut_threshold < 1):
        ji_shortcut = True
        if query_db['ji_sketchsize64'] != reference_db['ji_sketchsize64']:
            print(f"{prog}: error: query and reference sketch databases were built using different number of sketches", file=sys.stderr)
            sys.exit(1)
        if query_db['ji_kmerlen'] != reference_db['ji_kmerlen']:
            print(f"{prog}: error: query and reference sketch databases were built using different kmer lengths", file=sys.stderr)
            sys.exit(1)
        if query_db['ji_bbits'] != reference_db['ji_bbits']:
            print(f"{prog}: error: query and reference sketch databases were built using different number of bbits", file=sys.stderr)
            sys.exit(1)
        if query_db['ji_minhashtype'] != reference_db['ji_minhashtype']:
            print(f"{prog}: error: query and reference sketch databases were built using different minhash types", file=sys.stderr)
            sys.exit(1)
        if query_db['ji_dens'] != reference_db['ji_dens']:
            print(f"{prog}: error: query and reference sketch databases were built using different densification strategies", file=sys.stderr)
            sys.exit(1)
        if query_db['ji_randseed'] != reference_db['ji_randseed']:
            print(f"{prog}: error: query and reference sketch databases were built using different random seeds", file=sys.stderr)
            sys.exit(1)

        if not args.quiet:
            print(f"{prog}: computing Jaccard distances with BinDash", file=sys.stderr, flush=True)
        if args.reference_db is None:
            ji_comparisons = run_bindash_dist(query_db['ji_sketch_db'], '', args.mut_threshold, args.n_jobs)
        else:
            ji_comparisons = run_bindash_dist(query_db['ji_sketch_db'], reference_db['ji_sketch_db'], args.mut_threshold, args.n_jobs)
    else:
        ji_shortcut = False

    def has_fragments(db, seq_id: str) -> bool:
        return (
            seq_id in db['seqs']
            and db['seqs'][seq_id].get('fragments', 0) > 0
            and db['seqs'][seq_id].get('fwa', '') != ''
        )

    # Open the output file only once the databases have been checked
    try:
        output = sys.stdout if args.output == '-' else open(args.output, 'wt')
    except OSError as exc:
        print(f"{prog}: error: output file {args.output} could not be opened: {exc}", file=sys.stderr)
        sys.exit(1)

    # Parallel computation of pairwise ANI comparisons
    with multiprocessing.Pool(args.n_jobs) as pool:
        task_params = []
        if args.header:
            if ji_shortcut:
                output.write("#Qry_ID\tRef_ID\tQry_Size\tRef_Size\tQry_Fragments\tRef_Fragments\tMin_Fragments\tRBH_Fragments\tANI\tANI_SD\tMin_Smaller_Genome_Fragments_Coverage\tBD_Mut_Distance\tBD_P-value\tBD_JI\n")
            else:
                output.write("#Qry_ID\tRef_ID\tQry_Size\tRef_Size\tQry_Fragments\tRef_Fragments\tMin_Fragments\tRBH_Fragments\tANI\tANI_SD\tMin_Smaller_Genome_Fragments_Coverage\n")
        if ji_shortcut:
            # Only perform the comparisons indicated by BinDash
            for comp in ji_comparisons:
                (qry_id, ref_id, ji_output) = comp.split("\t", 2)
                if not has_fragments(query_db, qry_id):
                    continue
                if not has_fragments(reference_db, ref_id):
                    continue
                if (args.size_ratio == 0) or ((query_db['seqs'][qry_id]['size'] <= args.size_ratio * reference_db['seqs'][ref_id]['size']) and (reference_db['seqs'][ref_id]['size'] <= args.size_ratio * query_db['seqs'][qry_id]['size'])):
                    task_params.append((query_db['seqs'][qry_id]['fwa'], reference_db['seqs'][ref_id]['fwa'], \
                                        qry_id, ref_id, \
                                        query_db['seqs'][qry_id]['size'], reference_db['seqs'][ref_id]['size'], \
                                        query_db['seqs'][qry_id]['fragments'], reference_db['seqs'][ref_id]['fragments'], \
                                        window, step, fragment_pident, fragment_len, args.correct_ani, args.coverage, \
                                        ji_output))
        elif args.reference_db is not None:
            # Perform all comparisons of the matrix
            for qry_id in query_db['seqs']:
                if not has_fragments(query_db, qry_id):
                    continue
                for ref_id in reference_db['seqs']:
                    if not has_fragments(reference_db, ref_id):
                        continue
                    if (args.size_ratio == 0) or ((query_db['seqs'][qry_id]['size'] <= args.size_ratio * reference_db['seqs'][ref_id]['size']) and (reference_db['seqs'][ref_id]['size'] <= args.size_ratio * query_db['seqs'][qry_id]['size'])):
                        task_params.append((query_db['seqs'][qry_id]['fwa'], reference_db['seqs'][ref_id]['fwa'], \
                                            qry_id, ref_id, \
                                            query_db['seqs'][qry_id]['size'], reference_db['seqs'][ref_id]['size'], \
                                            query_db['seqs'][qry_id]['fragments'], reference_db['seqs'][ref_id]['fragments'], \
                                            window, step, fragment_pident, fragment_len, args.correct_ani, args.coverage, \
                                            ''))
        else:
            # Do not perform redundant comparisons if no reference databases are provided
            for i, qry_id in enumerate(query_db['seqs']):
                if not has_fragments(query_db, qry_id):
                    continue
                for ref_id in list(query_db['seqs'])[i:]:
                    if not has_fragments(query_db, ref_id):
                        continue
                    if (args.size_ratio == 0) or ((query_db['seqs'][qry_id]['size'] <= args.size_ratio * reference_db['seqs'][ref_id]['size']) and (reference_db['seqs'][ref_id]['size'] <= args.size_ratio * query_db['seqs'][qry_id]['size'])):
                        task_params.append((query_db['seqs'][qry_id]['fwa'], reference_db['seqs'][ref_id]['fwa'], \
                                            qry_id, ref_id, \
                                            query_db['seqs'][qry_id]['size'], reference_db['seqs'][ref_id]['size'], \
                                            query_db['seqs'][qry_id]['fragments'], reference_db['seqs'][ref_id]['fragments'], \
                                            window, step, fragment_pident, fragment_len, args.correct_ani, args.coverage, \
                                            ''))
        report = progress_reporter('ANI comparisons', len(task_params), args.quiet)
        failed = 0
        for done, task_output in enumerate(pool.imap_unordered(dist_task, task_params, chunksize=args.chunksize), 1):
            report(done)
            if task_output is None:
                failed += 1
            elif task_output != '':
                output.write(str(task_output))
                # TODO: Check if we need to flush the output
                output.flush()

    if output is not sys.stdout:
        output.close()

    if failed > 0:
        print(f"{prog}: warning: {failed} of {len(task_params)} comparisons were not reported "
              "because BLASTn failed (see the warnings above)", file=sys.stderr)


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    args.func(args)


if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
