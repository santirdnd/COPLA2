#!/usr/bin/env python3

"""
Prepare a COPLA 2 reference database for local use.

The database directory (e.g. databases/Copla_RS224_v1.0) must contain a manifest.json
(format version 1) describing the database files and the parameters used to build
it. This script performs the following steps:

  - Optionally downloads the files listed in the manifest that are missing, or whose
    MD5 checksum does not match the manifest (e.g. files of a previous version). Files
    marked as optional in the manifest (not needed by the COPLA2 scripts, such as the
    ANI table of the network) are downloaded only with --download-optional.

  - Verifies the MD5 checksum of each file (skipped with a warning when the
    manifest does not provide one).

  - Extracts into fna/ the per-plasmid compressed FASTA files (fna/<ID>.fna.gz)
    from the sequences ZIP archive, keeping only the curated plasmids listed in
    the curated file (one ID per line). The list is kept out of the archive, so a
    change in the curation does not change the archive.

  - Runs "psani.py makedb" with the parameters of the "makedb" section of the
    manifest "psani" entry to create the sequence fragment (fwa/) and sketch (kmer/)
    databases and the database JSON used by "psani.py dist".

  - Checks that the database contains every curated plasmid and removes the
    per-plasmid FASTA files, which are no longer needed after makedb.

Paths written to the database JSON are relative to the database directory, so the
database can be used from any working directory and moved without rebuilding it. With
--absolute-paths they are absolute; if the database directory is then moved, run this
script again.

The "dist" section of the "psani" entry (the parameters of the ANI between the
reference plasmids) and the "copla2" entry (the size ratio of the network links)
document how the reference network was built; they are not used by this script.
"""

import sys
import json
import shutil
import hashlib
import argparse
import subprocess
import zipfile
import urllib.request
from pathlib import Path, PurePosixPath


# Exit codes
EXIT_ERROR = 1
EXIT_INVALID_MANIFEST = 2
EXIT_INPUT_FILE = 3
EXIT_OUTPUT_EXISTS = 4

# Format version of the manifest.json read by this script (and by copla2.py)
MANIFEST_FORMAT_VERSION = 1

# Seconds without data from the server before a download is aborted
DOWNLOAD_TIMEOUT = 60

# Manifest "psani" "makedb" keys passed to "psani.py makedb", mapped to its arguments.
# Other keys are ignored.
PSANI_ARGS = {
    "db_name": "--db-name",
    "window": "--window",
    "step": "--step",
    "ji_sketchsize64": "--ji-sketchsize64",
    "ji_kmerlen": "--ji-kmerlen",
    "ji_minhashtype": "--ji-minhashtype",
    "ji_dens": "--ji-dens",
    "ji_bbits": "--ji-bbits",
    "ji_randseed": "--ji-randseed",
}

BIN_DIR = Path(__file__).resolve().parent
PROG = Path(sys.argv[0]).name


def die(msg: str, code: int = EXIT_ERROR) -> None:
    print(f"{PROG}: error: {msg}", file=sys.stderr)
    sys.exit(code)


def warn(msg: str) -> None:
    print(f"{PROG}: warning: {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"{PROG}: {msg}", file=sys.stderr)


def md5sum(path: Path) -> str:
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
    return md5.hexdigest()


def download(url: str, path: Path) -> None:
    # Download to a temporary file so an interrupted download is not taken as valid
    tmp_path = path.with_name(path.name + ".part")
    info(f"downloading {url}")
    try:
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT) as response, \
                open(tmp_path, "wb") as fh:
            # Report the progress every 10% when the server gives the file size
            total = int(response.headers.get("Content-Length") or 0)
            done = 0
            next_report = 10
            for chunk in iter(lambda: response.read(1 << 20), b""):
                fh.write(chunk)
                done += len(chunk)
                if total and (100 * done // total >= next_report):
                    info(f"  {done / 1e6:.0f}/{total / 1e6:.0f} MB ({100 * done // total}%)")
                    next_report = 100 * done // total // 10 * 10 + 10
    except OSError as exc:
        tmp_path.unlink(missing_ok=True)
        die(f"{url} could not be downloaded: {exc}", EXIT_INPUT_FILE)
    tmp_path.rename(path)


def run(cmd: list) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        die(f"command failed with exit code {exc.returncode}: {' '.join(cmd)}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a COPLA 2 reference database described by a manifest.json")
    parser.add_argument("dbdir", type=str,
        help="database directory containing manifest.json (e.g. databases/Copla_RS224_v1.0)")
    parser.add_argument("--download", action="store_true",
        help="download the files listed in the manifest that are missing or whose MD5 "
             "checksum does not match, except the optional ones")
    parser.add_argument("--download-optional", action="store_true",
        help="as --download, and also download the optional files (e.g. the ANI table of "
             "the reference network), which the COPLA2 scripts do not need")
    parser.add_argument("--absolute-paths", action="store_true",
        help="write absolute paths in the database JSON. Default: paths relative to the "
             "database directory, so it can be moved")
    parser.add_argument("--keep-fna", action="store_true",
        help="keep the per-plasmid FASTA files (fna/) after creating the database")
    parser.add_argument("-n", "--n-jobs", type=int, default=0,
        help="number of processes for psani makedb. Default: 0 (use all available CPUs)")
    parser.add_argument("--force", action="store_true",
        help="overwrite an existing database (database JSON, fwa/, kmer/ and fna/)")
    parser.add_argument("--version", action="version", version="%(prog)s 2.0.0")
    args = parser.parse_args()

    dbdir = Path(args.dbdir).resolve()
    manifest_path = dbdir / "manifest.json"
    if not manifest_path.is_file():
        die(f"{manifest_path} not found", EXIT_INPUT_FILE)
    try:
        with open(manifest_path, "rt", encoding="utf-8") as fh:
            manifest = json.load(fh)
        format_version = manifest.get("format_version")
    except (json.JSONDecodeError, AttributeError) as exc:
        die(f"invalid manifest {manifest_path}: {exc}", EXIT_INVALID_MANIFEST)
    # Check the version before the keys, whose layout depends on it
    if format_version != MANIFEST_FORMAT_VERSION:
        found = ("has no format version" if format_version is None
                 else f"has format version {format_version}")
        die(f"{manifest_path} {found}, but this script reads version "
            f"{MANIFEST_FORMAT_VERSION}. Use the manifest.json of this COPLA2 release",
            EXIT_INVALID_MANIFEST)
    try:
        files = manifest["files"]
        psani_params = manifest["psani"]["makedb"]
        db_name = psani_params["db_name"]
    except (KeyError, TypeError) as exc:
        die(f"invalid manifest {manifest_path}: {exc}", EXIT_INVALID_MANIFEST)
    # The sequences and the curated list build the psani database; the graph (the
    # reference network pickle) is needed by copla2.py
    for key in ("sequences", "curated", "graph"):
        if key not in files:
            die(f"manifest has no '{key}' file", EXIT_INVALID_MANIFEST)
        if files[key].get("optional", False):
            die(f"manifest file '{key}' cannot be optional", EXIT_INVALID_MANIFEST)
    # --download-optional also downloads the files that are not optional
    if args.download_optional:
        args.download = True

    info(f"database {manifest.get('id', '')} version {manifest.get('version', '')} in {dbdir}")

    # Do not overwrite a previous database unless requested
    db_json = dbdir / f"{db_name}.json"
    fna_dir = dbdir / "fna"
    outputs = [db_json, dbdir / "fwa", dbdir / "kmer", fna_dir]
    existing = [p for p in outputs if p.exists()]
    if existing and not args.force:
        die(f"{', '.join(str(p) for p in existing)} already exist(s). Use --force to overwrite",
            EXIT_OUTPUT_EXISTS)

    # Download (if requested) and verify the database files. With --download, a file whose
    # MD5 checksum does not match (e.g. from a previous version) is downloaded again.
    # Optional files (e.g. the ANI table of the network) are not used by the scripts: they
    # are downloaded only with --download-optional, but verified whenever present
    paths = {}
    for key, desc in files.items():
        path = dbdir / desc["path"]
        paths[key] = path
        url = desc.get("url")
        md5 = desc.get("md5")
        optional = bool(desc.get("optional", False))
        fetch = args.download_optional if optional else args.download
        downloaded = False
        if not path.is_file():
            if optional and not fetch:
                info(f"optional file {path.name} not present (use --download-optional to "
                     "download it)")
                continue
            if not fetch:
                die(f"{path} not found. Download it (see manifest.json) or use --download",
                    EXIT_INPUT_FILE)
            if not url:
                die(f"{path} not found and the manifest gives no URL for it: download it "
                    "manually (see manifest.json)", EXIT_INPUT_FILE)
            download(url, path)
            downloaded = True
        if not md5:
            warn(f"no MD5 checksum in the manifest for {path}; not verified")
            continue
        matches = (md5sum(path) == md5)
        if (not matches) and fetch and url and (not downloaded):
            warn(f"MD5 checksum of {path} does not match the manifest; downloading it again")
            download(url, path)
            matches = (md5sum(path) == md5)
        if not matches:
            option = "--download-optional" if optional else "--download"
            hint = (f"remove it and download it again ({option})" if not fetch
                    else "check the manifest URL or download it manually")
            die(f"MD5 checksum of {path} does not match the manifest. The file may be "
                f"incomplete, damaged or from another version: {hint}", EXIT_INPUT_FILE)

    # Remove the previous database only once the input files have been verified
    for p in existing:
        if p.is_dir():
            shutil.rmtree(p)
        else:
            p.unlink()

    # Extract the per-plasmid FASTA files of the curated plasmids.
    # Only the file name of each ZIP member is used, so nothing is written outside fna/
    with open(paths["curated"], "rt", encoding="utf-8") as fh:
        kept = set(fh.read().split())
    suffix = ".fna.gz"
    extracted = set()
    with zipfile.ZipFile(paths["sequences"]) as zf:
        info(f"extracting the {len(kept)} curated plasmids listed in {paths['curated'].name}")
        fna_dir.mkdir()
        for member in zf.infolist():
            fname = PurePosixPath(member.filename).name
            if member.is_dir() or not fname.endswith(suffix):
                continue
            seq_id = fname[:-len(suffix)]
            if seq_id not in kept:
                continue
            if seq_id in extracted:
                die(f"duplicated sequence '{seq_id}' in {paths['sequences']}", EXIT_INPUT_FILE)
            extracted.add(seq_id)
            with zf.open(member) as src, open(fna_dir / fname, "wb") as dst:
                shutil.copyfileobj(src, dst)
    not_found = kept - extracted
    if not_found:
        die(f"{len(not_found)} plasmids listed in {paths['curated'].name} have no sequence in "
            f"{paths['sequences']} (e.g. {', '.join(sorted(not_found)[:5])})", EXIT_INPUT_FILE)

    # FOFN with absolute paths and explicit IDs (fname - ID) for psani makedb
    fofn_path = dbdir / "fna.fofn"
    with open(fofn_path, "wt", encoding="utf-8") as fh:
        for seq_id in sorted(extracted):
            fh.write(f"{fna_dir / (seq_id + suffix)}\t{seq_id}\n")

    # Create the psani database
    cmd = [sys.executable, str(BIN_DIR / "psani.py"), "makedb", str(fofn_path), str(dbdir),
           "--n-jobs", str(args.n_jobs)]
    for key, value in psani_params.items():
        if key in PSANI_ARGS:
            cmd += [PSANI_ARGS[key], str(value)]
    if args.absolute_paths:
        cmd.append("--absolute-paths")
    info("creating the psani database")
    run(cmd)
    fofn_path.unlink()

    # Check that every curated plasmid is in the database (psani makedb skips the files
    # that are not FASTA)
    with open(db_json, "rt", encoding="utf-8") as fh:
        db_ids = set(json.load(fh)["seqs"])
    missing = kept - db_ids
    if missing:
        warn(f"{len(missing)} curated plasmids are not in the database "
             f"(e.g. {', '.join(sorted(missing)[:5])})")

    # The per-plasmid FASTA files are only read by psani makedb. The database JSON keeps
    # their paths, so a future psani that reads them again would need --keep-fna
    if not args.keep_fna:
        shutil.rmtree(fna_dir)

    info(f"database created: {db_json} ({len(db_ids)} sequences)")


if __name__ == "__main__":
    main()
