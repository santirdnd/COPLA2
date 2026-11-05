#!/usr/bin/env python3

"""
Prepare the query plasmids of a nucleotide multi-FASTA file for the COPLA 2 pipeline.

This script performs the following steps:

  - Splits a nucleotide multi-FASTA file into separate FASTA files grouped by
    the identifier prefix before a delimiter (default: '|'), following the
    pattern >genome_id|sequence_id. This allows representing multi-genome FASTA
    files in which sequences from different genomes are stored together.

  - Validates FASTA header consistency, ensuring that either all headers use
    the delimiter pattern or none of them do.

  - Normalizes sequences so that each FASTA record contains a single-line
    nucleotide sequence.

  - Computes basic sequence statistics per group, including sequence counts,
    length statistics, and the fractions of A/C/G/T/N characters and of N.

  - Reports potential issues such as sequences shorter than a configurable
    threshold, unusually long sequences, low ACGTN fraction (which may
    indicate amino acid sequences) or high N fraction (unknown bases, which
    reduce the fragments available for ANI).

  - Generates two auxiliary files used by the COPLA pipeline:
        * <input_basename>.fofn          (file of FASTA file paths)
        * <input_basename>.manifest.tsv  (sequence statistics per group)
"""

import io
import sys
import gzip
import argparse
from pathlib import Path
from typing import Optional, TextIO


# Exit codes
EXIT_OK = 0
EXIT_SUCCESS = 0
EXIT_ERROR = 1
EXIT_INVALID_ARGS = 2
EXIT_INVALID_FASTA = 3
EXIT_OUTPUT_COLLISION = 4
EXIT_INPUT_FILE = 5


def die(msg: str, code: int = EXIT_ERROR) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def strip_known_extensions(filename: str) -> str:
    """
    Remove common FASTA and gzip extensions from a filename.
    Examples:
      input.fasta.gz -> input
      input.fa -> input
      input.fna.gz -> input
    """
    base = filename

    if base.endswith(".gz"):
        base = base[:-3]

    for ext in (".fasta", ".fas", ".fa", ".fna"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break

    return base


def open_maybe_gzip(path: str, mode: str) -> TextIO:
    """
    Open a text file that may or may not be gzip-compressed.

    Detection is based on the gzip magic numbers instead of the filename.

    Supported modes
      rt, wt, at
    """
    allowed_modes = {"rt", "wt", "at"}
    if mode not in allowed_modes:
        die(
            f"Unsupported mode '{mode}'. Allowed modes: {sorted(allowed_modes)}",
            EXIT_INVALID_ARGS,
        )

    if mode == "rt":
        try:
            fh = open(path, "rb")
        except OSError as e:
            die(f"Cannot open input file '{path}': {e}", EXIT_INPUT_FILE)
        try:
            magic = fh.read(3)
            fh.seek(0)

            if magic == b"\x1f\x8b\x08":
                # Wrap existing handle without reopening file
                gz = gzip.GzipFile(fileobj=fh)
                return io.TextIOWrapper(gz, encoding="utf-8")

            if path.endswith(".gz"):
                die(
                    f"File '{path}' ends with '.gz' but is not a valid gzip file.",
                    EXIT_INPUT_FILE,
                )

            return io.TextIOWrapper(fh, encoding="utf-8")

        except Exception:
            fh.close()
            raise

    if path.endswith(".gz"):
        try:
            return gzip.open(path, mode, encoding="utf-8")
        except OSError as e:
            die(f"Cannot open output file '{path}': {e}", EXIT_OUTPUT_COLLISION)

    try:
        return open(path, mode, encoding="utf-8")
    except OSError as e:
        die(f"Cannot open output file '{path}': {e}", EXIT_OUTPUT_COLLISION)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare and validate the query plasmids of a nucleotide multi-FASTA file for the COPLA 2 "
            "pipeline by splitting sequences into grouped FASTA files based on the identifier prefix "
            "before a delimiter."
        )
    )
    parser.add_argument(
        "file",
        type=str,
        help="input nucleotide multi-FASTA file (optionally .gz)",
    )
    parser.add_argument(
        "outdir",
        type=str,
        help="destination directory for output FASTA files",
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default="|",
        help="delimiter character used to define the group (default: '|')",
    )
    parser.add_argument(
        "--name-mode",
        choices=["prefix", "enumerate"],
        default="prefix",
        help=(
            "output filename mode: 'prefix' uses the group label, "
            "'enumerate' uses <input_basename>_N.fasta (default: prefix)"
        ),
    )
    parser.add_argument(
        "--ext",
        type=str,
        default="fasta",
        help="output FASTA extension without dot (default: fasta)",
    )
    parser.add_argument(
        "--min-seq-len",
        type=int,
        default=1000,
        help="warning threshold for short sequences in bp (default: 1000)",
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=0,
        help=(
            "warning threshold for unusually long sequences in bp; "
            "0 disables the warning (default: 0)"
        ),
    )
    parser.add_argument(
        "--min-acgtn-fraction",
        type=float,
        default=0.85,
        help="warning threshold for ACGTN fraction per group (default: 0.85)",
    )
    parser.add_argument(
        "--max-n-fraction",
        type=float,
        default=0.2,
        help=(
            "warning threshold for N fraction per group, checked only when the ACGTN "
            "fraction is not low (default: 0.2)"
        ),
    )
    parser.add_argument(
        "--gzip-output",
        action="store_true",
        help="compress output FASTA files with gzip",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite existing output files (the input file is never overwritten)",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 2.0.0")
    args = parser.parse_args()

    if len(args.delimiter) != 1:
        die("--delimiter must be a single character.", EXIT_INVALID_ARGS)

    if not args.ext.lstrip("."):
        die("--ext must not be empty.", EXIT_INVALID_ARGS)

    if args.min_seq_len < 0:
        die("--min-seq-len must be >= 0.", EXIT_INVALID_ARGS)

    if args.max_seq_len < 0:
        die("--max-seq-len must be >= 0.", EXIT_INVALID_ARGS)

    if not 0.0 <= args.min_acgtn_fraction <= 1.0:
        die("--min-acgtn-fraction must be between 0 and 1.", EXIT_INVALID_ARGS)

    if not 0.0 <= args.max_n_fraction <= 1.0:
        die("--max-n-fraction must be between 0 and 1.", EXIT_INVALID_ARGS)

    input_path = Path(args.file)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_basename = strip_known_extensions(input_path.name)
    if not input_basename:
        input_basename = "input"

    group_to_filename: dict[str, str] = {}
    group_to_count: dict[str, int] = {}
    group_to_index: dict[str, int] = {}
    output_handles: dict[str, TextIO] = {}
    group_to_min_len: dict[str, int] = {}
    group_to_max_len: dict[str, int] = {}
    group_to_total_len: dict[str, int] = {}
    group_to_below_min_len: dict[str, int] = {}
    group_to_acgtn_count: dict[str, int] = {}
    group_to_n_count: dict[str, int] = {}

    seen_output_paths: set[Path] = set()

    next_index = 1
    total_headers = 0
    saw_with_delimiter = False
    saw_without_delimiter = False

    current_header: Optional[str] = None
    current_group: Optional[str] = None
    current_seq_chunks: list[str] = []

    def validate_group_name(group: str) -> None:
        if not group:
            die(
                f"Empty group detected in header pattern "
                f"'>group{args.delimiter}sequence_id'.",
                EXIT_INVALID_FASTA,
            )

        if group in {".", ".."}:
            die(f"Invalid group name '{group}'.", EXIT_INVALID_FASTA)

        invalid_chars = set('<>:"/\\|?*')
        if any(ch in invalid_chars for ch in group):
            die(
                f"Invalid group name '{group}': unsupported filename characters detected.",
                EXIT_INVALID_FASTA,
            )

    def reserve_output_path(filename: str) -> None:
        out_path = (outdir / filename).resolve()

        # The input is read while the outputs are written, so it can never be overwritten
        if out_path == input_path.resolve():
            die(f"Output file would overwrite the input file: '{out_path}'", EXIT_OUTPUT_COLLISION)

        if out_path.exists() and not args.force:
            die(
                f"Output file already exists: '{out_path}' (use --force to overwrite)",
                EXIT_OUTPUT_COLLISION,
            )

        if out_path in seen_output_paths:
            die(f"Output file collision detected: '{out_path}'", EXIT_OUTPUT_COLLISION)

        seen_output_paths.add(out_path)

    def output_filename_for_group(group: str) -> str:
        nonlocal next_index

        if group in group_to_filename:
            return group_to_filename[group]

        ext = args.ext.lstrip(".")

        if group == input_basename and not saw_with_delimiter:
            stem = input_basename
        elif args.name_mode == "prefix":
            validate_group_name(group)
            stem = group
        else:
            index = next_index
            group_to_index[group] = index
            stem = f"{input_basename}_{index}"
            next_index += 1

        filename = f"{stem}.{ext}"
        if args.gzip_output:
            filename += ".gz"

        reserve_output_path(filename)
        group_to_filename[group] = filename
        return filename

    def get_output_handle(group: str) -> TextIO:
        if group not in output_handles:
            filename = output_filename_for_group(group)
            out_path = outdir / filename
            output_handles[group] = open_maybe_gzip(str(out_path), "wt")
        return output_handles[group]

    def write_current_record() -> None:
        if current_header is None or current_group is None:
            return

        out_fh = get_output_handle(current_group)
        sequence = "".join(current_seq_chunks)

        seq_len = len(sequence)
        acgtn_count = sum(1 for ch in sequence if ch.upper() in {"A", "C", "G", "T", "N"})
        n_count = sequence.count("N") + sequence.count("n")

        out_fh.write(current_header + "\n")
        out_fh.write(sequence + "\n")

        group_to_count[current_group] = group_to_count.get(current_group, 0) + 1

        if current_group in group_to_min_len:
            group_to_min_len[current_group] = min(group_to_min_len[current_group], seq_len)
            group_to_max_len[current_group] = max(group_to_max_len[current_group], seq_len)
        else:
            group_to_min_len[current_group] = seq_len
            group_to_max_len[current_group] = seq_len

        group_to_total_len[current_group] = group_to_total_len.get(current_group, 0) + seq_len
        group_to_below_min_len[current_group] = (
            group_to_below_min_len.get(current_group, 0)
            + (1 if seq_len < args.min_seq_len else 0)
        )
        group_to_acgtn_count[current_group] = (
            group_to_acgtn_count.get(current_group, 0) + acgtn_count
        )
        group_to_n_count[current_group] = group_to_n_count.get(current_group, 0) + n_count

    try:
        with open_maybe_gzip(args.file, "rt") as fh_in:
            for raw in fh_in:
                if not raw:
                    continue

                line = raw.rstrip("\r\n")
                if not line:
                    continue

                if line.startswith(">"):
                    if current_header is not None:
                        write_current_record()
                        current_seq_chunks = []

                    total_headers += 1

                    seq_id = line[1:].partition(" ")[0]

                    if args.delimiter in seq_id:
                        saw_with_delimiter = True
                        group = seq_id.partition(args.delimiter)[0]
                    else:
                        saw_without_delimiter = True
                        group = input_basename

                    if saw_with_delimiter and saw_without_delimiter:
                        die(
                            "Mixed FASTA identifier patterns detected.\n"
                            f"  Some headers contain delimiter '{args.delimiter}' and others do not.\n"
                            f"  The pattern >genome_id{args.delimiter}sequence_id is used to distinguish contigs\n"
                            "  from different genomes within the same multi-FASTA file.\n"
                            "  All FASTA headers must consistently use this convention, or omit it.",
                            EXIT_INVALID_FASTA,
                        )

                    current_header = line
                    current_group = group

                else:
                    if current_header is None:
                        die("FASTA sequence data found before the first header.", EXIT_INVALID_FASTA)

                    current_seq_chunks.append(line.replace(" ", "").replace("\t", ""))

        if current_header is not None:
            write_current_record()

    finally:
        for fh in output_handles.values():
            fh.close()

    if total_headers == 0:
        die("No FASTA headers found.", EXIT_INVALID_FASTA)

    if args.name_mode == "enumerate" and saw_with_delimiter:
        groups = sorted(group_to_filename.keys(), key=lambda g: group_to_index[g])
    else:
        groups = sorted(group_to_filename.keys())

    fofn_path = outdir / f"{input_basename}.fofn"
    manifest_path = outdir / f"{input_basename}.manifest.tsv"

    def get_fofn_group(group: str, filename: str) -> str:
        if group == input_basename and not saw_with_delimiter:
            return input_basename
        if args.name_mode == "enumerate":
            fofn_group = Path(filename).name
            if fofn_group.endswith(".gz"):
                fofn_group = Path(fofn_group).stem
            return Path(fofn_group).stem
        return group

    with open(fofn_path, "wt", encoding="utf-8") as fofn_fh:
        for group in groups:
            filename = group_to_filename[group]
            full_path = (outdir / filename).resolve()
            n_records = group_to_count.get(group, 0)
            fofn_group = get_fofn_group(group, filename)
            fofn_fh.write(f"{full_path}\t{fofn_group}\t{n_records}\n")

    with open(manifest_path, "wt", encoding="utf-8") as manifest_fh:
        manifest_fh.write(
            "#QueryID\tGenome\tContigs\tTotal_Length\tMin_Contig\t"
            "Max_Contig\tContigs_Below_Min\tACGTN_Fraction\tN_Fraction\tNotes\n"
        )
        for group in groups:
            filename = group_to_filename[group]
            fofn_group = get_fofn_group(group, filename)
            n_records = group_to_count.get(group, 0)
            total_len = group_to_total_len.get(group, 0)
            min_len = group_to_min_len.get(group, 0)
            max_len = group_to_max_len.get(group, 0)
            n_below_min_len = group_to_below_min_len.get(group, 0)
            acgtn_count = group_to_acgtn_count.get(group, 0)
            acgtn_fraction = (acgtn_count / total_len) if total_len > 0 else 0.0
            n_count = group_to_n_count.get(group, 0)
            n_fraction = (n_count / total_len) if total_len > 0 else 0.0

            # Notes are separated by "; ", and the parts of a note by ", "
            notes: list[str] = []

            if n_below_min_len > 0:
                label = "contig" if n_below_min_len == 1 else "contigs"
                notes.append(
                    f"{n_below_min_len} {label} shorter than {args.min_seq_len} bp, "
                    "too short for ANI fragments"
                )

            if args.max_seq_len > 0 and max_len > args.max_seq_len:
                notes.append(f"max length exceeds {args.max_seq_len} bp")

            if acgtn_fraction < args.min_acgtn_fraction:
                notes.append(
                    f"low ACGTN fraction ({acgtn_fraction:.2f}), possible amino acid sequences"
                )
            elif n_fraction > args.max_n_fraction:
                # In amino acid sequences N is asparagine, so the N fraction is only
                # checked for nucleotide sequences
                notes.append(
                    f"too many unknown bases (N fraction {n_fraction:.2f}), "
                    "ANI links may be missed"
                )

            note = "; ".join(notes)
            manifest_fh.write(
                f"{fofn_group}\t{group}\t{n_records}\t{total_len}\t{min_len}\t"
                f"{max_len}\t{n_below_min_len}\t{acgtn_fraction:.6f}\t{n_fraction:.6f}\t{note}\n"
            )


if __name__ == "__main__":
    main()
