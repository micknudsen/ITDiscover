"""Command-line interface for ITDiscover."""

import argparse
from pathlib import Path
from typing import TextIO

from . import __version__
from .fastq import read_paired_fastq
from .pipeline import call_exact_itds_from_fragments


class OptionalDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults only for options with concrete defaults."""

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help
        if help_text is None:
            help_text = ""
        if (
            "%(default)" not in help_text
            and action.default is not argparse.SUPPRESS
            and action.default is not None
            and not action.required
        ):
            help_text += " (default: %(default)s)"
        return help_text


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="itdiscover",
        description="Discover FLT3 ITDs from amplicon sequencing of AML samples.",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command")

    call_parser = subparsers.add_parser(
        "call",
        help="Call exact-match ITDs from paired FASTQ files.",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    call_parser.add_argument(
        "--reference",
        required=True,
        help="Reference amplicon FASTA file containing exactly one sequence.",
    )
    call_parser.add_argument(
        "--r1",
        required=True,
        help="Forward-read FASTQ file.",
    )
    call_parser.add_argument(
        "--r2",
        required=True,
        help="Reverse-read FASTQ file.",
    )
    call_parser.add_argument(
        "--min-read-length",
        type=int,
        default=100,
        help="Minimum read length after trimming terminal Ns.",
    )
    call_parser.add_argument(
        "--min-mean-quality",
        type=float,
        default=30,
        help="Minimum mean Phred quality score per read.",
    )
    call_parser.add_argument(
        "--min-insert-length",
        type=int,
        default=6,
        help="Minimum insertion length to consider.",
    )
    call_parser.set_defaults(handler=_run_call_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the ITDiscover CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = getattr(args, "handler", None)
    if handler is None:
        return 0
    return handler(args)


def _run_call_command(args: argparse.Namespace) -> int:
    reference = _read_single_sequence_fasta(Path(args.reference))
    fragments = read_paired_fastq(args.r1, args.r2)
    calls = call_exact_itds_from_fragments(
        fragments,
        reference,
        min_read_length=args.min_read_length,
        min_mean_quality=args.min_mean_quality,
        min_insert_length=args.min_insert_length,
    )
    print("tandem_start\tinsertion_start\tsequence\tsupport_count\tcoverage\tvaf")
    for call in calls:
        print(
            "\t".join(
                [
                    str(call.itd.tandem_start),
                    str(call.itd.insertion.start),
                    call.itd.tandem_sequence,
                    str(call.support_count),
                    str(call.coverage),
                    f"{call.vaf:.6f}",
                ]
            )
        )
    return 0


def _read_single_sequence_fasta(path: Path) -> str:
    with path.open(mode="rt", encoding="utf-8") as handle:
        sequences = list(_iter_fasta_sequences(handle))
    if not sequences:
        raise ValueError("reference FASTA does not contain a sequence")
    if len(sequences) > 1:
        raise ValueError("reference FASTA must contain exactly one sequence")
    return sequences[0]


def _iter_fasta_sequences(handle: TextIO) -> list[str]:
    sequences: list[str] = []
    current_parts: list[str] = []

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_parts:
                sequences.append("".join(current_parts))
                current_parts = []
            continue
        current_parts.append(line)

    if current_parts:
        sequences.append("".join(current_parts))
    return sequences
