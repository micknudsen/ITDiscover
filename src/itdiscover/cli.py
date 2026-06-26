"""Command-line interface for ITDiscover."""

import argparse
import html
from pathlib import Path
from typing import TextIO

from . import __version__
from .alignment import align_read_to_reference
from .calls import (
    ITDCall,
    ITDFilter,
    UniqueSupportRepresentative,
    call_exact_itds_with_representatives,
    call_fuzzy_itds_with_representatives,
)
from .fastq import read_paired_fastq
from .insertions import Alignment
from .reads import preprocess_fragments


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
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference amplicon FASTA file containing exactly one sequence.",
    )
    parser.add_argument(
        "--r1",
        required=True,
        help="Forward-read FASTQ file.",
    )
    parser.add_argument(
        "--r2",
        required=True,
        help="Reverse-read FASTQ file.",
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=100,
        help="Minimum read length after trimming terminal Ns.",
    )
    parser.add_argument(
        "--min-mean-quality",
        type=float,
        default=30,
        help="Minimum mean Phred quality score per read.",
    )
    parser.add_argument(
        "--min-insert-length",
        type=int,
        default=6,
        help="Minimum insertion length to consider.",
    )
    parser.add_argument(
        "--max-mismatches",
        type=_non_negative_int,
        help="Allow fuzzy ITD calls with at most this many mismatches.",
    )
    parser.add_argument(
        "--min-support-count",
        type=int,
        default=1,
        help="Minimum fragment support count required for a call to pass filtering.",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=0,
        help="Minimum coverage required for a call to pass filtering.",
    )
    parser.add_argument(
        "--min-vaf",
        type=float,
        default=0.0,
        help="Minimum VAF required for a call to pass filtering.",
    )
    parser.add_argument(
        "--output",
        type=_html_output_path,
        help="Optional path for an HTML report with one representative alignment per unique support pattern.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the ITDiscover CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return _run_call_command(args)


def _run_call_command(args: argparse.Namespace) -> int:
    reference = _read_single_sequence_fasta(Path(args.reference))
    fragments = read_paired_fastq(args.r1, args.r2)
    processed_reads = preprocess_fragments(
        fragments,
        min_length=args.min_read_length,
        min_mean_quality=args.min_mean_quality,
    )
    alignments = [
        align_read_to_reference(read, reference)
        for read in processed_reads
    ]
    filters = ITDFilter(
        min_support_count=args.min_support_count,
        min_coverage=args.min_coverage,
        min_vaf=args.min_vaf,
    )
    if args.max_mismatches is None:
        calls, representatives = call_exact_itds_with_representatives(
            alignments,
            reference,
            min_insert_length=args.min_insert_length,
            filters=filters,
        )
    else:
        calls, representatives = call_fuzzy_itds_with_representatives(
            alignments,
            reference,
            max_mismatches=args.max_mismatches,
            min_insert_length=args.min_insert_length,
            filters=filters,
        )
    print(
        "tandem_start\tinsertion_start\tsequence\t"
        "support_count\tcoverage\tvaf\tstatus\tfilter_reasons"
    )
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
                    call.status,
                    _format_filter_reasons(call),
                ]
            )
        )
    if args.output:
        _write_unique_support_alignment_html_report(
            args.output,
            calls,
            representatives,
        )
    return 0


def _html_output_path(value: str) -> Path:
    path = Path(value)
    if path.suffix.lower() != ".html":
        raise argparse.ArgumentTypeError("output path must end with .html")
    return path


def _non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must not be negative")
    return parsed


def _format_filter_reasons(call: ITDCall) -> str:
    return "." if not call.filter_reasons else ";".join(call.filter_reasons)


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


def _write_unique_support_alignment_html_report(
    path: Path,
    calls: list[ITDCall],
    representatives: list[UniqueSupportRepresentative],
) -> None:
    representatives_by_key: dict[tuple[int, str, bool], list[UniqueSupportRepresentative]] = {}
    for representative in representatives:
        key = (
            representative.itd.tandem_start,
            representative.itd.tandem_sequence,
            representative.itd.insertion.trailing,
        )
        representatives_by_key.setdefault(key, []).append(representative)

    sections: list[str] = []
    for index, call in enumerate(calls, start=1):
        key = (
            call.itd.tandem_start,
            call.itd.tandem_sequence,
            call.itd.insertion.trailing,
        )
        call_representatives = representatives_by_key.get(key, [])
        sections.append(
            _render_html_call_section(index, call, call_representatives)
        )

    document = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>ITDiscover Unique Support Alignments</title>
  <style>
    :root {
      color-scheme: light;
      --ink: #18232f;
      --muted: #5a6875;
      --line: #cfd8df;
      --panel: #f7fafc;
      --diff-bg: #ffe7a3;
      --diff-fg: #6e4a00;
      --insert-bg: #d8f3dc;
      --insert-fg: #164b23;
    }
    body {
      margin: 24px;
      color: var(--ink);
      font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
      line-height: 1.4;
    }
    h1, h2, h3 {
      margin: 0 0 12px;
    }
    .itd {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 16px;
      margin-bottom: 20px;
      background: white;
    }
    .summary {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
      gap: 10px 16px;
      margin-bottom: 16px;
    }
    .summary div {
      background: var(--panel);
      border-radius: 4px;
      padding: 8px 10px;
    }
    .summary dt {
      margin: 0 0 4px;
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }
    .summary dd {
      margin: 0;
      font-size: 16px;
      font-weight: 600;
    }
    .support {
      border-top: 1px solid var(--line);
      padding-top: 14px;
      margin-top: 14px;
    }
    .support-header {
      display: flex;
      gap: 16px;
      flex-wrap: wrap;
      margin-bottom: 10px;
      align-items: baseline;
    }
    .support-header strong {
      font-size: 16px;
    }
    .support-meta {
      color: var(--muted);
      font-size: 14px;
    }
    .signature {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: var(--panel);
      border-radius: 4px;
      padding: 6px 8px;
      margin-bottom: 10px;
      display: inline-block;
    }
    .alignment-block {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
      overflow-x: auto;
      background: #fbfdff;
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 10px;
    }
    .alignment-row {
      white-space: pre;
    }
    .label {
      color: var(--muted);
    }
    .diff {
      background: var(--diff-bg);
      color: var(--diff-fg);
      font-weight: 700;
    }
    .insert {
      background: var(--insert-bg);
      color: var(--insert-fg);
      font-weight: 700;
    }
  </style>
</head>
<body>
  <h1>Unique Support Alignments</h1>
  __SECTIONS__
</body>
</html>
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        document.replace("__SECTIONS__", "\n".join(sections)),
        encoding="utf-8",
    )


def _render_html_call_section(
    index: int,
    call: ITDCall,
    representatives: list[UniqueSupportRepresentative],
) -> str:
    summary = (
        ('Tandem Start', str(call.itd.tandem_start)),
        ('Insertion Start', str(call.itd.insertion.start)),
        ('Sequence', call.itd.tandem_sequence),
        ('Support Count', str(call.support_count)),
        ('Coverage', str(call.coverage)),
        ('VAF', f"{call.vaf:.6f}"),
        ('Status', call.status),
        ('Filter Reasons', _format_filter_reasons(call)),
    )
    summary_html = "".join(
        f"<div><dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd></div>"
        for label, value in summary
    )
    baseline_alignment = representatives[0].alignment if representatives else None
    support_html = "".join(
        _render_html_support_block(
            support_index,
            representative,
            baseline_alignment,
        )
        for support_index, representative in enumerate(representatives, start=1)
    )
    return (
        f'<section class="itd">'
        f"<h2>ITD {index}</h2>"
        f'<dl class="summary">{summary_html}</dl>'
        f"{support_html}"
        f"</section>"
    )


def _render_html_support_block(
    support_index: int,
    representative: UniqueSupportRepresentative,
    baseline_alignment: Alignment | None,
) -> str:
    alignment = representative.alignment
    reference_html = _highlight_alignment_differences(
        alignment.aligned_reference,
        comparison_classes=None,
    )
    comparison_classes = (
        _alignment_comparison_classes(alignment, baseline_alignment)
        if baseline_alignment is not None
        else None
    )
    read_html = _highlight_alignment_differences(
        alignment.aligned_read,
        comparison_classes=comparison_classes,
    )
    return (
        '<section class="support">'
        f'<div class="support-header"><strong>Support {support_index}</strong>'
        f'<span class="support-meta">count {representative.support_count}</span>'
        f'<span class="support-meta">read {html.escape(alignment.read_id)}</span>'
        f'<span class="support-meta">fragment {html.escape(alignment.fragment_id)}</span>'
        '</div>'
        f'<div class="signature">{html.escape(representative.signature)}</div>'
        '<div class="alignment-block">'
        f'<div class="alignment-row"><span class="label">reference  </span>{reference_html}</div>'
        f'<div class="alignment-row"><span class="label">read       </span>{read_html}</div>'
        '</div>'
        '</section>'
    )


def _highlight_alignment_differences(
    sequence: str,
    *,
    comparison_classes: list[str | None] | None,
) -> str:
    fragments: list[str] = []
    for index, base in enumerate(sequence):
        escaped_base = html.escape(base)
        css_class = None
        if comparison_classes is not None and index < len(comparison_classes):
            css_class = comparison_classes[index]
        if css_class is not None:
            fragments.append(f'<span class="{css_class}">{escaped_base}</span>')
            continue
        fragments.append(escaped_base)
    return "".join(fragments)


def _alignment_comparison_classes(
    alignment,
    baseline_alignment,
) -> list[str | None]:
    baseline_reference_bases, baseline_insertions = _alignment_features(
        baseline_alignment.aligned_reference,
        baseline_alignment.aligned_read,
    )
    classes: list[str | None] = []
    ref_pos = -1
    insertion_offsets: dict[int, int] = {}

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        css_class = None
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                baseline_base = baseline_reference_bases.get(ref_pos)
                if baseline_base is not None and baseline_base != read_base:
                    css_class = "diff"
        elif read_base != "-":
            insertion_site = ref_pos
            offset = insertion_offsets.get(insertion_site, 0)
            baseline_insertion = baseline_insertions.get(insertion_site, "")
            if offset < len(baseline_insertion):
                if baseline_insertion[offset] != read_base:
                    css_class = "diff"
            else:
                css_class = "insert"
            insertion_offsets[insertion_site] = offset + 1

        classes.append(css_class)

    return classes


def _alignment_features(
    aligned_reference: str,
    aligned_read: str,
) -> tuple[dict[int, str], dict[int, str]]:
    reference_bases: dict[int, str] = {}
    insertions: dict[int, list[str]] = {}
    ref_pos = -1

    for read_base, ref_base in zip(aligned_read, aligned_reference, strict=True):
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                reference_bases[ref_pos] = read_base
            continue
        if read_base != "-":
            insertions.setdefault(ref_pos, []).append(read_base)

    return (
        reference_bases,
        {site: "".join(bases) for site, bases in insertions.items()},
    )
