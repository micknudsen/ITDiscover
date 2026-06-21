"""Inter-base coverage and variant allele frequency calculations."""

from collections import Counter
from collections.abc import Iterable

from .insertions import Alignment


def covered_reference_positions(alignment: Alignment) -> set[int]:
    """Return zero-based reference positions covered by read bases."""
    covered: set[int] = set()
    ref_pos = -1

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        if ref_base == "-":
            continue
        ref_pos += 1
        if read_base != "-":
            covered.add(ref_pos)

    return covered


def reference_length(alignment: Alignment) -> int:
    """Return the ungapped reference length represented by an alignment."""
    return sum(1 for base in alignment.aligned_reference if base != "-")


def spans_insertion_site(alignment: Alignment, site: int) -> bool:
    """Return whether an alignment spans an inter-base insertion site.

    `site` uses the same convention as `Insertion.start`: an insertion at
    `site` occurs after reference base `site`. Leading insertions use `-1`.
    Internal sites require read bases on both sides of the inter-base position.
    """
    ref_length = reference_length(alignment)
    if site < -1 or site >= ref_length:
        raise ValueError("site is outside the aligned reference")

    covered = covered_reference_positions(alignment)
    if site == -1:
        return 0 in covered
    if site == ref_length - 1:
        return site in covered
    return site in covered and site + 1 in covered


def interbase_coverage(alignments: Iterable[Alignment]) -> dict[int, int]:
    """Return read-count weighted coverage for every spanned insertion site."""
    coverage: Counter[int] = Counter()

    for alignment in alignments:
        ref_length = reference_length(alignment)
        for site in range(-1, ref_length):
            if spans_insertion_site(alignment, site):
                coverage[site] += alignment.count

    return dict(coverage)


def variant_allele_frequency(supporting_count: int, spanning_count: int) -> float:
    """Return VAF as a fraction of supporting reads among spanning reads."""
    if supporting_count < 0:
        raise ValueError("supporting_count must not be negative")
    if spanning_count < 0:
        raise ValueError("spanning_count must not be negative")
    if spanning_count == 0:
        return 0.0
    if supporting_count > spanning_count:
        raise ValueError("supporting_count must not exceed spanning_count")
    return supporting_count / spanning_count
