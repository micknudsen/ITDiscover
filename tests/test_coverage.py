import pytest

from itdiscover.coverage import (
    covered_reference_positions,
    interbase_coverage,
    spans_insertion_site,
    variant_allele_frequency,
)
from itdiscover.insertions import Alignment


def test_covered_reference_positions_excludes_insertions_and_deletions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AACCCT",
        aligned_read="AA---C-T",
        aligned_reference="AAGGG-CT",
        direction="forward",
    )

    assert covered_reference_positions(alignment) == {0, 1, 6}


def test_internal_site_requires_coverage_on_both_sides() -> None:
    spanning = Alignment(
        read_id="spanning",
        read_sequence="AAACCC",
        aligned_read="AAACCC",
        aligned_reference="AAACCC",
        direction="forward",
    )
    deleted_right_side = Alignment(
        read_id="deleted",
        read_sequence="AAACC",
        aligned_read="AAA-CC",
        aligned_reference="AAACCC",
        direction="forward",
    )

    assert spans_insertion_site(spanning, 2)
    assert not spans_insertion_site(deleted_right_side, 2)


def test_edge_sites_require_terminal_reference_base_coverage() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AACCC",
        aligned_read="-AACCC",
        aligned_reference="AAACCC",
        direction="forward",
    )

    assert not spans_insertion_site(alignment, -1)
    assert spans_insertion_site(alignment, 5)


def test_interbase_coverage_is_weighted_by_alignment_count() -> None:
    alignments = [
        Alignment(
            read_id="read-1",
            read_sequence="AAACCC",
            aligned_read="AAACCC",
            aligned_reference="AAACCC",
            direction="forward",
            count=3,
        ),
        Alignment(
            read_id="read-2",
            read_sequence="CCC",
            aligned_read="---CCC",
            aligned_reference="AAACCC",
            direction="reverse",
            count=2,
        ),
    ]

    assert interbase_coverage(alignments) == {
        -1: 3,
        0: 3,
        1: 3,
        2: 3,
        3: 5,
        4: 5,
        5: 5,
    }


def test_spans_insertion_site_rejects_out_of_range_site() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAA",
        aligned_read="AAA",
        aligned_reference="AAA",
        direction="forward",
    )

    with pytest.raises(ValueError, match="outside"):
        spans_insertion_site(alignment, 3)


def test_variant_allele_frequency_returns_fraction() -> None:
    assert variant_allele_frequency(2, 8) == 0.25
    assert variant_allele_frequency(0, 8) == 0.0
    assert variant_allele_frequency(0, 0) == 0.0


def test_variant_allele_frequency_rejects_impossible_counts() -> None:
    with pytest.raises(ValueError, match="must not be negative"):
        variant_allele_frequency(-1, 10)
    with pytest.raises(ValueError, match="must not be negative"):
        variant_allele_frequency(1, -10)
    with pytest.raises(ValueError, match="must not exceed"):
        variant_allele_frequency(11, 10)
