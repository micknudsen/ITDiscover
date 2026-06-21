import pytest

from itdiscover.insertions import Alignment, Insertion, extract_insertions


def test_alignment_requires_equal_aligned_lengths() -> None:
    with pytest.raises(ValueError, match="equal length"):
        Alignment(
            read_id="read-1",
            read_sequence="ACGT",
            aligned_read="ACGT",
            aligned_reference="ACG",
            direction="forward",
        )


def test_extracts_internal_in_frame_insertion() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCCCGGGTTT",
        aligned_read="AAACCCCCCGGGTTT",
        aligned_reference="AAA------GGGTTT",
        direction="forward",
        count=4,
    )

    assert extract_insertions(alignment) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCCCC",
            direction="forward",
            count=4,
            trailing=False,
        )
    ]


def test_extracts_multiple_insertions_from_one_alignment() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCCCGGGTTTTTTAAA",
        aligned_read="AAACCCCCCGGGTTTTTTAAA",
        aligned_reference="AAA------GGG------AAA",
        direction="forward",
    )

    assert extract_insertions(alignment) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCCCC",
            direction="forward",
        ),
        Insertion(
            read_id="read-1",
            start=5,
            sequence="TTTTTT",
            direction="forward",
        ),
    ]


def test_filters_short_insertions_by_default() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCGGGTTT",
        aligned_read="AAACCGGGTTT",
        aligned_reference="AAA--GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment) == []
    assert extract_insertions(alignment, min_length=2, require_in_frame=False) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CC",
            direction="forward",
        )
    ]


def test_filters_internal_out_of_frame_insertions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACCCCGGGTTT",
        aligned_read="AAACCCCGGGTTT",
        aligned_reference="AAA----GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment, min_length=4) == []
    assert extract_insertions(alignment, min_length=4, require_in_frame=False) == [
        Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCC",
            direction="forward",
        )
    ]


def test_allows_out_of_frame_trailing_insertions() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="CCCCAAA",
        aligned_read="CCCCAAA",
        aligned_reference="----AAA",
        direction="reverse",
    )

    assert extract_insertions(alignment, min_length=4) == [
        Insertion(
            read_id="read-1",
            start=-1,
            sequence="CCCC",
            direction="reverse",
            trailing=True,
        )
    ]


def test_filters_insertions_with_ambiguous_bases() -> None:
    alignment = Alignment(
        read_id="read-1",
        read_sequence="AAACNCGGGTTT",
        aligned_read="AAACNCGGGTTT",
        aligned_reference="AAA---GGGTTT",
        direction="forward",
    )

    assert extract_insertions(alignment, min_length=3) == []
