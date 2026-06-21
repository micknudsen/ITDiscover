import pytest

from itdiscover.reads import (
    SequencingRead,
    collapse_identical_reads,
    orient_read,
    passes_read_filters,
    preprocess_reads,
    trim_terminal_ns,
)
from itdiscover.sequences import reverse_complement


def test_sequencing_read_requires_quality_for_each_base() -> None:
    with pytest.raises(ValueError, match="equal length"):
        SequencingRead(
            read_id="read-1",
            sequence="ACGT",
            qualities=(30, 30, 30),
            direction="forward",
        )


def test_sequencing_read_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        SequencingRead(
            read_id="read-1",
            sequence="ACgT",
            qualities=(30, 30, 30, 30),
            direction="forward",
        )


def test_reverse_complement_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        reverse_complement("ACgT")


def test_orient_read_keeps_forward_reads_as_is() -> None:
    read = orient_read(
        read_id="r1",
        sequence="ACGTN",
        qualities=(10, 20, 30, 40, 35),
        direction="forward",
    )

    assert read == SequencingRead(
        read_id="r1",
        sequence="ACGTN",
        qualities=(10, 20, 30, 40, 35),
        direction="forward",
    )


def test_orient_read_rejects_lowercase_forward_reads() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        orient_read(
            read_id="r1",
            sequence="ACgTN",
            qualities=(10, 20, 30, 40, 35),
            direction="forward",
        )


def test_orient_read_reverse_complements_reverse_reads() -> None:
    read = orient_read(
        read_id="r2",
        sequence="ACGTT",
        qualities=(10, 20, 30, 40, 35),
        direction="reverse",
    )

    assert read == SequencingRead(
        read_id="r2",
        sequence="AACGT",
        qualities=(35, 40, 30, 20, 10),
        direction="reverse",
    )


def test_orient_read_rejects_lowercase_reverse_reads() -> None:
    with pytest.raises(ValueError, match="invalid bases"):
        orient_read(
            read_id="r2",
            sequence="ACgTT",
            qualities=(10, 20, 30, 40, 35),
            direction="reverse",
        )


def test_trim_terminal_ns_preserves_internal_ns_and_quality_alignment() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="NNACNGTN",
        qualities=(5, 6, 30, 31, 32, 33, 34, 7),
        direction="forward",
    )

    assert trim_terminal_ns(read) == SequencingRead(
        read_id="read-1",
        sequence="ACNGT",
        qualities=(30, 31, 32, 33, 34),
        direction="forward",
    )


def test_passes_read_filters_uses_trimmed_length_and_mean_quality() -> None:
    passing = SequencingRead(
        read_id="passing",
        sequence="ACGT",
        qualities=(30, 31, 32, 33),
        direction="forward",
    )
    low_quality = SequencingRead(
        read_id="low-quality",
        sequence="ACGT",
        qualities=(10, 10, 10, 10),
        direction="forward",
    )

    assert passes_read_filters(passing, min_length=4, min_mean_quality=30)
    assert not passes_read_filters(passing, min_length=5, min_mean_quality=30)
    assert not passes_read_filters(low_quality, min_length=4, min_mean_quality=30)


def test_collapse_identical_reads_keeps_directions_separate() -> None:
    reads = [
        SequencingRead("read-1", "ACGT", (30, 30, 30, 30), "forward", count=2),
        SequencingRead("read-2", "ACGT", (31, 31, 31, 31), "forward", count=3),
        SequencingRead("read-3", "ACGT", (32, 32, 32, 32), "reverse", count=4),
    ]

    assert collapse_identical_reads(reads) == [
        SequencingRead("read-1", "ACGT", (30, 30, 30, 30), "forward", count=5),
        SequencingRead("read-3", "ACGT", (32, 32, 32, 32), "reverse", count=4),
    ]


def test_preprocess_reads_trims_filters_and_collapses() -> None:
    reads = [
        SequencingRead("read-1", "NACGTN", (5, 30, 30, 30, 30, 5), "forward"),
        SequencingRead("read-2", "ACGT", (31, 31, 31, 31), "forward"),
        SequencingRead("read-3", "NNNN", (40, 40, 40, 40), "forward"),
        SequencingRead("read-4", "ACGT", (10, 10, 10, 10), "forward"),
    ]

    assert preprocess_reads(reads, min_length=4, min_mean_quality=30) == [
        SequencingRead("read-1", "ACGT", (30, 30, 30, 30), "forward", count=2),
    ]
