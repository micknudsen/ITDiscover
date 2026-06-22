import pytest

from itdiscover.alignment import AlignmentScoring, align_read_to_reference
from itdiscover.insertions import Alignment
from itdiscover.reads import SequencingRead

_Alignment = Alignment
_SequencingRead = SequencingRead


def Alignment(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Alignment(**kwargs)


def SequencingRead(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _SequencingRead(**kwargs)


def test_align_read_to_reference_aligns_perfect_full_length_read() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="AAACCCGGG",
        qualities=(40,) * 9,
        direction="forward",
    )

    assert align_read_to_reference(read, "AAACCCGGG") == Alignment(
        read_id="read-1",
        read_sequence="AAACCCGGG",
        aligned_read="AAACCCGGG",
        aligned_reference="AAACCCGGG",
        direction="forward",
    )


def test_align_read_to_reference_uses_unpenalized_end_gaps() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="AAACCC",
        qualities=(40,) * 6,
        direction="forward",
    )

    assert align_read_to_reference(read, "AAACCCGGG") == Alignment(
        read_id="read-1",
        read_sequence="AAACCC",
        aligned_read="AAACCC---",
        aligned_reference="AAACCCGGG",
        direction="forward",
    )


def test_align_read_to_reference_represents_insertions_as_reference_gaps() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="AAACCCGGGTTT",
        qualities=(40,) * 12,
        direction="reverse",
    )

    assert align_read_to_reference(read, "AAAGGGTTT") == Alignment(
        read_id="read-1",
        read_sequence="AAACCCGGGTTT",
        aligned_read="AAACCCGGGTTT",
        aligned_reference="AAA---GGGTTT",
        direction="reverse",
    )


def test_align_read_to_reference_accepts_custom_scoring() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="AAATTT",
        qualities=(40,) * 6,
        direction="forward",
    )

    alignment = align_read_to_reference(
        read,
        "AAACCC",
        scoring=AlignmentScoring(match=2, mismatch=-1, gap_open=-5, gap_extend=-1),
    )

    assert alignment.aligned_read == "AAATTT"
    assert alignment.aligned_reference == "AAACCC"


def test_align_read_to_reference_rejects_lowercase_reference() -> None:
    read = SequencingRead(
        read_id="read-1",
        sequence="AAACCC",
        qualities=(40,) * 6,
        direction="forward",
    )

    with pytest.raises(ValueError, match="reference contains invalid bases"):
        align_read_to_reference(read, "AAAccc")
