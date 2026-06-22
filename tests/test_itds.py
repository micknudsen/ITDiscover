import pytest

from itdiscover.insertions import Insertion
from itdiscover.itds import ITD, classify_exact_itd

_Insertion = Insertion


def Insertion(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Insertion(**kwargs)


def test_classifies_exact_upstream_tandem_duplication() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=8,
        sequence="CCCGGG",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="upstream",
    )


def test_classifies_exact_downstream_tandem_duplication() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="reverse",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
    )


def test_prefers_upstream_match_when_both_adjacent_segments_match() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAAAAA") == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="AAA",
        orientation="upstream",
    )


def test_does_not_classify_non_adjacent_reference_match() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="TTT",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") is None


def test_does_not_classify_sequence_absent_from_reference() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="GGGGGG",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") is None


def test_classify_exact_itd_rejects_lowercase_reference() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="forward",
    )

    with pytest.raises(ValueError, match="reference contains invalid bases"):
        classify_exact_itd(insertion, "AAAcccGGGTTT")


def test_itd_reports_inclusive_tandem_end_and_length() -> None:
    itd = ITD(
        insertion=Insertion(
            read_id="read-1",
            start=8,
            sequence="CCCGGG",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="upstream",
    )

    assert itd.tandem_end == 8
    assert itd.length == 6
