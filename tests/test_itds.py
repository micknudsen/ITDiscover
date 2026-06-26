import pytest

from itdiscover.insertions import Insertion
from itdiscover.itds import (
    ITD,
    TandemSimilarity,
    classify_exact_itd,
    classify_fuzzy_itd,
    score_tandem_similarity,
)

_Insertion = Insertion


def Insertion(**kwargs):
    kwargs.setdefault("fragment_id", kwargs["read_id"])
    return _Insertion(**kwargs)


def test_classifies_exact_upstream_tandem_duplication_as_canonical_downstream() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=8,
        sequence="CCCGGG",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAACCCGGGTTT") == ITD(
        insertion=Insertion(
            read_id="read-1",
            start=2,
            sequence="CCCGGG",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
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


def test_uses_right_most_canonical_downstream_break_when_both_sides_match() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAAAAA") == ITD(
        insertion=Insertion(
            read_id="read-1",
            start=2,
            sequence="AAA",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="AAA",
        orientation="downstream",
    )


def test_slides_repetitive_exact_itds_to_right_most_equivalent_breakpoint() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=0,
        sequence="AAA",
        direction="forward",
    )

    assert classify_exact_itd(insertion, "AAAAAA") == ITD(
        insertion=Insertion(
            read_id="read-1",
            start=2,
            sequence="AAA",
            direction="forward",
        ),
        tandem_start=3,
        tandem_sequence="AAA",
        orientation="downstream",
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


def test_scores_exact_tandem_similarity() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGG",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAACCCGGGTTT")
    
    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        mismatches=0,
    )


def test_scores_best_tandem_similarity_with_one_mismatch() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAACCCGGGTTT")

    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        mismatches=1,
    )
    assert similarity.matches == 5
    assert similarity.identity == 5 / 6


def test_scores_right_most_tandem_on_tie() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=5,
        sequence="AAA",
        direction="forward",
    )

    similarity = score_tandem_similarity(insertion, "AAAAAA")

    assert similarity == TandemSimilarity(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="AAA",
        mismatches=0,
    )


def test_scores_tandem_similarity_returns_none_for_long_insertions() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGGTTTAAAA",
        direction="forward",
    )

    assert score_tandem_similarity(insertion, "AAACCCGGGTTT") is None


def test_classifies_fuzzy_downstream_tandem_duplication_with_one_mismatch() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=1,
    ) == ITD(
        insertion=insertion,
        tandem_start=3,
        tandem_sequence="CCCGGG",
        orientation="downstream",
    )


def test_does_not_classify_fuzzy_itd_when_mismatches_exceed_threshold() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=0,
    ) is None


def test_does_not_classify_fuzzy_itd_when_best_window_is_non_adjacent() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="TTT",
        direction="forward",
    )

    assert classify_fuzzy_itd(
        insertion,
        "AAACCCGGGTTT",
        max_mismatches=3,
    ) is None


def test_does_not_classify_fuzzy_itd_when_adjacent_windows_tie() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="AAT",
        direction="forward",
    )

    assert classify_fuzzy_itd(insertion, "AAAAAA", max_mismatches=1) is None


def test_classify_fuzzy_itd_rejects_negative_mismatch_threshold() -> None:
    insertion = Insertion(
        read_id="read-1",
        start=2,
        sequence="CCCGGA",
        direction="forward",
    )

    with pytest.raises(ValueError, match="max_mismatches must not be negative"):
        classify_fuzzy_itd(insertion, "AAACCCGGGTTT", max_mismatches=-1)
