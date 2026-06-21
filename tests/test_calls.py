from itdiscover.calls import ITDCall, call_exact_itds
from itdiscover.insertions import Alignment, Insertion
from itdiscover.itds import ITD


def test_call_exact_itds_reports_support_coverage_and_vaf() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="itd-read",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAACCCGGG------TTT",
            direction="forward",
            count=3,
        ),
        Alignment(
            read_id="wt-read",
            read_sequence=reference,
            aligned_read=reference,
            aligned_reference=reference,
            direction="reverse",
            count=7,
        ),
    ]

    assert call_exact_itds(alignments, reference) == [
        ITDCall(
            itd=ITD(
                insertion=Insertion(
                    read_id="itd-read",
                    start=8,
                    sequence="CCCGGG",
                    direction="forward",
                    count=3,
                ),
                tandem_start=3,
                tandem_sequence="CCCGGG",
                orientation="upstream",
            ),
            support_count=3,
            coverage=10,
            vaf=0.3,
        )
    ]


def test_call_exact_itds_merges_upstream_and_downstream_representations() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="upstream-representation",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAACCCGGG------TTT",
            direction="forward",
            count=2,
        ),
        Alignment(
            read_id="downstream-representation",
            read_sequence="AAACCCGGGCCCGGGTTT",
            aligned_read="AAACCCGGGCCCGGGTTT",
            aligned_reference="AAA------CCCGGGTTT",
            direction="reverse",
            count=3,
        ),
        Alignment(
            read_id="wt-read",
            read_sequence=reference,
            aligned_read=reference,
            aligned_reference=reference,
            direction="forward",
            count=7,
        ),
    ]

    calls = call_exact_itds(alignments, reference)

    assert len(calls) == 1
    assert calls[0].support_count == 5
    assert calls[0].coverage == 12
    assert calls[0].vaf == 5 / 12
    assert calls[0].itd.insertion.read_id == "downstream-representation"
    assert calls[0].itd.orientation == "downstream"


def test_call_exact_itds_ignores_non_itd_insertions() -> None:
    reference = "AAACCCGGGTTT"
    alignments = [
        Alignment(
            read_id="insertion-read",
            read_sequence="AAACCCGGGAAAAAATTT",
            aligned_read="AAACCCGGGAAAAAATTT",
            aligned_reference="AAACCCGGG------TTT",
            direction="forward",
            count=4,
        ),
        Alignment(
            read_id="wt-read",
            read_sequence=reference,
            aligned_read=reference,
            aligned_reference=reference,
            direction="forward",
            count=6,
        ),
    ]

    assert call_exact_itds(alignments, reference) == []


def test_call_exact_itds_returns_sorted_calls() -> None:
    reference = "AAACCCGGGTTTAAACCC"
    alignments = [
        Alignment(
            read_id="later-itd",
            read_sequence="AAACCCGGGTTTAAACCCCCC",
            aligned_read="AAACCCGGGTTTAAACCCCCC",
            aligned_reference="AAACCCGGGTTTAAACCC---",
            direction="forward",
        ),
        Alignment(
            read_id="earlier-itd",
            read_sequence="AAACCCGGGCCCGGGTTTAAACCC",
            aligned_read="AAACCCGGGCCCGGGTTTAAACCC",
            aligned_reference="AAACCCGGG------TTTAAACCC",
            direction="forward",
        ),
    ]

    calls = call_exact_itds(alignments, reference, min_insert_length=3)

    assert [call.itd.insertion.read_id for call in calls] == [
        "earlier-itd",
        "later-itd",
    ]
