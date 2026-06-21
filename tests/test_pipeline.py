import pytest

from itdiscover.calls import ITDCall
from itdiscover.insertions import Insertion
from itdiscover.itds import ITD
from itdiscover.pipeline import call_exact_itds_from_reads
from itdiscover.reads import SequencingRead


def test_call_exact_itds_from_reads_reports_support_coverage_and_vaf() -> None:
    reference = "AAACCCGGGTTT"
    reads = [
        SequencingRead(
            read_id="itd-read-1",
            sequence="AAACCCGGGCCCGGGTTT",
            qualities=(40,) * 18,
            direction="forward",
        ),
        SequencingRead(
            read_id="itd-read-2",
            sequence="AAACCCGGGCCCGGGTTT",
            qualities=(40,) * 18,
            direction="forward",
            count=2,
        ),
        SequencingRead(
            read_id="wt-read",
            sequence=reference,
            qualities=(40,) * 12,
            direction="reverse",
            count=7,
        ),
    ]

    assert call_exact_itds_from_reads(
        reads,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    ) == [
        ITDCall(
            itd=ITD(
                insertion=Insertion(
                    read_id="itd-read-1",
                    start=2,
                    sequence="CCCGGG",
                    direction="forward",
                    count=3,
                ),
                tandem_start=3,
                tandem_sequence="CCCGGG",
                orientation="downstream",
            ),
            support_count=3,
            coverage=10,
            vaf=0.3,
        )
    ]


def test_call_exact_itds_from_reads_filters_low_quality_reads() -> None:
    reference = "AAACCCGGGTTT"
    reads = [
        SequencingRead(
            read_id="low-quality-itd",
            sequence="AAACCCGGGCCCGGGTTT",
            qualities=(10,) * 18,
            direction="forward",
        ),
        SequencingRead(
            read_id="wt-read",
            sequence=reference,
            qualities=(40,) * 12,
            direction="forward",
            count=5,
        ),
    ]

    assert (
        call_exact_itds_from_reads(
            reads,
            reference,
            min_read_length=12,
            min_mean_quality=30,
        )
        == []
    )


def test_call_exact_itds_from_reads_trims_terminal_ns_before_calling() -> None:
    reference = "AAACCCGGGTTT"
    reads = [
        SequencingRead(
            read_id="itd-read",
            sequence="NAAACCCGGGCCCGGGTTTN",
            qualities=(5,) + (40,) * 18 + (5,),
            direction="forward",
        ),
        SequencingRead(
            read_id="wt-read",
            sequence=reference,
            qualities=(40,) * 12,
            direction="forward",
        ),
    ]

    calls = call_exact_itds_from_reads(
        reads,
        reference,
        min_read_length=12,
        min_mean_quality=30,
    )

    assert len(calls) == 1
    assert calls[0].support_count == 1
    assert calls[0].coverage == 2
    assert calls[0].vaf == 0.5


def test_call_exact_itds_from_reads_rejects_lowercase_reference() -> None:
    with pytest.raises(ValueError, match="reference contains invalid bases"):
        call_exact_itds_from_reads([], "AAAccc")
