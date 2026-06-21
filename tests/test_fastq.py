import gzip

import pytest

from itdiscover.fastq import phred33_to_quality, read_fastq
from itdiscover.reads import SequencingRead


def test_phred33_to_quality_decodes_ascii_qualities() -> None:
    assert phred33_to_quality("!?I") == (0, 30, 40)


def test_read_fastq_reads_plain_forward_records(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text(
        "@read-1 extra fields\n"
        "ACGTN\n"
        "+\n"
        "!?I55\n"
        "@read-2\n"
        "TTAA\n"
        "+\n"
        "IIII\n",
        encoding="utf-8",
    )

    assert list(read_fastq(fastq_path, direction="forward")) == [
        SequencingRead(
            read_id="read-1",
            sequence="ACGTN",
            qualities=(0, 30, 40, 20, 20),
            direction="forward",
        ),
        SequencingRead(
            read_id="read-2",
            sequence="TTAA",
            qualities=(40, 40, 40, 40),
            direction="forward",
        ),
    ]


def test_read_fastq_reads_gzipped_reverse_records(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq.gz"
    with gzip.open(fastq_path, mode="wt", encoding="utf-8") as handle:
        handle.write("@read-1\nACGTT\n+\n!5?II\n")

    assert list(read_fastq(fastq_path, direction="reverse")) == [
        SequencingRead(
            read_id="read-1",
            sequence="AACGT",
            qualities=(40, 40, 30, 20, 0),
            direction="reverse",
        )
    ]


def test_read_fastq_rejects_incomplete_record(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text("@read-1\nACGT\n+\n", encoding="utf-8")

    with pytest.raises(ValueError, match="record 1 is incomplete"):
        list(read_fastq(fastq_path, direction="forward"))


def test_read_fastq_rejects_invalid_header(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text("read-1\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="header must start with @"):
        list(read_fastq(fastq_path, direction="forward"))


def test_read_fastq_rejects_invalid_separator(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text("@read-1\nACGT\n-\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="separator must start with \\+"):
        list(read_fastq(fastq_path, direction="forward"))


def test_read_fastq_rejects_mismatched_sequence_and_quality_lengths(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text("@read-1\nACGT\n+\nIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must have equal length"):
        list(read_fastq(fastq_path, direction="forward"))


def test_read_fastq_rejects_lowercase_bases(tmp_path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text("@read-1\nACgT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid bases"):
        list(read_fastq(fastq_path, direction="forward"))
