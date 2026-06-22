import gzip

import pytest

from itdiscover.fastq import phred33_to_quality, read_paired_fastq
from itdiscover.reads import Fragment, SequencingRead


def test_phred33_to_quality_decodes_ascii_qualities() -> None:
    assert phred33_to_quality("!?I") == (0, 30, 40)


def test_read_paired_fastq_reads_plain_records(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text(
        "@fragment-1/1 extra fields\n"
        "ACGTN\n"
        "+\n"
        "!?I55\n",
        encoding="utf-8",
    )
    r2_path.write_text(
        "@fragment-1/2 extra fields\n"
        "ACGTT\n"
        "+\n"
        "!5?II\n",
        encoding="utf-8",
    )

    assert list(read_paired_fastq(r1_path, r2_path)) == [
        Fragment(
            fragment_id="fragment-1",
            forward_read=SequencingRead(
                read_id="fragment-1/1",
                fragment_id="fragment-1",
                sequence="ACGTN",
                qualities=(0, 30, 40, 20, 20),
                direction="forward",
            ),
            reverse_read=SequencingRead(
                read_id="fragment-1/2",
                fragment_id="fragment-1",
                sequence="AACGT",
                qualities=(40, 40, 30, 20, 0),
                direction="reverse",
            ),
        )
    ]


def test_read_paired_fastq_reads_gzipped_records(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq.gz"
    r2_path = tmp_path / "reads_R2.fastq.gz"
    with gzip.open(r1_path, mode="wt", encoding="utf-8") as handle:
        handle.write("@fragment-1/1\nACGT\n+\nIIII\n")
    with gzip.open(r2_path, mode="wt", encoding="utf-8") as handle:
        handle.write("@fragment-1/2\nACGT\n+\nIIII\n")

    fragments = list(read_paired_fastq(r1_path, r2_path))

    assert len(fragments) == 1
    assert fragments[0].fragment_id == "fragment-1"
    assert fragments[0].forward_read.sequence == "ACGT"
    assert fragments[0].reverse_read.sequence == "ACGT"


def test_read_paired_fastq_rejects_incomplete_record(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("@fragment-1/1\nACGT\n+\n", encoding="utf-8")
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="record 1 is incomplete"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_mismatched_fragments(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("@fragment-1/1\nACGT\n+\nIIII\n", encoding="utf-8")
    r2_path.write_text("@fragment-2/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="mismatched fragment IDs"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_incomplete_pairs(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text(
        "@fragment-1/1\nACGT\n+\nIIII\n"
        "@fragment-2/1\nACGT\n+\nIIII\n",
        encoding="utf-8",
    )
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="pair 2 is incomplete"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_invalid_header(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("fragment-1/1\nACGT\n+\nIIII\n", encoding="utf-8")
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="header must start with @"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_invalid_separator(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("@fragment-1/1\nACGT\n-\nIIII\n", encoding="utf-8")
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="separator must start with \\+"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_mismatched_sequence_and_quality_lengths(
    tmp_path,
) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("@fragment-1/1\nACGT\n+\nIII\n", encoding="utf-8")
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must have equal length"):
        list(read_paired_fastq(r1_path, r2_path))


def test_read_paired_fastq_rejects_lowercase_bases(tmp_path) -> None:
    r1_path = tmp_path / "reads_R1.fastq"
    r2_path = tmp_path / "reads_R2.fastq"
    r1_path.write_text("@fragment-1/1\nACgT\n+\nIIII\n", encoding="utf-8")
    r2_path.write_text("@fragment-1/2\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid bases"):
        list(read_paired_fastq(r1_path, r2_path))
