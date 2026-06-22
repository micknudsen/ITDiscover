"""FASTQ file parsing."""

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .reads import Fragment, orient_read


@dataclass(frozen=True)
class FastqRecord:
    """A raw FASTQ record."""

    read_id: str
    fragment_id: str
    sequence: str
    qualities: tuple[int, ...]


def read_paired_fastq(
    forward_path: str | Path,
    reverse_path: str | Path,
) -> Iterator[Fragment]:
    """Yield paired-end fragments from R1 and R2 FASTQ files."""
    forward_records = _read_fastq_records(Path(forward_path))
    reverse_records = _read_fastq_records(Path(reverse_path))
    pair_number = 0

    while True:
        forward_record = next(forward_records, None)
        reverse_record = next(reverse_records, None)
        if forward_record is None and reverse_record is None:
            return

        pair_number += 1
        if forward_record is None or reverse_record is None:
            raise ValueError(f"FASTQ pair {pair_number} is incomplete")
        if forward_record.fragment_id != reverse_record.fragment_id:
            raise ValueError(
                f"FASTQ pair {pair_number} has mismatched fragment IDs: "
                f"{forward_record.fragment_id!r} != {reverse_record.fragment_id!r}"
            )

        yield Fragment(
            fragment_id=forward_record.fragment_id,
            forward_read=orient_read(
                read_id=forward_record.read_id,
                fragment_id=forward_record.fragment_id,
                sequence=forward_record.sequence,
                qualities=forward_record.qualities,
                direction="forward",
            ),
            reverse_read=orient_read(
                read_id=reverse_record.read_id,
                fragment_id=reverse_record.fragment_id,
                sequence=reverse_record.sequence,
                qualities=reverse_record.qualities,
                direction="reverse",
            ),
        )


def phred33_to_quality(quality_string: str) -> tuple[int, ...]:
    """Decode an ASCII Phred+33 quality string."""
    return tuple(ord(char) - 33 for char in quality_string)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(mode="rt", encoding="utf-8")


def _read_fastq_records(path: Path) -> Iterator[FastqRecord]:
    with _open_text(path) as handle:
        yield from _iter_fastq_records(handle)


def _iter_fastq_records(handle: TextIO) -> Iterator[FastqRecord]:
    record_number = 0
    while True:
        header = handle.readline()
        if header == "":
            return

        record_number += 1
        sequence = handle.readline()
        separator = handle.readline()
        qualities = handle.readline()
        if not sequence or not separator or not qualities:
            raise ValueError(f"FASTQ record {record_number} is incomplete")

        header = header.rstrip("\n\r")
        sequence = sequence.rstrip("\n\r")
        separator = separator.rstrip("\n\r")
        qualities = qualities.rstrip("\n\r")

        if not header.startswith("@"):
            raise ValueError(f"FASTQ record {record_number} header must start with @")
        if not separator.startswith("+"):
            raise ValueError(f"FASTQ record {record_number} separator must start with +")
        if len(sequence) != len(qualities):
            raise ValueError(
                f"FASTQ record {record_number} sequence and qualities "
                "must have equal length"
            )

        read_id = _read_id_from_header(header)
        yield FastqRecord(
            read_id=read_id,
            fragment_id=_fragment_id_from_read_id(read_id),
            sequence=sequence,
            qualities=phred33_to_quality(qualities),
        )


def _read_id_from_header(header: str) -> str:
    return header[1:].split()[0]


def _fragment_id_from_read_id(read_id: str) -> str:
    if read_id.endswith(("/1", "/2")):
        return read_id[:-2]
    return read_id
