"""FASTQ file parsing."""

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

from .insertions import Direction
from .reads import SequencingRead, orient_read


def read_fastq(path: str | Path, *, direction: Direction) -> Iterator[SequencingRead]:
    """Yield oriented reads from a plain or gzipped four-line FASTQ file."""
    path = Path(path)
    with _open_text(path) as handle:
        yield from _iter_fastq_records(handle, direction=direction)


def phred33_to_quality(quality_string: str) -> tuple[int, ...]:
    """Decode an ASCII Phred+33 quality string."""
    return tuple(ord(char) - 33 for char in quality_string)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(mode="rt", encoding="utf-8")


def _iter_fastq_records(
    handle: TextIO,
    *,
    direction: Direction,
) -> Iterator[SequencingRead]:
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

        yield orient_read(
            read_id=header[1:].split()[0],
            sequence=sequence,
            qualities=phred33_to_quality(qualities),
            direction=direction,
        )
