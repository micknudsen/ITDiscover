"""Sequencing read preprocessing utilities."""

from collections.abc import Iterable
from dataclasses import dataclass

from .insertions import Direction
from .sequences import reverse_complement, validate_sequence


@dataclass(frozen=True)
class SequencingRead:
    """A sequencing read oriented to the forward WT reference."""

    read_id: str
    sequence: str
    qualities: tuple[int, ...]
    direction: Direction
    count: int = 1

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.qualities):
            raise ValueError("sequence and qualities must have equal length")
        validate_sequence(self.sequence)
        if self.count < 1:
            raise ValueError("count must be at least 1")

    @property
    def length(self) -> int:
        """Return the read sequence length."""
        return len(self.sequence)

    @property
    def mean_quality(self) -> float:
        """Return the mean Phred quality score."""
        if not self.qualities:
            return 0.0
        return sum(self.qualities) / len(self.qualities)


def orient_read(
    *,
    read_id: str,
    sequence: str,
    qualities: Iterable[int],
    direction: Direction,
) -> SequencingRead:
    """Return a read oriented to the forward WT reference.

    Reverse-direction reads are reverse-complemented, and their qualities are
    reversed to keep base qualities aligned with the oriented sequence.
    """
    qualities = tuple(qualities)
    if direction == "forward":
        return SequencingRead(
            read_id=read_id,
            sequence=sequence,
            qualities=qualities,
            direction=direction,
        )
    return SequencingRead(
        read_id=read_id,
        sequence=reverse_complement(sequence),
        qualities=tuple(reversed(qualities)),
        direction=direction,
    )


def trim_terminal_ns(read: SequencingRead) -> SequencingRead:
    """Trim ambiguous N bases from both ends of a read."""
    start = 0
    end = len(read.sequence)

    while start < end and read.sequence[start] == "N":
        start += 1
    while end > start and read.sequence[end - 1] == "N":
        end -= 1

    return SequencingRead(
        read_id=read.read_id,
        sequence=read.sequence[start:end],
        qualities=read.qualities[start:end],
        direction=read.direction,
        count=read.count,
    )


def passes_read_filters(
    read: SequencingRead,
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
) -> bool:
    """Return whether a read passes length and mean-quality filters."""
    if min_length < 0:
        raise ValueError("min_length must not be negative")
    if min_mean_quality < 0:
        raise ValueError("min_mean_quality must not be negative")
    return read.length >= min_length and read.mean_quality >= min_mean_quality


def collapse_identical_reads(reads: Iterable[SequencingRead]) -> list[SequencingRead]:
    """Collapse reads with identical oriented sequence and direction."""
    collapsed: dict[tuple[str, Direction], SequencingRead] = {}
    counts: dict[tuple[str, Direction], int] = {}

    for read in reads:
        key = (read.sequence, read.direction)
        if key not in collapsed:
            collapsed[key] = read
            counts[key] = 0
        counts[key] += read.count

    return [
        SequencingRead(
            read_id=read.read_id,
            sequence=read.sequence,
            qualities=read.qualities,
            direction=read.direction,
            count=counts[key],
        )
        for key, read in collapsed.items()
    ]


def preprocess_reads(
    reads: Iterable[SequencingRead],
    *,
    min_length: int = 100,
    min_mean_quality: float = 30,
) -> list[SequencingRead]:
    """Trim, filter, and collapse reads."""
    passing_reads = []
    for read in reads:
        trimmed = trim_terminal_ns(read)
        if passes_read_filters(
            trimmed,
            min_length=min_length,
            min_mean_quality=min_mean_quality,
        ):
            passing_reads.append(trimmed)
    return collapse_identical_reads(passing_reads)
