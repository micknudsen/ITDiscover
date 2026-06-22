"""FLT3 internal tandem duplication classification."""

from dataclasses import replace
from dataclasses import dataclass
from typing import Literal

from .insertions import Insertion
from .sequences import validate_sequence

TandemOrientation = Literal["upstream", "downstream"]


@dataclass(frozen=True)
class ITD:
    """An insertion classified as an internal tandem duplication."""

    insertion: Insertion
    tandem_start: int
    tandem_sequence: str
    orientation: TandemOrientation

    @property
    def tandem_end(self) -> int:
        """Return the inclusive end coordinate of the duplicated WT segment."""
        return self.tandem_start + len(self.tandem_sequence) - 1

    @property
    def length(self) -> int:
        """Return the duplicated sequence length."""
        return len(self.tandem_sequence)


def classify_exact_itd(insertion: Insertion, reference: str) -> ITD | None:
    """Classify an insertion as an exact-match tandem duplication.

    This first classifier only accepts insertions that exactly duplicate the
    adjacent upstream or downstream WT reference segment. Fuzzy matching can be
    layered on later without changing the call model. Exact ITDs are reported
    in a canonical downstream representation, choosing the right-most
    equivalent breakpoint in repetitive sequence contexts.
    """
    _validate_reference(reference)
    sequence = insertion.sequence
    if not sequence:
        return None

    if not _has_adjacent_exact_match(insertion, reference):
        return None

    return _canonical_downstream_itd(insertion, reference)


def _validate_reference(reference: str) -> None:
    validate_sequence(reference, field_name="reference")


def _has_adjacent_exact_match(insertion: Insertion, reference: str) -> bool:
    sequence = insertion.sequence

    upstream_start = insertion.start - len(sequence) + 1
    if upstream_start >= 0:
        upstream_sequence = reference[upstream_start : insertion.start + 1]
        if upstream_sequence == sequence:
            return True

    downstream_start = insertion.start + 1
    downstream_end = downstream_start + len(sequence)
    if downstream_end <= len(reference):
        downstream_sequence = reference[downstream_start:downstream_end]
        if downstream_sequence == sequence:
            return True

    return False


def _canonical_downstream_itd(insertion: Insertion, reference: str) -> ITD:
    sequence = insertion.sequence
    mutant_sequence = _mutant_sequence(reference, insertion.start, sequence)
    tandem_starts = [
        tandem_start
        for tandem_start in range(len(reference) - len(sequence) + 1)
        if reference[tandem_start : tandem_start + len(sequence)] == sequence
        and _mutant_sequence(reference, tandem_start - 1, sequence) == mutant_sequence
    ]

    tandem_start = max(tandem_starts)
    return ITD(
        insertion=replace(insertion, start=tandem_start - 1),
        tandem_start=tandem_start,
        tandem_sequence=sequence,
        orientation="downstream",
    )


def _mutant_sequence(reference: str, insertion_start: int, sequence: str) -> str:
    insertion_index = insertion_start + 1
    return reference[:insertion_index] + sequence + reference[insertion_index:]
