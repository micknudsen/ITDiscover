"""FLT3 internal tandem duplication classification."""

from dataclasses import dataclass
from typing import Literal

from .insertions import Insertion

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
    layered on later without changing the call model.
    """
    reference = reference.upper()
    sequence = insertion.sequence.upper()
    if not sequence:
        return None

    upstream_start = insertion.start - len(sequence) + 1
    if upstream_start >= 0:
        upstream_sequence = reference[upstream_start : insertion.start + 1]
        if upstream_sequence == sequence:
            return ITD(
                insertion=insertion,
                tandem_start=upstream_start,
                tandem_sequence=upstream_sequence,
                orientation="upstream",
            )

    downstream_start = insertion.start + 1
    downstream_end = downstream_start + len(sequence)
    if downstream_end <= len(reference):
        downstream_sequence = reference[downstream_start:downstream_end]
        if downstream_sequence == sequence:
            return ITD(
                insertion=insertion,
                tandem_start=downstream_start,
                tandem_sequence=downstream_sequence,
                orientation="downstream",
            )

    return None
