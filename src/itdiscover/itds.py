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


@dataclass(frozen=True)
class TandemSimilarity:
    """Similarity between an inserted sequence and a WT tandem window."""

    insertion: Insertion
    tandem_start: int
    tandem_sequence: str
    mismatches: int

    @property
    def matches(self) -> int:
        """Return the number of matching bases in the best window."""
        return len(self.tandem_sequence) - self.mismatches

    @property
    def identity(self) -> float:
        """Return the fraction of matching bases in the best window."""
        if not self.tandem_sequence:
            return 0.0
        return self.matches / len(self.tandem_sequence)


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


def score_tandem_similarity(
    insertion: Insertion,
    reference: str,
) -> TandemSimilarity | None:
    """Score how well an insertion matches any WT tandem window.

    This is a minimal fuzzy-matching primitive for later ITD classification.
    It searches all WT windows of the same length as the insertion and returns
    the best match, preferring the right-most window when there is a tie.
    """
    _validate_reference(reference)
    sequence = insertion.sequence
    if not sequence:
        return None
    if len(sequence) > len(reference):
        return None

    best_start = None
    best_mismatches = None
    for tandem_start in range(len(reference) - len(sequence) + 1):
        tandem_sequence = reference[tandem_start : tandem_start + len(sequence)]
        mismatches = sum(
            1 for observed, expected in zip(sequence, tandem_sequence, strict=True)
            if observed != expected
        )
        key = (mismatches, -tandem_start)
        if best_start is None or key < (best_mismatches, -best_start):
            best_start = tandem_start
            best_mismatches = mismatches

    assert best_start is not None
    assert best_mismatches is not None
    return TandemSimilarity(
        insertion=insertion,
        tandem_start=best_start,
        tandem_sequence=reference[best_start : best_start + len(sequence)],
        mismatches=best_mismatches,
    )


def classify_fuzzy_itd(
    insertion: Insertion,
    reference: str,
    *,
    max_mismatches: int,
) -> ITD | None:
    """Classify an insertion as an adjacent fuzzy-match tandem duplication.

    The best same-length WT window must be adjacent to the insertion and must
    be within the allowed mismatch threshold. If both adjacent windows are tied
    for best mismatch score, classification is rejected as ambiguous.
    """
    if max_mismatches < 0:
        raise ValueError("max_mismatches must not be negative")

    similarity = score_tandem_similarity(insertion, reference)
    if similarity is None or similarity.mismatches > max_mismatches:
        return None

    if not _is_adjacent_tandem_start(insertion, similarity.tandem_start):
        return None

    adjacent_candidates = _adjacent_tandem_candidates(insertion, reference)
    best_adjacent_mismatches = min(candidate[2] for candidate in adjacent_candidates)
    if sum(
        1
        for _, _, mismatches in adjacent_candidates
        if mismatches == best_adjacent_mismatches
    ) > 1:
        return None

    orientation: TandemOrientation
    if similarity.tandem_start == insertion.start + 1:
        orientation = "downstream"
    else:
        orientation = "upstream"

    return ITD(
        insertion=insertion,
        tandem_start=similarity.tandem_start,
        tandem_sequence=similarity.tandem_sequence,
        orientation=orientation,
    )


def _validate_reference(reference: str) -> None:
    validate_sequence(reference, field_name="reference")


def _is_adjacent_tandem_start(insertion: Insertion, tandem_start: int) -> bool:
    sequence_length = len(insertion.sequence)
    return tandem_start in (
        insertion.start - sequence_length + 1,
        insertion.start + 1,
    )


def _adjacent_tandem_candidates(
    insertion: Insertion,
    reference: str,
) -> list[tuple[TandemOrientation, str, int]]:
    sequence = insertion.sequence
    candidates: list[tuple[TandemOrientation, str, int]] = []

    upstream_start = insertion.start - len(sequence) + 1
    if upstream_start >= 0:
        upstream_sequence = reference[upstream_start : insertion.start + 1]
        candidates.append(
            (
                "upstream",
                upstream_sequence,
                _mismatch_count(sequence, upstream_sequence),
            )
        )

    downstream_start = insertion.start + 1
    downstream_end = downstream_start + len(sequence)
    if downstream_end <= len(reference):
        downstream_sequence = reference[downstream_start:downstream_end]
        candidates.append(
            (
                "downstream",
                downstream_sequence,
                _mismatch_count(sequence, downstream_sequence),
            )
        )

    return candidates


def _mismatch_count(observed: str, expected: str) -> int:
    return sum(
        1 for observed_base, expected_base in zip(observed, expected, strict=True)
        if observed_base != expected_base
    )


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
