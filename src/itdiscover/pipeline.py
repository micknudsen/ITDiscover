"""High-level in-memory ITD calling pipeline."""

from collections.abc import Iterable

from .alignment import AlignmentScoring, align_read_to_reference
from .calls import ITDCall, ITDFilter, call_exact_itds
from .reads import Fragment, preprocess_fragments
from .sequences import validate_sequence


def call_exact_itds_from_fragments(
    fragments: Iterable[Fragment],
    reference: str,
    *,
    min_read_length: int = 100,
    min_mean_quality: float = 30,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
    scoring: AlignmentScoring = AlignmentScoring(),
) -> list[ITDCall]:
    """Call exact-match ITDs from paired-end sequencing fragments."""
    validate_sequence(reference, field_name="reference")
    processed_reads = preprocess_fragments(
        fragments,
        min_length=min_read_length,
        min_mean_quality=min_mean_quality,
    )
    alignments = [
        align_read_to_reference(read, reference, scoring=scoring)
        for read in processed_reads
    ]
    return call_exact_itds(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        filters=filters,
    )
