"""Build reportable ITD calls from aligned reads."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .coverage import interbase_coverage, variant_allele_frequency
from .insertions import Alignment, Insertion, extract_insertions
from .itds import ITD, classify_exact_itd


@dataclass(frozen=True)
class ITDCall:
    """An ITD with support, coverage, and allele frequency."""

    itd: ITD
    support_count: int
    coverage: int
    vaf: float


ITDCallKey = tuple[int, str, bool]


def call_exact_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
) -> list[ITDCall]:
    """Call exact-match ITDs and attach support, coverage, and VAF."""
    alignments = list(alignments)
    coverage_by_site = interbase_coverage(alignments)
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)

    for alignment in alignments:
        insertions = extract_insertions(alignment, min_length=min_insert_length)
        for insertion in insertions:
            itd = classify_exact_itd(insertion, reference)
            if itd is not None:
                grouped_itds[_itd_call_key(itd)].append(itd)

    calls: list[ITDCall] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        support_count = sum(itd.insertion.count for itd in itds)
        coverage = coverage_by_site.get(representative.insertion.start, 0)
        calls.append(
            ITDCall(
                itd=representative,
                support_count=support_count,
                coverage=coverage,
                vaf=variant_allele_frequency(support_count, coverage),
            )
        )

    return sorted(calls, key=_sort_key)


def _itd_call_key(itd: ITD) -> ITDCallKey:
    return (
        itd.tandem_start,
        itd.tandem_sequence,
        itd.insertion.trailing,
    )


def _representative_itd(itds: list[ITD]) -> ITD:
    return max(itds, key=lambda itd: itd.insertion.count)


def _sort_key(call: ITDCall) -> tuple[int, int, str]:
    return (
        call.itd.insertion.start,
        call.itd.tandem_start,
        call.itd.insertion.sequence,
    )
