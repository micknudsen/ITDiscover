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
    unique_support_count: int
    coverage: int
    vaf: float
    status: str = "PASS"
    filter_reasons: tuple[str, ...] = ()

    @property
    def passes_filters(self) -> bool:
        """Return whether the call passes the configured thresholds."""
        return self.status == "PASS"


@dataclass(frozen=True)
class ITDFilter:
    """Thresholds used to label exact-match ITD calls."""

    min_support_count: int = 1
    min_unique_support_count: int = 1
    min_coverage: int = 0
    min_vaf: float = 0.0

    def __post_init__(self) -> None:
        if self.min_support_count < 1:
            raise ValueError("min_support_count must be at least 1")
        if self.min_unique_support_count < 1:
            raise ValueError("min_unique_support_count must be at least 1")
        if self.min_coverage < 0:
            raise ValueError("min_coverage must not be negative")
        if self.min_vaf < 0:
            raise ValueError("min_vaf must not be negative")


@dataclass(frozen=True)
class UniqueSupportRepresentative:
    """One representative alignment for a unique local ITD support pattern."""

    itd: ITD
    signature: str
    alignment: Alignment
    support_count: int


ITDCallKey = tuple[int, str, bool]
SupportRepresentativeMap = dict[ITDCallKey, dict[str, UniqueSupportRepresentative]]


def call_exact_itds(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
) -> list[ITDCall]:
    """Call exact-match ITDs and attach support, coverage, and VAF."""
    calls, _ = call_exact_itds_with_representatives(
        alignments,
        reference,
        min_insert_length=min_insert_length,
        filters=filters,
    )
    return calls


def call_exact_itds_with_representatives(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int = 6,
    filters: ITDFilter = ITDFilter(),
) -> tuple[list[ITDCall], list[UniqueSupportRepresentative]]:
    """Call exact-match ITDs and retain one alignment per unique support pattern."""
    alignments = list(alignments)
    coverage_by_site = interbase_coverage(alignments)
    grouped_itds, representative_map = _collect_exact_itd_support(
        alignments,
        reference,
        min_insert_length=min_insert_length,
    )

    calls: list[ITDCall] = []
    representatives: list[UniqueSupportRepresentative] = []
    for itds in grouped_itds.values():
        representative = _representative_itd(itds)
        key = _itd_call_key(representative)
        support_count = len({itd.insertion.fragment_id for itd in itds})
        unique_support_count = len(representative_map[key])
        coverage = coverage_by_site.get(representative.insertion.start, 0)
        vaf = variant_allele_frequency(support_count, coverage)
        filter_reasons = _call_filter_reasons(
            support_count=support_count,
            unique_support_count=unique_support_count,
            coverage=coverage,
            vaf=vaf,
            filters=filters,
        )
        call = ITDCall(
            itd=representative,
            support_count=support_count,
            unique_support_count=unique_support_count,
            coverage=coverage,
            vaf=vaf,
            status="PASS" if not filter_reasons else "FAIL",
            filter_reasons=filter_reasons,
        )
        calls.append(call)
        representatives.extend(
            _sorted_representatives(representative_map[key].values())
        )

    calls.sort(key=_sort_key)
    representatives.sort(key=_representative_sort_key)
    return calls, representatives


def _itd_call_key(itd: ITD) -> ITDCallKey:
    return (
        itd.tandem_start,
        itd.tandem_sequence,
        itd.insertion.trailing,
    )


def _representative_itd(itds: list[ITD]) -> ITD:
    return itds[0]


def _sort_key(call: ITDCall) -> tuple[int, int, str]:
    return (
        call.itd.insertion.start,
        call.itd.tandem_start,
        call.itd.insertion.sequence,
    )


def _support_signature(
    alignment: Alignment,
    itd: ITD,
    reference: str,
    *,
    flank_size: int,
) -> str:
    observed_bases = _observed_bases_by_reference_position(alignment)
    left_start = max(0, itd.insertion.start - flank_size + 1)
    left = "".join(
        observed_bases.get(position, "-")
        for position in range(left_start, itd.insertion.start + 1)
    )
    right_end = min(len(reference), itd.tandem_start + flank_size)
    right = "".join(
        observed_bases.get(position, "-")
        for position in range(itd.tandem_start, right_end)
    )
    return f"{left}[{itd.tandem_sequence}]{right}"


def _observed_bases_by_reference_position(
    alignment: Alignment,
) -> dict[int, str]:
    ref_pos = -1
    observed_bases: dict[int, str] = {}

    for read_base, ref_base in zip(
        alignment.aligned_read,
        alignment.aligned_reference,
        strict=True,
    ):
        if ref_base != "-":
            ref_pos += 1
            if read_base != "-":
                observed_bases[ref_pos] = read_base

    return observed_bases


def _collect_exact_itd_support(
    alignments: Iterable[Alignment],
    reference: str,
    *,
    min_insert_length: int,
) -> tuple[dict[ITDCallKey, list[ITD]], SupportRepresentativeMap]:
    grouped_itds: dict[ITDCallKey, list[ITD]] = defaultdict(list)
    representative_map: dict[ITDCallKey, dict[str, Alignment]] = defaultdict(dict)
    fragment_ids_by_signature: dict[ITDCallKey, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )

    for alignment in alignments:
        insertions = extract_insertions(alignment, min_length=min_insert_length)
        for insertion in insertions:
            itd = classify_exact_itd(insertion, reference)
            if itd is None:
                continue

            key = _itd_call_key(itd)
            grouped_itds[key].append(itd)
            signature = _support_signature(
                alignment,
                itd,
                reference,
                flank_size=itd.length,
            )
            representative_map[key].setdefault(signature, alignment)
            fragment_ids_by_signature[key][signature].add(alignment.fragment_id)

    finalized_map: SupportRepresentativeMap = defaultdict(dict)
    for key, alignments_by_signature in representative_map.items():
        representative_itd = _representative_itd(grouped_itds[key])
        for signature, alignment in alignments_by_signature.items():
            finalized_map[key][signature] = UniqueSupportRepresentative(
                itd=representative_itd,
                signature=signature,
                alignment=alignment,
                support_count=len(fragment_ids_by_signature[key][signature]),
            )

    return grouped_itds, finalized_map


def _sorted_representatives(
    representatives: Iterable[UniqueSupportRepresentative],
) -> list[UniqueSupportRepresentative]:
    return sorted(representatives, key=_representative_sort_key)


def _representative_sort_key(
    representative: UniqueSupportRepresentative,
) -> tuple[int, int, str, int, str, str]:
    return (
        representative.itd.insertion.start,
        representative.itd.tandem_start,
        representative.itd.tandem_sequence,
        -representative.support_count,
        representative.signature,
        representative.alignment.read_id,
    )


def _call_filter_reasons(
    *,
    support_count: int,
    unique_support_count: int,
    coverage: int,
    vaf: float,
    filters: ITDFilter,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if support_count < filters.min_support_count:
        reasons.append("LOW_SUPPORT")
    if unique_support_count < filters.min_unique_support_count:
        reasons.append("LOW_UNIQUE_SUPPORT")
    if coverage < filters.min_coverage:
        reasons.append("LOW_COVERAGE")
    if vaf < filters.min_vaf:
        reasons.append("LOW_VAF")
    return tuple(reasons)
