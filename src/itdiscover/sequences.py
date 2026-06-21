"""Shared DNA sequence validation and transformation helpers."""

VALID_BASES = frozenset("ACGTN")
VALID_ALIGNMENT_CHARS = VALID_BASES | frozenset("-")
DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def validate_sequence(
    sequence: str,
    *,
    valid_chars: frozenset[str] = VALID_BASES,
    field_name: str = "sequence",
) -> None:
    """Raise ValueError if a sequence contains characters outside an alphabet."""
    if invalid_chars := sorted(set(sequence) - valid_chars):
        raise ValueError(f"{field_name} contains invalid bases: {''.join(invalid_chars)}")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an uppercase DNA sequence."""
    validate_sequence(sequence)
    return sequence.translate(DNA_COMPLEMENT)[::-1]
