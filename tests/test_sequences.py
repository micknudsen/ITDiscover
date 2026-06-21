import pytest

from itdiscover.sequences import (
    VALID_ALIGNMENT_CHARS,
    VALID_BASES,
    reverse_complement,
    validate_sequence,
)


def test_valid_base_alphabet() -> None:
    assert VALID_BASES == frozenset("ACGTN")
    assert VALID_ALIGNMENT_CHARS == frozenset("ACGTN-")


def test_validate_sequence_accepts_uppercase_dna() -> None:
    validate_sequence("ACGTN")


def test_validate_sequence_rejects_lowercase_bases() -> None:
    with pytest.raises(ValueError, match="read contains invalid bases"):
        validate_sequence("ACgTN", field_name="read")


def test_reverse_complement_preserves_ambiguous_bases() -> None:
    assert reverse_complement("ACGTN") == "NACGT"
