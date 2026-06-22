import pytest

import itdiscover
import itdiscover.cli as cli


def test_main_version(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"itdiscover {itdiscover.__version__}\n"


def test_main_requires_command(capsys) -> None:
    assert cli.main([]) == 0
    assert capsys.readouterr().err == ""


def test_call_command_reports_exact_itd_from_paired_fastq(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment/1\n"
            "AAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "AAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
                "call",
                "--reference",
                str(reference_path),
                "--r1",
                str(r1_path),
                "--r2",
                str(r2_path),
                "--min-read-length",
                "12",
                "--min-mean-quality",
                "30",
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == (
        "tandem_start\tinsertion_start\tsequence\tsupport_count\tcoverage\tvaf\n"
        "3\t2\tCCCGGG\t1\t2\t0.500000\n"
    )


def test_call_command_rejects_multi_sequence_reference(tmp_path) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">ref1\nAAAA\n>ref2\nCCCC\n", encoding="utf-8")

    with pytest.raises(ValueError, match="exactly one sequence"):
        cli.main(
            [
                "call",
                "--reference",
                str(reference_path),
                "--r1",
                "unused_R1.fastq",
                "--r2",
                "unused_R2.fastq",
            ]
        )
