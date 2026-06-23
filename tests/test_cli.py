import pytest

import itdiscover
import itdiscover.cli as cli
from itdiscover.insertions import Alignment


def test_main_version(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"itdiscover {itdiscover.__version__}\n"


def test_main_requires_arguments(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main([])

    assert exc_info.value.code == 2
    assert "required" in capsys.readouterr().err


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
        "tandem_start\tinsertion_start\tsequence\tsupport_count\tunique_support_count\tcoverage\tvaf\n"
        "3\t2\tCCCGGG\t1\t1\t2\t0.500000\n"
    )


def test_call_command_rejects_multi_sequence_reference(tmp_path) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">ref1\nAAAA\n>ref2\nCCCC\n", encoding="utf-8")

    with pytest.raises(ValueError, match="exactly one sequence"):
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                "unused_R1.fastq",
                "--r2",
                "unused_R2.fastq",
            ]
        )


def test_call_command_writes_unique_support_alignment_html_report(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nTTTAAACCCGGGTTT\n", encoding="utf-8")
    report_path = tmp_path / "reports" / "unique-support.html"

    r1_path = tmp_path / "sample_R1.fastq"
    r1_path.write_text(
        (
            "@itd-fragment-1/1\n"
            "TTTAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@itd-fragment-2/1\n"
            "TTCAAACCCGGGCCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIIIIIIIII\n"
            "@wt-fragment/1\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    r2_path = tmp_path / "sample_R2.fastq"
    r2_path.write_text(
        (
            "@itd-fragment-1/2\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
            "@itd-fragment-2/2\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
            "@wt-fragment/2\n"
            "TTTAAACCCGGGTTT\n"
            "+\n"
            "IIIIIIIIIIIIIII\n"
        ),
        encoding="utf-8",
    )

    assert (
        cli.main(
            [
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
                "--output",
                str(report_path),
            ]
        )
        == 0
    )

    assert "unique_support_count" in capsys.readouterr().out
    report = report_path.read_text(encoding="utf-8")
    assert "<h1>Unique Support Alignments</h1>" in report
    assert "<h2>ITD 1</h2>" in report
    assert "Support Count" in report
    assert "Unique Support Count" in report
    assert "count 1" in report
    assert "TTCAAA[CCCGGG]CCCGGG" in report
    assert "TTTAAA[CCCGGG]CCCGGG" in report
    assert '<span class="diff">T</span>' in report


def test_call_command_rejects_non_html_output_path(tmp_path, capsys) -> None:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">FLT3\nAAACCCGGGTTT\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "--reference",
                str(reference_path),
                "--r1",
                "unused_R1.fastq",
                "--r2",
                "unused_R2.fastq",
                "--output",
                "report.txt",
            ]
        )

    assert exc_info.value.code == 2
    assert "must end with .html" in capsys.readouterr().err


def test_alignment_comparison_classes_ignore_leading_gap_shift() -> None:
    baseline = Alignment(
        read_id="baseline",
        fragment_id="baseline",
        read_sequence="GCAATTTAGGT",
        aligned_reference="GCAATTTAGGT",
        aligned_read="GCAATTTAGGT",
        direction="forward",
    )
    shifted = Alignment(
        read_id="shifted",
        fragment_id="shifted",
        read_sequence="AGCAATTTAGGT",
        aligned_reference="-GCAATTTAGGT",
        aligned_read="AGCAATTTAGGT",
        direction="forward",
    )

    classes = cli._alignment_comparison_classes(shifted, baseline)

    assert classes[0] == "insert"
    assert all(css_class is None for css_class in classes[1:])
