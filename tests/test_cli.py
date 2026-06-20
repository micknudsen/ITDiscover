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
