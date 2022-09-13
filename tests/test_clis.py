"""Test spycipdb client interfaces."""
import argparse
import inspect
import pytest

from spycipdb import __main__ as main_module
from spycipdb.clis import(
    # add new clis here and below
    cli,
    cli_cs,
    cli_jc,
    cli_noe,
    cli_pre,
    cli_rdc,
    cli_rh,
    cli_saxs,
    cli_smfret,    
    )

from spycipdb.libs import libcli as lc

subclients = [
    cli_cs,
    cli_jc,
    cli_noe,
    cli_pre,
    cli_rdc,
    cli_rh,
    cli_saxs,
    cli_smfret,
]

@pytest.fixture(params=subclients)
def client(request):
        """Loop all individual clients"""
        return request.param
    

@pytest.fixture(params=[cli] + subclients)
def all_clients(request):
    """Consider individual clients plus main cli"""
    return request.param


def test_main_module():
    """Test main module entry point."""
    assert hasattr(main_module, 'maincli')


def test_main_module_if():
    """Test cli main entry point."""
    lines = inspect.getsource(main_module).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''


def test_cli__version(all_clients):
    """Test clients have version parameter."""
    with pytest.raises(SystemExit) as err:
        all_clients.ap.parse_args(['-v'])
    assert err.value.code == 0


def test_cli_script_1():
    """Test cli main entry point."""
    assert hasattr(cli, 'maincli')


def test_cli_script_2():
    """Test cli has load_args function."""
    assert hasattr(cli, 'load_args')


def test_cli_script_3():
    """Test cli has _ap function."""
    assert hasattr(cli, '_ap')


def test_cli__ap_returns():
    """Test cli _ap returns ap."""
    assert isinstance(cli._ap(), argparse.ArgumentParser)


def test_cli_script_4():
    """Test cli main if __name__ part."""
    lines = inspect.getsource(cli).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''
