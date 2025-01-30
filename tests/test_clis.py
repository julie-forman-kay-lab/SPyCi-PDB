"""Test internal spycipdb client interfaces."""
import inspect
from pathlib import Path

import pytest

# import bellow by alphabetical order the cli interfaces implemented
from spycipdb.clis import cli_jc, cli_noe, cli_smfret

from . import (
    asyn_test_tar,
    drk_test_tar,
    fret_exp_expected,
    jc_exp_expected,
    noe_exp_expected,
    )


subclients = [
    # add your new client to this list.
    cli_jc,
    cli_noe,
    cli_smfret,
    ]


@pytest.fixture(params=subclients)
def client(request):
    """Loop all individual clients."""
    return request.param


def test_clients_have_main(client):
    """Test all modules have main."""
    assert hasattr(client, 'main')


def test_clients_have_ap(client):
    """Test all modules have ap."""
    assert hasattr(client, 'ap')


def test_clients_have_help(client):
    """Test all modules have help."""
    assert hasattr(client, '_help')


def test_clients_if_main_code(client):
    """Test if __name__ == '__main__' code lines."""
    lines = inspect.getsource(client).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    libcli.maincli(ap, main)"
    assert lines[-1] == ''


@pytest.mark.parametrize(
    'module,name',
    [
        # add here the name of your client, this ensures that changes
        # in the command interface get noticed by the tests.
        # (cli_NAME, 'NAME'),
        (cli_jc, 'jc'),
        (cli_noe, 'noe'),
        (cli_smfret, 'smfret'),
        ],
    )
def test_cli__name(module, name):
    """Test name messages."""
    assert module._name == name


# tox -e test hanging up after trying to read paths...
def test_cli_jc():
    """Test jc module."""
    cli_jc.main(
        str(drk_test_tar),
        str(jc_exp_expected),
        output='jc_output.json',
        )
    o = Path('jc_output.json')
    debug = Path('.spycipdb_jc.debug')
    error = Path('.spycipdb_jc.error')
    log = Path('.spycipdb_jc.log')
    assert o.exists()
    assert debug.exists()
    assert error.exists()
    assert log.exists()
    o.unlink()
    debug.unlink()
    error.unlink()
    log.unlink()


def test_cli_noe():
    """Test noe module."""
    cli_noe.main(
        str(drk_test_tar),
        str(noe_exp_expected),
        output='noe_output.json',
        )
    o = Path('noe_output.json')
    debug = Path('.spycipdb_noe.debug')
    error = Path('.spycipdb_noe.error')
    log = Path('.spycipdb_noe.log')
    assert o.exists()
    assert debug.exists()
    assert error.exists()
    assert log.exists()
    o.unlink()
    debug.unlink()
    error.unlink()
    log.unlink()


def test_cli_smfret():
    """Test smfret module."""
    cli_smfret.main(
        str(asyn_test_tar),
        str(fret_exp_expected),
        output='fret_output.json',
        )
    o = Path('fret_output.json')
    debug = Path('.spycipdb_smfret.debug')
    error = Path('.spycipdb_smfret.error')
    log = Path('.spycipdb_smfret.log')
    assert o.exists()
    assert debug.exists()
    assert error.exists()
    assert log.exists()
    o.unlink()
    debug.unlink()
    error.unlink()
    log.unlink()
