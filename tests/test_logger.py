"""Test logging functions of spycipdb."""
from functools import partial

import pytest

from spycipdb import Path, log
from spycipdb.core.exceptions import ReportOnCrashError
from spycipdb.logger import S, T, init_files, report_on_crash


def test_init_files():
    """Test init log files."""
    init_files(log, '.dummy')
    paths = [Path('.dummy').with_suffix(p) for p in ['.log', '.error', '.debug']]
    assert all(p.exists() for p in paths)
    try:
        for p in paths:
            p.unlink()
    # Should only be in the case for windows test cases
    except PermissionError:
        assert True


def test_T():
    """Test T formatter."""
    logmsg = T('my title {}', 'IDP')
    assert str(logmsg) == 'My Title IDP:'


@pytest.mark.parametrize(
    'msg,args,spacer,indent,expected',
    [
        (
            'a log message with param {}',
            'IDP',
            '+',
            8,
            '++++++++a log message with param IDP',
            ),
        ('string {}', (), '', 1, 'string ()'),
        ('string', (), ' ', 4, '    string'),
        ],
    )
def test_S(msg, args, spacer, indent, expected):
    """Test S formatter."""
    sobj = S(msg, args, spacer=spacer, indent=indent)
    assert str(sobj) == expected


def test_report_on_crash():
    """Test record func on error to file."""
    def funca(a, b, c=1, d=2):
        raise TypeError

    ext = 'testing_ROC'
    with pytest.raises(ReportOnCrashError):
        report_on_crash(
            funca,
            'spycipdb', c=range(10), d=dict.fromkeys('qwerty'),
            ROC_exception=TypeError,
            ROC_ext=ext,
            )

    errfiles = list(Path.cwd().glob(f'*.{ext}'))
    assert len(errfiles) > 0
    for p in errfiles:
        p.unlink()


def test_report_on_crash_partial():
    """Test record func on error to file."""
    def funca(a, b, c=1, d=2):
        raise TypeError

    funcb = partial(funca, 58)

    ext = 'testing_ROC'
    with pytest.raises(ReportOnCrashError):
        report_on_crash(
            funcb,
            'spycipdb', c=range(10), d=dict.fromkeys('qwerty'),
            ROC_exception=TypeError,
            ROC_ext=ext,
            )

    errfiles = list(Path.cwd().glob(f'*.{ext}'))
    assert len(errfiles) > 0
    for p in errfiles:
        p.unlink()
