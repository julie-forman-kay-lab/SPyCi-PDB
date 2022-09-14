"""Test custom spycipdb exceptions."""
import inspect

from hypothesis import given
from hypothesis import strategies as st

from spycipdb.core import exceptions as EXCPTNS


EXCPT_classes = inspect.getmembers(EXCPTNS, predicate=inspect.isclass)
error_classes = [
    t[1] for t in EXCPT_classes
    if issubclass(t[1], EXCPTNS.SPyCiPDBException)
    and t[0].endswith('Error')
    ]


def test_all_errors_names_end_in_error():
    """Test whether all custom error classes end with 'Error'."""
    endswitherror = [t[1] for t in EXCPT_classes if t[0].endswith('Error')]
    subclss = [
        t[1] for t in EXCPT_classes
        if issubclass(t[1], EXCPTNS.SPyCiPDBException)
        ]
    assert len(endswitherror) == len(subclss) - 1


def test_SPyCiPDBException_type():
    """Test IDPConfGenException is Exception."""
    assert issubclass(EXCPTNS.SPyCiPDBException, Exception)


def test_SPyCiPDBException_no_error_mg_0():
    """Test clean instation gives error msg."""
    errmsg = EXCPTNS.SPyCiPDBException.errmsg
    assert str(EXCPTNS.SPyCiPDBException()) == errmsg


@given(st.none())
def test_SPyCiPDB_errmsg_None(errmsg):
    """Test init with errmsg=None, should be ignored."""
    err = EXCPTNS.SPyCiPDBException(errmsg=errmsg)
    assert str(err) == EXCPTNS.SPyCiPDBException.errmsg
