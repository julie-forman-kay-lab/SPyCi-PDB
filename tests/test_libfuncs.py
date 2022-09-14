"""Test libfuncs."""
from math import isclose

from spycipdb.libs import libfuncs


def test_get_scalar():
    """Test get scalar function."""
    scalar_pve = libfuncs.get_scalar(0.1, 0.2, 0.3)
    scalar_nve = libfuncs.get_scalar(-0.1, -0.2, -0.3)
    
    assert isclose(scalar_pve, 0.37416, abs_tol=1e-5)
    assert isclose(scalar_nve, 0.37416, abs_tol=1e-5)
