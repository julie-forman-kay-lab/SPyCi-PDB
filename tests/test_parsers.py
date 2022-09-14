"""Test parsing functions in SPyCi-PDB."""
from spycipdb.core.parsers import (
    get_exp_format_noe,
    get_exp_format_pre,
    get_exp_format_smfret,
    )

from . import fret_exp_expected, noe_exp_expected, pre_exp_expected


def test_get_exp_format_noe():
    """Test getting the format from noe exp file."""
    format = get_exp_format_noe(noe_exp_expected)
    expected = {
        'res1': [],
        'atom1': [],
        'atom1_multiple_assignments': [],
        'res2': [],
        'atom2': [],
        'atom2_multiple_assignments': [],
        }
    for key in format:
        # Purposefully not using `isinstance()` as I'd like to
        # check for the same keys as well ;)
        assert type(format[key]) == type(expected[key])  # noqa: E721


def test_get_exp_format_pre():
    """Test getting format from pre exp file."""
    format = get_exp_format_pre(pre_exp_expected)
    expected = {
        'res1': [],
        'atom1': [],
        'res2': [],
        'atom2': [],
        }
    for key in format:
        assert type(format[key]) == type(expected[key])  # noqa: E721


def test_get_exp_format_smfret():
    """Test getting format from smfret exp file."""
    format = get_exp_format_smfret(fret_exp_expected)
    expected = {
        'res1': [],
        'res2': [],
        'scale': [],
        }
    for key in format:
        assert type(format[key]) == type(expected[key])  # noqa: E721
