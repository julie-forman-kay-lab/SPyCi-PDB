"""Test internal calculator functions."""
import json

from spycipdb.core.calculators import calc_jc, calc_noe, calc_pre, calc_smfret

from . import (
    asyn_test,
    drk_test,
    fret_exp_expected,
    fret_output,
    jc_exp_expected,
    jc_output,
    noe_exp_expected,
    noe_output,
    pre_exp_expected,
    pre_output,
    )


def test_calc_jc():
    """Test the internal `calc_jc` module."""
    _pdb, jc_bc = calc_jc(jc_exp_expected, drk_test)
    with open(jc_output, 'r') as f:
        loadedf = json.load(f)
        expected = list(loadedf.values())[0]
        assert expected == jc_bc


def test_calc_noe():
    """Test the internal `calc_noe` module."""
    _pdb, noe_bc = calc_noe(noe_exp_expected, drk_test)
    with open(noe_output, 'r') as f:
        loadedf = json.load(f)
        expected = list(loadedf.values())[0]
        assert expected == noe_bc


def test_calc_pre():
    """Test the internal `calc_pre` module."""
    _pdb, pre_bc = calc_pre(pre_exp_expected, drk_test)
    with open(pre_output, 'r') as f:
        loadedf = json.load(f)
        expected = list(loadedf.values())[0]
        assert expected == pre_bc


def test_calc_smfret():
    """Test the internal `calc_smfret` module."""
    _pdb, fret_bc = calc_smfret(fret_exp_expected, asyn_test)
    with open(fret_output, 'r') as f:
        loadedf = json.load(f)
        expected = list(loadedf.values())[0]
        assert expected == fret_bc
