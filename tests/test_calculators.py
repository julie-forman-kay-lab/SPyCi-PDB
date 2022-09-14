"""Test internal calculator functions."""
import json

from spycipdb.core.calculators import (
    calc_jc,
    calc_noe,
    calc_pre,
    calc_smfret,
    )

from . import (
    asyn_test,
    drk_test,
    jc_exp_expected,
    noe_exp_expected,
    pre_exp_expected,
    fret_exp_expected,
    jc_output,
    noe_output,
    pre_output,
    fret_output,
    )


# Not working right now due to idpconfgen parsing path?
# def test_calc_jc():
#     pdb, jc_bc = calc_jc(jc_exp_expected, drk_test)
#     with open(jc_output, 'r') as f:
#         expected = json.load(f)

        