"""Test internal calculator functions."""
# `noqa` temporary for now
import json  # noqa: F401

from spycipdb.core.calculators import (  # noqa: F401
    calc_jc,
    calc_noe,
    calc_pre,
    calc_smfret,
    )

from . import (  # noqa: F401
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


# Not working right now due to idpconfgen parsing path?
# def test_calc_jc():
#     pdb, jc_bc = calc_jc(jc_exp_expected, drk_test)
#     with open(jc_output, 'r') as f:
#         expected = json.load(f)
