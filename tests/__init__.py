"""Contains general variables for spycipdb unit tests."""

from spycipdb.libs import Path

datap = Path(
    Path(__file__).myparents(),
    'data',
    )

jc_exp_expected = Path(datap, 'drksh3_JC.txt')
noe_exp_expected = Path(datap, 'drksh3_NOE.txt')
pre_exp_expected = Path(datap, 'drksh3_PRE.txt')
fret_exp_expected = Path(datap, 'asyn_FRET.txt')

jc_output = Path(datap, 'drksh3_jc_out.json')
noe_output = Path(datap, 'drksh3_noe_out.json')
pre_output = Path(datap, 'drksh3_pre_out.json')
fret_output = Path(datap, 'asyn_fret_out.json')

drk_test = Path(datap, 'drksh3_conf.pdb')
asyn_test = Path(datap, 'asyn_conf.pdb')
