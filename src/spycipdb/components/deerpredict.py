# pylint: skip-file
# flake8: noqa
"""Helper script to run DEERPREdict calculations for the `pre` module."""
import pandas as pd
from DEERPREdict.PRE import PREpredict
from MDAnalysis import Universe


def calc_pre_predict(fexp, pdb):
    """Back calculates PRE data intensity ratios based on DEERPREdict."""
    ratios = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    atom1_name = exp.atom1.values
    res2 = exp.res2.values.astype(int)
    atom2_name = exp.atom2.values
    
    u = Universe(pdb)
    
    return pdb, ratios
    