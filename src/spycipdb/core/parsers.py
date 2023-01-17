"""
Functions to retrieve formats of experimental files for SPyCi-PDB.

These functions will also double as validators to check against a
sample PDB input and the data in question.

I.e. Sees if the experimental file has any negative values for residue
numbers and if there are inconsistencies between atom names.
"""
import pandas as pd

from spycipdb.core.exceptions import SPyCiPDBException


def get_exp_format_noe(fexp):
    """Get format from experimental template."""
    format = {}
    exp = pd.read_csv(fexp)
    
    try:
        format['res1'] = exp.res1.values.astype(int).tolist()
        format['atom1'] = exp.atom1.values.tolist()
        format['atom1_multiple_assignments'] = exp.atom1_multiple_assignments.values.tolist()  # noqa: E501
        format['res2'] = exp.res2.values.astype(int).tolist()
        format['atom2'] = exp.atom2.values.tolist()
        format['atom2_multiple_assignments'] = exp.atom2_multiple_assignments.values.tolist()  # noqa: E501
    except AttributeError as err:
        errmsg = (
            'Incorrect experimental file format for NOE subclient. '
            'Text file must have the following columns: '
            'res1,atom1,atom1_multiple_assignments,res2,atom2,atom2_multiple_assignments'  # noqa: E501
            )
        raise SPyCiPDBException(errmsg) from err
    
    return format


def get_exp_format_pre(fexp):
    """Get format based on experimental file."""
    format = {}
    exp = pd.read_csv(fexp)
    
    try:
        format['res1'] = exp.res1.values.astype(int).tolist()
        format['atom1'] = exp.atom1.values.tolist()
        format['res2'] = exp.res2.values.astype(int).tolist()
        format['atom2'] = exp.atom2.values.tolist()
    except AttributeError as err:
        errmsg = (
            'Incorrect experimental file format for PRE subclient. '
            'Text file must have the following columns: '
            'res1,atom1,res2,atom2'
            )
        raise SPyCiPDBException(errmsg) from err
    
    return format


def get_exp_format_smfret(fexp):
    """Get format from experimental template."""
    format = {}
    exp = pd.read_csv(fexp)
    
    try:
        format['res1'] = exp.res1.values.astype(int).tolist()
        format['res2'] = exp.res2.values.astype(int).tolist()
        format['scale'] = exp.scale.values.tolist()
    except AttributeError as err:
        errmsg = (
            'Incorrect experimental file format for smFRET subclient. '
            'Text file must have the following columns: '
            'res1,res2,scaler'
            )
        raise SPyCiPDBException(errmsg) from err
    
    return format
