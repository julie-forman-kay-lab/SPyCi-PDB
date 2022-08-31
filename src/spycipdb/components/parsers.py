"""Contains functions to retrieve formats of experimental files for SPyCi-PDB"""
import pandas as pd


def get_exp_format_noe(fexp):
    """Get format from experimental template."""
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['atom1'] = exp.atom1.values.tolist()
    format['atom1_multiple_assignments'] = exp.atom1_multiple_assignments.values.tolist()  # noqa: E501
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['atom2'] = exp.atom2.values.tolist()
    format['atom2_multiple_assignments'] = exp.atom2_multiple_assignments.values.tolist()  # noqa: E501
    
    return format


def get_exp_format_pre(fexp):
    """Get format based on experimental file."""
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['atom1'] = exp.atom1.values.tolist()
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['atom2'] = exp.atom2.values.tolist()
    
    return format


def get_exp_format_smfret(fexp):
    """Get format from experimental template."""
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['scale'] = exp.scale.values.tolist()
    
    return format
