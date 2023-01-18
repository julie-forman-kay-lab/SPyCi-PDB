"""
Functions to retrieve formats of experimental files for SPyCi-PDB.

These functions will also double as validators to check against a
sample PDB input and the data in question.

I.e. Sees if the experimental file has any negative values for residue
numbers and if there are inconsistencies between atom names.
"""
import numpy as np
import pandas as pd
from idpconfgen.libs.libstructure import Structure, col_name, col_resSeq

from spycipdb.core.exceptions import SPyCiPDBException


def get_exp_format_noe(fexp, fpdb):
    """Get format from experimental template."""
    errmsgs = []
    format = {}
    exp = pd.read_csv(fexp)
    
    struc = Structure(fpdb)
    struc.build()
    struc_residues = struc.data_array[:, col_resSeq]
    struc_atomnames = struc.data_array[:, col_name]
    last_residue = int(struc.data_array[:, col_resSeq][-1])

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
    
    for i, res1 in enumerate(format['res1']):
        errmsg = False
        atom1 = format['atom1'][i]
        res2 = format['res2'][i]
        atom2 = format['atom2'][i]
        
        if res1 <= 0 or res2 <= 0:
            errmsg = 'res1/res2 cannot contain 0 or negative values'
        elif res1 > last_residue or res2 > last_residue:
            errmsg = 'res1/res2 cannot be greater than the maximum number of residues in your PDB structure.'  # noqa: E501
        
        try:
            res1_idxs = np.where(struc_residues == str(res1))
            res2_idxs = np.where(struc_residues == str(res2))
            atom1_idxs = np.where(struc_atomnames == atom1)
            atom2_idxs = np.where(struc_atomnames == atom2)
            
            one = np.isin(res1_idxs, atom1_idxs)
            two = np.isin(res2_idxs, atom2_idxs)
            
            if True not in one:
                errmsgs.append(f'{atom1} in PDB does not exist for res1. Skipping {res1} {atom1}')  # noqa: E501
            if True not in two:
                errmsgs.append(f'{atom2} in PDB does not exist for res2. Skipping {res2} {atom2}')  # noqa: E501
            
        except Exception:
            errmsg = 'Residues and/or atoms are not compatible with PDB file(s)'  # noqa: E501
        
        if errmsg is not False:
            raise SPyCiPDBException(errmsg)
    
    errmsgs = [*set(errmsgs)]
    
    return format, errmsgs


def get_exp_format_pre(fexp, fpdb):
    """Get format based on experimental file."""
    errmsgs = []
    format = {}
    exp = pd.read_csv(fexp)
    
    struc = Structure(fpdb)
    struc.build()
    struc_residues = struc.data_array[:, col_resSeq]
    struc_atomnames = struc.data_array[:, col_name]
    last_residue = int(struc.data_array[:, col_resSeq][-1])
    
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
    
    for i, res1 in enumerate(format['res1']):
        errmsg = False
        atom1 = format['atom1'][i]
        res2 = format['res2'][i]
        atom2 = format['atom2'][i]
        
        if res1 <= 0 or res2 <= 0:
            errmsg = 'res1/res2 cannot contain 0 or negative values'
        elif res1 > last_residue or res2 > last_residue:
            errmsg = 'res1/res2 cannot be greater than the maximum number of residues in your PDB structure.'  # noqa: E501
        
        try:
            res1_idxs = np.where(struc_residues == str(res1))
            res2_idxs = np.where(struc_residues == str(res2))
            atom1_idxs = np.where(struc_atomnames == atom1)
            atom2_idxs = np.where(struc_atomnames == atom2)
            
            one = np.isin(res1_idxs, atom1_idxs)
            two = np.isin(res2_idxs, atom2_idxs)
            
            if True not in one:
                errmsgs.append(f'{atom1} in PDB does not exist for res1. Skipping {res1} {atom1}')  # noqa: E501
            if True not in two:
                errmsgs.append(f'{atom2} in PDB does not exist for res2. Skipping {res2} {atom2}')  # noqa: E501
            
        except Exception:
            errmsg = 'Residues and/or atoms are not compatible with PDB file(s)'  # noqa: E501
        
        if errmsg is not False:
            raise SPyCiPDBException(errmsg)
    
    errmsgs = [*set(errmsgs)]
    
    return format, errmsgs


def get_exp_format_smfret(fexp, fpdb):
    """Get format from experimental template."""
    format = {}
    errmsgs = []
    exp = pd.read_csv(fexp)
    
    struc = Structure(fpdb)
    struc.build()
    struc_residues = struc.data_array[:, col_resSeq]
    last_residue = int(struc.data_array[:, col_resSeq][-1])
    
    try:
        format['res1'] = exp.res1.values.astype(int).tolist()
        format['res2'] = exp.res2.values.astype(int).tolist()
        format['scale'] = exp.scale.values.tolist()
    except AttributeError as err:
        errmsg = (
            'Incorrect experimental file format for smFRET subclient. '
            'Text file must have the following columns: '
            'res1,res2,scale'
            )
        raise SPyCiPDBException(errmsg) from err
    
    for i, res1 in enumerate(format['res1']):
        errmsg = False
        res2 = format['res2'][i]
        
        if res1 <= 0 or res2 <= 0:
            errmsg = 'res1/res2 cannot contain 0 or negative values'
        elif res1 > last_residue or res2 > last_residue:
            errmsg = 'res1/res2 cannot be greater than the maximum number of residues in your PDB structure.'  # noqa: E501
        
        try:
            res1_idxs = np.where(struc_residues == str(res1))
            res2_idxs = np.where(struc_residues == str(res2))
            
            if len(res1_idxs) == 0:
                errmsgs.append(f'{res1} in PDB does not exist for res1. Skipping {res1}')  # noqa: E501
            if len(res2_idxs) == 0:
                errmsgs.append(f'{res2} in PDB does not exist for res2. Skipping {res2}')  # noqa: E501
            
        except Exception:
            errmsg = 'Residues and/or atoms are not compatible with PDB file(s)'
        
        if errmsg is not False:
            raise SPyCiPDBException(errmsg)
    
    errmsgs = [*set(errmsgs)]
    
    return format, errmsgs
