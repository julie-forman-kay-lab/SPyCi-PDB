"""Houses main internal back-calculators for SPyCi-PDB."""
import numpy as np
import pandas as pd
from idpconfgen.libs.libhigherlevel import get_torsions
from idpconfgen.libs.libstructure import Structure, col_name, col_resSeq

from spycipdb.libs.libfuncs import get_scalar


def calc_jc(fexp, pdb):
    """
    Back-calculate JC data.
    
    Requires residues of interest derived from experimental template.
    
    Parameters
    ----------
    fexp : str or Path
        To the experimental file template
    
    pdb : Path
        To the PDB file
    
    Returns
    -------
    jc_bc : list
        Of JC values
    
    pdb : Path
        Of the PDB calculated
    """
    exp = pd.read_csv(fexp)
    # align torsion index as the first residue doesn't have phi torsion
    resn = exp.resnum.values - 2
    jc_bc = np.cos(get_torsions(pdb)[2::3][resn] - np.radians(60))
    
    return pdb, jc_bc.tolist()


def calc_noe(fexp, pdb):
    """
    Back-calculate NOE data.
    
    Atom-pairs and multi-assigns derived from experimental template.
    """
    dist = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    atom1_name = exp.atom1.values
    res2 = exp.res2.values.astype(int)
    atom2_name = exp.atom2.values
    multi1 = exp.atom1_multiple_assignments.values
    multi2 = exp.atom2_multiple_assignments.values
    
    s = Structure(pdb)
    s.build()
    
    for i in range(exp.shape[0]):
        r1 = int(res1[i])
        r2 = int(res2[i])
        atom1_list = []
        atom2_list = []
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r1:
                if atom1_name[i] == 'H':
                    atom1_list.append(s.coords[j, :])
                    break
                if atom1_name[i] in s.data_array[j, col_name]:
                    atom1_list.append(s.coords[j, :])
                if len(atom1_list) == 2:
                    break
                if not multi1[i] and len(atom1_list) == 1:
                    break
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r2:
                if atom2_name[i] == 'H':
                    atom2_list.append(s.coords[j, :])
                    break
                if atom2_name[i] in s.data_array[j, col_name]:
                    atom2_list.append(s.coords[j, :])
                if len(atom2_list) == 2:
                    break
                if not multi2[i] and len(atom2_list) == 1:
                    break

        combos = 0.0
        num_combos = 0

        for first_atom in atom1_list:
            for second_atom in atom2_list:
                dv = first_atom - second_atom
                assert dv.shape == (3,)
                combos += (get_scalar(dv[0], dv[1], dv[2])) ** (-6.)
                num_combos += 1

        dist.append((combos / float(num_combos)) ** (-1 / 6))
    
    return pdb, dist


def calc_pre(fexp, pdb):
    """Back calculates PRE data based on atom-pairs from template."""
    dist = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    atom1_name = exp.atom1.values
    res2 = exp.res2.values.astype(int)
    atom2_name = exp.atom2.values
    
    s = Structure(pdb)
    s.build()
    
    for i in range(exp.shape[0]):
        r1 = int(res1[i])
        r2 = int(res2[i])
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r1:
                if atom1_name[i] == 'H':
                    atom1 = s.coords[j, :]
                    break
                if atom1_name[i] in s.data_array[j, col_name]:
                    atom1 = s.coords[j, :]
                    break
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r2:
                if atom2_name[i] == 'H':
                    atom2 = s.coords[j, :]
                    break
                if atom2_name[i] in s.data_array[j, col_name]:
                    atom2 = s.coords[j, :]
                    break
                
        dv = atom1 - atom2
        assert dv.shape == (3,)
        dist.append(get_scalar(dv[0], dv[1], dv[2]))
    
    return pdb, dist


def calc_smfret(fexp, pdb):
    """
    Back calculate smFRET values.
    
    Take into consideration residue pairs and
    scale from experimental data.
    """
    fret_bc = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    res2 = exp.res2.values.astype(int)
    scale = exp.scale.values
    
    s = Structure(pdb)
    s.build()
    
    for i in range(exp.shape[0]):
        r1 = int(res1[i])
        r2 = int(res2[i])
        r1p = 0
        r2p = 0
        
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r1 and s.data_array[j, col_name] == 'CA':
                r1p = j
                break
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r2 and s.data_array[j, col_name] == 'CA':
                r2p = j
                break
            
        # distance between 2 CA atoms
        dv = s.coords[r1p, :] - s.coords[r2p, :]
        assert dv.shape == (3,)
        d = get_scalar(dv[0], dv[1], dv[2])

        # scale_factor to adjust for dye size and CA to label distances
        # TODO: change r1 and r2 to r1p and r2p respectively
        scale_factor = ((np.abs(r1 - r2) + 7) / np.abs(r1 - r2)) ** 0.5
        d = d * scale_factor
        eff = 1.0 / (1.0 + (d / scale[i]) ** 6.0)
        fret_bc.append(eff)
    
    fret_bc = np.reshape(fret_bc, (-1, exp.shape[0])).tolist()[0]
    return pdb, fret_bc
