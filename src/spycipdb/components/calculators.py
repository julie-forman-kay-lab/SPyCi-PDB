"""Contains all the different calculators used throughout SPyCi-PDB"""
import os
import sys
import subprocess
import pandas as pd
import numpy as np

from idpconfgen.libs.libhigherlevel import get_torsions
from idpconfgen.libs.libstructure import Structure, col_name, col_resSeq

from spycipdb.libs.libfuncs import get_scalar
from spycipdb.components.hullrad import Sved, model_from_pdb

# Interesting way to import from repository that cannot be
# installed as a module ;-)
# https://www.geeksforgeeks.org/python-import-module-from-different-directory/
current_file_path = os.path.realpath(__file__)
curr_fp_split = current_file_path.split('/')
cspred_fp = ""
for item in curr_fp_split:
    if item == "SPyCi-PDB":
        cspred_fp += item + "/" + "CSpred"
        break
    else:
        cspred_fp += item + "/"
sys.path.insert(0, cspred_fp)

from CSpred import calc_sing_pdb  # noqa: E402


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


# obtaining absolute path of pales executable from current file path
current_file_path = os.path.realpath(__file__)
curr_fp_split = current_file_path.split('/')
PALES_FP = ""
for item in curr_fp_split:
    if item == "SPyCi-PDB":
        PALES_FP += item + "/thirdparty/pales/linux/pales"
        break
    else:
        PALES_FP += item + "/"
        

def pales_helper(exp, pdb_path):
    """
    Handle external PALES shell command.

    Parameters
    ----------
    exp : str
        Absolute path of experimental file formatted per PALES
        standard.
    
    pdb_path : str
        Absolute path of PDB file.

    Returns
    -------
    format : dict
        Format of what atoms from which residues are
        back-calculated
    
    pdb_name_ext : str
        PDB file name with extension.
    
    rdc_bc : list
        List of RDC values. Formatting is given already.
    """
    rdc_bc = []
    format = {
        "resnum1": [],
        "resname1": [],
        "atomname1": [],
        "resnum2": [],
        "resname2": [],
        "atomname2": [],
        }
    pdb_name_ext = pdb_path.rsplit('/', 1)[-1]
    outpath = pdb_name_ext + ".txt"
    
    subprocess.run(
        f"{PALES_FP} -inD {exp} -pdb {pdb_path} -outD {outpath}",
        shell=True,
        capture_output=True,
        )
    
    with open(outpath, 'r') as pales_out:
        for line in pales_out:
            linesplit = line.split()
            try:
                if linesplit[0].isdigit():
                    format['resnum1'].append(int(linesplit[0]))
                    format['resname1'].append(linesplit[1])
                    format['atomname1'].append(linesplit[2])
                    format['resnum2'].append(int(linesplit[3]))
                    format['resname2'].append(linesplit[4])
                    format['atomname2'].append(linesplit[5])
                    rdc_bc.append(float(linesplit[8]))
            except IndexError:
                continue
    
    os.remove(outpath)
    
    return format, pdb_name_ext, rdc_bc


def hullrad_helper(pdb_path):
    """Return translational hydrodynamic radius given PDB."""
    pdb_name_ext = pdb_path.rsplit('/', 1)[-1]
    
    all_atm_rec, num_MG, num_MN, model_array = model_from_pdb(pdb_path)
    
    s, Dt, Dr, vbar_prot, Rht, ffo_hyd_P, M, Ro, Rhr, int_vis, a_b_ratio, \
        Ft, Rg, Dmax, tauC, asphr, AA, NA, GL, DT, useNumpy \
        = Sved(all_atm_rec, num_MG, num_MN, model_array)

    return pdb_name_ext, Rht


def crysol_helper(pdb_path, lm):
    """
    Handle external crysol shell command.

    Parameters
    ----------
    pdb_path : str
        Absolute path of PDB file.
    
    lm : int
        Maximum order of harmonics used for CRYSOL

    Returns
    -------
    pdb_name_ext : str
        PDB file name with extension.
    
    saxs_bc : dict
        Dictionary of index and values for each back-calculation.
    """
    saxs_bc = {}
    index = []
    value = []
    
    wrkdir = os.getcwd()
    pdb_name_ext = pdb_path.rsplit('/', 1)[-1]
    pdb_name = pdb_name_ext[0: pdb_name_ext.index('.')]
    paths = wrkdir + "/" + pdb_name
    
    p = subprocess.Popen(
        f"crysol {pdb_path} --lm={lm} --shell=water",
        stdout=subprocess.PIPE,
        shell=True,
        )
    p.communicate()  # waits for subprocess to stop running
    
    with open(paths + ".abs", mode='r') as crysol_out:
        data = crysol_out.readlines()
        data.pop(0)
        for line in data:
            splitted = line.split()
            index.append(float(splitted[0]))
            value.append(float(splitted[1]))
        
    saxs_bc['index'] = index
    saxs_bc['value'] = value
    
    # removing crysol generated files
    os.remove(paths + ".abs")
    os.remove(paths + ".alm")
    os.remove(paths + ".log")
    os.remove(paths + ".int")
    
    return pdb_name_ext, saxs_bc


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
        scale_factor = ((np.abs(r1 - r2) + 7) / np.abs(r1 - r2)) ** 0.5
        d = d * scale_factor
        eff = 1.0 / (1.0 + (d / scale[i]) ** 6.0)
        fret_bc.append(eff)
    
    fret_bc = np.reshape(fret_bc, (-1, exp.shape[0])).tolist()
    
    return pdb, fret_bc
