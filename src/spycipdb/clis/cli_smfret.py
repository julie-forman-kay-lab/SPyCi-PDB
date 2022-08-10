"""
Back-calculates scaled smFRET distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD.
Error = 0.0074 as reported in Lincoff et al. 2020.

USAGE:
    $ spycipdb smfret <PDB-FILES> [--exp-file]
    $ spycipdb smfret <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with at least the following columns:
    
    res1,res2,scaler
    
    Where res1/res2 is the residue number for the first and second residue respectively.
    Scaler is the r0 Foster radius of the dye pair.

OUTPUT:
    Output is in standard .JSON format as follows, with the first
    key-value pair being the reference formatting for residues and
    scaler values:
    
    {
        'format': { 'res1': [],
                    'res2': [],
                    'scale': [],
                    },
        'pdb1': [values],
        'pdb2': [values],
        ...
    }
    
TODO: currently assumes 'CA' as atom labeled
TODO: provide alternative strategies of back-calculating and
    interpreting smFRET data.
"""
import json
import argparse
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from functools import partial

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.logger import S, T, init_files, report_on_crash
from spycipdb.libs.libfuncs import get_pdb_paths, get_scalar

from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import(
    Structure,
    col_name,
    col_resSeq,
    )

LOGFILESNAME = '.spycipdb_noe'
_name = 'noe'
_help = 'NOE back-calculator given experimental data template.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)
libcli.add_argument_exp_file(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

TMPDIR = '__tmpsmfret__'
ap.add_argument(
    '--tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def get_exp_format_pre(fexp):
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['scale'] = exp.scale.values.tolist()
    
    return format


def calc_smfret(fexp, pdb):
    """
    Main logic for back-calculating smFRET values
    taking into consideration residue pairs and
    scale from experimental data.
    """
    fret_bc = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    res2 = exp.res2.values.astype(int)
    scale = exp.scale.values
    
    s = Structure(pdb)
    s.build()
    # assumes CA as atom labeled
    s.add_filter(lambda x: x[col_name] == 'CA')
    
    for j in range(exp.shape[0]):
        r1 = int(res1[j])
        r2 = int(res2[j])
        # distance between 2 CA atoms
        dv = s.coords[r1, :] - s.coords[r2, :]
        assert dv.shape == (3,)
        d = get_scalar(dv[0], dv[1], dv[2])
        
        # scale_factor to adjust for dye size and CA to label distances
        scale_factor = ((np.abs(r1 - r2) + 7) / np.abs(r1 - r2)) ** 0.5
        d = d * scale_factor
        eff = 1.0 / (1.0 + (d / scale[j]) ** 6.0)

        fret_bc.append(eff)
    
    fret_bc = np.reshape(fret_bc, (-1, exp.shape[0]))
    
    return pdb, fret_bc


def main(
        pdb_files,
        exp_file,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Main logic for back-calculating smFRET values from PDB structures
    given experimental file template.

    Parameters
    ----------
    pdb_files : str or Path, required
        Path to a .TAR or folder of PDB files.
        
    exp_file : str or Path, required
        Path to experimental file template.
        Required to know which distances to calculate.
    
    output : str or Path, optional
        Where to store the back-calculated data.
        Defaults to working directory.
        
    ncores : int, optional
        The number of cores to use.
        Defaults to 1.
    
    tmpdir : str or Path, optional
        Path to the temporary directory if working with .TAR files.
        Defaults to TMPDIR.
    """
    init_files(log, LOGFILESNAME)
    
    log.info(T('reading input paths'))
    pdbs2operate, _istarfile = get_pdb_paths(pdb_files, tmpdir)
    log.info(S('done'))
    
        
    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
