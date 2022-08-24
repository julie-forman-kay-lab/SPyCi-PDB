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
    Experimental data must be comma-delimited with the following columns:
    
    res1,res2,scaler
    
    Where res1/res2 is the residue number for the
    first and second residue respectively.
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
        'pdb1': [smfret_values],
        'pdb2': [smfret_values],
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
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resSeq,
    )

LOGFILESNAME = '.spycipdb_smfret'
_name = 'smfret'
_help = 'smFRET back-calculator given experimental data template.'

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


def get_exp_format_smfret(fexp):
    """Get format from experimental template."""
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['scale'] = exp.scale.values.tolist()
    
    return format


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


def main(
        pdb_files,
        exp_file,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Back calculate smFRET values from PDB structures.
    
    Requires experimental file template.

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
    
    log.info(T(f'back calculating using {ncores} workers'))
    execute = partial(
        report_on_crash,
        calc_smfret,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    
    _output = {}
    _output['format'] = get_exp_format_smfret(exp_file)
    for results in execute_pool:
        _output[results[0].stem] = results[1]
    log.info(S('done'))
    
    log.info(T('Writing output onto disk'))
    with open(output, mode="w") as fout:
        fout.write(json.dumps(_output, indent=4))
    log.info(S('done'))
        
    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
