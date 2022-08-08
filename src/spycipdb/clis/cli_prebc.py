"""
Back-calculates the PRE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from Oufan Zhang @Oufan75

USAGE:
    $ spycipdb prebc <PDB-FILES> [--exp-file]
    $ spycipdb prebc <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with at least the following columns:
    
    res1,atom1,res2,atom2
    
    Where res1/atom1 is the atom number and name respectively for the first residue
    and res2/atom2 is the atom number and name respectively for the second residue.
"""
import argparse
import shutil
import numpy as np
import pandas as pd

from pathlib import Path

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.logger import S, T, init_files, report_on_crash
from spycipdb.libs.libfuncs import get_scalar

from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import Structure, col_name

LOGFILESNAME = '.spycipdb_prebc'
_name = 'prebc'
_help = 'PRE back-calculator given experimental data template.'

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

TMPDIR = '__tmpprebc__'
ap.add_argument(
    '--tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def calc_pre(fexp, pdb):
    dist = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(np.int)
    atom1_name = exp.atom1.values
    res2 = exp.res2.values.astype(np.int)
    atom2_name = exp.atom2.values
    
    s = Structure(pdb)
    s.build()
    # TODO: see how idpconfgen struct processes atom names and coords
    
    for idx in range (exp.shape[0]):
        r1 = np.int(res1[idx])
        r2 = np.int(res2[idx])
        
    
    return dist


def main(
        pdb_files,
        exp_file,
        output="prebc.json",
        ncores=1,
        tmpdir=TMPDIR,
        ):
    """
    Main logic for back-calculating PRE values from PDB structures
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
        Defaults to "prebc.json".
        
    ncores : int, optional
        The number of cores to use.
        Defaults to 1.
    
    tmpdir : str or Path, optional
        Path to the temporary directory if working with .TAR files.
        Defaults to TMPDIR.
    """
    init_files(log, LOGFILESNAME)
    
    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))


    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
