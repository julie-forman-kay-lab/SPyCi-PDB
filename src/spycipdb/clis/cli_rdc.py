"""
Back-calculates RDC values from PDB structure file.

Uses third-party software to do so.
Default is PALES v6.0 for Red Hat Linux.

Publication:
"Prediction of Sterically Induced Alignment in a Dilute Liquid Crystalline Phase:
Aid to Protein Structure Determination by NMR"
Markus Zweckstetter and Ad Bax
Journal of the American Chemical Society 2000 122 (15), 3791-3792
DOI: 10.1021/ja0000908

Back-calculator error 0.88 from Lincoff, J. et al, 2020.
Units for values are in Hz.

USAGE:
    $ spycipdb rdc <PDB-FILES>
    $ spycipdb rdc <PDB-FILES> [--output] [--ncores]

TODO: figure out how the output is supposed to look like
OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'pdb1': value,
        'pdb2': value,
        ...
    }
"""
from multiprocessing import current_process
import os
import json
import argparse
import shutil
from pathlib import Path
from functools import partial

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash
from spycipdb.components.hullrad import model_from_pdb, Sved

from idpconfgen.libs.libmulticore import pool_function

LOGFILESNAME = '.spycipdb_rdc'
_name = 'rdc'
_help = 'RDC back-calculator using PALES v6.0.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

TMPDIR = '__tmprdc__'
ap.add_argument(
    '--tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )

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
        

def pales_helper(pdb_path):
    # TODO: incomplete, need to figure out output format
    rdc_bc = []
    
    return pdb_path, rdc_bc



def main(
        pdb_files,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Main logic for using UCBShift to predict chemical shift
    values for PDB structures and output.

    Parameters
    ----------
    pdb_files : str or Path, required
        Path to a .TAR or folder of PDB files.
    
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
    str_pdbpaths = [str(path) for path in pdbs2operate]
    log.info(S('done'))
    
    
    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
