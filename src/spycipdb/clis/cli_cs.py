"""
Back-calculates Chemical Shift  distances from PDB structure file.

Uses third-party software to do so. Default is UCBShift:
https://github.com/THGLab/CSpred

Back-calculator errors as reported for H, HA, C, CA, CB, N:
0.38, 0.22, 1.31, 0.97, 1.29 and 2.16 ppm respectively.

USAGE:
    $ spycipdb cs <PDB-FILES> [--ph]
    $ spycipdb cs <PDB-FILES> [--ph] [--output] [--ncores]

REQUIREMENTS:
    Installation of UCBShift, please refer to documentation.

OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'pdb1': {
                    'H': [],
                    'HA': [],
                    'C': [],
                    'CA': [],
                    'CB': [],
                    'N': [],
                },
        'pdb2': {
                    'H': [],
                    'HA': [],
                    'C': [],
                    'CA': [],
                    'CB': [],
                    'N': [],
                },
        ...
    }
"""
import pandas as pd
import json
import argparse
import shutil
from pathlib import Path
from functools import partial

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash

from idpconfgen.libs.libmulticore import pool_function

from CSpred.CSpred import calc_sing_pdb

LOGFILESNAME = '.spycipdb_cs'
_name = 'cs'
_help = 'CS back-calculator using UCBShift.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '--ph',
    help='pH value to be considered. Defaults to 5',
    type=int,
    default=5,
    )

libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

TMPDIR = '__tmpcs__'
ap.add_argument(
    '--tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def main(
        pdb_files,
        output,
        ph=5,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Main logic for using UCBShift to predict chemical shift
    values for PDB structures and output.

    Parameters
    ----------
        pdb_files (_type_): _description_
        output (_type_): _description_
        ph (int, optional): _description_. Defaults to 5.
        ncores (int, optional): _description_. Defaults to 1.
        tmpdir (_type_, optional): _description_. Defaults to TMPDIR.
    """
    init_files(log, LOGFILESNAME)
    
    if ph < 2 or ph > 12:
        log.info(S('WARNING: Predictions for proteins in extreme pH '
                   'conditions are likely to be erroneous. '
                   'Take prediction results at your own risk!'
                   ))
    
    log.info(T('reading input paths'))
    pdbs2operate, _istarfile = get_pdb_paths(pdb_files, tmpdir)
    log.info(S('done'))
    
    log.info(T(f'back calculaing using {ncores} workers'))

    execute = partial(
        report_on_crash,
        calc_sing_pdb,
        pH=ph,
        )
    execute_pool = pool_function(execute, pdb_file_name=pdbs2operate, ncores=ncores)
    
    _output = {}
    for results in execute_pool:
        per_struct = {'H': [], 'HA': [], 'C': [], 'CA': [], 'CB': [], 'N': []}
        preds = results[1]
        
        _output[results[0].stem] = per_struct
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
