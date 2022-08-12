"""
Back-calculates SAXS distances from PDB structure file.

Uses third-party software to do so. Default is CRYSOL 3.0 from ATSAS 3.1.1:
https://www.embl-hamburg.de/biosaxs/crysol.html

Publication:
Franke, D., Petoukhov, M.V., Konarev, P.V., Panjkovich, A., Tuukkanen, A.,
Mertens, H.D.T., Kikhney, A.G., Hajizadeh, N.R., Franklin, J.M., Jeffries,
C.M. and Svergun, D.I. (2017) ATSAS 2.8: a comprehensive data analysis suite
for small-angle scattering from macromolecular solutions. J. Appl. Cryst.
50(4), 1212-1225.

Back-calculator error ? (0.006 from Lincoff et al. 2020)

USAGE:
    $ spycipdb saxs <PDB-FILES> [--lm]
    $ spycipdb saxs <PDB-FILES> [--lm] [--output] [--ncores]

REQUIREMENTS:
    Installation of CRYSOL 3.0 from ATSAS 3.1.1, please refer to documentation.

OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'pdb1': {
                },
        'pdb2': {
                },
        ...
    }
"""
import sys
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

from idpconfgen.libs.libmulticore import pool_function

LOGFILESNAME = '.spycipdb_saxs'
_name = 'saxs'
_help = 'SAXS back-calculator using CRYSOLv3.0.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '--lm',
    help=('Maximum number of harmonics (1-100). '
          'For large/extended structures, higher = better results '
          'at a cost of computational time.'
          'Defaults to 20.'
          ),
    type=int,
    default=20,
    )

libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

TMPDIR = '__tmpsaxs__'
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
        lm=20,
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
        
    lm : int, optional
        Maximum order of harmonics for CRYSOL back-calc.
        Defaults to 20.
    
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
    
    if lm < 1 or lm > 100:
        log.info(S('WARNING: maximum order of harmonics '
                   'is not within the range of 1-100 inclusive. '
                   ))
        return
    
    log.info(T('reading input paths'))
    pdbs2operate, _istarfile = get_pdb_paths(pdb_files, tmpdir)
    log.info(S('done'))




    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
