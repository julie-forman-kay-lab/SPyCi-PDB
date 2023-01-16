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

Back-calculator error ? (0.006 from Lincoff et al. 2020).
Units for values are: I_abs(s)[cm^-1]/c[mg/ml]

USAGE:
    $ spycipdb saxs <PDB-FILES> [--lm]
    $ spycipdb saxs <PDB-FILES> [--lm] [--output] [--ncores]

REQUIREMENTS:
    Installation of CRYSOL 3.0 from ATSAS 3.1.1, please refer to documentation.

OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'format': [index],
        'pdb1': [value],
        'pdb2': [value],
        ...
    }
"""
import argparse
import json
import shutil
import subprocess
from functools import partial
from pathlib import Path

from idpconfgen.libs.libmulticore import pool_function

from spycipdb import log
from spycipdb.components.helpers import crysol_helper
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.spycipdb_saxs'
_name = 'saxs'
_help = 'SAXS back-calculator using CRYSOL v3.0.'

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
    Use CRYSOLv3 to predict SAXS values from PDBs.

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
    
    # check to see if crysol is installed
    try:
        subprocess.check_output(["crysol", "-h"])
    except FileNotFoundError:
        log.info(S('WARNING: ATSAS v3.1.1 is not installed. Exiting...'))
        return
    
    if lm < 1 or lm > 100:
        log.info(S('WARNING: maximum order of harmonics is not within the range of 1-100 inclusive.'))  # noqa: E501
        
        return
    
    log.info(T('reading input paths'))
    pdbs2operate, _istarfile = get_pdb_paths(pdb_files, tmpdir)
    str_pdbpaths = [str(path) for path in pdbs2operate]
    if len(pdbs2operate) == 0:
        log.info(
            'No .pdb files were found based on the input. Make sure the '
            'folder/tarball contains .pdb files. Only .tar, .tar.xz, .tar.gz '
            'tarballs are accepted.'
            )
        return
    log.info(S('done'))
    
    log.info(T(f'back calculaing using {ncores} workers'))
    execute = partial(
        report_on_crash,
        crysol_helper,
        lm=lm,
        )
    execute_pool = pool_function(execute, str_pdbpaths, ncores=ncores)
    
    _output = {}
    for result in execute_pool:
        _output['format'] = result[1]['index']
        _output[result[0]] = result[1]['value']
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
