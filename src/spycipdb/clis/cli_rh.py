"""
Back-calculates translational Rh values from PDB structure file.

Uses third-party software to do so.
Default is HullRad 8.1.

Publication:
Fleming, P.J. and Fleming, K.G.
"HullRad: Fast Calculations of Folded and Disordered Protein
and Nucleic Acid Hydrodynamic Properties"
Biophysical Journal, 114:856-869, February 27, 2018
DOI: 10.1016/j.bpj.2018.01.002

Back-calculator error 0.1 from Fleming, P.J. and Fleming, K.G. 2018.
Units for values are in Angstroms.

USAGE:
    $ spycipdb rh <PDB-FILES>
    $ spycipdb rh <PDB-FILES> [--output] [--ncores]

OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'pdb1': value,
        'pdb2': value,
        ...
    }
"""
import argparse
import json
import shutil
from functools import partial
from pathlib import Path

from idpconfgen.libs.libmulticore import pool_function

from spycipdb import log
from spycipdb.components.helpers import hullrad_helper
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.spycipdb_rh'
_name = 'rh'
_help = 'Rh back-calculator using HullRad v8.1.'

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

TMPDIR = '__tmprh__'
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
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Use HullRad to predict Rh values.

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
        hullrad_helper,
        )
    execute_pool = pool_function(execute, str_pdbpaths, ncores=ncores)
    
    _output = {}
    for result in execute_pool:
        _output[result[0]] = result[1]
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
