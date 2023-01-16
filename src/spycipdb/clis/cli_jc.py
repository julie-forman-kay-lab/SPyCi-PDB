"""
Back-calculates the J-Couping angles (radians) from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

USAGE:
    $ spycipdb jc <PDB-FILES>
    $ spycipdb jc <PDB-FILES> [--output] [--ncores]
    
REQUIREMENTS:
    Experimental data must be comma-delimited with the following column:
    
    resnum
    
    Where `resnum` indicates the JC for a specific residue.

OUTPUT:
    Output is in standard .JSON format as follows, with `jc_values`
    being for each residue aligned with `resnum` format.
    
    {
        'format': [resnum],
        'pdb1': [jc_values],
        'pdb2': [jc_values],
        ...
    }
"""
import argparse
import json
import shutil
from functools import partial
from pathlib import Path

import pandas as pd
from idpconfgen.libs.libmulticore import pool_function

from spycipdb import log
from spycipdb.core.calculators import calc_jc
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.spycipdb_jc'
_name = 'jc'
_help = 'J-Coupling back-calculator given optional experimental data template.'

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

TMPDIR = '__tmpjc__'
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
        exp_file,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Process PDB structures and return back-calculated JC values.
    
    Parameters
    ----------
    pdb_files : str or Path, required
        Path to a .TAR or folder of PDB files.
        
    exp_file : str or Path, required
        Path to experimental file template.
        Required to know for which residues to calculate.
    
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
        calc_jc,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    
    _output = {}
    exp = pd.read_csv(exp_file)
    _output['format'] = exp.resnum.values.tolist()
    for result in execute_pool:
        _output[result[0].stem] = result[1]
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
