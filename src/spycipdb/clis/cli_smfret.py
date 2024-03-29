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
import argparse
import json
import shutil
from functools import partial
from pathlib import Path

from idpconfgen.libs.libmulticore import pool_function
from natsort import os_sorted

from spycipdb import log
from spycipdb.core.calculators import calc_smfret
from spycipdb.core.parsers import get_exp_format_smfret
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


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
    if len(pdbs2operate) == 0 or pdbs2operate is None:
        log.info(
            'No .pdb files were found based on the input. Make sure the '
            'folder/tarball contains .pdb files. Only .tar, .tar.xz, .tar.gz '
            'tarballs are accepted.'
            )
        return
    log.info(S('done'))
    pdbs2operate = os_sorted(pdbs2operate)
    log.info(T(f'back calculating using {ncores} workers'))
    execute = partial(
        report_on_crash,
        calc_smfret,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, method='imap', ncores=ncores)  # noqa: E501
    
    _output = {}
    _output['format'], _ = get_exp_format_smfret(exp_file, pdbs2operate[0])
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
