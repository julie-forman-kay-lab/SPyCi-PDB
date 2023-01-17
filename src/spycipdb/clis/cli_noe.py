"""
Back-calculates the NOE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD.
Error = 0.0001 as reported in Lincoff et al. 2020.

USAGE:
    $ spycipdb noe <PDB-FILES> [--exp-file]
    $ spycipdb noe <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with the following columns:
    
    res1,atom1,atom1_multiple_assignments,res2,atom2,atom2_multiple_assignments
    
    Where res1/atom1 is the residue number and atom name respectively for
    the first residue and res2/atom2 is the residue number and atom name
    respectively for the second residue.
    Multiple assignments are either 0 or 1 (no/yes respectively).

OUTPUT:
    Output is in standard .JSON format as follows, with the first
    key-value pair being the reference formatting for residues and
    atom-names:
    
    {
        'format': { 'res1': [],
                    'atom1': [],
                    'atom1_multiple_assignments': [],
                    'res2': [],
                    'atom2': [],
                    'atom2_multiple_assignments': []
                    },
        'pdb1': [dist_values],
        'pdb2': [dist_values],
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
from spycipdb.core.calculators import calc_noe
from spycipdb.core.parsers import get_exp_format_noe
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


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

TMPDIR = '__tmpnoe__'
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
    Process PDB structures and output back-calculated NOE values.
    
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
    
    log.info(T(f'back calculating using {ncores} workers'))
    execute = partial(
        report_on_crash,
        calc_noe,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    
    _output = {}
    _output['format'] = get_exp_format_noe(exp_file)
    if get_exp_format_noe(exp_file) is False:
        log.info(
            'Incorrect experimental file format for NOE subclient. '
            'Text file must have the following columns: '
            )
        log.info('res1,atom1,atom1_multiple_assignments,res2,atom2,atom2_multiple_assignments')  # noqa: E501
        return
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
