"""
Back-calculates the NOE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD.
Error = 0.0001 as reported in Lincoff et al. 2020.

USAGE:
    $ spycipdb noe <PDB-FILES> [--exp-file]
    $ spycipdb noe <PDB-FILES> [--exp-file] [--output] [--ncores] [--plot]

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

import pandas as pd
from idpconfgen.libs.libmulticore import pool_function
from natsort import os_sorted

from spycipdb import log
from spycipdb.core.calculators import calc_noe
from spycipdb.core.parsers import get_exp_format_noe
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths, plot_data_and_ranges
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
libcli.add_argument_plot(ap)

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
        plot=False,
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
    
    plot : Bool, optional
        Whether to plot the back-calculated results or not.
        Defaults to False.
    
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
        calc_noe,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, method='imap', ncores=ncores)  # noqa: E501
    
    _output = {}
    _output['format'], _ = get_exp_format_noe(exp_file, pdbs2operate[0])
    
    for result in execute_pool:
        _output[result[0].stem] = result[1]
    log.info(S('done'))
    
    log.info(T('Writing output onto disk'))
    with open(output, mode="w") as fout:
        fout.write(json.dumps(_output, indent=4))
    log.info(S('done'))
    
    if _istarfile:
        shutil.rmtree(tmpdir)
    
    if plot:
        log.info(T('Plotting back-calculated data'))
        log.info(S('Please be patient if you have a large ensemble...'))
        noe_exp_df = pd.read_csv(exp_file)
        noe_ranges = []
        noe_indices = []
        for row in noe_exp_df.itertuples():
            val = row.dist_value
            upper = row.upper
            lower = row.lower
            ubound = val + upper
            lbound = val - lower
            noe_ranges.append((lbound, ubound))
            noe_indices.append(row.Index)
        
        noe_vals = []
        _output.pop("format")
        for i in noe_indices:
            temp_vals = []
            for conf in _output:
                temp_vals.append(_output[conf][i])
            noe_vals.append(temp_vals)

        plot_data_and_ranges(
            noe_vals,
            noe_ranges,
            noe_indices,
            "NOE",
            "noe_plot.png"
            )
        
        log.info(S('done'))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
