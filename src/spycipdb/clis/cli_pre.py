"""
Back-calculates the PRE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD and DEERPREdict.
Error for the default calculator = 0.0001 as reported in Lincoff et al. 2020.
Error for PREpredict = 0.1 as reported in Tesei and Martins et al. 2021.

Use `-m deerpredict` for PREpredict. Otherwise, leave out the `-m` flag.

USAGE:
    $ spycipdb pre <PDB-FILES> [--exp-file]
    $ spycipdb pre <PDB-FILES> [--exp-file] [--method] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with the following columns:
    
    res1,atom1,res2,atom2
    
    Where res1/atom1 is the residue number and atom name respectively
    for the first residue and res2/atom2 is the residue number and
    atom name respectively for the second residue.

OUTPUT:
    Output is in standard .JSON format as follows, with the first
    key-value pair being the reference formatting for residues and
    atom-names:
    
    {
        'format': {'res1': [], 'atom1': [], 'res2': [], 'atom2': []},
        'pdb1': [dist_values],
        'pdb2': [dist_values],
        ...
    }
    
    or when using DEERPREdict:
    
    {
        'format': {'res1': [], 'atom1': [], 'res2': [], 'atom2': []},
        'pdb1': [intensity_ratios],
        'pdb2': [intensity_ratios],
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
from spycipdb.components.deerpredict import (
    calc_pre_predict,
    default_deerpredict,
    )
from spycipdb.core.calculators import calc_pre
from spycipdb.core.parsers import get_exp_format_pre
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths, plot_data_and_ranges
from spycipdb.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.spycipdb_pre'
_name = 'pre'
_help = 'PRE back-calculator given experimental data template.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)
libcli.add_argument_parameters(ap)
libcli.add_argument_method(ap)
libcli.add_argument_exp_file(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)
libcli.add_argument_plot(ap)

TMPDIR = '__tmppre__'
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
        parameters,
        method="default",
        ncores=1,
        plot=False,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Back-calculate PRE values from PDB structures.

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
    
    parameters : str or Path, optional
        Path to parameters file for methods of back-calculation
        other than the default distance-based calculator.
    
    method : str, optional
        Either uses default back-calculator logic or
        DEERPREdict.
        Defaults to default.
    
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
    log.info(S('done'))
    pdbs2operate = os_sorted(pdbs2operate)
    
    log.info(T(f'back calculating using {ncores} workers'))
    
    if method.lower() == "default":
        execute = partial(
            report_on_crash,
            calc_pre,
            exp_file,
            )
        execute_pool = pool_function(
            execute,
            pdbs2operate,
            method='imap',
            ncores=ncores
            )

        _output = {}
        _output['format'], _ = get_exp_format_pre(exp_file, pdbs2operate[0])
        for result in execute_pool:
            _output[result[0].stem] = result[1]
        log.info(S('done'))
    elif method.lower() == "deerpredict":
        log.info(T('reading parameters'))
        try:
            cust_params = {}
            with open(parameters, 'r') as p_in:
                lines = p_in.readlines()
                for line in lines:
                    p = line.split("=")
                    val = p[1].strip()
                    
                    if p[0] == "temp":
                        val = int(eval(val))
                    elif p[0] == "atom":
                        # Do nothing since it's already a string
                        pass
                    else:
                        val = float(eval(val))
                    cust_params[p[0]] = val
        except FileNotFoundError:
            log.info('Error with parameters file. Using default parameters.')
            cust_params = default_deerpredict
        
        execute = partial(
            report_on_crash,
            calc_pre_predict,
            fexp=exp_file,
            parameters=cust_params,
            )
        execute_pool = pool_function(
            execute,
            pdbs2operate,
            method='imap',
            ncores=ncores
            )

        _output = {}
        _output['format'], _ = get_exp_format_pre(exp_file, pdbs2operate[0])
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
        if method.lower() == "deerpredict":
            log.info(S('plotting is incompatible with DEERPREdict yet...'))
            log.info(S('done'))
        else:
            log.info(T('Plotting back-calculated data'))
            log.info(S('Please be patient if you have a large ensemble...'))
            pre_exp_df = pd.read_csv(exp_file)
            pre_ranges = []
            pre_indices = []
            for row in pre_exp_df.itertuples():
                val = row.dist_value
                upper = row.upper
                lower = row.lower
                ubound = val + upper
                lbound = val - lower
                pre_ranges.append((lbound, ubound))
                pre_indices.append(row.Index)

            pre_vals = []
            _output.pop("format")
            for i in pre_indices:
                temp_vals = []
                for conf in _output:
                    temp_vals.append(_output[conf][i])
                pre_vals.append(temp_vals)

            plot_data_and_ranges(
                pre_vals,
                pre_ranges,
                pre_indices,
                "PRE",
                "pre_plot.png"
                )

            log.info(S('done'))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
