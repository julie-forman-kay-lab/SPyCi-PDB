"""
Back-calculates the J-Couping angles (radians) from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

USAGE:
    $ spycipdb jc <PDB-FILES>
    $ spycipdb jc <PDB-FILES> [--output] [--ncores] [--plot]
    
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

import matplotlib.pyplot as plt
import pandas as pd
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import Structure, col_resSeq
from natsort import os_sorted

from spycipdb import log
from spycipdb.core.calculators import calc_jc
from spycipdb.core.exceptions import SPyCiPDBException
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
libcli.add_argument_plot(ap)

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


def karplus_j(x):
    """
    Convert back-calculated units to hertz using the Karplus equation.
    
    Constants defined from Lincoff et al. 2020
    (https://doi.org/10.1038/s42004-020-0323-0)

    Parameters
    ----------
    x : float
        Back calculated value.

    Returns
    -------
        Value converted to Hertz.
    """
    a = 6.51
    b = -1.76
    c = 1.6
    
    return (a * (x ** 2)) + b * x + c


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
    _output = {}
    exp = pd.read_csv(exp_file)
    try:
        _output['format'] = exp.resnum.values.tolist()
        struc = Structure(pdbs2operate[0])
        struc.build()
        last_residue = int(struc.data_array[:, col_resSeq][-1])
        resnums = exp.resnum.values.astype(int).tolist()
        
        for res in resnums:
            if res <= 0 or res > last_residue:
                errmsg = (
                    'resnum cannot contain 0 or negative values '
                    'and cannot be greater than the maximum number of '
                    'residues in your PDB structure.'
                    )
                raise SPyCiPDBException(errmsg)
        
    except AttributeError as err:
        errmsg = (
            'Incorrect experimental file format for JC subclient. '
            'Text file must have the following columns: '
            'resnum'
            )
        raise SPyCiPDBException(errmsg) from err
    
    log.info(T(f'back calculaing using {ncores} workers'))
    execute = partial(
        report_on_crash,
        calc_jc,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, method='imap', ncores=ncores)  # noqa: E501

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
        
        rawdata = _output
        RESIDUES = rawdata['format']
        del rawdata['format']

        aligned_data = {}
        for res in RESIDUES:
            aligned_data[res] = []

        for conf in rawdata:
            jc = rawdata[conf]
            for i, res in enumerate(RESIDUES):
                j = karplus_j(jc[i])
                aligned_data[res].append(j)

        dataframe = pd.DataFrame(aligned_data)

        fig, ax = plt.subplots()
        fig.set_size_inches(12, 4)
        ax.boxplot(dataframe, flierprops={'marker': 'o', 'markersize': 3})
        ax.set_xticklabels(RESIDUES)

        plt.xlabel('Residue Number', fontsize=14)
        plt.ylabel(r'$^3$J-HNHA Coupling (Hz)', fontsize=14)
        plt.savefig("jc_plot.png", dpi=300, bbox_inches='tight')
        
        log.info(S('done'))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
