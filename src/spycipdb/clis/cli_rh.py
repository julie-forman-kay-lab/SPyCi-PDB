"""
Back-calculates translational Rh values from PDB structure file.

Uses third-party software to do so.
Default is HullRadSAS 3.1.

Publication:
Fleming, P.J., Corriea, J.J. and Fleming, K.G.
"Revisiting Macromolecular Hydration with HullRadSAS"
https://www.biorxiv.org/content/10.1101/2022.10.20.513022v1

USAGE:
    $ spycipdb rh <PDB-FILES>
    $ spycipdb rh <PDB-FILES> [--output] [--ncores] [--plot]

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

import matplotlib.pyplot as plt
import pandas as pd
from idpconfgen.libs.libmulticore import pool_function
from natsort import os_sorted

from spycipdb import log
from spycipdb.components.helpers import hullrad_helper
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.spycipdb_rh'
_name = 'rh'
_help = 'Rh back-calculator using HullRadSAS v3.1.'

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
libcli.add_argument_plot(ap)

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
        plot=False,
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
    str_pdbpaths = [str(path) for path in pdbs2operate]
    str_pdbpaths = os_sorted(str_pdbpaths)
    if len(pdbs2operate) == 0 or pdbs2operate is None:
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
    execute_pool = pool_function(execute, str_pdbpaths, method='imap', ncores=ncores)  # noqa: E501
    
    _output = {}
    rh = []
    for result in execute_pool:
        _output[result[0]] = result[1]
        rh.append(result[1])
    log.info(S('done'))
    
    log.info(T('Writing output onto disk'))
    with open(output, mode="w") as fout:
        fout.write(json.dumps(_output, indent=4))
    log.info(S('done'))

    if _istarfile:
        shutil.rmtree(tmpdir)
    
    if plot:
        log.info(T('Plotting back-calculated data'))
        
        dataframe = pd.DataFrame(rh)
        fig, ax = plt.subplots()
        fig.set_size_inches(2, 4)
        ax.boxplot(dataframe, flierprops={'marker': 'o', 'markersize': 5})

        ax.xaxis.set_ticklabels([])

        plt.xlabel('Protein', fontsize=14)
        plt.ylabel('Rh (Å)', fontsize=14)
        plt.savefig("rh_plot.png", dpi=300, bbox_inches='tight')
        
        log.info(S('done'))
        
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
