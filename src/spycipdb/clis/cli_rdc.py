"""
Back-calculates RDC values from PDB structure file.

Uses third-party software to do so.
Default is PALES v6.0 for Red Hat Linux.

Publication:
"Prediction of Sterically Induced Alignment in a Dilute
Liquid Crystalline Phase: Aid to Protein Structure Determination by NMR"
Markus Zweckstetter and Ad Bax
Journal of the American Chemical Society 2000 122 (15), 3791-3792
DOI: 10.1021/ja0000908

Back-calculator error 0.88 from Lincoff, J. et al, 2020.
Units for values are in Hz.

USAGE:
    $ spycipdb rdc <PDB-FILES> [--exp-file]
    $ spycipdb rdc --pales <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    To use PALES, experimental data must be formatted to PALES specifications:
    DATA SEQUENCE MTYKLILNGK TLKGETTTEA VDAATAEKVF KQYANDNGVD
    DATA SEQUENCE GEWTYDDATK TFTVTE

    VARS   RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D DD W
    FORMAT %5d     %6s       %6s        %5d     %6s    %6s %9.3f %9.3f %.2f

        2    THR     HN     2    THR      N      1.464      1.000  1.00
        3    TYR     HN     3    TYR      N      7.057      1.000  1.00
        4    LYS     HN     4    LYS      N      8.654      1.000  1.00
        5    LEU     HN     5    LEU      N     12.180      1.000  1.00
        7    LEU     HN     7    LEU      N     12.691      1.000  1.00
        8    ASN     HN     8    ASN      N      5.202      1.000  1.00
       12    LEU     HN    12    LEU      N     11.677      1.000  1.00
       14    GLY     HN    14    GLY      N     10.553      1.000  1.00
       15    GLU     HN    15    GLU      N     11.154      1.000  1.00
    
    For more information, please see the "Data Format" section of:
    https://spin.niddk.nih.gov/bax/software/PALES/

OUTPUT:
    Output is in standard .JSON format as follows:
    {
        'format': {resnum1: [],
                   resname1: [],
                   atomname1: [],
                   resnum2: [],
                   resname2: [],
                   atomname2: [],
                   },
        'pdb1': [rdc_values],
        'pdb2': [rdc_values],
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
from spycipdb.libs import libcli
from spycipdb.libs.libfuncs import get_pdb_paths
from spycipdb.logger import S, T, init_files, report_on_crash
from spycipdb.components.calculators import pales_helper


LOGFILESNAME = '.spycipdb_rdc'
_name = 'rdc'
_help = 'RDC back-calculator using PALES v6.0.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '--pales',
    help='Uses PALES for back-calculator, on by default',
    action='store_true',
    default=True,
    )

libcli.add_argument_pdb_files(ap)
libcli.add_argument_exp_file(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

TMPDIR = '__tmprdc__'
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
        pales=True,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Back-calculate RDCs from PDB and experimental file.

    Parameters
    ----------
    pdb_files : str or Path, required
        Path to a .TAR or folder of PDB files.
        
    exp_file : str or Path, required
        Path to experimental file formatted per PALES requirements.
    
    output : str or Path, optional
        Where to store the back-calculated data.
        Defaults to working directory.
        
    pales : Bool, optional
        Whether to use PALES or another back-calculator.
        Defaults to True.
        
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
    log.info(S('done'))
    
    log.info(T(f'back calculaing using {ncores} workers'))
    if pales:
        execute = partial(
            report_on_crash,
            pales_helper,
            exp_file,
            )
        execute_pool = pool_function(execute, str_pdbpaths, ncores=ncores)
        
        _output = {}
        for result in execute_pool:
            _output['format'] = result[0]
            _output[result[1]] = result[2]
            log.info(S('done'))
    
    # TODO: future PR, do LRDC module
        
    # process output
    log.info(T('Writing output onto disk'))
    with open(output, mode="w") as fout:
        fout.write(json.dumps(_output, indent=4))
    log.info(S('done'))
    
    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
