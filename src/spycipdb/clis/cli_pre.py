"""
Back-calculates the PRE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD.
Error = 0.0001 as reported in Lincoff et al. 2020.

USAGE:
    $ spycipdb pre <PDB-FILES> [--exp-file]
    $ spycipdb pre <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with at least the following columns:
    
    res1,atom1,res2,atom2
    
    Where res1/atom1 is the residue number and atom name respectively for the first residue
    and res2/atom2 is the residue number and atom name respectively for the second residue.

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
"""
import json
import argparse
import shutil
import pandas as pd
from pathlib import Path
from functools import partial

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.logger import S, T, init_files, report_on_crash
from spycipdb.libs.libfuncs import get_scalar, get_pdb_paths

from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import(
    Structure,
    col_name,
    col_resSeq,
    )

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
libcli.add_argument_exp_file(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

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


def get_exp_format_pre(fexp):
    format = {}
    exp = pd.read_csv(fexp)
    
    format['res1'] = exp.res1.values.astype(int).tolist()
    format['atom1'] = exp.atom1.values.tolist()
    format['res2'] = exp.res2.values.astype(int).tolist()
    format['atom2'] = exp.atom2.values.tolist()
    
    return format


def calc_pre(fexp, pdb):
    """
    Main logic for back-calculating PRE data
    with atom-pairs derived from experimental template.
    """
    dist = []
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    atom1_name = exp.atom1.values
    res2 = exp.res2.values.astype(int)
    atom2_name = exp.atom2.values
    
    s = Structure(pdb)
    s.build()
    
    for i in range (exp.shape[0]):
        r1 = int(res1[i])
        r2 = int(res2[i])
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r1:
                if atom1_name[i] == 'H':
                    atom1 = s.coords[j, :]
                    break
                if atom1_name[i] in s.data_array[j, col_name]:
                    atom1 = s.coords[j, :]
                    break
        for j, r in enumerate(s.data_array[:, col_resSeq].astype(int)):
            if r == r2:
                if atom2_name[i] == 'H':
                    atom2 = s.coords[j, :]
                    break
                if atom2_name[i] in s.data_array[j, col_name]:
                    atom2 = s.coords[j, :]
                    break
                
        dv = atom1 - atom2
        assert dv.shape == (3,)
        dist.append(get_scalar(dv[0], dv[1], dv[2]))
    
    return pdb, dist


def main(
        pdb_files,
        exp_file,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    """
    Main logic for back-calculating PRE values from PDB structures
    given experimental file template.

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
    log.info(S('done'))
    
    log.info(T(f'back calculating using {ncores} workers'))
    execute = partial(
        report_on_crash,
        calc_pre,
        exp_file,
        )
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    
    _output = {}
    _output['format'] = get_exp_format_pre(exp_file)
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
