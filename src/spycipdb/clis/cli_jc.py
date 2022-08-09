"""
Back-calculates the PRE distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

USAGE:
    $ spycipdb jcbc <PDB-FILES>
    $ spycipdb jcbc <PDB-FILES> [--output] [--ncores]
    
REQUIREMENTS:
    Experimental data must be comma-delimited with at least the following columns:
    
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
import numpy as np
import json
import argparse
import shutil
from pathlib import Path
from functools import partial

from spycipdb import log
from spycipdb.libs import libcli
from spycipdb.logger import S, T, init_files, report_on_crash

from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libcalc import calc_torsion_angles

LOGFILESNAME = '.spycipdb_jc'
_name = 'jc'
_help = 'J-Coupling back-calculator given optional experimental data template.'

