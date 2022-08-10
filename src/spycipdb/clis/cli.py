"""
SPyCi-PDB.

Command-line-interface for various back-calculators for
experimental NMR, SAXS, and smFRET data.

Back-calculator logic for PRE, NOE, and JC can be found here:
https://github.com/Oufan75/X-EISD/blob/79df80b3ee76dd37a30d784fed3d2e46ee5cd144/EISD_back_calc.ipynb

USAGE:
    For help:
    >>> spycipdb -h
"""
import argparse
import sys

from spycipdb.libs import libcli
from spycipdb.logger import S
from spycipdb import __version__, log
from spycipdb.clis import(
    cli_pre,
    cli_noe,
    cli_jc,
    cli_smfret,
    )

_prog, _description, _usageage = libcli.parse_doc_params(__doc__)

description = f"""
{_description}

Core back-calculator functions:

    * {cli_pre._name}
    * {cli_noe._name}
    * {cli_jc._name}
    * {cli_smfret._name}
"""

ap = libcli.CustomParser(
    prog='spycipdb',  # _prog,
    description=libcli.detailed.format(description),
    usage=_usageage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version(ap)

subparsers = ap.add_subparsers(
    title='SPyCi-PDB routines',
    help='Short description:',
    )

libcli.add_subparser(subparsers, cli_pre)
libcli.add_subparser(subparsers, cli_noe)
libcli.add_subparser(subparsers, cli_jc)
libcli.add_subparser(subparsers, cli_smfret)

def load_args():
    """Load user input arguments."""
    return ap.parse_args()

def maincli():
    """
    Execute subroutine.

    Arguments are read from user command line input.
    """
    # prints help if not arguments are passed
    # if >2 prints help of subroutines.
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = load_args()

    with open('spycipdb.version', 'w') as fout:
        fout.write(f'version: {__version__}')

    cmd.func(**vars(cmd))
    log.info(S('finished properly'))


if __name__ == '__main__':
    maincli()
