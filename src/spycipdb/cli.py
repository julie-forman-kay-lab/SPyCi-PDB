"""
SPyCi-PDB.

Command-line-interface for various back-calculators for
    experimental NMR, SAXS, and smFRET data.

USAGE:
    For help:
    >>> spycipdb -h
"""
import argparse
import sys

from spycipdb.libs import libcli
from spycipdb.logger import S
from spycipdb import(
    __version__,
    log,
    # Rest of the CLI modules go here
    )

_prog, _description, _usageage = libcli.parse_doc_params(__doc__)

description = f"""
{_description}

    * Name goes here

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

    with open('idpconfgen.version', 'w') as fout:
        fout.write(f'version: {__version__}')

    cmd.func(**vars(cmd))
    log.info(S('finished properly'))


if __name__ == '__main__':
    maincli()
