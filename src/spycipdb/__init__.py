"""
SPyCi-PDB.

Command-line-interface for various back-calculators for
experimental NMR, SAXS, and smFRET data.

Logic imported from IDPConformerGenerator:
https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/3aef6b085ec09eeebc5812639a5eb6832c0215cd/src/idpconfgen/__init__.py
"""
import logging
import string
from os import fspath, get_terminal_size
from pathlib import Path as _Path

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

try:
    get_terminal_size()
except OSError:
    has_terminal = False
    log.addHandler(logging.NullHandler())
else:
    _ch = logging.StreamHandler()
    _ch.setLevel(logging.INFO)
    _ch.setFormatter(logging.Formatter('[%(asctime)s]%(message)s'))
    log.addHandler(_ch)
    has_terminal = True


class Path(type(_Path())):
    """
    Path object dedicated to this software.

    Inherits from pathlib.Path.

    This creates an interface so that if new methods are required
    the Path interface does not need to be refactored across.
    """

    def str(self):
        """
        Return string version of Path.

        Avoids using os.fspath around libs.
        """
        return fspath(self)

    def myparents(self):
        """Return the Path to the parent folder resolved to absolute."""
        return self.resolve().parent

    @property
    def absparent(self):
        """Return the Path to the parent folder resolved to absolute."""
        return self.resolve().parent


def count_string_formatters(s):
    """
    Count string formatters: ``{}``.

    Returns
    -------
    int
        The number of string formatters.
    """
    assert isinstance(s, str), f'`s` of wrong type: {type(s)}'
    return sum(1 for f in list(string.Formatter().parse(s)) if f[1] is not None)


source_folder = Path(__file__).absparent

__version__ = '0.3.5'
