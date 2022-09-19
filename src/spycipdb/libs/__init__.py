"""General Libraries for the project."""
import os
from pathlib import Path as _Path


class Path(type(_Path())):
    """
    Extends Python's `Path object`_ interface.
    
    Derived from:
    https://github.com/joaomcteixeira/taurenmd/blob/4e087f0bdb7ea7b7b5ae70c589c73a72437f3de6/src/taurenmd/core.py
    
    .. _Path object: https://docs.python.org/3/library/pathlib.html
    """
    
    def str(self):
        """
        Represent path as string.
        
        Alias to ``os.fspath(self)``.
        
        Returns
        -------
        str
           ``os.fspath(self)``.
        """
        return os.fspath(self)
    
    def myparents(self):
        """
        List of the path parent folders.
        
        Alias to ``pathlib.Path.resolve().parents[0]``.
        
        Returns
        -------
        list
            Parent paths. Name file or folder are excluded.
        """
        return self.resolve().parents[0]
