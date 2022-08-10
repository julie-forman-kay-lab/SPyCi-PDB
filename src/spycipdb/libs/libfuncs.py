"""Useful functions required throughout."""
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle


def get_scalar(x, y, z):
    return (x**2 + y**2 + z**2)**0.5


def get_pdb_paths(pdb_files, tmpdir):
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    
    return pdbs2operate, _istarfile
