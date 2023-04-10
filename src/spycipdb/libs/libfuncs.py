"""Useful functions required throughout."""
import matplotlib.pyplot as plt
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle


def get_scalar(x, y, z):
    """Find scalar quantity using Pythagorean theorem."""
    return (x**2 + y**2 + z**2)**0.5


def get_pdb_paths(pdb_files, tmpdir):
    """Get paths of PDB files from a tarball or folder."""
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    
    return pdbs2operate, _istarfile


def plot_data_and_ranges(
        values,
        ranges,
        indices,
        header,
        output
        ):
    """General plotting function for overlaying bc to exp."""
    fig, ax = plt.subplots()
    positions = list(range(1, len(ranges) + 1))
    for i, r in enumerate(ranges):
        ax.fill_between([i + 0.8, i + 1.2], r[0], r[1], color='gray', alpha=0.5)

    for i, val_group in enumerate(values):
        for val in val_group:
            ax.scatter(
                i % len(values) + 1,
                val,
                color='r',
                zorder=2,
                alpha=0.2,
                s=6
                )

    ax.set_xlabel("Index Number")
    ax.set_ylabel(f"{header} Distance (Å)")

    ax.set_xticks(positions)
    ax.set_xticklabels(indices)
    plt.xticks(fontsize=4)
    
    fig.set_size_inches(15, 5)
    fig.savefig(output, dpi=300)
