"""
Back-calculates scaled smFRET distances from PDB structure file.

Uses idpconfgen libraries for coordinate parsing as it's proven
to be faster than BioPython.

Back-calculator logic inspired from X-EISD.
Error = 0.0074 as reported in Lincoff et al. 2020.

USAGE:
    $ spycipdb smfret <PDB-FILES> [--exp-file]
    $ spycipdb smfret <PDB-FILES> [--exp-file] [--output] [--ncores]

REQUIREMENTS:
    Experimental data must be comma-delimited with at least the following columns:
    
    res1,res2,scaler
    
    Where res1/res2 is the residue number for the first and second residue respectively.
    Scaler is the r0 Foster radius of the dye pair.

OUTPUT:
    Output is in standard .JSON format as follows, with the first
    key-value pair being the reference formatting for residues and
    scaler values:
    
    {
        'format': { 'res1': [],
                    'res2': [],
                    'scaler': [],
                    },
        'pdb1': [values],
        'pdb2': [values],
        ...
    }
"""