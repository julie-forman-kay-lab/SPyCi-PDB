Usage
=====

SPyCi-PDB runs entirely through command-lines. Theoretically, it is compatible
with any machine that can un Python. However, it's only been thoroughly tested
on WSL 2.0 on Windows 11, Linux Ubuntu 20.04.X LTS and 18.04.X LTS.

Please follow the explanations in this page plus the documentation on the command-line
themselves via::

    spycipdb <MODULE> -h


Please note that SPyCi-PDB can also be used as a library as back-calculator functions
are completely modular::

    import spycipdb

Command-Lines
-------------

To execute :code:`spycipdb` command-line, run :code:`spycipdb` in your
terminal window, after :ref:`installation <Installation>`::

    spycipdb

or::

    spycipdb -h

Both will output the help menu.

.. note::
    All subclients have the :code:`-h` option to show help information.

Formatting for Input and Output Files
-------------------------------------

Conformers will be required in the PDB format with the ``.pdb`` file extension.
For tarballs and folders, only the ``.pdb`` files in the folder/tarball will
be used. Accepted tarballs include ``.tar``, ``.tar.gz``, and ``.tar.xz`` file
extensions.

Most modules will require a sample experimental results template to base
the back-calculations off of. Please note that experimental result ``values``
are not required. For the mathematical equations for the default internal calculators
(``pre``, ``noe``, ``jc``, ``smfret``) please refer to the paper.

All input files can be in the ``.txt`` file format with comma-delimitation
for values. Examples can be found in the ``example/drksh3_exp_data/`` folder.
The required header formatting for the ``.txt`` file for each experimental
module is highlighted below.

Output files are saved in a standard human-readable ``.JSON`` format.
In most cases, the first key-value pair gives the format for each of the values
in subsequent key-value pairs.

Paramagnetic Resonance Entropy (PRE) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Input:** ``res1,atom1,res2,atom2``
| Where res1/atom1 is the residue number and atom name respectively for
  the first residue and res2/atom2 is the residue number and atom name
  respectively for the second residue.
| **Output:** ``'format': {'res1': [], 'atom1': [], 'res2': [], 'atom2': []}``
| Subsequent keys are the names of the PDB file with a value of: ``[dist_values]``.
| **About:** The default back-calculator is a simple distance based interpretation of PRE data.
  Scalar distances between pairs of atoms according to the pairs derived from the experimental template are used.

Nuclear Overhauser Effect (NOE) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Input:** ``res1,atom1,atom1_multiple_assignments,res2,atom2,atom2_multiple_assignments``
| Where res1/atom1 is the residue number and atom name respectively for the first residue and res2/atom2 is 
  the residue number and atom name respectively for the second residue. 
  Multiple assignments are either 0 or 1 (no/yes respectively).
| **Output:** ``'format': {'res1': [], 'atom1': [], 'atom1_multiple_assignments': [], 'res2': [], 'atom2': [], 'atom2_multiple_assignments': []}``
| Subsequent keys are the names of the PDB file with a value of: ``[dist_values]``.
| **About:** The default back-calculator is a simple distance based interpretation of NOE data.
  Scalar distances between pairs of atoms according to the pairs derived from the experimental template are used.

3J-HNHA coupling (JC) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Input:** ``resnum``
| Where ``resnum`` indicates the JC for a specific residue.
| **Output:** ``'format': [resnum]``
| Subsequent keys are the names of the PDB file with a value of: ``[jc_values]``.
| **About:** The default back-calculator uses the Karplus curve, a cosine function, to back-calculate desired
  J-couplings according to the residue number as provided by the experimental template file.

single-molecule Fluoresence Resonance Energy Transfer (smFRET) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Input:** ``res1,res2,scaler``
| Where res1/res2 is the residue number for the first and second residue respectively.
  Scaler is the r0 Foster radius of the dye pair.
| **Output:** ``'format': { 'res1': [], 'res2': [], 'scale': []}``
| Subsequent keys are the names of the PDB file with a value of: ``[smfret_values]``.
| **About:** The default back-calculator, abliet distance based, takes into consideration residue pairs
  and a scale factor to adjust for dye size from the experimental setup to back-calculate distances between
  two CA backbone atoms.

Residual Dipolar Coupling (RDC) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Input:** The default back-calculator for RDCs is PALES. To run PALES, however, you
will need the provide an experimental file template that is the same format as for PALES.
For more information, please see the `Data Format <https://spin.niddk.nih.gov/bax/software/PALES/>`_
section.
| **Output:** ``'format': {resnum1: [], resname1: [], atomname1: [], resnum2: [], resname2: [], atomname2: []}``
| Subsequent keys are the names of the PDB file with a value of: ``[rh_values]``.
| **About:** The default back-calculator uses the third-party program PALES, which uses the steric obstruction model to derive the RDC. It has also been chosen due to its popularity in the field for RDC back-calculations.

.. note::
    SAXS, Rh, and CS modules do not require input files. The following is formatting for the output.

Small Angle X-ray Scattering (SAXS) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Output:** Each PDB name is assigned a list of ``index`` and ``value``
  where ``index`` represents the X-axis and ``value`` represents the Y-axis with
  units ``I_abs(s)[cm^-1]/c[mg/ml]``.
| **About:** The default back-calculator uses the third-party program CRYSOLv3, from the ATSAS suite of biological softwares.
  it has been chosen due to its popularity as well as its ability to evaluate the hydration shell using dummy water parameters.

Hydrodynamic Radius (Rh) module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Output:** Each PDB name is assigned a singular hydrodynamic radius value.
| **About:** The default back-calculator uses the third-party program HullRad.
  It has been chosen for its speed and usage of the convex hull model to estimate hydrodynamic properties.

Chemical Shift (CS) module
^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Output:** ``'format': {'res': [], 'resname': []}``.
  Subsequent keys are the names of the PDB file with a value of:
  ``{'H': [], 'HA': [], 'C': [], 'CA': [], 'CB': [], 'N': []}``.
| **About:** The default back-calculator uses the third-party program UCBShift.
  It has been chosen for its two-pronged machine learning approach for both feature and sequence alignment
  in order to provide an accurate chemical shift prediction.

Basic Usage Examples
--------------------

The ``example/`` folder contains instructions to test native back-calculators
on 100 test conformers of the unfolded state of the drkN SH3 domain.

.. include:: ../example/README.rst
    :start-after: .. start-description
    :end-before: .. end-description
