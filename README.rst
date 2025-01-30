SPyCi-PDB
=========
.. image:: https://github.com/julie-forman-kay-lab/SPyCi-PDB/blob/0e7b60d50a021d741dd7db4e2b9fbb9605fec95b/docs/spycipdb_ascii1.png

.. start-description

.. image:: https://github.com/julie-forman-kay-lab/SPyCi-PDB/actions/workflows/ci.yml/badge.svg?branch=main
    :target: https://github.com/julie-forman-kay-lab/SPyCi-PDB/actions/workflows/ci.yml
    :alt: Test Status

.. image:: https://readthedocs.org/projects/spyci-pdb/badge/?version=stable
    :target: https://spyci-pdb.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status

.. image:: https://joss.theoj.org/papers/10.21105/joss.04861/status.svg
   :target: https://doi.org/10.21105/joss.04861
   :alt: JOSS Paper

.. image:: https://zenodo.org/badge/518984240.svg
   :target: https://zenodo.org/badge/latestdoi/518984240
   :alt: Zenodo Archive

**Structural Python (Back) Calculator Interface for PDBs**

**Goal:** User friendly Python3 based interface to generate back-calculated experimental data for singular and multiple (ensembles of) PDB structures.

**Back-Calculators:** The current back calculators integrated are:

#. NOE, J-Coupling, smFRET: Internal
#. Paramagnetic Relaxation Enhancement (PRE): distance-based and `DEERPREdict <https://github.com/KULL-Centre/DEERpredict>`_
#. Chemical shifts (CS): `UCBShift <https://github.com/THGLab/CSpred>`_
#. SAXS: `CRYSOLv3 <https://www.embl-hamburg.de/biosaxs/crysol.html>`_
#. Hydrodynamic Radius (Rh): `HullRadSAS <http://52.14.70.9/HullRadSAS_descrip.html>`_
#. Residual Dipolar Couplings (RDC): `PALES <https://spin.niddk.nih.gov/bax/>`_

Please note for third-party software, installation instructions have been fully
documented and tested for the following Linux Ubuntu versions:
24.04 LTS, 22.04 LTS, 20.04 LTS, and 18.04 LTS.

To make new requests and/or additions, please see ``docs/contributing.rst``.

**Developer Notes:** project CI based on `@joaomcteixeira <https://github.com/joaomcteixeira>`_'s `Python-Project-Skeleton template <https://github.com/joaomcteixeira/python-project-skeleton>`_.
Developed as a standalone program with integration into the IDPConformerGenerator platform in mind.

.. end-description

Documentation
=============

More detailed documentation can be found at: https://spyci-pdb.readthedocs.io/en/stable/

Within the repository you can find:

#. Installation instructions in ``docs/installation.rst``.
#. Usage instructions in ``docs/usage.rst``.
#. See also the examples in the ``example/`` folder.

How to Cite
-----------

.. start-citing

If you use SPyCi-PDB, please cite::

    Liu et al., (2023). SPyCi-PDB: A modular command-line interface for back-calculating experimental datatypes of protein structures.. Journal of Open Source Software, 8(85), 4861, https://doi.org/10.21105/joss.04861

.. end-citing

Version
=======

v0.5.1
