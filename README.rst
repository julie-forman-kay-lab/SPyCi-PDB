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

**Structural Python (Back) Calculator Interface for PDBs**

**Goal:** User friendly Python3 based interface to generate back-calculated experimental data for singular and multiple (ensembles of) PDB structures.

**Back-Calculators:** The current back calculators integrated are:

#. PRE, NOE, J-Coupling, smFRET: Internal
#. Chemical shifts (CS): `UCBShift <https://github.com/THGLab/CSpred>`_
#. SAXS: `CRYSOLv3 <https://www.embl-hamburg.de/biosaxs/crysol.html>`_
#. Hydrodynamic Radius (Rh): `HullRad <http://52.14.70.9/>`_
#. Residual Dipolar Couplings (RDC): `PALES <https://spin.niddk.nih.gov/bax/>`_

Please note for third-party software, installation instructions have been fully
documented and tested for Linux Ubuntu Focal Fossa (20.04.X LTS) and Bionic Beaver (18.04.X LTS).

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

    TBD

.. end-citing

Version
=======

v0.3.5
