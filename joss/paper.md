---
title: 'SPyCi-PDB: A modular command-line interface for back-calculating experimental datatypes of protein structures.'
tags:
  - Python
  - NMR
  - SAXS
  - smFRET
  - Back-calculator
  - PDB
  - Structural Biology
  - Proteins
  - Biochemistry
authors:
  - name: Zi Hao Liu
    orcid: 0000-0002-8357-8507
    affiliation: "1, 2"
    corresponding: true
  - name: Oufan Zhang
    affiliation: "3, 4"
  - name: Jie Li
    affiliation: "3, 4"
  - name: Teresa Head-Gordon
    orcid: 0000-0003-0025-8987
    affiliation: "3, 4, 5, 6"
  - name: Julie D. Forman-Kay
    orcid: 0000-0001-8265-972X
    affiliation: "1, 2"
    corresponding: true
affiliations:
  - name: "Molecular Medicine Program, Hospital for Sick Children, Toronto, Ontario M5G 0A4, Canada"
    index: 1
  - name: "Department of Biochemistry, University of Toronto, Toronto, Ontario, M5S 1A8, Canada"
    index: 2
  - name: "Pitzer Center for Theoretical Chemistry, University of California, Berkeley, California 94720-1460, USA"
    index: 3
  - name: "Department of Chemistry, University of California, Berkeley, California 94720-1460, USA"
    index: 4
  - name: "Department of Chemical and Biomolecular Engineering, University of California, Berkeley, California 94720-1462, USA"
    index: 5
  - name: "Department of Bioengineering, University of California, Berkeley, California 94720-1762, USA"
    index: 6
date: 25 August 2022
bibliography: paper.bib
---

# Summary

The protein folding problem has been sought out since the early 1960s [@Dill2008] and although
recent technological advances has made it possible to solve structures of stable, folded proteins: such as Nuclear Magnetic
Resonance spectroscopy [@Kanelis2001], X-ray crystallography [@Smyth2000], and more recently Cryo-electron microscopy 
[@Malhotra2019]. Modeling intrinsically disordered proteins and disordered regions (IDPs/IDRs) remain challenging due to
their highly dynamic nature and low-propensity to conform to an energy-minima folded structure [@Mittag2007].

Currently, the approaches to model IDPs/IDRs can be generalized into two groups of thought. The first, and more traditional
method is to generate conformational ensembles of IDPs/IDRs *de novo* using sampling techniques presented in TraDES [@Feldman2000; @Feldman2001],
Flexible-meccano [@Ozenne2012], FastFloppyTail [@Ferrie2020], IDPConformerGenerator [@Teixeira2022], and others [@Estaa2019] that primarily uses the torsion angle distributions
found in high-resolution folded protein structures deposited in the RCSB Protein Data Bank [@Berman2000]. Another popular, but computationally
expensive approach to generate conformational ensembles *ab initio* is the use of different force-fields within an Molecular Dynamics (MD)
simulation.

After generating the initial pool of structures, back-calculations to experimental data and reweighting using Monte-Carlo [@Krzeminski2012]
or Bayesian statistics [@Lincoff2020; @Bottaro2020] must be performed to obtain high quality structures that have a better fitment
with solution NMR, SAXS, smFRET, and other experimentally obtained data from these IDPs/IDRs.
An emerging method to generate conformations of IDPs/IDRs uses machine learning generative models that use ensembles generated 
from sampling or MD techniques as training data and reinforces learning with experimental and back-calculated datatypes.

**SPyCi-PDB** focuses on streamlining the back-calculation stage by acting as a platform for internal back-calculator functions as well as
published third-party software. The development of **SPyCi-PDB** was inspired as a way to minimize the existing issues with different
data-formats from softwares and scripts within the IDP/IDR research community and improves accessibility to non-computational researchers.
In this release, **SPyCi-PDB** can back-calculate paramagnetic resonance entropy (`pre`), nuclear overhauser effect (`noe`), 3J-HNHA coupling (`jc`),
chemical shifts (`cs`), small angle X-ray scattering (`saxs`), hydrodynamic radius (`rh`), residual dipolar couplings (`rdc`), and single-molecule 
fluorescence resonance energy transfer (`smfret`) values from all-atom PDB structures of IDP/IDR conformations.

# Statement of Need

As new software packages and *in silico* methodologies emerge to better model IDP/IDR structures, back-calculations to
multiple experimental datatypes are required to quantitatively assess the conformers generated. **SPyCi-PDB** hopes to
assist by providing a user-friendly, all-in-one package to reduce the time and confusion this step may bring as well as
open opportunities for future collaborations and integrations of new experimental datatypes.

Furthermore, **SPyCi-PDB** aims to unify different input and output data formats from different experimental datatypes
to increase productivity and accelerate research. As stated in the documentation hosted by ReadTheDocs, input formats
are conventional comma-delimited tables (e.g. `.CSV, .TXT`), while the output format is human-readable `.JSON` files that
can be easily manipulated using Python or other software based on the user's ultimate needs.

Lastly, **SPyCi-PDB** was developed with integration into the IDPConformerGenerator platform; created with
modularity and best practices in mind, we hope to use this example as what we hope other researchers could do to contribute
towards our ultimate goal of modelling IDPs and IDRs.

# Implementation

As of the current release of 0.1.X, four out of eight modules of  **SPyCi-PDB**'s back-calculators (`pre`, `noe`, `jc`, `smfret`)
use internal mathematical algorithms and PDB structure processing using IDPConformerGenerator libraries [@Teixeira2022].

The `pre` and `noe` module calculates scalar distances between pairs of atoms according to the pairs
derived from the experimental template. It utilizes an algorithm that matches atom names per residue
with allowance for multiple assignments for `noe`. The `jc` module uses a simple cosine function
to back-calculate the desired J-couplings according to residue number as provided by the experimental
template file. Finally, the `smfret` module takes into consideration residue pairs and a scale factor
to adjust for dye size from experimental data to back-calculate distances between two CA atoms.

The remaining 4 modules (`cs`, `saxs`, `rh`, `rdc`) call upon third-party academic software per the following:
UCBShift for chemical shifts [@Li2020], CRYSOL v3 for SAXS [@Franke2017], HullRad for Rh [@Fleming2018],
and PALES for RDC [@Zweckstetter2000]. 

Thorough testing have been performed to ensure smooth installation and troubleshooting as well as retaining
the multiprocessing features that may not have been implemented in their standalone form. Future additions
to the **SPyCi-PDB** interface suite are welcome and easy to perform given its design with modularity in mind.

# Installation

To install **spycipdb**, follow the instructions provided on the corresponding documentation
web page or `docs/installation.rst` in the repository. As **spycipdb** is written completely in Python,
it's comptable with any platform able to execute Python (>=3.8 but <4.0). However, certain third-party
extensions to perform back-calculations (SAXS and RDC) have only been tested on 64-bit Ubuntu 
18.04.X LTS and 20.04.X LTS, as well as WSL 2.0 on 64-bit Windows 11. It is recommended to read the
detailed installation and troubleshooting documentation before using Anaconda or virtualenv
to install and configure **spycipdb**.

# Example Usage

The current version of **spycipdb** has eight command-line modules that execute different
back-calculation routines for PDB formatted protein structure files. Usage examples, along
with input/output formats are provided both in the project's documentation hosted on ReadTheDocs and by the command:

```bash
$ spycipdb -h
```

To provide an example, we showcase the built-in PRE module with `pre` used to back-calculate
PRE distances given an experimental data template:

```bash
$ spycipdb pre /path/pdbs/ -e pdbs_PRE.txt -n 16 -o bc_pdbs_PRE.json
```

Users can specify a directory or tarball (.TAR) with `N` number of PDB files and an experimental template file
`--exp-file`,`-e`, along with how many cores/workers to use `--ncores`,`-n` and the output file `--output`,`-o`.

# Acknowledgements

The motivation behind this project was to create a modular-yet-standalone software package to back-calculate
experimental datatypes for conformers generated by the IDPConformerGenerator [@Teixeira2022] platform. The
author thanks João M.C. Teixeira (ORCID: 0000-0002-9113-0622) for his mentorship regarding Best-Practices and
Python continuous integration programming practices [@pyskeleton].

# References