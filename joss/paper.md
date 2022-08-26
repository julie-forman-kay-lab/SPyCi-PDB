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
  - name: Julie D. Forman-Kay
    orcid: 0000-0001-8265-972X
    affiliation: "1, 2"
    corresponding: true
affiliations:
  - name: "Molecular Medicine Program, Hospital for Sick Children, Toronto, Ontario M5G 0A4, Canada"
    index: 1
  - name: "Department of Biochemistry, University of Toronto, Toronto, Ontario, M5S 1A8, Canada"
    index: 2
date: 25 August 2022
bibliography: paper.bib
---

# Summary

The protein folding problem has been sought out since the early 1960s [10.1146/annurev.biophys.37.092707.153558] and although
recent technological advances has made it possible to solve structures of stable, folded proteins: such as Nuclear Magnetic
Resonance spectroscopy [@nmrstructure], X-ray crystallography [@xraystructure], and more recently Cryo-electron microscopy 
[@cryostructure]. Modeling intrinsically disordered proteins and disorered regions (IDPs/IDRs) remain challenging due to
their highly dynamic nature and low-propensity to conform to an energy-minima folded structure [@idpproblem1].

Currently, the approaches to model IDPs/IDRs can be generalized into two groups of thought. The first, and more traditional
method is to generate conformational ensembles of IDPs/IDRs *de novo* using sampling techniques presented in TraDES [https://pubmed.ncbi.nlm.nih.gov/10737933/; https://pubmed.ncbi.nlm.nih.gov/11746699/],
Flexible-meccano [10.1093/bioinformatics/bts172], FastFloppyTail [@FFT], IDPConformerGenerator [], and others [@ESTANA2019381] that primarily uses the torsion angle distributions
found in high-resolution folded protein structures deposited in the RCSB Protein Data Bank [@RCSB]. Another popular, but computationally
expensive approach to generate conformational ensembles *ab initio* is the use of different force-fields within an Molecular Dyanmics (MD)
simulation.

After generating the initial pool of structures, back-calculations to experimental data and reweighting using Monte-Carlo [@ENSEMBLE]
or Bayesian statistics [https://doi.org/10.1038/s42004-020-0323-0; 10.1007/978-1-0716-0270-6_15] must be performed to obtain
high quality structures that have a better fitment with solution NMR, SAXS, smFRET, and other experimentally obtained data from these IDPs/IDRs.
An emerging method to generate conformations of IDPs/IDRs uses machine learning generative models that use ensembles generated 
from sampling or MD techniques as training data and reinforces learning with experimental and back-calculated datatypes.

**SPyCi-PDB** focuses on streamlining the back-calculation stage by acting as a platform for internal back-calculator functions as well as
published third-party software. The development of **SPyCi-PDB** was inspired as a way to minimize the existing issues with different
data-formats from softwares and scripts within the IDP/IDR research community and improves accessibility to non-computational researchers.
In this release, **SPyCi-PDB** can back-calculate paramagnetic resonance entropy (`pre`), nuclear overhauser effect (`noe`), 3J-HNHA coupling (`jc`),
chemical shifts (`cs`), small angle X-ray scattering (`saxs`), hydrodynamic radius (`rh`), residual dipolar couplings (`rdc`), and single-molecule 
fluoresence resonance energy transfer (`smfret`) values from all-atom PDB structures of IDP/IDR conformations.

# Statement of Need




# Installation

To install **spycipdb**, follow the instructions provided on the corresponding documentation
webpage or `docs/installation.rst` in the repository. As **spycipdb** is written completely in Python,
it's comptable with any platform able to execute Python (>=3.8 but <4.0). However, certain third-party
extensions to perform back-calculations (SAXS and RDC) have only been tested on 64-bit Ubuntu 
18.04.X LTS and 20.04.X LTS, as well as WSL 2.0 on 64-bit Windows 11. It is recommended to read the
detailed installation and troubleshooting documentation before using Anaconda or virtualenv
to install and configure **spycipdb**.

# Example Usage

The current version of **spycipdb** has eight command-line modules that execute different
back-calculation routines for .PDB formatted protein structure files. Usage examples, along
with input/output formats are provided both in the project's documentation hosted on ReadTheDocs
as well as by the command:

```bash
$ spycipdb -h
```

To showcase the built-in PRE module, `pre` is used to back-calculate PRE distances
given an experimental data template:

```bash
$ spycipdb pre /path/pdbs/ -e pdbs_PRE.txt -n 16 -o bc_pdbs_PRE.json
```

Where users could specify a directory or tarball (.TAR) with `N` number of .PDBs and an experimental file
template `--exp-file`,`-e`, along with how many cores/workers to use `--ncores`,`-n` and the output file
`--output`,`-o`.

# Acknowledgements

The motivation behind this project was to create a modular-yet-standalone software package to back-calculate
experimental datatypes for conformers generated by the IDPConformerGenerator [@idpconfgen] platform. The author
thanks João M.C. Teixeira (ORCID: 0000-0002-9113-0622) for his mentorship regarding Best-Practices and Python
continous integration programming practices [@pyskeleton].

# References