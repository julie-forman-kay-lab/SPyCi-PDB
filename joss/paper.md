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
affiliations:
  - name: "Molecular Medicine Program, Hospital for Sick Children, Toronto, Ontario M5G 0A4, Canada"
    index: 1
  - name: "Department of Biochemistry, University of Toronto, Toronto, Ontario, M5S 1A8, Canada"
    index: 2
date: 25 August 2022
bibliography: paper.bib
---

# Summary


# Implementation


# Installation


# Use cases

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


# References