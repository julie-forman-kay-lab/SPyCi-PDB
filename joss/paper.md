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

Structural determination of proteins has been a central scientific focus since the early 1960s [@Dill2008] with technological advances facilitating experimental structures of stable, folded proteins by nuclear magnetic resonance (NMR) spectroscopy [@Kanelis2001], X-ray crystallography [@Smyth2000], and cryo-electron microscopy [@Malhotra2019], as well as the recent computational prediction of structures [@Jumper2021; @Baek2021]. Modeling intrinsically disordered proteins (IDPs) and intrinsically disordered regions (IDRs), however, remains challenging due to their highly dynamic nature and low propensity to form low energy folded structures [@Mittag2007]. 

Currently, approaches to model IDPs/IDRs generally start with initial pools of structures that sample potentially accessible conformations and then utilize experimental data to narrow the pool. One method to generate initial conformational ensembles of IDPs/IDRs uses sampling techniques such as in TraDES [@Feldman2000; @Feldman2001], Flexible-meccano [@Ozenne2012], FastFloppyTail [@Ferrie2020], IDPConformerGenerator [@Teixeira2022], and others [@Estaa2019] that rely on the torsion angle distributions found in high-resolution folded protein structures deposited in the RCSB Protein Data Bank [@Berman2000]. Another more computationally expensive approach generates conformational ensembles using molecular dynamics (MD) simulations with different force-fields [@Salvi2016; @Robustelli2018]. 

After generating the initial pool of structures, back-calculations to experimental data and reweighting using Monte-Carlo [@Krzeminski2012] or Bayesian statistics [@Lincoff2020; @Bottaro2020] can be performed to define structural ensembles that better match solution NMR, small-angle X-ray scattering (SAXS), single molecule fluorescence (SMF), and other experimentally obtained data from these IDPs/IDRs. An emerging method to generate conformations of IDPs/IDRs uses machine learning generative models based on ensembles generated from sampling or MD techniques as training data and reinforces learning with experimental data [@Zhang2022]. Both of these general approaches rely on back-calculation of "experimental observables" from coordinates of conformers within the ensembles, a task that is increasingly complex due to the various models for interpretation of experimental data and the numerous tools available. 

Here we present **SPyCi-PDB**, designed to facilitate and streamline this back-calculation stage by acting as a platform for internal back-calculator functions as well as published third-party software, utilizing PDB structures of disordered protein conformations. One goal of **SPyCi-PDB** is to minimize the existing issues with different data-formats from software and scripts within the IDP/IDR research community and improve accessibility to researchers with a range of computational expertise. In this release, **SPyCi-PDB** can back-calculate NMR chemical shift (`cs`), paramagnetic resonance enhancement (`pre`), nuclear Overhauser effect (`noe`), <sup>3</sup>J-HNHA coupling (`jc`), and residual dipolar coupling (`rdc`) data; hydrodynamic radius (`rh`) data from NMR, light scattering or size exclusion chromatography; SAXS (`saxs`); and single-molecule fluorescence resonance energy transfer (`smfret`) values from all-atom PDB structures of IDP/IDR conformations. 

# Statement of Need

As new software packages and *in silico* methodologies emerge to better model IDP/IDR structures,
back-calculations to multiple experimental datatypes are required to quantitatively assess the conformers generated.
The current state of back-calculating PRE and NOE distances for dynamic protein systems, are not robust as the correlation
of two interacting groups of protons, in the case of NOEs, cannot be described by a simple distance. Furthermore, the data
from PRE back-calculations are improperly used in IDPs as a wider range of distances are sampled for disordered compared
to folded proteins. Although advancements in back-calculation of chemical shifts has been made, as seen in UCBShift [@Li2020],
improvements can still be made as the error in some cases is larger than the deviation of the back-calculated values. Given
the rapid developing nature of different software tools to perform back-calculations, **SPyCi-PDB** should assist by providing
a user-friendly, all-in-one package to reduce time and confusion in this back-calculation step as well as
open opportunities for future collaborations and integration of new experimental datatypes. Furthermore,
**SPyCi-PDB** aims to unify different input and output data formats from different experimental datatypes
to increase productivity and accelerate research. As stated in the documentation hosted by ReadTheDocs,
input formats are conventional comma-delimited tables (e.g. `.CSV, .TXT`), while the output format is human-readable `.JSON` files that
can be easily manipulated using Python or other software based on the user's ultimate needs. **SPyCi-PDB**
was also developed to integrate into the IDPConformerGenerator platform.

Ultimately, given the complicated and dynamic exchanging nature of IDPs, back-calculators will still need to be developed
with these features in mind. By creating a tool with modularity and best practices, we aim to encourage other
researchers to contribute towards this platform to further the goal of improved modelling of IDPs and IDRs.

# Implementation

As **spycipdb** is written completely in Python, it is compatible with any platform able to execute Python
(>=3.8, <4.0). However, certain third-party extensions to perform back-calculations (SAXS and RDC) have
only been tested on 64-bit Ubuntu 18.04.X LTS and 20.04.X LTS, as well as WSL 2.0 on 64-bit Windows 11.

In the production version 0.1.X, four out of eight modules of  **SPyCi-PDB**'s back-calculators
(`pre`, `noe`, `jc`, `smfret`) use internal mathematical algorithms and PDB structure processing using
IDPConformerGenerator libraries [@Teixeira2022]. The `pre` (1) and `noe` (2) module calculates scalar distances
between pairs of atoms according to the pairs derived from the experimental template. It utilizes an algorithm
that matches atom names of each residue with allowance for multiple assignments for `noe`. The `jc` (3) module uses
the Karplus curve, a simple cosine function, to back-calculate the desired J-couplings according to residue
number as provided by the experimental template file [@Perez2001]. Finally, the `smfret` (4) module takes into
consideration residue pairs and a scale factor to adjust for dye size from the experimental setup to back-calculate
distances between two alpha-Carbon (CA) atoms. The equations mentioned above are as follows:

$$
\begin{aligned}
\sqrt{\delta x^2 + \delta y^2 + \delta z^2}\ (1)\\\\
\sqrt[6]{(\frac{(\delta x^2 + \delta y^2 + \delta z^2)^3}{N}}\ (2)\\\\
\cos(\varphi - \frac{\pi}{3})\ (3)\\\\
\frac{1}{1 + (\frac{D \cdot \sqrt{\frac{|R_1 - R_2|+7}{R_1 - R_2}}}{S})^6}\ (4)
\end{aligned}
$$

Where δx, δy, δz are the cartesian differences between two atoms of interest (1, 2), N represents the number of combinations
for NOE atom pairs (2), φ is the Phi torsion angle of interest (3), D is the scalar distance between the residues of interest
with R_1 and R_2 being the vector cartesian co-ordinates for the residues and S being the scale factor according to experimental information.

The remaining 4 modules (`cs`, `saxs`, `rh`, `rdc`) call upon third-party academic software: UCBShift, a machine learning algorithm that
uses structural alignment for experimental chemical shift replication and employs a random forest regression on curated data to most accurately
predict protein chemical shifts [@Li2020]; CRYSOL v3, an updated version of the well established SAXS back-calculator from ATSAS that can now
evaluate the hydration shell by populating the protein structure with dummy water for SAXS calculations [@Franke2017];
HullRad, to calculate hydrodynamic radis (Rh) by using a convex hull model to estimate the hydrodynamic properties of a macromolecule [@Fleming2018];
PALES, using the obstruction model to derive dipolar coupling (RDC) information from the average orientation of the a 3-D structural model [@Zweckstetter2000].
Thorough testing of each module has been performed to ensure smooth installation and troubleshooting as well as retaining or providing
multiprocessing capabilities that may not have been implemented in their standalone form. When choosing third-party software, we prioritized
recent published endeavors written in Python for ease of integration such as UCBShift, and HullRad.

We are also open to integrate other options for back-calculators such as DEER-PREdict [@Tesei2021] for PREs, and alternative
methods to calculating experimental datatypes internally such as using a parameterizable fluoresence liftime and the Förster
distance as used in the Naudi-Fabra et al. study of describing intrinsically disordered proteins using
smFRET, NMR, and SAXS [@Naudi-Fabra2021]. Future additions to the **SPyCi-PDB** interface suite are welcome and easy to perform
given its design with modularity in mind.

# Example Usage

The current version of **spycipdb** has eight command-line modules that execute different
back-calculation routines for PDB formatted protein structure files. Usage examples, along
with input/output formats are provided both in the project's documentation hosted on ReadTheDocs and by the command:

```bash
$ spycipdb -h
```

To provide an example, we showcase the built-in PRE module with `pre` used to back-calculate
PRE distances given an experimental data template as seen in the `example/` folder:

```bash
$ spycipdb pre /drksh3_csss_100/ -e /example/drksh3_exp_data/drksh3_PRE.txt -n 16 -o bc_pdbs_PRE.json
```

Users can specify a directory or tarball (.TAR) with `N` number of PDB files and an experimental template file
`--exp-file`,`-e`, along with how many cores/workers to use `--ncores`,`-n` and the output file `--output`,`-o`.

# Acknowledgements

The motivation behind this project was to create a modular-yet-standalone software package to back-calculate
experimental datatypes for conformers generated by the IDPConformerGenerator [@Teixeira2022] platform. The
author thanks João M.C. Teixeira (ORCID: 0000-0002-9113-0622) for his mentorship regarding Best-Practices and
Python continuous integration programming practices [@pyskeleton].

# References