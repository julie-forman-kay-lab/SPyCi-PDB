
Changelog
=========

v0.5.3 (2025-02-04)
------------------------------------------------------------

* Fixed DEERPREdict parsing bug for ``pre`` module

v0.5.2 (2025-01-30)
------------------------------------------------------------

* Update README.rst to account for DEERPREdict and future requests

v0.5.1 (2025-01-30)
------------------------------------------------------------

* Correct documentation for DEERPREdict versioning
* Fixed tests for ``pre`` module

v0.5.0 (2025-01-30)
------------------------------------------------------------

* Add DEERPREdict as one of the calculators for ``pre`` data.
* Updated documentation for installation and usage of DEERPREdict within SPyCi-PDB

v0.4.3 (2024-11-20)
------------------------------------------------------------

* Fix a bug where ``cs`` module was duplicating ``CB`` shifts into ``N`` during the output.

v0.4.2 (2024-11-14)
------------------------------------------------------------

* Update CLI and documentation for Rh calculation as seen in v0.4.0

v0.4.1 (2024-11-13)
------------------------------------------------------------

* Linting

v0.4.0 (2024-11-13)
------------------------------------------------------------

* Update ``rh`` module to use HullRadSAS v3.1
* Update RtD documentation compilation
* Update installation instructions for UCBShift

v0.3.15 (2024-09-18)
------------------------------------------------------------

* Correct multiprocessing bug in the ``rdc`` module

v0.3.14 (2024-01-12)
------------------------------------------------------------

* Update UCBShift requirements

v0.3.13 (2023-07-27)
------------------------------------------------------------

* Correct formatting error in documentation

v0.3.12 (2023-07-27)
------------------------------------------------------------

* Update MANIFEST to include existing_requirements.txt

v0.3.11 (2023-07-27)
------------------------------------------------------------

* Update requirements for existing IDPConformerGenerator environment
* Update documentations for installation for clarity

v0.3.10 (2023-05-11)
------------------------------------------------------------

* Update README to include citation and JOSS button

v0.3.9 (2023-05-10)
------------------------------------------------------------

* Correct repository for versioning and Git actions

v0.3.8 (2023-05-10)
------------------------------------------------------------

* Update minor "co-ordinates" to "coordinates" in docs

v0.3.7 (2023-05-10)
------------------------------------------------------------

* Update changelog and minor formatting for README

v0.3.6 (2023-05-10)
------------------------------------------------------------

* Update paper affiliations

v0.3.5 (2023-05-02)
------------------------------------------------------------

* Implement minor editorial suggestions from issue #51

v0.3.4 (2023-04-26)
------------------------------------------------------------

* Minor corrections to documentation

v0.3.3 (2023-04-15)
------------------------------------------------------------

* Update acknowledgements section of manuscript

v0.3.2 (2023-04-13)
------------------------------------------------------------

* Update and make references uniform for manuscript

v0.3.1 (2023-04-12)
------------------------------------------------------------

* Correct 'ATSAS 2.8' formatting in references of paper
* Include in RtD documentation explanations of back-calculators
* Combine formatting for input/output in usage section of RtD documentation

v0.3.0 (2023-04-10)
------------------------------------------------------------

* Update results in the figure to include PRE and NOE as per #43
* Add plotting feature to ``rh``, ``jc``, ``pre``, ``noe`` modeules

v0.2.3 (2023-03-27)
------------------------------------------------------------

* Update example usage documentation

v0.2.2 (2023-03-07)
------------------------------------------------------------

* Update installation documentation
* Temporarily change requirements for default IDPConformerGenerator repository link

v0.2.1 (2023-03-06)
------------------------------------------------------------

* Update client setup
* Add warning for updating the environment when pulling

v0.2.0 (2023-01-27)
------------------------------------------------------------

* Adds ``natsort`` as dependency to yield ordered results
* Results will be sorted as how they will appear on your OS

v0.1.16 (2023-01-18)
------------------------------------------------------------

* Fixes input file extension validataion for PDB files
* Implements data validation for input experimental files
* Clarifies documentation for these changes
* Addresses issues in #35

v0.1.15 (2023-01-16)
------------------------------------------------------------

* Install IDPConformerGenerator while resolving the conda env
* Addresses bug in issue #33

v0.1.14 (2022-11-28)
------------------------------------------------------------

* Automatically catch ``CSpred`` missing import issue
* Addresses issue #31

v0.1.13 (2022-09-30)
------------------------------------------------------------

* Modify SAXS module to follow "format" in output

v0.1.12 (2022-09-29)
------------------------------------------------------------

* Minor edits for the paper

v0.1.11 (2022-09-28)
------------------------------------------------------------

* Bugfix for smFRET module with output formatting error

v0.1.10 (2022-09-20)
------------------------------------------------------------

* Re-test for PR

v0.1.9 (2022-09-20)
------------------------------------------------------------

* Add buttons on README

v0.1.8 (2022-09-19)
------------------------------------------------------------

* Update authors
* Add ASCII art for SPyCi-PDB

v0.1.7 (2022-09-19)
------------------------------------------------------------

* Create unit tests for internal calculators and parsers

v0.1.6 (2022-09-14)
------------------------------------------------------------

* Edits to manuscript per Sept 7 comments from Dr. Julie Forman-Kay
* Add figure 1 to paper

v0.1.5 (2022-09-08)
------------------------------------------------------------

* Edits to the manuscript per Aug 29 comments from Dr. Julie Forman-Kay
* Fix documentation error for RDC module
* Fix small issue of CS module output

v0.1.4 (2022-09-01)
------------------------------------------------------------

* Update failing tests
* Upload manuscript, bibliography, and tests for JOSS (#21)

v0.1.3 (2022-09-01)
------------------------------------------------------------

* Update RtD link in README.rst

v0.1.2 (2022-08-31)
------------------------------------------------------------

* Minor fix to gitworkflows for tests

v0.1.1 (2022-08-31)
------------------------------------------------------------

* Modularize all calculator components
* Remove Python 3.7 from requirements

v0.1.0 (2022-08-24)
------------------------------------------------------------

* Lint everything

v0.0.15 (2022-08-24)
------------------------------------------------------------

* Update README documentation
* Update ReadTheDocs format and associated docs

v0.0.14 (2022-08-23)
------------------------------------------------------------

* Upgrade CS module for multiprocessing with UCBShift
* Update installation instructions for UCBShift

v0.0.13 (2022-08-22)
------------------------------------------------------------

* Logic/module to link PALES v6.0 for RDC back-calculator (#14)
* Documentation for installing dependencies for PALES v6.0 for Ubuntu 20.04 LTS

v0.0.12 (2022-08-12)
------------------------------------------------------------

* Logic/module to link HullRad for Rh back-calculator (#13)

v0.0.11 (2022-08-12)
------------------------------------------------------------

* Logic/module to link CRYSOL 3.0 for SAXS back-calculator (#12)
* Documentation for installing CRYSOL 3.0 on top of SPyCi-PDB

v0.0.10 (2022-08-12)
------------------------------------------------------------

* Logic/module to link UCBShift for CS back-calculator (#10)
* Documentation for installing UCBShift on top of SPyCi-PDB

v0.0.9 (2022-08-10)
------------------------------------------------------------

* Logic/module for smFRET back-calculator (#9)

v0.0.8 (2022-08-10)
------------------------------------------------------------

* Logic/module for NOE back-calculator (#8)
* Refractor get_pdb_paths

v0.0.7 (2022-08-10)
------------------------------------------------------------

* Examples folder and some usage documentation (#7)

v0.0.6 (2022-08-10)
------------------------------------------------------------

* Logic/module for JC back-calculator (#6)

v0.0.5 (2022-08-09)
------------------------------------------------------------

* Logic/module for PRE back-calculator (#5)

v0.0.4 (2022-08-08)
------------------------------------------------------------

* Documentation for installing IDPConformerGenerator as a library (#4)

v0.0.3 (2022-08-08)
------------------------------------------------------------

* Core CLI backbone and base libs required (#2)
* Basic documentation for installation and updates

v0.0.2 (2022-08-08)
------------------------------------------------------------

* Fix reference to python-project-skeleton (#3)

v0.0.1 (2022-07-28)
------------------------------------------------------------

* Housekeeping items (#1)
* Building based on python-project-skeleton
* Renaming and changing base structure
