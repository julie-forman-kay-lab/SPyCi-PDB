Formatting and Usage Examples
=============================

.. start-description

To get started with using the different modules for back-calculating
experimental data, we will be using the unfolded state of drkN SH3, an
intensively studied intrinsically disordered protein.

The goal of these examples is to walk you through expected experimental data
formatting styles, as well as how to get started with each module. The
experimental files could be found in the :code:`example/drksh3_exp_data`
folder. A set of 100 conformers generated using `IDPConformerGenerator
<https://github.com/julie-forman-kay-lab/IDPConformerGenerator>`_ has
been provided as a tarball: :code:`example/drksh3_csss_100.tar.xz`.

Note that these experimental data files are comma-delimited per CSV
formatting for ease of use in :code:`pandas` dataframe as well as
Microsoft Excel usage.

To use the bare-bones version of :code:`spycipdb`, the PRE, NOE, JC, and
smFRET modules do not require third-party applications.

For example, in the :code:`pre` module, view its help using the following::

    spycipdb pre -h

To perform the back-calculation, give the tarball of provided structures
as the first argument and location of the PRE experimental data file::

    spycipdb pre ./drksh3_csss_100.tar.xz -e ./drksh3_exp_data/drksh3_PRE.txt -n

Please note for large number of PDB ensembles, it is recommended to specify
the number of CPU cores you are comfortable with using to maximize speed.
Having just :code:`-n` will utilize all but one CPU thread.

Every module is equipped with :code:`--help` sections with detailed usage
examples and documentation for each functionality.

.. end-description
