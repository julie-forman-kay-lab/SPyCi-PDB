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

To use the bare-bones version of :code:`spycipdb`, the PRE, NOE, JC, Rh,
and smFRET modules do not require third-party installation instructions.

Every module is equipped with :code:`--help` sections with detailed usage
examples and documentation for each functionality. For customized :code:`--output`
file names, the flag ``-o`` can be used with every module.

Please note for large number of PDB ensembles, it is recommended to specify
the number of CPU cores you are comfortable with using to maximize speed.
Having just :code:`-n` will utilize all but one CPU thread.

NOE, and JC Modules
-------------------
To perform the back-calculation using the default internal calculators,
give the tarball of provided structures as the first argument as well as the
path to the experimental data file to use as a template. For example, using the
NOE module::

    spycipdb noe ./drksh3_csss_100.tar.xz -e ./drksh3_exp_data/drksh3_NOE.txt -n

Likewise to the NOE module, the tarball of the provided structures and sample
JC experimental data file are required::

    spycipdb jc ./drksh3_csss_100.tar.xz -e ./drksh3_exp_data/drksh3_JC.txt -n

PRE Module
----------
There are two method of back-calculating PRE values within SPyCi-PDB. The first
uses a distance-based approach as seen in the original X-EISD `paper <https://www.nature.com/articles/s42004-020-0323-0>`_.
The second method uses `DEERPREdict <https://github.com/KULL-Centre/DEERpredict>`_ with
adjustable experimental parameters and yields intensity ratios.

Using the default method, the usage is akin to the NOE and JC modules above.::

    spycipdb pre ./drksh3_csss_100.tar.xz -e ./drksh3_exp_data/drksh3_PRE.txt -n

To use DEERPREdict within SPyCi-PDB, you must specify ``--method deerpredict``
as well as a parameters text file with the following default parameters.::
    
    atom=H
    temp=298
    tau_c=2*1e-9
    tau_t=0.5*1e-9
    delay=10e-3
    r_2=10
    wh=750

Where ``temp`` is the integer temperature value of the experiment in Kelvin, 
``tau_c`` and ``tau_t`` are the rotational tumbling time and internal
correlation time respectively, ``delay`` is the indept delay within the pulse
sequence, ``r_2`` is the diamagnetic transverse relaxation rate in the diamagnetic
molecule, and ``wh`` is the strength of the magnetic field. These parameters can be
changed and specified by using ``--parameters ./deerpredict_params.txt``. An example
of usage could be.::

    spycipdb pre ./PATH_TO_CONFORMERS -m deerpredict --parameters ./PATH_TO_TXT -e ./PATH_TO_EXP_DATA -n

The output will represent intensity ratios as ``Ipara/Idia``.

CS Module - Using UCBShift
--------------------------
After ensuring UCBShift is installed, the CS module does not require experimental
file samples. Furthermore, you could also adjust the pH value to be considered.
The default pH value is 5.0. A sample command is as follows::

    spycipdb cs ./drksh3_csss_100.tar.xz --ph 7 -n

The above command sets a custom pH of 7.0. Please also note that UCBShift is fairly
RAM intensive, it's recommended to run with less than 10 CPUs (can be changed with
the flag ``-n 10`` for example).

SAXS Module - Using CRYSOLv3
----------------------------
After ensuring the proper ATSAS/CRYSOL version is installed, the following command
can be used to run the SAXS module. Please note again that an experimental template
is not required::

    spycipdb saxs ./drksh3_csss_100.tar.xz -n

The SAXS module is equipped with a ``--lm`` flag, as CRYSOLv3 uses an adjustable
number of harmonics to perform the appropriate SAXS back-calculation. The default
value is 20 (from 1-100). Increasing the number of harmonics will also increase
the cost of computational time.

Rh Module - Using HullRadSAS v3.1
-----------------------------
HullRadSAS should be working out of the box with the basic installation instructions.
An experimental file template is not required::

    spycipdb rh ./drksh3_csss_100.tar.xz -n

.. end-description
