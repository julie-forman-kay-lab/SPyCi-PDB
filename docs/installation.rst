============
Installation
============

SPyCi-PDB uses only Python-based APIs for which we expect it to run
native on any system Python can run, as long as third-party back-calculator
installation requirements are met.

Please note that the IDPConformerGenerator library is required for parsing
PDB structure co-ordinates. Details for IDPConformerGenerator
can be found `here <https://github.com/julie-forman-kay-lab/IDPConformerGenerator>`_.

Full installation instructures are highlighted below.

From Source
-----------

Clone from the official repository::

    git clone https://github.com/julie-forman-kay-lab/SPyCi-PDB

Navigate to the new ``SPyCi-PDB`` folder::

    cd SPyCi-PDB

Run the following commands to install ``spycipdb`` dependencies if
Anaconda is used as your Python package manager::

    conda env create -f requirements.yml
    conda activate spycipdb

.. note::
    If you don't use Anaconda to manage your Python installations, you can use
    ``virtualenv`` and the ``requirements.txt`` file following the commands:

    | ``virtualenv spycienv --python=3.9``
    | ``source spycienv/bin/activate``
    | ``pip install -r requirements.txt``

    If you have difficulties installing ``spycipdb``, raise an issue in the
    main GitHub repository, and we will help you.

Install ``spycipdb`` in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps

.. note::
    The above applies also if you used ``virtualenv`` instead of ``conda``.

**Remember** to active the ``spycipdb`` environment every time you open a new
terminal window, from within the repository folder, choose yours::

    # Installation with Anaconda
    conda activate spycipdb

    # Installation with virtualenv
    source spycienv/bin/activate

To update to the latest version, navigate to the repository folder, activate the
``spycipdb`` python environment as described above, and run the commands::

    git pull

    # if you used anaconda to create the python environment, run:
    conda env update -f requirements.yml

    # if you used venv to create the python environment, run:
    pip install -r requirements.txt  --upgrade

    python setup.py develop --no-deps

Your installation will become up to date with the latest developments.


Installing IDPConformerGenerator
--------------------------------

.. note::
    You should be in the parent directory where ``SPyCi-PDB`` was cloned to.
    The ``spycipdb`` conda environment should not be active. Deactivate using::
        
        conda deactivate
    
    if you're using ``virtualenv``, remain in the environment.

Clone from the official repository::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator

Navigate to the new ``IDPConformerGenerator`` folder::

    cd IDPConformerGenerator

Run the following commands to install ``idpconfgen`` dependencies if
Anaconda is used as your Python package manager::

    conda env update --name spycipdb --file requirements.yml --prune
    conda activate spycipdb
    python setup.py develop --no-deps
    
Run the following commands to install ``idpconfgen`` dependencies if
virtualenv was used to install SPyCi-PDB::

    pip install -r requirements.txt
    python setup.py develop --no-deps

Go back to the ``SPyCi-PDB`` directory and reinstall ``spycipdb``::

    cd ..
    cd SPyCi-PDB
    python setup.py develop --no-deps


Installing Third-party Software
---------------------------------------

Some functionalities of ``SPyCi-PDB`` require third-party software.
These are not mandatory to install unless you want to use such operations.

UCBShift
````````

.. note::
    Module is only required if you wish to perform Chemical Shift back-calculations.
    
    You should be in the parent directory where ``SPyCi-PDB`` was cloned to
    with the Python environment is deactivated.

Clone the UCBShift repository from `my fork <https://github.com/menoliu/CSpred>`_ on
GitHub.::

    git clone https://github.com/menoliu/CSpred

Enter the ``CSpred`` folder and make a ``models`` directory, then download and
extract the latest `trained models <https://datadryad.org/stash/dataset/doi:10.6078/D1B974>`_
to ``CSpred/models`` directory.

Move back into the ``SPyCi-PDB`` directory and enter the ``thirdparty/ucbshift_reqs`` folder::

    cd ..
    cd SPyCi-PDB
    cd ./thirdparty/ucbshift_reqs

Run the following commands to install ``UCBShift`` dependencies if
Anaconda is used as your Python package manager::

    conda env update --name spycipdb --file ucbshift_requirements.yml --prune

Run the following commands to install ``UCBShift`` dependencies if
virtualenv was used to install SPyCi-PDB::

    pip install -r ucbshift_requirements.txt

Go back to the ``SPyCi-PDB`` directory and reinstall ``spycipdb``::
    
    conda activate spycipdb
    cd ../..
    python setup.py develop --no-deps

Again with virtualenv::

    source spycienv/bin/activate
    cd ../..
    python setup.py develop --no-deps

.. note::
    ``idpconfgen`` may need to be reinstalled while the ``spycipdb`` or
    ``spycienv`` is active as well after installing UCBShift using::
        ``python setup.py develop --no-deps``

ATSAS v3.1.1 - CRYSOL v3.0
``````````````````````````

.. note::
    ATSAS installation is only required for the :code:`saxs` module.

Please visit the `ATSAS website <https://www.embl-hamburg.de/biosaxs/download.html>`_
to download v3.1.1 of ATSAS. Theoretically, SPyCi-PDB will work if you already
have ATSAS v3.X installed.

Test your installation via::

    crysol -h

PALES v6.0
``````````

.. note::
    PALES installation is only required for the :code:`rdc` module.

A package of PALES v6.0 for Linux is already included in the ``thirdparty/`` folder.
Downloaded and extracted from the `Bax website <https://spin.niddk.nih.gov/bax/software/PALES/index.html>`_.

For use with x64 bit Linux Ubuntu 20.04.X LTS and 18.04.X, you must install the i386 architecture
along with required package libraries::

    sudo dpkg --add-architecture i386
    sudo apt update
    sudo apt install libc6:i386 libncurses5:i386 libstdc++6:i386 libx11-6:i386

If the last command above fails, run the following instead::

    sudo apt install multiarch-support
