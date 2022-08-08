============
Installation
============

SPyCi-PDB uses only Python-based APIs for which we expect it to run
native on any system Python can run, as long as third-party back-calculator
installation requirements are met.

Please note that the IDPConformerGenerator library is required for parsing
PDB structure co-ordinates.

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

    If you have difficulties installing ``spycipdb``, raise an Issue in the
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


Installing IDPConformerGenerator On Top of SPyCi-PDB
----------------------------------------------------

.. note::
    You should already be inside the ``SPyCi-PDB`` directory.
    The ``spycipdb`` conda environment should not be active. Deactivate using:
    
    | conda deactivate
    
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
    python setup.py develop --no-deps
