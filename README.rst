.. image:: https://github.com/jdtibochab/coralme/blob/main/docs/logo.png

The **COmprehensive Reconstruction ALgorithm for ME-models (coralME)** is an automatic pipeline for the reconstruction of ME-models. coralME integrates existing ME-modeling packages `COBRAme`_, `ECOLIme`_, and `solveME`_, generalizes their functions for implementation on any prokaryote, and processes readily available organism-specific inputs for the automatic generation of a working ME-model.

coralME has four main objectives:

1. **Synchronize** input files to remove contradictory entries.
2. **Complement** input files from homology with a template organism to complete the E-matrix.
3. **Reconstruct** a ME-model.
4. **Troubleshoot** the ME-model to make it functional.

For more information on coralME's inputs, architecture and usage, see coralME's documentation.

Installation
------------

Install using pip (available upon publication)
==============================================
1. ``pip install coralme``

Install locally
===============
1. Clone repository and navigate to coralme/.
2. ``pip install -r requirements.txt``
3. ``python3 setup.py clean build install``

Install using docker
====================
1. Clone repository and navigate to coralme/
2. ``docker build --file "./Dockerfile-Python3.10" . -t "python3.10-coralme"``
3. ``docker run --detach -p 10000:8888 -v USER/PATH/TO/coralme/:/opt/notebooks/ python3.10-coralme``
4. In your browser, go to ``localhost:10000``

Requirements
------------

coralME was tested with the following package versions:

- Python3, versions 3.7, 3.8, 3.9, and 3.10
- COBRApy version 0.26.3
- GUROBIpy version 9.5.2 (license is required)
- Ubuntu 22.04 is recommended (libgfortran.so.5 is required to execute MINOS and quad MINOS)
- Windows and MacOS users need to install `Gurobi`_ or `IBM CPLEX Optimizer <cplex_>`_. Alternatively, Windows users can install `WSL <wsl_>`_ and Ubuntu. Windows and MacOS users can use as well Docker Desktop to install it. We recommend the installation of Jupyter in the guest and its access through a browser from the host.

Compiled MINOS and quad MINOS are provided here as ``*.so`` files under ``coralme/solver``, and have been compiled using:

- Python3, versions 3.7.17, 3.8.17, 3.9.17, and 3.10.12
- wheel 0.38.4
- numpy 1.21.6
- scipy 1.7.3
- cython 0.29.32
- cpython 0.0.6

Documentation
-------------

You can find the documentation as a combined PDF called coralME_Documentation.pdf

.. refs
.. _COBRAme: https://github.com/SBRG/cobrame
.. _ECOLIme: https://github.com/SBRG/ecolime
.. _solveME: https://github.com/SBRG/solvemepy
.. _readthedocs: https://coralme.readthedocs.io/
.. _Gurobi: https://www.gurobi.com/
.. _cplex: https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer
.. _wsl: https://learn.microsoft.com/en-us/windows/wsl/install
