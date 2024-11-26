Welcome to coralME documentation!
======================================

Description
~~~~~~~~~~~

The **COmprehensive Reconstruction ALgorithm for ME-models (coralME)** is an automatic pipeline for the reconstruction of ME-models. coralME integrates existing ME-modeling packages `COBRAme`_, `ECOLIme`_, and `solveME`_, generalizes their functions for implementation on any prokaryote, and processes readily available organism-specific inputs for the automatic generation of a working ME-model.

coralME has four main objectives:

1. **Synchronize** input files to remove contradictory entries.
2. **Complement** input files from homology with a template organism to complete the E-matrix.
3. **Build** a working ME-model
4. **Inform** the user about necessary steps to curate the ME-model.

This resource is intended to:

1. Describe basic inputs required for ME-model reconstruction.
2. Describe the architecture of coralME.
3. Demonstrate how to build a ME-model with coralME.
4. Describe how to perform manual curation guided by coralME's curation notes.

Content
~~~~~~~~~~~~~~~~~~

.. toctree::
   :numbered:
   :maxdepth: 2
    
   GettingStarted
   BasicInputs
   coralMEArchitecture
   ManualCuration
   FAQ

Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. refs
.. _COBRAme: https://github.com/SBRG/cobrame
.. _ECOLIme: https://github.com/SBRG/ecolime
.. _solveME: https://github.com/SBRG/solvemepy