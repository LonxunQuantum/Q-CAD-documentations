PWmat
=====

Download
--------
You can download PWmat from `LongXun  <http://www.pwmat.com/>`_ offcial website. Currently, we only provide 
the executable binary files, which are compiled according to different mpi and cuda versions. Before download,
it is necessary to confirm the mpi and cuda versions on your machines.

1. confirm the mpi version:

.. code-block::

   $ mpirun --version

2. confirm the cuda version:

.. code-block::
   
   $ cuda --version

3. according to mpi and cuda versions, download different installation packages:

a. mpi \--> Open MPI 2.1.0,  cuda \--> 8.0:

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-beta.zip

b. mpi \--> Open MPI 2.1.0,  cuda \--> 10.1:

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-10.1-beta.zip

c. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 8.0:

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-betaintel.zip

d. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 10.1:

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-10.1-betaintel.zip


Installation
------------

After downloading the corresponding installation package, you can install it in the following way:

a. mpi \--> Open MPI 2.1.0,  cuda \--> 8.0:

.. code-block:: 

   $ unzip cuda-8.0-beta.zip
   $ sh install.beta.run
   $ bash

b. mpi \--> Open MPI 2.1.0,  cuda \--> 10.1:

.. code-block:: 

   $ unzip cuda-10.1-beta.zip
   $ sh install.beta.run
   $ bash

c. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 8.0:

.. code-block:: 

   $ unzip cuda-8.0-betaintel.zip
   $ sh install.betaintel.run
   $ bash

d. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 10.1:

.. code-block:: 

   $ unzip cuda-10.1-betaintel.zip
   $ sh install.betaintel.run
   $ bash

Manual
------

The `manual <http://www.pwmat.com/pwmat-resource/Manual.pdf>`_ of PWmat is available on the LongXun offical website, please check it!

Tutorials and examples
----------------------

.. toctree::
   /tutorials_and_examples/Si_SCF_Calculation
   /tutorials_and_examples/Si_DOS_Calculation
   /tutorials_and_examples/Si_Band_Calculation
   /tutorials_and_examples/GaAs_HSE_SCF_Calculation
   /tutorials_and_examples/GaAs_HSE_Band_Calculation

Modules
-------

Utilities
---------

Cases
-----

