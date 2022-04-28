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
   
   $ nvcc --version

3. according to mpi and cuda versions, download different installation packages:

a. mpi \--> Open MPI 2.1.0,  cuda \--> 8.0:

> download `cuda-8.0-beta.zip <http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-beta.zip>`_, or use wget

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-beta.zip

b. mpi \--> Open MPI 2.1.0,  cuda \--> 10.1:

> download `cuda-10.1-beta.zip <http://www.pwmat.com/pwmat-resource/mstation-download/cuda-10.1-beta.zip>`_, or use wget

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-10.1-beta.zip

c. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 8.0:

> download `cuda-8.0-betaintel.zip <http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-betaintel.zip>`_, or use wget

.. code-block::
   
   $ wget http://www.pwmat.com/pwmat-resource/mstation-download/cuda-8.0-betaintel.zip

d. mpi \--> Intel(R) MPI Version 5.1.3,  cuda \--> 10.1:

> download `cuda-10.1-betaintel.zip <http://www.pwmat.com/pwmat-resource/mstation-download/cuda-10.1-betaintel.zip>`_, or use wget

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

Input Files
------------

As a minimal setup, PWmat requires the user to prepare the following input files: 
**etot.input, atom.config, pseudopotential file (for example: Si.SG15.PBE.UPF)**.

If present (in the directory where the calculation runs) the following output files 
of previous runs may be read as restart information **OUT.WG** and/or **OUT.RHO** 
**(Note: you must rename IN.WG and or IN.RHO, and add the tag 
"IN.WG = T" and/or "IN.RHO = T" in etot.input)**.

Some specific features of PWmat require additional output files.

Output Files
------------

he main output file of PWmat is the **REPORT**. The **RELAXSTEPS** file contains the total energies of the electronic and ionic SCF steps, and it is useful for the monitoring of the calculations.

Here is a comprehensive list of all important output files:

Directly readable files:
~~~~~~~~~~~~~~~~~~~~~~~~

============ ======================================================================
File name    Simple introduction
============ ======================================================================  
REPORT       Main output file  
RELAXSTEPS   Information about each electronic and ionic SCF step
MDSTEPS      Information about each electronic and molecular dynamics step
MOVEMENT     Contains the atomic position, atomic force, et al for each ionic step
NEB.BARRIER  Contains the energies along the images for each relaxation steps
final.config Is the updated atom.config file after each calculation
OUT.KPT      Contains the k-point vectors and their weights
OUT.SYMM     Contains symmetry operation information
OUT.OCC      The occupation of eigen states
OUT.VATOM    The atom center potential for SCF or MD simulation
OUT.FERMI    Contains the fermi energy 
OUT.FORCE    Forces on the atoms
OUT.STRESS   Stess tensor
OUT.QDIV     The atomic charge on each atom
OUT.ENDIV    The decomposed atomic energy on each atom
OUT.ATOMSPIN Contains local charge and magnetic moment when SPIN = 222
============ ======================================================================  

Binary files:
~~~~~~~~~~~~~

=================== ==============================================================================================
File name           Simple introduction
=================== ==============================================================================================
OUT.WG              Contains wave function
OUT.HSEWR(i)        Real space wave functions for the Fork exchange kernel for all the extended k-points on GPU(i)
OUT.REAL.RHOWF_SP   The charge density or wave function in real space
OUT.RHO             Charge density output file
OUT.RHO_2           Charge density for spin down components
OUT.RHO_SOM         A 2x2 complex spin matrix density
OUT.RHO_4DIELECTRIC The rho_e to be used to generate the dielectric function
OUT.RHO_POLARIZE    The solvent induced polarization charge
OUT.V_POLARIZE      The polarization potential generated by the polarization charge OUT.RHO_POLARIZE
OUT.RHOP_VHION      The polarization charge multiplied by the electric static potentail of the solute molecule
OUT.VR              Total potential output file
OUT.VR_hion         Hartree + Vion, the electrostatic potential without XC potential
OUT.VR_2            Potential for spin down components
OUT.VR_SOM          A 2x2 complex spin matrix potential
OUT.VR_DELTA        A real up-down potential
OUT.SPIN_X/Y/Z      Spin charge density in x/y/z direction at every r point
OUT.EIGEN           The eigen energies output file
bpsiiofil10000x     Wave function to atomic orbital projection file for kpoint x           
=================== ==============================================================================================

Examples
---------

.. toctree::
   /Examples/Si_SCF_Calculation
   /Examples/Si_DOS_Calculation
   /Examples/Si_Band_Calculation
   /Examples/GaAs_HSE_SCF_Calculation
   /Examples/GaAs_HSE_Band_Calculation
   /Examples/Si_RELAX_Calculation
   /Examples/Graphene_CELL_REL_Calculation
   /Examples/GaAs_HSE_REL_Calculation
   /Examples/BN_C_Charge_Relaxation_Cal
   /Examples/Cu-Au-Alloy-MD-Calculation
   /Examples/Ni_SCF_Spin_Calculation
   /Examples/CdSe_SOC_Band_Calculation
   /Examples/Fe_Spin222_SCF_Calculation
   /Examples/Ni_SCF_DFTU_Calculation
   /Examples/Trans_dipole_moment_exclude_Calculation
   /Examples/Trans_dipole_moment_include_Calculation
   /Examples/NEB_Calculation
   /Examples/Dimer_Calculation

Tutorials
---------

Modules
-------
.. toctree::
   /phonon/phonon
   /optical/optical
   /defect/defect
   /ultrafast/ultrafast

Utilities
---------

Cases
-----

