===============
Getting Started
===============


We need to get a copy of ChiDG up and running on your machine. There are two options 
to accomplish this:

    - :ref:`Install from Python Package`
    - :ref:`Build from source`

|
|
|
|


.. _Install from Python Package:

Install from Python Package
===========================

**WARNING: Currently in pre-release. The package-based install is currently being
developed and not yet well-supported. We are working towards supporting this capability.**


ChiDG can be installed using the Python package manager ``pip``.

:: 

    pip install chidg


|
|
|
|

.. _Build from source:

Build from source
=================


Step 1: Get ChiDG
-----------------
ChiDG is hosted on GitHub and may be cloned using the https protocol as

::

    git clone https://github.com/nwukie/ChiDG.git ChiDG



Step 2: Install Dependencies
----------------------------

CMake
~~~~~
ChiDG uses the CMake build-sytem to generate the ``make`` configuration
files for the project. If you do not already have CMake installed on your
system, you should be able to download a package/binary from the CMake Website.


HDF5
~~~~
HDF5 is used as a container for the ChiDG-formatted grid/solution files and
is a required dependency of ChiDG. The HDF5 library gives the option to build 
various portions of the library itself, depending on what the user needs available. 
ChiDG uses the HDF5 Fortran module interfaces and the high-level API. Be sure to 
configure with the appropriate flags when building HDF5. Refer to the HDF5 documentation 
and installation guide for complete directions. Likely configuration parameters 
that are needed are:

::

    --enable-fortran
    --enable-fortran2003

The ChiDG CMake build will try to automatically locate the HDF5 libraries and
modules. To ensure they are found, **be sure to set the HDF5_ROOT environment 
variable** such that the directory includes the ``bin``, ``lib``, and ``include``
directories where the libraries and modules can be located.


``export HDF5_ROOT='HDF5-Install-Directory'``







Step 3: Configure the build
---------------------------

.. note:: We should probably double-check a few things...

    - **Check** the ``HDF5_ROOT`` environment variable is exported and is set to the install
      directory of the HDF5 library. The following should be valid: ``$HDF5_ROOT/bin``, ``$HDF5_ROOT/lib``.





+----------------------------------------------------+----------------------------------------------------------+
| **Configure steps**                                                                                           |
+----------------------------------------------------+----------------------------------------------------------+
|                                                    |                                                          |
| **#1** Change to the ChiDG root directory          | ``cd ChiDG``                                             |
|                                                    |                                                          |
+----------------------------------------------------+----------------------------------------------------------+
|                                                    |                                                          |
| **#2** Create a new build directory                | ``mkdir build``                                          |
|                                                    |                                                          |
+----------------------------------------------------+----------------------------------------------------------+
|                                                    |                                                          |
| **#3** Change to the build directory               | ``cd build``                                             |
|                                                    |                                                          |
+----------------------------------------------------+----------------------------------------------------------+
|                                                    |                                                          |
| **#4** Configure with CMake + user options         | ``cmake ..``  or ``cmake -DCMAKE_BUILD_TYPE=Release ..`` |
|                                                    |                                                          |
+----------------------------------------------------+----------------------------------------------------------+




The following parameters can be passed when invoking ``cmake`` in order to configure the build with 
different options. They are passed with the ``-D`` flag as:

::

    cmake -DParameter=Option ..


=============================== ======================================================= ================
Build parameter                 Description                                             Options
=============================== ======================================================= ================
CMAKE_Fortran_COMPILER          Specify a Fortran Compiler        
                                                                                        
                                                                                        
CMAKE_BUILD_TYPE                Specify what type of compiler settings to build with    ``Release``
                                                                                        ``Debug``
USER_MPI                        Use User environment variables for finding MPI          ``True``
                                                                                        ``False``
BUILD_TESTS                     Build the test binaries. Requires pFUnit Library        ``True``
                                                                                        ``False``
=============================== ======================================================= ================

**Up-to-date Fortran(F2008) compiler required. Tested with gfortran 5.3**













Step 4: Build ChiDG
-------------------


+-------------------------------------------+--------------------------------------------+
| **Build steps**                                                                        |
+-------------------------------------------+--------------------------------------------+
|                                           |                                            |
| Run ``make`` to build the ChiDG library   |   ``make`` or ``make -j 4``                |
|                                           |                                            |
+-------------------------------------------+--------------------------------------------+




|
|
|
|


.. _Usage:

Usage
=====

So, how does one actually use ChiDG? The interface for using the ChiDG library 
is still being designed and changed. One thing to keep in mind, is that 
the interface provided is just an interface. One could create their own driver 
interface by linking to the library and compiling an executable. For now, we 
will just focus on detailing how ChiDG currently gets used.


ChiDG Executable: ``chidg``
---------------------------

When you build and install ChiDG, it builds an executable file. This is simply
a driver file(``src/interfaces/driver.f90``) that uses routines from the ChiDG library. 
The driver file gets compiled, linked to the ChiDG library, and put in an executable, ``chidg``.
The ChiDG executable, ``chidg``, works as:

    - a utility for converting Plot3D grid files 
    - a utility for editing boundary conditions in the ChiDG-format HDF file
    - a utility for processing a solution for viewing in Tecplot or Paraview
    - a driver to run ChiDG simulations


.. attribute:: chidg convert file.x

    - Convert a Plot3D grid file(example in this case: 'file.x') to a ChiDG-formatted HDF5 file.


.. attribute:: chidg edit file.h5

    - Edit a ChiDG-formatted HDF5 grid file(example in this case: 'file.h5').


.. attribute:: chidg post file.h5

    - Post-process a ChiDG-formatted HDF5 grid/solution file. 
    - Creates Tecplot/Paraview files for visualization

.. attribute:: chidg

    - Run a ChiDG simulation
    - Assumes chidg.nml namelist file is present in current directory

