===============
Getting Started
===============


We need to get a copy of ChiDG up and running on your machine. Let's do it! The current option
we have to do this is to build ChiDG from source on your machine:

    - :ref:`Build from source`


..     - :ref:`Install from Python Package`


.. |
.. |
.. |
.. |
.. 
.. 
.. .. _Install from Python Package:
.. 
.. Install from Python Package
.. ===========================
.. 
.. **WARNING: Currently in pre-release. The package-based install is currently being
.. developed and not yet well-supported. We are working towards supporting this capability.**
.. 
.. 
.. ChiDG can be installed using the Python package manager ``pip``.
.. 
.. :: 
.. 
..     pip install chidg
.. 

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



Step 2: Check the compilers
---------------------------
To build ChiDG, you will need a Fortran compiler that is relatively up-to-date. Here are some
requirements and recommendations regarding compilers:

.. note:: 

    - The GNU Fortran compiler, ``gfortran``, is recommended. **Version >= 5.0 is required**
    - **Do not** use the Intel Fortran compiler, ``ifort``. There are currently several bugs in all versions of the compiler that make it unusable. One bug in the compiler creates a bad memory leak in the compiled code. Intel has been notified and was able to reproduce the issue. 
    - Other Fortran compilers have not been tested and are not supported.



Step 3: Install the dependencies
--------------------------------

.. note::  

    **Make sure** the dependencies listed here are available on your machine or that you
    build them and install them in an appropriate location.

    When installing dependencies, please place them somewhere other than the ChiDG project directory.
           
        - `CMake`_
        - `HDF5`_
        - `METIS`_
        - `MPI`_
        - `BLAS and LAPACK`_
   
   





CMake
~~~~~
ChiDG uses the CMake build-sytem to generate the ``Makefile`` configuration
files that are used when you run ``make`` to build the project. If you do not already have 
CMake installed on your system, you should be able to download a package/binary from 
the `CMake Website <https://cmake.org>`_. You will run the ``cmake`` executable in 
`Step 4: Configure the build`_.


HDF5
~~~~

.. note:: 

    | Version: ``>= 1.10``
    | Configure with: ``--enable-fortran``

HDF5 is used as a container for the ChiDG-formatted grid/solution files. The HDF5 library 
gives the option to build various portions of the library itself, depending on what the 
user needs available. ChiDG uses the HDF5 Fortran module interfaces and the high-level API. 
HDF5 documentation and installation guide can be found on
`HDF5 Website <https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL>`_.
Be sure to configure with the appropriate flags when building HDF5. You likely need to add the 
following flags to the configuration:

::

    ./configure --enable-fortran


The ChiDG CMake build will try to automatically locate the HDF5 libraries and
modules. If the HDF5 libraries were installed in a non-standard location, you can
help CMake find them by **setting an environment variable**, ``HDF5_ROOT``, such that 
the directory includes the ``bin``, ``lib``, and ``include``
directories where the libraries and modules can be located.

So, if you installed HDF5 in the directory ``/home/sparky/hdf5``, you might place 
the line 
``export HDF5_ROOT=/home/sparky/hdf5``
in your ``.bashrc`` or ``.bash_profile`` files if you are using the Bash shell. 



METIS
~~~~~
METIS is a partitioning library that is used to facilitate domain decomposition
for ChiDG and can be obtained from the 
`METIS Website <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_.

If METIS is installed in a non-standard location and CMake can't find it, 
you can help CMake find it by **setting an environment variable**: ``METIS_ROOT``.

So, if you installed METIS in the directory ``/home/sparky/metis``, you might place
the line 
``export METIS_ROOT=/home/sparky/metis``
in your ``.bashrc`` or ``.bash_profile`` files if you are using the Bash shell.


MPI
~~~

.. note::

    ChiDG uses the MPI Fortran 2008 module interfaces. Check that the MPI implementation
    being used supports and is configured to build these interfaces. OpenMPI v2.0.1 tries 
    to build these by default.

    Also:
        - MPI should be built with the same compiler that is being used to build ChiDG.
        - Once MPI is built, it is recommended to use the MPI Fortran compiler wrapper
          ``mpifort`` as the compiler for ChiDG in `Step 4: Configure the build`_

    Tested:
       - `OpenMPI <https://www.open-mpi.org>`_, Version 2.0.1.
       - Intel MPI  :  Found bugs in Intel implementation. Do not use.
    
MPI is used to facilitate the parallel communication in ChiDG. If your MPI installation is 
in a non-standard location, CMake will not be able to find it and it will tell you so. 
You may also want to force CMake to use a particular MPI installation if you have multiple 
installations on your machine. Either way, you can direct CMake find your MPI library 
by **setting the environment variable**: ``MPI_Fortran_COMPILER``.

So, if your MPI installation is located in the directory ``/home/sparky/mpi``, you might
place the line
``export MPI_Fortran_COMPILER=/home/sparky/my_mpi_install/bin/mpif90``
in your ``.bashrc`` or ``.bash_profile`` files if you are using the Bash shell.






BLAS and LAPACK
~~~~~~~~~~~~~~~
BLAS and LAPACK are libraries that contain routines for standard linear algebra operations
and algorithms. There are also many implementations of the BLAS and LAPACK libraries
by different vendors. Often times, a particular installation will have been optimized
for running on a given machine. Some common implementations of BLAS and LAPACK are the
Intel MKL, ATLAS, and the Apple Accelerate Framework.

If you do not have the BLAS and LAPACK libraries installed on your machine, reference
implementations can be downloaded at:

| `Reference BLAS implementation <http://www.netlib.org/blas>`_
| `Reference LAPACK implementation <http://www.netlib.org/lapack>`_
|

These will give the correct answers and can be used to get things up and running, 
however they are not optimized and so will degrade performance for ChiDG.

If your BLAS/LAPACK installations are in a non-standard location, CMake will not be able
to find it and it will tell you so. You can help CMake find them by **appending the location
of the libraries to the Linker path**.

.. note::

    On machines running LINUX:
        - export LD_LIBRARY_PATH=/my/path/to/blas:$LD_LIBRARY_PATH
        - export LD_LIBRARY_PATH=/my/path/to/lapack:$LD_LIBRARY_PATH

    On machines running Apple's OS X or macOS operating system:
        - export DYLD_LIBRARY_PATH=/my/path/to/blas:$DYLD_LIBARARY_PATH
        - export DYLD_LIBRARY_PATH=/my/path/to/lapack:$DYLD_LIBARARY_PATH




Step 4: Configure the build
---------------------------

.. note:: We should probably double-check a few things...

    - **Check** that all environment variables that were set for the dependencies are initialized in your environment.
      You may consider opening a new shell or running ``source ~/.bashrc`` or ``source ~/.bash_profile``.




+----------------------------------------------------+------------------------------------------------------------------+
| **Configure steps**                                                                                                   |
+----------------------------------------------------+------------------------------------------------------------------+
|                                                    |                                                                  |
| **#1** Change to the ChiDG root directory          | ``cd ChiDG``                                                     |
|                                                    |                                                                  |
+----------------------------------------------------+------------------------------------------------------------------+
|                                                    |                                                                  |
| **#2** Create a new build directory                | ``mkdir build``                                                  |
|                                                    |                                                                  |
+----------------------------------------------------+------------------------------------------------------------------+
|                                                    |                                                                  |
| **#3** Change to the build directory               | ``cd build``                                                     |
|                                                    |                                                                  |
+----------------------------------------------------+------------------------------------------------------------------+
|                                                    |                                                                  |
| **#4** Configure with CMake + user options         | ``cmake ..``  or ``cmake -DCMAKE_Fortran_COMPILER=mpifort ..``   |
|                                                    |                                                                  |
+----------------------------------------------------+------------------------------------------------------------------+




Regarding configure stage **#4**, configuration options can be passed when invoking ``cmake`` in order to influence the build
process. They are passed with the ``-D`` flag as:
 
:: 

    cmake -DParameter=Option ..


A typical build configure looks like:

::

    cmake -DCMAKE_Fortran_COMPILER=mpifort ..


A developer might configure the build using the following option:
 
::

    cmake -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_BUILD_TYPE=Debug ..


.. =============================== ======================================================= ================
.. Parameter                       Description                                             Options
.. =============================== ======================================================= ================
.. CMAKE_Fortran_COMPILER          Specify a Fortran Compiler       			``gfortran`` 
..                                                                                         
..                                                                                         
.. CMAKE_BUILD_TYPE                Specify what type of compiler settings to build with    ``Release``
..                                                                                         ``Debug``
.. =============================== ======================================================= ================


+------------------------------+-------------------------------------------------------+------------------------+
| Parameter                    |  Description                                          |  Options               |
+------------------------------+-------------------------------------------------------+------------------------+
|                              |                                                       |                        |
| CMAKE_Fortran_COMPILER       | Specify a Fortran Compiler                            | ``mpifort``            |
|                              |                                                       |                        |
+------------------------------+-------------------------------------------------------+------------------------+
|                              |                                                       |                        |
| CMAKE_BUILD_TYPE             | Specify what type of compiler settings to build with  | ``Release``            |
|                              |                                                       | ``Debug``              |
+------------------------------+-------------------------------------------------------+------------------------+










Step 5: Build ChiDG
-------------------


+-------------------------------------------+--------------------------------------------+
| **Build steps**                                                                        |
+-------------------------------------------+--------------------------------------------+
|                                           |                                            |
| Run ``make`` to build the ChiDG library   |   ``make`` or ``make -j 4``                |
|                                           |                                            |
+-------------------------------------------+--------------------------------------------+





Step 6: Test ChiDG
-------------------
ChiDG uses `pFUnit <http://pfunit.sourceforge.net>`_ to support unit and integration 
testing of the compiled ChiDG library.

.. note::

    Requires web access to retrieve pFUnit source.

+--------------------------------------------+------------------------------------------+
| **Testing steps**                                                                     |
+--------------------------------------------+------------------------------------------+
|                                            |                                          |
| Run ``make check`` to:                     |   ``make check`` or ``make check -j 4``  |
|   - download/build pFUnit                  |                                          |
|   - build ChiDG tests                      |                                          |
+--------------------------------------------+------------------------------------------+
|                                            |                                          |
| Run ``make test`` to:                      |                                          |
|   - run tests on ChiDG                     |   ``make test``                          |
+--------------------------------------------+------------------------------------------+








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
    - Assumes chidg.nml namelist and grid files are present in the working directory




|
|
|
|
|




.. _Running ChiDG:

Running ChiDG
=============

Step 1: Create a ChiDG-formatted grid file
------------------------------------------
To create a ChiDG-formatted grid, generate your grid as a **multi-block, 
unformatted, double-precision Plot3D file**. As an example, we will just 
assume you called this file, ``awesome_grid.x``.

    - Run ``chidg convert awesome_grid.x`` to create a ChiDG-formatted grid, ``awesome_grid.h5``.
    - Run ``chidg edit awesome_grid.h5`` and edit the boundary conditions for your grid.

|
|

Step 2: Create a ``chidg.nml`` input file
-----------------------------------------

    - Download a default :download:`chidg.nml <../examples/chidg.nml>` and place it in the 
      directory in which you plan to run ``chidg``.

    - Edit the ``chidg.nml`` entries as:
        
        +---------------------------------------------+
        | | ``gridfile         = 'awesome_grid.h5'``  |
        | | ``solutionfile_in  = 'none'``             |
        | | ``solutionfile_out = 'awesome_grid.h5'``  |
        +---------------------------------------------+
        
        
|
|
        
Step 3: Run ``chidg``
---------------------
In the working directory, call the ``chidg`` executable:

.. attribute:: Serial calculation

    ``chidg``


.. attribute:: Parallel calculation

    ``mpirun -np 3 chidg``



.. note:: 
    
    Once you get the hang of how the process works, you may consider playing around
    with the other entries in ``chidg.nml`` to understand how they affect the behavior
    and performance of ChiDG.









