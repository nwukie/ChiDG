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

The ChIDG CMake build will try to automatically locate the HDF5 libraries and
modules. To ensure they are found, **be sure to set the HDF5_ROOT environment 
variable** such that the directory includes the ``bin``, ``lib``, and ``include``
directories where the libraries and modules can be located







Step 3: Configure the build
---------------------------



+-------------------------------------------------------+---------------------------------------------------------------+
| **Out-of-source build configuration**                                                                                 |
+-------------------------------------------------------+---------------------------------------------------------------+
|                                                       |                                                               |
| #1 Change to the ChiDG root directory                 |   ``cd ChiDG``                                                |
|                                                       |                                                               |
+-------------------------------------------------------+---------------------------------------------------------------+
|                                                       |                                                               |
| #2 Create a new build directory                       |   ``mkdir build``                                             |
|                                                       |                                                               |
+-------------------------------------------------------+---------------------------------------------------------------+
|                                                       |                                                               |
| #3 Change to the build directory                      |   ``cd build``                                                |
|                                                       |                                                               |
+-------------------------------------------------------+---------------------------------------------------------------+
|                                                       |                                                               |
| #4 Configure with CMake + user options                |   ``cmake ..``  or ``cmake -DCMAKE_BUILD_TYPE=Release ..``    |
|                                                       |                                                               |
+-------------------------------------------------------+---------------------------------------------------------------+



The following parameters can be passed when invoking ``cmake`` in order to configure the build with different options. They are passed with the ``-D`` flag as:

::

    cmake -DParameter=Option ..


================================= ========================================================= ===========================================================================
Build parameter                   Description                                               Options
================================= ========================================================= ===========================================================================
CMAKE_Fortran_COMPILER            Specify a Fortran Compiler        
                                                                                            
                                                                                            
CMAKE_BUILD_TYPE                  Specify what type of compiler settings to build with      ``Release``
                                                                                            ``Debug``
USER_MPI                          Use User environment variables for finding MPI            ``True``
                                                                                            ``False``
BUILD_TESTS                       Build the test binaries. Requires pFUnit Library          ``True``
                                                                                            ``False``
================================= ========================================================= ===========================================================================

**Up-to-date Fortran(F2008) compiler required. Tested with gfortran 5.3**













Step 4: Build ChiDG
-------------------


+-------------------------------------------------------+---------------------------------------------------------------+
|                                                       |                                                               |
| Run ``make`` to build the ChiDG library               |   ``make`` or ``make -j 4``                                   |
|                                                       |                                                               |
+-------------------------------------------------------+---------------------------------------------------------------+




