=============================
Finding your way around ChiDG
=============================


ChiDG is a fairly large code and exists primarily as a composition of small concepts that
create the entire simulation environment. Keep in mind that the ChiDG environment 
exists as the ``chidg_t`` datatype, which resides in ``src/chidg/type_chidg.f90``.
To understand ChiDG, it is recommended to start with the ``chidg_t`` class, look
at what it contains. Start there and work your way into the code from the top down.
You can see what modules are being used by looking at ``use`` statements at the top 
of the file. 

.. seealso::

    :ref:`chidg_t`








---------------------------
Source Layout ``chidg/src``
---------------------------


ChiDG Modules
=============
Almost *everything* in ChiDG exists in Fortran Modules. Modules are nice because they allow 
a program to be broken up into logical units and used as necessary. They can also store 
data as module variables that could then be ``use``'d throughout the program as global data.
Additionally, they eliminate the need to write and maintain interface files like ``.h`` files
as is common in C and C++. In Fortran, the interfaces for routines that are located in a module
are automatically generated. If you need to make use of a ``mesh_t`` somewhere in the program, 
then ``use type_mesh, only: mesh_t`` will import the ``mesh_t`` data type from the ``type_mesh`` 
module and make it available to use in the current file.


Convention on ChiDG Modules:

 - Modules that define new data-types begin their name with ``type_``, as in ``type_mesh.f90``
 - Modules that mostly contain procedures begin their name with ``mod_``, as in ``mod_chimera.f90``



ChiDG Datatypes ``type_datatype.f90``
=====================================
ChiDG uses a lot of datatypes and classes for handling data and algorithms. Any user-defined
datatype in ChiDG exists in its own module. This helps keep ideas separate and avoids clashes
between routines of the same name. 

The convention for user-defined types in ChiDG is as follows:

    - A user-defined datatype(ex. ``newtype``) in ChiDG exists in its own module, ``type_newtype.f90``.
    - In the module ``type_newtype.f90``, the naming convention is the type with an appended ``_t`` as in ``newtype_t``
      to indicate throughout the program that is is a datatype.
    - The datatype can then be used throughout a program by importing its definition as, ``use type_newtype, only: newtype_t``.




.. code-block:: Fortran


    module type_newtype
        implicit none

        !>  This is a description of the datatype, newtype_t
        !!
        !!  @author Bob the Builder
        !!  @date   9/4/2016
        !!
        !-------------------------------------------------------
        type, public :: newtype_t 

            real :: my_real

        end type newtype_t
        !********************************************************

    end module







ChiDG Topics ``src/topic``
==========================
The ChiDG source directory is broken up into topic area directories. This is not out of 
necessity, rather for convenience and organization of ideas. Source files that pertain
to the representation of computational grids are probably located in ``src/grid``.
The base class for a linear system solver is probably located in ``src/linear_solvers``.





Algorithm Topics
----------------

Some topic directories refer to *algorithms* that may be dynamically switched out at run-time.
Some examples of these are ``src/linear_solvers``, ``src/nonlinear_solvers``, ``src/time_integrators``,
``src/preconditioners``. Inside an *algorithm* topic directory, convention is that there is a
base class with the same name. For example, inside ``src/linear_solvers`` there is a file 
``type_linear_solver.f90`` that implements the class ``linear_solver_t``. There are also concrete
implementations of the ``linear_solver_t``, such as ``fgmres_t`` inside of the ``src/linear_solvers``
topic directory.








