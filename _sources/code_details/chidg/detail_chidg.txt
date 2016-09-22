==============
ChiDG Instance
==============


ChiDG exists primarily as a Fortran Class(``chidg_t``) that provides an interface to the simulation environment. 
A ``chidg_t`` instance contains five primary components. These are:

.. code-block:: fortran
    
    type :: chidg_t
        type(chidg_data_t)          data
        type(time_scheme_t)         time_scheme
        type(nonlinear_solver_t)    nonlinear_solver
        type(linear_solver_t)       linear_solver
        type(preconditioner_t)      preconditioner
        ...
    end type chidg_t



|
|
|

----------
Components
----------


Data
----

::

    type(chidg_data_t)  data


The ``chidg_data_t`` data type contains all of the working information to run a ChiDG 
calculation. This includes:

::

    type(mesh_t)                    mesh(:)
    type(bcset_t)                   bcset(:)
    type(equationset_wrapper_t)     eqnset(:)
    type(solverdata_t)              sdata


