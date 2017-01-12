.. _chidg_t:

==============
ChiDG Instance
==============


ChiDG exists primarily as a Fortran Class(``chidg_t``) that provides an interface to 
the simulation environment. A ``chidg_t`` instance contains five primary components. 
These are:

.. code-block:: fortran
    
    type :: chidg_t
        type(chidg_data_t)          data
        type(time_integrator_t)     time_integrator
        type(nonlinear_solver_t)    nonlinear_solver
        type(linear_solver_t)       linear_solver
        type(preconditioner_t)      preconditioner
        ...
    end type chidg_t





.. function:: chidg%start_up('option')

    Performs start-up routines on a ChiDG instance. Options are:

    +---------------+-------------------------------------------------------------------------------+
    |  ``mpi``      | | Calls ``MPI_Init``, initializes ``IRANK``, ``NRANK``                        |
    +---------------+-------------------------------------------------------------------------------+
    | ``core``      | | Register functions, equations, models, operators, bcs, compute reference    |
    |               | | matrices. Initialize ``ChiDG_COMM`` MPI communicator. Default communicator  |
    |               | | is ``MPI_COMM_WORLD``.                                                      |
    +---------------+-------------------------------------------------------------------------------+
    | ``namelist``  | | Read namelist input file, chidg.nml                                         |
    +---------------+-------------------------------------------------------------------------------+

.. function:: chidg%shut_down('option')

    Performs shut-down routines on a ChiDG instance. Options are:

    +--------------+--------------------------------------------------------------------------------+
    |  ``mpi``     | | Calls ``MPI_Finalize``.                                                      |
    +--------------+--------------------------------------------------------------------------------+
    | ``core``     | | Close log files, close HDF interface.                                        |
    +--------------+--------------------------------------------------------------------------------+
    | ``log``      | | Close log files.                                                             |
    +--------------+--------------------------------------------------------------------------------+


.. function:: chidg%init('option')


    Performs shut-down routines on a ChiDG instance. Options are:

    +-------------------+--------------------------------------------------------------------------------+
    |  ``all``          | | Calls ``MPI_Finalize``.                                                      |
    +-------------------+--------------------------------------------------------------------------------+
    | ``communication`` | | Close log files, close HDF interface.                                        |
    +-------------------+--------------------------------------------------------------------------------+
    | ``chimera``       | | Close log files.                                                             |
    +-------------------+--------------------------------------------------------------------------------+
    | ``finalize``      | | Close log files.                                                             |
    +-------------------+--------------------------------------------------------------------------------+

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

    type(mesh_t)            mesh(:)
    type(bcset_t)           bcset(:)
    type(equation_set_t)    eqnset(:)
    type(solverdata_t)      sdata


