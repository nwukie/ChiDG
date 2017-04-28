.. _chidg_t:

=========
ChiDG API
=========


The ChiDG API is a Fortran Class(``chidg_t``) that provides an interface to 
the ChiDG library. 

.. code-block:: fortran
    
    type :: chidg_t

        type(chidg_data_t)                      :: data
        class(time_integrator_t),   allocatable :: time_integrator
        class(nonlinear_solver_t),  allocatable :: nonlinear_solver
        class(linear_solver_t),     allocatable :: linear_solver
        class(preconditioner_t),    allocatable :: preconditioner

    contains

        ! Open/Close
        procedure   :: start_up
        procedure   :: shut_down

        ! IO
        procedure   :: read_grid
        procedure   :: read_solution
        procedure   :: write_grid
        procedure   :: write_solution

        ! Initialization
        procedure   :: set
        procedure   :: init

        ! Run
        procedure   :: run
        
    end type chidg_t



.. .. note:: **Example:**

.. rubric:: **Example:**


.. code-block:: fortran

    ! This is an example of how a driver might use the ChiDG API to run an
    ! analysis. This is similar to the file src/driver.f90 that is used
    ! to create the ChiDG executable.

    ! Import the ChiDG API class
    use type_chidg, only: chidg_t

    ! Create an instance of the ChiDG API
    type(chidg_t)   :: chidg

    ! Initialize MPI and ChiDG infrastructure
    call chidg%start_up('mpi')
    call chidg%start_up('core')

    ! Set algorithms
    call chidg%set('Time Integrator' , algorithm='Steady'      )
    call chidg%set('Nonlinear Solver', algorithm='Quasi-Newton')
    call chidg%set('Linear Solver'   , algorithm='FGMRES'      )
    call chidg%set('Preconditioner'  , algorithm='ILU0'        )
    call chidg%set('Solution Order'  , integer_input=2         )

    ! Read grid (solution if available)
    call chidg%read_grid('grid.h5')
    call chidg%read_solution('grid.h5')

    ! Run ChiDG analysis
    call chidg%run()




.. rubric:: **API procedures:**


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


.. function:: chidg%set('String Selector', algorithm='string', integer_index=value)

    Set algorithms and parameters for the ChiDG instance. 

    +-----------------------+-----------------------------------------------------------------------+
    |  ``String Selector``  |  Expected input                                                       |
    +=======================+=======================================================================+
    |  ``Time Integrator``  |  algorithm='string'                                                   |
    +-----------------------+-----------------------------------------------------------------------+
    |  ``Nonlinear Solver`` |  algorithm='string'                                                   |
    +-----------------------+-----------------------------------------------------------------------+
    |  ``Linear Solver``    |  algorithm='string'                                                   |
    +-----------------------+-----------------------------------------------------------------------+
    |  ``Preconditioner``   |  algorithm='string'                                                   |
    +-----------------------+-----------------------------------------------------------------------+
    |  ``Solution Order``   |  integer_input=value                                                  |
    +-----------------------+-----------------------------------------------------------------------+

.. function:: chidg%read_grid('file name')

    Read domains and boundary conditions from the ChiDG-formatted HDF file indicated by the 
    string parameter passed in. 
    
    .. note:: All calls to ``chidg%set()`` shall occur before this.


.. function:: chidg%read_solution('file name')

    Read solution from the ChiDG-formatted HDF file indicated by the 
    string parameter passed in. 
    
    .. note:: A grid shall already have been read/initialized before a call to ``chidg%read_solution('file')``.



.. function:: chidg%run()

    Run the ChiDG analysis. This calls a certain number of 'steps' on the time integrator.

    .. note:: All reading/initializing/setting shall occur before this.


.. .. function:: chidg%init('option')
.. 
.. 
..     Performs initialization activities on data. This should be executed after add data has been
..     read in. In an operational setting, using the option ``all`` is recommended. The more granular
..     initialization options are sometimes useful for testing purposes, where not all initialization
..     is necessary. Options are:
.. 
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     |  ``all``          | | Calls all initialization routines.                                                          |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``domains``       | | Initialize domain data that depend on solution expansion.                                   |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``bc``            | | Initialize boundary condition coupling, bc parallel comm.                                   |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``communication`` | | Initialize interior communication. Local/Global interior faces.                             |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``chimera``       | | Initialize chimera communication.                                                           |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``solvers``       | | Allocate matrix/vector storage used by solvers.                                             |
..     +-------------------+-----------------------------------------------------------------------------------------------+
..     | ``finalize``      | | Check necessary algorithms have been set. Initialize time integrator and preconditioner.    |
..     +-------------------+-----------------------------------------------------------------------------------------------+
.. 












