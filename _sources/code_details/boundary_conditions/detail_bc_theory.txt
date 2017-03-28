=====================
Theory + Abstractions
=====================

Boundary conditions in ChiDG are abstracted and represented as a boundary condition 
data type, ``bc_t``. The boundary condition data type itself is separated into two
distinct concepts. These are **geometry** and **state**. The concepts of geometry and 
boundary condition state functions exist separetely from each other. Together, they 
form a ``bc_t``:


.. code-block:: Fortran

    type, public :: bc_t

        integer                                 :: BC_ID
        character(:),               allocatable :: bc_family
        type(bc_patch_t)                        :: bc_patch
        type(bc_state_wrapper_t),   allocatable :: bc_state(:)

    contains

        ...

    end type bc_t


``bc_t`` objects are responsible for defining an external solution state over 
a defined geometric boundary. That is, ``bc_state_t`` objects provide
an external solution state :math:`Q_{BC}` over a geometry region defined by 
``bc_patch_t``. The concept is that boundary conditions provide an external
solution state :math:`Q_{BC}` but do not actually compute a flux contribution
:math:`F(Q_{BC})` to any elements. Flux contributions are left to ``operator_t``
objects that are defined and added to ``equation_set_t`` objects that compute
:math:`F(Q_{BC})` assuming :math:`Q_{BC}` has already been defined by some 
boundary condition state function.

.. note:: 
    
    Boundary conditions that do not have any allocated ``bc_state_t`` function
    objects are flagged as Chimera boundaries. The ChiDG infrastructure will
    then automatically attempt to find Chimera donors for the faces associated
    with the ``bc_patch_t`` geometry representation.

    So, to initiate Chimera block-to-block communication, one is required only
    to **not** associate any boundary condition state functions with the block
    patches on the interior grid regions.


------------------------
Boundary condition patch
------------------------

The geometry over which a boundary condition is to be applied is represented 
by a boundary condition patch, ``bc_patch_t``. This consists of sets of triple 
indices(Domain,Element,Face) defining particular faces associated with the 
boundary condition patch. The ``bc_patch_t`` contains vectors for storing
these sets of integers.


.. code-block:: Fortran

    type, public :: bc_patch_t

        type(ivector_t)     :: idomain_l_
        type(ivector_t)     :: ielement_l_
        type(ivector_t)     :: iface_

    contains

        ...

    end type bc_patch_t


.. note::

    Boundary condition patches are currently created in the conversion process from 
    Plot3D-formatted grid files to ChiDG-formatted grid files. In the conversion
    process, the block boundaries in the Plot3D grid file are taken as boundary condition
    patches.




----------------------------------
Boundary condition state functions
----------------------------------

Boundary condition state functions, ``bc_state_t``, are responsible for defining 
the solution state exterior to a computational domain, :math:`Q_{BC}`. This can 
be computed as a function of the interior state and some user-specified parameters, 
:math:`Q_{BC}(Q_-,U_{input})`.



.. code-block:: Fortran

    type, public, abstract :: bc_state_t

        character(:),   allocatable :: name
        character(:),   allocatable :: family

        type(bcproperty_set_t)      :: bcproperties

    contains

        procedure(bc_state_init),       deferred    :: init
        procedure(bc_state_compute),    deferred    :: compute_bc_state

        ...

    end type bc_state_t













