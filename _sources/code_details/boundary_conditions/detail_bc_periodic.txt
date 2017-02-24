===========================
Periodic boundary condition
===========================

Set a ``Periodic`` condition.

============   =======================================================================
Parameter      Description
============   =======================================================================
``Offset-1``      Offset in coordinate direction 1 of the opposite periodic boundary.
``Offset-2``      Offset in coordinate direction 2 of the opposite periodic boundary.
``Offset-3``      Offset in coordinate direction 3 of the opposite periodic boundary.
============   =======================================================================



The ``Periodic`` boundary condition is actually a little bit special in ChiDG in that 
it does not actually compute an exterior state :math:`Q_{BC}` like all other boundary
condition state functions are expected to do. Rather, any faces in a ``bc_patch_t`` 
that are associated with a ``Periodic`` boundary condition are reset from BOUNDARY
faces to CHIMERA faces. Additionally, the offset parameters from the ``Periodic``
``bc_state_t`` object are stored in the ``face_t`` object to inform the Chimera
infrastructure to search for a donor, not at the original location, but at a location
defined by the offset coordinates. In this way, the ``Periodic`` boundary supports
fully non-matching boundaries through the Chimera infrstructure.



.. figure:: periodic_chimera.png
    :width: 270 pt
    :align: center
    :figclass: align-center








.. note::

    In setting up ``Periodic`` boundaries, recognize that opposite boundaries will 
    require opposite offset coordinates. To accomplish this, create two different
    boundary condition groups. Add a ``Periodic`` state function to each. When
    setting the coordinate offset parameters, make sure these parameters are 
    opposite sign between the two groups.



.. warning::

    A current limitation to the ``Periodic`` boundary condition implementation 
    is that it shall not be used to specify ``Periodic`` boundaries that 
    couple an element with itself. The Chimera infrastructure automatically
    throws the receiver element out as a potential donor candidate. This
    situation might be encountered when trying to run a '2D' calculation
    with 3D elements, but only using 1 element in the 3rd dimension. One
    might consider applying ``Periodic`` boundaries in the 3rd dimension,
    however the infrastructure will not set up this coupling correctly.
    An alternative would be to apply ``Symmetry``-like boundary conditions.





