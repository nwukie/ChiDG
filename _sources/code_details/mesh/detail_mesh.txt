====
Mesh
====

The ``mesh`` data structure contains an entire geometry description for a single
domain. This exists as an array of ``element`` types, and array of ``face`` types,
and a ``chimera`` instance. An ``element`` exists for every element in the ``mesh``
domain. For a given ``element``, a ``face`` instance exists for each face.



.. image:: d_mesh_exploded.png
    :width: 90 %
    :align: center

.. image:: d_mesh_arrays.png
    :width: 90 %
    :align: center






-------------
Elements
-------------

An ``element`` instance contains information needed by the framework and also
general information that could be useful to users. This includes:

::

    elem_pts(:)     An array of points defining the element in real space(cartesian, cylindrical, etc.)
    quad_pts(:)     An array of points defining the location of each volume quadrature node in real space.
    metric(3,3,:)   An array, defining for each quadrature point, a matrix of element metric values.
    jinv(:)         An array of inverse element jacobian values at each volume quadrature node.
    ddx(:,:)        An array of derivatives of the basis functions with respect to real coordinates at volume quadrature nodes.


.. image:: d__element.png
    :width: 80 %
    :align: center


Metric terms
------------

The metric terms are defined at each quadrature point in the ``metric(:,:,:)`` component 
of a given ``element``. To access the matrix of metric components for a given quadrature 
node 'igq', the component can be used as

::

    metric(:,:,igq)

This returns the metric components at the quadrature node in a 3x3 matrix as:

.. math::

    \begin{pmatrix}
      \xi_x   \quad \xi_y   \quad   \xi_z \\
      \eta_x  \quad \eta_y  \quad   \eta_z \\
      \zeta_x \quad \zeta_y \quad   \zeta_z
    \end{pmatrix} 


Alternatively, a given metric term can be accessed for the set of quadrature nodes as

::

    metric(1,1,:)

which would return a 1D array of values for :math:`\xi_x` corresponding to each volume
quadrature node.

The inverse element jacobian terms ``jinv(:)`` are defined at each quadrature node as

.. math::

    J^{-1} = ( x_\xi \xi_x + x_\eta \eta_x + x_\zeta \zeta_x )


Derivatives
-----------

The derivatives of basis functions with respect to the computational coordinates on a 
reference element are already defined in a quadrature instance associated with an 
element in the component ``element%gq%vol``. For example, the component 
``element%gq%vol%ddxi`` gives:

.. math::

    \frac{\partial \psi_{igq, imode}}{\partial \xi} =
        \begin{pmatrix}
            \frac{\partial \psi_{1,1}}{\partial \xi} &  \frac{\partial \psi_{1,2}}{\partial \xi}  & \cdots  & \frac{\partial \psi_{1,N}}{\partial \xi} \\
            \frac{\partial \psi_{2,1}}{\partial \xi}  & \frac{\partial \psi_{2,2}}{\partial \xi}  & \cdots  & \frac{\partial \psi_{2,N}}{\partial \xi} \\
            \vdots & \vdots & \vdots & \vdots \\
            \frac{\partial \psi_{{ngq},1}}{\partial \xi} & \frac{\partial \psi_{{ngq},2}}{\partial \xi} &  \cdots &  \frac{\partial \psi_{{ngq},N}}{\partial \xi} \\
        \end{pmatrix}


Derivatives in real space coordinates in an ``element`` can be computed using 
``ddx(:,:)`` components. The derivatives of basis functions with respect to real 
coordinates(:math:`x,y,z` , :math:`r,\theta,z` ) are specific to each ``element`` 
and these derivatives can be accessed in the ``ddx``, ``ddy``, ``ddz`` components. 
The ``element%ddx`` component for example gives


.. math::

    \frac{\partial \psi_{igq, imode}}{\partial x} =
        \begin{pmatrix}
            \frac{\partial \psi_{1,1}}{\partial x} &  \frac{\partial \psi_{1,2}}{\partial x}  & \cdots  & \frac{\partial \psi_{1,N}}{\partial x} \\
            \frac{\partial \psi_{2,1}}{\partial x}  & \frac{\partial \psi_{2,2}}{\partial x}  & \cdots  & \frac{\partial \psi_{2,N}}{\partial x} \\
            \vdots & \vdots & \vdots & \vdots \\
            \frac{\partial \psi_{{ngq},1}}{\partial x} & \frac{\partial \psi_{{ngq},2}}{\partial x} &  \cdots &  \frac{\partial \psi_{{ngq},N}}{\partial x} \\
        \end{pmatrix}















-------------
Faces
-------------

.. image:: d__face.png
    :width: 90%
    :align: center


Face metrics
------------

Metric terms for the ``face`` data structure are defined exactly the same as for the 
``element`` data structure. The difference is that the ``metric`` and ``jinv`` components of 
``face`` return values for boundary quadrature nodes. This contrasts the ``element`` 
structure, which returns values for volume quadrature nodes.


Face normals
------------


.. math::

    \vec{x} = [x, y, z]     \quad   \vec{\xi} = [\xi, \eta, \zeta]

Face normal vectors are stored for each face quadrature node. The component ``norm`` is the
face normal vector with respect to computational coordinates on a reference element
(:math:`\xi`, :math:`\eta`, :math:`\zeta`) as

.. math::

    \vec{n}_{\xi_k} = \frac{\partial \vec{x}}{\partial \xi_i} \times \frac{\partial \vec{x}}{\partial \xi_j}


where :math:`\xi_i` and :math:`\xi_j` are the in-place coordinates of face :math:`\xi_k`.
In this was, the normal vectors for :math:`\xi`, :math:`\eta`, and :math:`\zeta` faces
are defined respectively as

.. math:: 

    \vec{n}_\xi = 
    \frac{\partial \vec{x}}{\partial \eta} \times \frac{\partial \vec{x}}{\partial \zeta} = 
    [ y_\eta z_\zeta - y_\zeta z_\eta, \quad x_\zeta z_\eta - x_\eta z_\zeta, \quad x_\eta y_\zeta - x_\zeta y_\eta] =
    [ \xi_x, \quad \xi_y, \quad \xi_z ]

    \vec{n}_\eta = 
    \frac{\partial \vec{x}}{\partial \zeta} \times \frac{\partial \vec{x}}{\partial \xi} = 
    [ y_\zeta z_\xi - y_\xi z_\zeta, \quad x_\xi z_\zeta - x_\zeta z_\xi, \quad x_\zeta y_\xi - x_\xi y_\zeta ] =
    [ \eta_x, \quad \eta_y, \quad \eta_z ]

    \vec{n}_\zeta = 
    \frac{\partial \vec{x}}{\partial \xi} \times \frac{\partial \vec{x}}{\partial \eta} = 
    [ y_\xi z_\eta - y_\eta z_\xi, \quad x_\eta z_\xi - x_\xi z_\eta, \quad x_\xi y_\eta - x_\eta y_\xi] =
    [ \zeta_x, \quad \zeta_y, \quad \zeta_z ]




Applying the above formula to element faces produces normal vectors that are inward 
facing for :math:`\xi = -1` faces and outward facing for :math:`\xi = 1` faces.
Inward facing vectors are negated so that all resultant normal vectors in the ``norm`` 
component are outward facing. This applied to :math:`\eta` and :math:`\zeta` faces as well.

Unit normal vectors can be accessed in the ``unorm`` component and are computed as

.. math::

    \hat{n}_{\xi_i} = \frac{\vec{n}_{\xi_i}}{||\vec{n}_{\xi_i}||_2}






------------------
Chimera Interfaces
------------------

Each ``mesh`` instance contains a ``mesh%chimera`` component that holds all information
regarding chimera communication for that particular mesh block. This takes the
form of ``chimera_receiver`` and ``chimera_donor`` components. Currently, only
the ``chimera_receiver`` is utilized. ``chimera_donor`` will be used to facilitate 
communication between processors for parallel code execution.


.. image:: d__chimera_receiver.png
    :width: 90 %
    :align: center



In a given ``mesh`` block, every face that gets information from a separate block is 
designated as a CHIMERA face, it is assigned an integer ID ``face%ChiID``, and it gets an 
entry in the ``mesh%chimera%recv%data`` components. It can be accessed as

::

    mesh%chimera%recv%data(ChiID)

Example
-------

Consider an example with two mesh domains, as shown below.
``mesh(1)`` contains four elements. ``mesh(2)`` contains eight elements.
``mesh(1)`` overlaps with ``mesh(2)``. In particular, the top faces of elements E3 and E4 lie 
inside ``mesh(2)``. These faces are designated as CHIMERA faces and are given a mesh-global
chimera ID. The top face of E3 is given the ID ChiID=1 and the top face of E4 is given
the ID ChiID=2.


.. image:: d__chimera_demo_a.png
    :width: 90 %
    :align: center


Each CHIMERA face has its own set of chimera information, which can be accessed via 
``mesh%chimera%recv%data(ChiID)``. This is shown below for the two faces in this example.



.. image:: d__chimera_demo_b.png
    :width: 90 %
    :align: center







































