====
Mesh
====



------------------
Coordinate Systems
------------------

ChiDG supports two coordinate systems: ``Cartesian`` and ``Cylindrical``. In the
code, many things are referred to only by coordinate index. The convention for these
indices is defined here as:


.. math:: 

    (x,y,z) \rightarrow (1,2,3)

    (r, \theta, z)  \rightarrow (1,2,3)


Additionally, effort has been put forth to represent data in a manner consistent
with vector-calculus. For example, test function gradient components are computed
as:

.. math::

    \nabla \psi = \frac{\partial \psi}{\partial x}\hat{x} + \frac{\partial \psi}{\partial y} \hat{y} + \frac{\partial \psi}{\partial z} \hat{z}

    \nabla \psi = \frac{\partial \psi}{\partial r}\hat{r} + \frac{1}{r}\frac{\partial \psi}{\partial \theta} \hat{\theta} + \frac{\partial \psi}{\partial z} \hat{z}









----
Mesh
----

The ``mesh_t`` data structure contains an entire geometry description for a single
domain. This exists as an array of ``element_t`` types, and array of ``face_t`` types,
and a ``chimera_t`` instance. An ``element_t`` exists for every element in the ``mesh_t``
domain. For a given ``element_t``, a ``face_t`` instance exists for each face.



.. image:: d_mesh_exploded.png
    :width: 90 %
    :align: center

.. image:: d_mesh_arrays.png
    :width: 90 %
    :align: center







Elements
-------------

An ``element_t`` instance contains information needed by the framework and also
general information that could be useful to users. This includes:

::

    elem_pts(:)     An array of points defining the element in real space(cartesian, cylindrical, etc.)
    quad_pts(:)     An array of points defining the location of each volume quadrature node in real space.
    metric(3,3,:)   An array, defining for each quadrature point, a matrix of element metric values.
    jinv(:)         An array of inverse element jacobian values at each volume quadrature node.
    grad1(:,:)      An array of gradients of the basis functions with respect to physical coordinates at volume quadrature nodes.


.. image:: d__element.png
    :width: 80 %
    :align: center


Metric terms
~~~~~~~~~~~~

Metric terms are computed by considering the transformation of a differential
volume in physical space to reference space as

.. math::

    \begin{pmatrix}
        dx \\ dy \\ dz
    \end{pmatrix}
    =
    \begin{pmatrix}
      x_\xi  \quad  x_\eta  \quad  x_\zeta \\
      y_\xi  \quad  y_\eta  \quad  y_\zeta \\
      z_\xi  \quad  z_\eta  \quad  z_\zeta
    \end{pmatrix} 
    \begin{pmatrix}
        d\xi \\ d\eta \\ d\zeta
    \end{pmatrix}
    \quad
    \quad
    \quad
    \begin{pmatrix}
        dr \\ rd\theta \\ dz
    \end{pmatrix}
    =
    \begin{pmatrix}
      r_\xi       \quad  r_\eta         \quad  r_\zeta \\
      r\theta_\xi \quad  r\theta_\eta   \quad  r\theta_\zeta \\
      z_\xi       \quad  z_\eta         \quad  z_\zeta
    \end{pmatrix} 
    \begin{pmatrix}
        d\xi \\ d\eta \\ d\zeta
    \end{pmatrix}

The terms :math:`\partial \vec{x}/\partial \vec{\xi}` are computed from the polynomial
expansion representing the element coordinates as

.. math::

    x = \sum \psi \hat{x}   \quad\rightarrow\quad  \frac{\partial x}{\partial \xi} = \sum \frac{\partial \psi}{\partial \xi} \hat{x}

The metric terms are computed by inverting the matrices 
:math:`\partial \vec{x}/\partial \vec{\xi}` to give

.. math::

    \begin{pmatrix}
      \xi_x   \quad \xi_y   \quad   \xi_z \\
      \eta_x  \quad \eta_y  \quad   \eta_z \\
      \zeta_x \quad \zeta_y \quad   \zeta_z
    \end{pmatrix} 
    \quad
    \quad
    \quad
    \begin{pmatrix}
       \xi_r   \quad   \xi_\theta   \quad   \xi_z  \\
       \eta_r  \quad   \eta_\theta  \quad   \eta_z \\
       \zeta_r \quad   \zeta_\theta \quad   \zeta_z
    \end{pmatrix} 

The metric terms are defined at each quadrature point in the ``metric(:,:,:)`` component 
of a given ``element_t``. To access the matrix of metric components for a given quadrature 
node ``igq``, the component can be used as

::

    metric(:,:,igq)

This returns the metric components(``Cartesian`` or ``Cylindrical``) at the quadrature node in a 3x3 matrix as

.. math::

    \begin{pmatrix}
      \xi_x   \quad \xi_y   \quad   \xi_z \\
      \eta_x  \quad \eta_y  \quad   \eta_z \\
      \zeta_x \quad \zeta_y \quad   \zeta_z
    \end{pmatrix} 
    \quad
    \quad
    \quad
    \begin{pmatrix}
       \xi_r   \quad   \xi_\theta   \quad   \xi_z  \\
       \eta_r  \quad   \eta_\theta  \quad   \eta_z \\
       \zeta_r \quad   \zeta_\theta \quad   \zeta_z
    \end{pmatrix} 

Alternatively, a given metric term can be accessed for the set of quadrature nodes as

::

    metric(1,1,:)

which would return a 1D array of values for (:math:`\xi_x` or :math:`\xi_r`) corresponding to each 
quadrature node in the set.

The inverse element jacobian terms(``Cartesian`` or ``Cylindrical``) ``jinv(:)`` are defined at each quadrature node as

.. math::

    J^{-1} = ( x_\xi y_\eta z_\zeta  -  x_\eta y_\xi z_\zeta  -  x_\xi y_\zeta z_\eta  +  x_\zeta y_\xi z_\eta  +  x_\eta y_\zeta z_\xi  -  x_\zeta y_\eta z_\xi )

    J^{-1} = r ( r_\xi \theta_\eta z_\zeta  -  r_\eta \theta_\xi z_\zeta  -  r_\xi \theta_\zeta z_\eta  +  r_\zeta \theta_\xi z_\eta  +  r_\eta \theta_\zeta z_\xi  -  r_\zeta \theta_\eta z_\xi )



Derivatives + Gradients
~~~~~~~~~~~~~~~~~~~~~~~

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


Gradients in pysical coordinates in an ``element_t`` can be computed using 
``grad1(:,:)``, ``grad2(:,:)``, and ``grad3(:,:)``  components. The gradients
of basis functions with respect to physical 
coordinates(:math:`x,y,z` , :math:`r,\theta,z` ) are specific to each 
``element_t`` and these derivatives can be accessed in the 
``grad1``, ``grad2``, ``grad3`` components. For example, the ``element%grad1`` component 
contains the gradient along the 1st physical coordinate  for all test functions at all
quadrature nodes as:


.. math::

    \nabla_1 \psi_{ngq, nmode} =
        \begin{pmatrix}
            \nabla_1 \psi_{1,1}  &  \nabla_1 \psi_{1,2}  & \cdots  & \nabla_1 \psi_{1,N} \\
            \nabla_1 \psi_{2,1}  &  \nabla_1 \psi_{2,2}  & \cdots  & \nabla_1 \psi_{2,N} \\
            \vdots & \vdots & \vdots & \vdots \\
            \nabla_1 \psi_{{ngq},1}  & \nabla_1 \psi_{{ngq},2}  &  \cdots &  \nabla_1 \psi_{{ngq},N} \\
        \end{pmatrix}
















Faces
-------------

.. image:: d__face.png
    :width: 90%
    :align: center


Face metrics
~~~~~~~~~~~~

Metric terms for the ``face_t`` data structure are defined exactly the same as for the 
``element_t`` data structure. The difference is that the ``metric`` and ``jinv`` components of 
``face_t`` return values for boundary quadrature nodes. This contrasts the ``element_t`` 
structure, which returns values for volume quadrature nodes.


Face normals
~~~~~~~~~~~~


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







Chimera Interfaces
------------------

Each ``mesh_t`` instance contains a ``mesh%chimera`` component that holds all information
regarding chimera communication for that particular mesh block. This takes the
form of ``chimera_receiver`` and ``chimera_donor`` components. Currently, only
the ``chimera_receiver`` is utilized. ``chimera_donor`` will be used to facilitate 
communication between processors for parallel code execution.


.. image:: d__chimera_receiver.png
    :width: 90 %
    :align: center



In a given ``mesh_t`` block, every face that gets information from a separate block is 
designated as a CHIMERA face, it is assigned an integer ID ``face%ChiID``, and it gets an 
entry in the ``mesh%chimera%recv%data`` components. It can be accessed as

::

    mesh%chimera%recv%data(ChiID)

Example
~~~~~~~

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




































