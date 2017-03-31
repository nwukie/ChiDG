===============
Fluid Equations
===============

Primary field variables for the fluid equations in ChiDG are ``Density``, 
``Momentum-1``, ``Momentum-2``, ``Momentum-3``, and ``Energy``


.. math::

      Q_{cart} = 
      \begin{pmatrix}
        \rho \\ \rho u_x \\ \rho u_y \\ \rho u_z \\ \rho E 
      \end{pmatrix}
      \quad\quad
      Q_{cyl} = 
      \begin{pmatrix}
        \rho \\ \rho u_r \\ r \rho u_\theta \\ \rho u_z \\ \rho E 
      \end{pmatrix}



.. note::

    In ``Cylindrical`` coordinates, angular momentum, :math:`r \rho u_\theta`, is being conserved. Not tangential momentum, :math:`\rho u_\theta`.


The transformation of the momentum vector between ``Cartesian`` and ``Cylindrical`` coordinates is

.. math::

    \begin{bmatrix}
        \rho u_x \\ \rho u_y \\ \rho u_z
    \end{bmatrix}
    =
    \begin{bmatrix}
        cos(\theta) &   -sin(\theta)  &  0 \\
        sin(\theta) &    cos(\theta)  &  0 \\
        0 & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
        \rho u_r \\ \rho u_\theta \\ \rho u_z
    \end{bmatrix}

The velocity vector in ``Cartesian`` and ``Cylindrical`` coordinates is

.. math::

    \vec{u} = [u_x, u_y, u_z]       \quad\quad  \vec{u} = [u_r, u_\theta, u_z]

At times, when a formula is independent of the coordinate system, generic 
components of velocity will be used

.. math::

    \vec{u} = [u,v,w]


Within the framework, effort has been invested in representing quantities
such that they are consistent with vector calculus. As such, representations
of a gradient in ``Cartesian`` and ``Cylindrical`` coordinates as

.. math:: 

    \nabla f &= \bigg[\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z}\bigg] 
    \\
    \\
    \nabla f &= \bigg[\frac{\partial f}{\partial r}, \frac{1}{r}\frac{\partial f}{\partial \theta}, \frac{\partial f}{\partial z}\bigg]

are represented in terms of the components of the gradient vector as

.. math::

    \nabla f = [\nabla_1 f, \nabla_2 f, \nabla_3 f]

    
|
|
|
|
|
|

---------------
Fluid Advection
---------------


The fluxes governing fluid advection are


Cartesian
---------

.. math:: 

      \vec{F}_a(\vec{x},Q) =
      \begin{bmatrix}
        \rho     u_x        \\ 
        \rho u_x u_x  +  p  \\
        \rho u_y u_x        \\
        \rho u_z u_x        \\
        \rho H   u_x
      \end{bmatrix}
      ,
      \begin{bmatrix}
        \rho     u_y        \\ 
        \rho u_x u_y        \\
        \rho u_y u_y  +  p  \\
        \rho u_z u_y        \\
        \rho H   u_y 
      \end{bmatrix}
      ,
      \begin{bmatrix}
        \rho     u_z       \\ 
        \rho u_x u_z        \\
        \rho u_y u_z        \\
        \rho u_z u_z  +  p  \\
        \rho H   u_z
      \end{bmatrix}
      \quad\quad
      S_a(\vec{x},Q) = 
      \begin{pmatrix}
          0 \\ 0 \\ 0 \\ 0 \\ 0
      \end{pmatrix}


Cylindrical
-----------

.. math:: 

      \vec{F}_a(\vec{r},Q) =
      \begin{bmatrix}
          \rho          u_r        \\ 
          \rho u_r      u_r  +  p  \\
        r \rho u_\theta u_r        \\
          \rho u_z      u_r        \\
          \rho H        u_r
      \end{bmatrix}
      ,
      \begin{bmatrix}
          \rho          u_\theta        \\ 
          \rho u_r      u_\theta        \\
        r \rho u_\theta u_\theta  +  r p  \\
          \rho u_z      u_\theta        \\
          \rho H        u_\theta
      \end{bmatrix}
      ,
      \begin{bmatrix}
          \rho          u_z         \\ 
          \rho u_r      u_z         \\
        r \rho u_\theta u_z         \\
          \rho u_z      u_z  +   p  \\
          \rho H        u_z
      \end{bmatrix}
      \quad\quad
      S_a(\vec{r},Q) = 
      \begin{pmatrix}
          0 \\ \frac{\rho u_\theta u_\theta - p}{r} \\ 0 \\ 0 \\ 0
      \end{pmatrix}


:math:`H = \frac{\rho E + p}{\rho}` is the total enthalpy and :math:`p` is the static pressure.

|
|
|
|
|
|

---------------
Fluid Diffusion
---------------


The fluxes governing fluid diffusion are

Cartesian
---------

.. math::

      \vec{F}_d(Q,\nabla Q) = 
      -
      \begin{pmatrix}
        0 \\ \overline{\overline{\tau}} \\ k \nabla T  +  \overline{\overline{\tau}} \cdot \vec{u}
      \end{pmatrix} 
      \quad\quad
      S_d(\vec{x},Q) = 
      \begin{pmatrix}
          0 \\ 0 \\ 0 \\ 0 \\ 0
      \end{pmatrix}


Cylindrical
-----------

.. math::

      \vec{F}_d(Q,\nabla Q) = 
      -
      \begin{pmatrix}
        0 \\ \overline{\overline{\tau}} \\ k \nabla T  +  \overline{\overline{\tau}} \cdot \vec{u}
      \end{pmatrix} 
      \quad\quad
      S_d(\vec{r},Q) = 
      \begin{pmatrix}
          0 \\ -\frac{\tau_{\theta\theta}}{r} \\ 0 \\ 0 \\ 0
      \end{pmatrix}


|
|
|
|
|
|

------
Models
------






Equations of State
------------------


Ideal Gas
~~~~~~~~~

+-------------------+-------------------------------------------------------------------------------------------+
| **Model name:**   | ``Ideal Gas``                                                                             |
+-------------------+-------------------------------------------------------------------------------------------+
| **Model fields:** | ``Pressure`` ``Temperature``                                                              |
+-------------------+-------------------------------------------------------------------------------------------+

The ideal gas equation of state computes ``Pressure`` and ``Temperature`` as

.. math::

    p &= (\gamma - 1)\bigg(\rho E - \frac{1}{2}\frac{\vec{\rho u} \cdot \vec{\rho u}}{\rho} \bigg)
    \\
    \\
    T &= \frac{p}{\rho R}

where :math:`R` is the specific gas constant

.. math::

    R = 287.15  \quad\quad \bigg[ \frac{\text{J}}{\text{kg} \cdot \text{K}} \bigg]



    




Viscosity
---------


Sutherland's Law
~~~~~~~~~~~~~~~~

+-------------------+-------------------------------------------------------------------------------------------+
| **Model name:**   | ``Sutherlands Law``                                                                       |
+-------------------+-------------------------------------------------------------------------------------------+
| **Model fields:** | ``Laminar Viscosity``                                                                     |
+-------------------+-------------------------------------------------------------------------------------------+

Sutherland's Law computes ``Laminar Viscosity`` as a function of temperature
using

.. math::

    \mu = \mu_0 \bigg(\frac{T}{T_0}\bigg)^{3/2} \frac{T_0 + S}{T + S}


where the model constants are

.. math::

    \mu_0 &= 1.7894e^{-5}    \quad\quad  &\bigg[\frac{\text{kg}}{\text{m} \cdot \text{s}}\bigg] \\
    T_0   &= 273.11          \quad\quad  &\big[K\big] \\
    S     &= 110.56          \quad\quad  &\big[K\big]




Constant Viscosity
~~~~~~~~~~~~~~~~~~


Velocity Gradients
------------------


Gradients of velocity are computed using the chain rule. From the ChiDG framework, we have gradients of the 
primary field variables. Here, we have gradients of the components of 
momentum: :math:`(\nabla \rho u,\nabla \rho v, \nabla \rho w)`. Gradients
of velocity, :math:`(\nabla u, \nabla v, \nabla w)` are computed by recognizing that in general

.. math::

    \nabla (\phi f) = \nabla(\phi) f  +  \phi \nabla (f)



So the gradient of the velocity components can be computed as

.. math::

    \nabla (u) &= \frac{\nabla(\rho u)}{\rho} - \frac{u \nabla(\rho)}{\rho}
    \\
    \nabla (v) &= \frac{\nabla(\rho v)}{\rho} - \frac{v \nabla(\rho)}{\rho}
    \\
    \nabla (w) &= \frac{\nabla(\rho w)}{\rho} - \frac{w \nabla(\rho)}{\rho}


.. note::

    In ``Cylindrical`` coordinates, we have :math:`\nabla(r \rho u_\theta)` instead 
    of :math:`\nabla(\rho u_\theta)`. The gradient of tangential momentum is 
    computed from the angular momentum gradient as

    .. math::

        \nabla(\rho u_\theta) = 
        \begin{bmatrix}
            \frac{\nabla_1(r \rho u_\theta)}{r} - \frac{u_\theta}{r}, &
            \frac{\nabla_2(r \rho u_\theta)}{r}, &
            \frac{\nabla_3(r \rho u_\theta)}{r}
        \end{bmatrix}




Shear Stress
------------

+-------------------+-------------------------------------------------------------------------------------------+
| **Model name:**   | ``Shear Stress``                                                                          |
+-------------------+-------------------------------------------------------------------------------------------+
| **Model fields:** | ``Shear-11``, ``Shear-22``, ``Shear-33``, ``Shear-12``, ``Shear-13``, ``Shear-23``        |
+-------------------+-------------------------------------------------------------------------------------------+

The shear stress tensor is defined as

.. math::

    \overline{\overline{\tau}} = \mu(\nabla \vec{u} + \nabla \vec{u}^T) +  \lambda \overline{\overline{I}} \nabla \cdot \vec{u}

The tensor compnents are

.. math::

    \overline{\overline{\tau}} = 
    \begin{bmatrix}
        \tau_{11} & \tau_{12} & \tau_{13} \\
        \tau_{21} & \tau_{22} & \tau_{23} \\
        \tau_{31} & \tau_{32} & \tau_{33} \\
    \end{bmatrix}

.. note::

    The stress tensor is symmetric. So, only the upper triangular components of the tensor 
    are computed.


The components of the stress tensor are computed as

Cartesian
~~~~~~~~~

.. math::

    \tau_{11} &= 2 \mu \bigg(\nabla_1 u \bigg)  + \lambda(\nabla \cdot \vec{u}) \\
    \tau_{22} &= 2 \mu \bigg(\nabla_2 v \bigg)  + \lambda(\nabla \cdot \vec{u}) \\
    \tau_{33} &= 2 \mu \bigg(\nabla_3 w\bigg)   + \lambda(\nabla \cdot \vec{u}) \\
 \\
    \tau_{12} &= \mu \bigg( \nabla_2 u + \nabla_1 v \bigg) \\
    \tau_{13} &= \mu \bigg( \nabla_3 u + \nabla_1 w \bigg) \\
    \tau_{23} &= \mu \bigg( \nabla_2 w + \nabla_3 v \bigg)

.. math::

    \nabla \cdot \vec{u} = \bigg( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z} \bigg) = \bigg[ \nabla_1 u + \nabla_2 v + \nabla_3 w \bigg]

Cylindrical
~~~~~~~~~~~

.. math::

    \tau_{11} &= 2 \mu \bigg(\nabla_1 u \quad\quad \bigg)     + \lambda(\nabla \cdot \vec{u}) \\
    \tau_{22} &= 2 \mu \bigg(\nabla_2 v + \frac{u}{r} \bigg)  + \lambda(\nabla \cdot \vec{u}) \\
    \tau_{33} &= 2 \mu \bigg(\nabla_3 w \quad\quad \bigg)     + \lambda(\nabla \cdot \vec{u}) \\
 \\
    \tau_{12} &= \mu \bigg( \nabla_2 u + \nabla_1 v - \frac{v}{r} \bigg) \\
    \tau_{13} &= \mu \bigg( \nabla_3 u + \nabla_1 w \quad\quad \bigg) \\
    \tau_{23} &= \mu \bigg( \nabla_2 w + \nabla_3 v \quad\quad \bigg)

.. math::

    \nabla \cdot \vec{u} = \bigg( \frac{1}{r}\frac{\partial r u}{\partial r} + \frac{1}{r}\frac{\partial v}{\partial \theta} + \frac{\partial w}{\partial z}\bigg)  =  \bigg( \frac{\partial u}{\partial r} + \frac{1}{r}\frac{\partial v}{\partial \theta} + \frac{\partial w}{\partial z}  +  \frac{u}{r} \bigg) = \bigg[\nabla_1 u + \nabla_2 v + \nabla_3 w + \frac{u}{r}\bigg]










Temperature Gradient
--------------------

+-------------------+-------------------------------------------------------------------------------------------+
| **Model name:**   | ``Temperature Gradient``                                                                  |
+-------------------+-------------------------------------------------------------------------------------------+
| **Model fields:** | ``Temperature Gradient - 1``  ``Temperature Gradient - 2``  ``Temperature Gradient - 3``  |
+-------------------+-------------------------------------------------------------------------------------------+


Gradients of temperature are computed using the chain rule. From the ChiDG framework, 
we have gradients of the primary field variables. The gradient of the scalar 
temperature field :math:`\nabla T(Q)` can be computed by expanding

.. math::

    \nabla T(Q) = \frac{\partial T}{\partial Q} \nabla Q

as

.. math::

    \nabla T(\rho, \rho u, \rho v, \rho w, \rho E) = 
    \frac{\partial T}{\partial \rho} \nabla \rho + 
    \frac{\partial T}{\partial \rho u} \nabla \rho u + 
    \frac{\partial T}{\partial \rho v} \nabla \rho v + 
    \frac{\partial T}{\partial \rho w} \nabla \rho w + 
    \frac{\partial T}{\partial \rho E} \nabla \rho E


.. note::

    In ``Cylindrical`` coordinates, we have :math:`\nabla(r \rho u_\theta)` instead 
    of :math:`\nabla(\rho u_\theta)`. The gradient of tangential momentum is 
    computed from the angular momentum gradient as

    .. math::

        \nabla(\rho u_\theta) = 
        \begin{bmatrix}
            \frac{\nabla_1(r \rho u_\theta)}{r} - \frac{u_\theta}{r}, &
            \frac{\nabla_2(r \rho u_\theta)}{r}, &
            \frac{\nabla_3(r \rho u_\theta)}{r}
        \end{bmatrix}



