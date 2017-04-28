===============
Fluid Equations
===============

Primary field variables for the fluid equations in ChiDG are ``Density``, 
``Momentum-1``, ``Momentum-2``, ``Momentum-3``, and ``Energy``


.. math::

      Q_{cart} = 
      \begin{pmatrix}
        \rho \\ \rho v_x \\ \rho v_y \\ \rho v_z \\ \rho E 
      \end{pmatrix}
      \quad\quad
      Q_{cyl} = 
      \begin{pmatrix}
        \rho \\ \rho v_r \\ r \rho v_\theta \\ \rho v_z \\ \rho E 
      \end{pmatrix}



.. note::

    In ``Cylindrical`` coordinates, angular momentum, :math:`r \rho v_\theta`, 
    is being conserved. Not tangential momentum, :math:`\rho v_\theta`.


The transformation of the momentum vector between ``Cartesian`` and 
``Cylindrical`` coordinates is

.. math::

    \begin{bmatrix}
        \rho v_x \\ \rho v_y \\ \rho v_z
    \end{bmatrix}
    =
    \begin{bmatrix}
        cos(\theta) &   -sin(\theta)  &  0 \\
        sin(\theta) &    cos(\theta)  &  0 \\
        0 & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
        \rho v_r \\ \rho v_\theta \\ \rho v_z
    \end{bmatrix}

The velocity vector in ``Cartesian`` and ``Cylindrical`` coordinates is

.. math::

    \vec{u} = [v_x, v_y, v_z]       \quad\quad  \vec{v} = [v_r, v_\theta, v_z]

Some distinctions are made here between absolute, relative, and frame velocities.
The absolute velocity(:math:`v`) is that which is perceived in the inertial frame of 
reference. The relative(or advection)(:math:`w`) velocity is the velocity, which is perceived
in the non-inertial frame of reference. The frame velocity(:math:`u`) is the velocity of 
the moving coordinate system. These are defined in the Cartesian and Cylindrical
coordinate systems as

.. math::

    \vec{v} &= [v_x, v_y, v_z]   \quad\quad\:\:\:\:  \vec{v} = [v_r, v_\theta, v_y] \\
    \vec{w} &= [w_x, w_y, w_z]   \quad\quad          \vec{w} = [w_r, w_\theta, w_y] \\
    \vec{u} &= [u_x, u_y, u_z]   \quad\quad\:\:\:    \vec{u} = [u_r, u_\theta, u_y]

For calculations using an inertial frame, the frame velocity is zero.

.. math::

    \vec{u} = [0, 0, 0]

For a rotating Cylindrical coordinate system, the frame velocity is

.. math::

    \vec{u}_{cyl} = [0, \omega r, 0]

where :math:`\omega` is the rotational rate in :math:`rad/s`. The relation
between absolute, relative, and frame velocities is

.. math::

    \vec{v} = \vec{w} + \vec{u}

At times, components of the velocity vectors are represented in the same manner 
for both coordinate systems

.. math::

    \vec{v}  = [v_1,v_2,v_3] 


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


The fluxes representing advective transport for the fluid equations
are given compactly as

.. math::

    \vec{F}_a(Q) = 
    \begin{pmatrix}
        \rho \vec{w} \\
        \mathcal{S} (\rho \vec{v} \otimes \vec{w} + \overline{\overline{I}}p) \\
        \rho H \vec{w} + \vec{u} p
    \end{pmatrix}
    \quad\quad
    \mathcal{S}_{cart} = 
    \begin{pmatrix}
        1   \\  1   \\  1
    \end{pmatrix}
    \quad\quad
    \mathcal{S}_{cyl} = 
    \begin{pmatrix}
        1   \\  r   \\  1
    \end{pmatrix}

where :math:`\mathcal{S}` is a scaling vector that is sued to provide the correct
dimensional scaling for the angular momentum equation. Expanded forms for the 
fluxes in Cartesian and Cylindrical coordinates are given as follows.


Cartesian
---------

.. math:: 

    \vec{F}_a(\vec{x},Q) =
    \begin{pmatrix}
    \begin{bmatrix}
        \rho     w_x       \\ 
        \rho v_x w_x  +  p \\
        \rho v_y w_x       \\
        \rho v_z w_x       \\
        \rho H   w_x  + u_x p
    \end{bmatrix}
    ,
    \begin{bmatrix}
        \rho     w_y        \\ 
        \rho v_x w_y        \\
        \rho v_y w_y  +  p  \\
        \rho v_z w_y        \\
        \rho H   w_y  + u_y p
    \end{bmatrix}
    ,
    \begin{bmatrix}
        \rho     w_z       \\ 
        \rho v_x w_z        \\
        \rho v_y w_z        \\
        \rho v_z w_z  +  p  \\
        \rho H   w_z  +  u_z p
    \end{bmatrix}
    \end{pmatrix}
    \quad\quad
    S_a(\vec{x},Q) = 
    \begin{pmatrix}
        0 \\ 0 \\ 0 \\ 0 \\ 0
    \end{pmatrix}


Cylindrical
-----------

.. math:: 

      \vec{F}_a(\vec{r},Q) =
      \begin{pmatrix}
      \begin{bmatrix}
          \rho          w_r              \\ 
          \rho v_r      w_r  +  p        \\
        r \rho v_\theta w_r              \\
          \rho v_z      w_r              \\
          \rho H        w_r  +  u_r p
      \end{bmatrix}
      ,
      \begin{bmatrix}
          \rho          w_\theta         \\ 
          \rho v_r      w_\theta         \\
        r \rho v_\theta w_\theta  +  r p \\
          \rho v_z      w_\theta         \\
          \rho H        w_\theta  + u_\theta p
      \end{bmatrix}
      ,
      \begin{bmatrix}
          \rho          w_z              \\ 
          \rho v_r      w_z              \\
        r \rho v_\theta w_z              \\
          \rho v_z      w_z  +   p       \\
          \rho H        w_z  +  u_z p
      \end{bmatrix}
      \end{pmatrix}
      \quad\quad
      S_a(\vec{r},Q) = 
      \begin{pmatrix}
          0 \\ \frac{\rho v_\theta w_\theta - p}{r} \\ 0 \\ 0 \\ 0
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


The fluxes representing diffusive transport for the fluid equations are given compactly as

.. math::

      \vec{F}_d(Q,\nabla Q) = 
      -
      \begin{pmatrix}
        0 \\ \mathcal{S} \overline{\overline{\tau}} \\ k \nabla T  +  \overline{\overline{\tau}} \cdot \vec{v}
      \end{pmatrix} 
      \quad\quad
      \mathcal{S}_{cart} = 
      \begin{pmatrix}
        1   \\  1   \\  1
      \end{pmatrix}
      \quad\quad
      \mathcal{S}_{cyl} = 
      \begin{pmatrix}
        1   \\  r   \\  1
      \end{pmatrix}

where :math:`\mathcal{S}` is a scaling vector that is used to provide the correct dimensional
scaling for the angular momentum equation. Expanded forms for the fluxes in Cartesian 
and Cylindrical coordinates are given as follows.


Cartesian
---------

.. math::

      \vec{F}_d(\vec{x},Q,\nabla Q) =
      -
      \begin{pmatrix}
      \begin{bmatrix}
                   0               \\ 
                  \tau_{11}        \\
                  \tau_{21}        \\
                  \tau_{31}        \\
             k \nabla_1 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      ,
      \begin{bmatrix}
                   0                    \\ 
                  \tau_{12}             \\
                  \tau_{22}             \\
                  \tau_{32}             \\
             k \nabla_2 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      ,
      \begin{bmatrix}
                   0                \\ 
                  \tau_{13}         \\
                  \tau_{23}         \\
                  \tau_{33}         \\
             k \nabla_3 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      \end{pmatrix}
      \quad
      S_d(\vec{x},Q) = 
      \begin{pmatrix}
        0 \\ 0 \\ 0 \\ 0 \\ 0
      \end{pmatrix}



Cylindrical
-----------


.. math:: 

      \vec{F}_d(\vec{r},Q,\nabla Q) =
      -
      \begin{pmatrix}
      \begin{bmatrix}
                   0               \\ 
                  \tau_{11}        \\
                 r\tau_{21}        \\
                  \tau_{31}        \\
             k \nabla_1 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      ,
      \begin{bmatrix}
                   0                    \\ 
                  \tau_{12}             \\
                 r\tau_{22}             \\
                  \tau_{32}             \\
             k \nabla_2 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      ,
      \begin{bmatrix}
                   0                \\ 
                  \tau_{13}         \\
                 r\tau_{23}         \\
                  \tau_{33}         \\
             k \nabla_3 T  + \overline{\overline{\tau}} \cdot \vec{v}
      \end{bmatrix}
      \end{pmatrix}
      \quad
      S_d(\vec{r},Q) = 
      \begin{pmatrix}
        0 \\ -\frac{\tau_{22}}{r} \\ 0 \\ 0 \\ 0
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

    p &= (\gamma - 1)\bigg(\rho E - \frac{1}{2}\frac{\vec{\rho v} \cdot \vec{\rho v}}{\rho} \bigg)
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
| **Model fields:** | ``Dynamic Viscosity``                                                                     |
+-------------------+-------------------------------------------------------------------------------------------+

Sutherland's Law computes ``Dynamic Viscosity`` as a function of temperature
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
momentum: :math:`(\nabla \rho v_1,\nabla \rho v_2, \nabla \rho v_3)`. Gradients
of velocity, :math:`(\nabla v_1, \nabla v_2, \nabla v_3)` are computed by recognizing that in general

.. math::

    \nabla (\phi f) = \nabla(\phi) f  +  \phi \nabla (f)



So the gradient of the velocity components can be computed as

.. math::

    \nabla (v_1) &= \frac{\nabla(\rho v_1)}{\rho} - \frac{u \nabla(\rho)}{\rho}
    \\
    \nabla (v_2) &= \frac{\nabla(\rho v_2)}{\rho} - \frac{v \nabla(\rho)}{\rho}
    \\
    \nabla (v_3) &= \frac{\nabla(\rho v_3)}{\rho} - \frac{w \nabla(\rho)}{\rho}


.. note::

    In ``Cylindrical`` coordinates, we have :math:`\nabla(r \rho v_\theta)` instead 
    of :math:`\nabla(\rho v_\theta)`. The gradient of tangential momentum is 
    computed from the angular momentum gradient as

    .. math::

        \nabla(\rho v_\theta) = 
        \begin{bmatrix}
            \frac{\nabla_1(r \rho v_\theta)}{r} - \frac{v_\theta}{r}, &
            \frac{\nabla_2(r \rho v_\theta)}{r}, &
            \frac{\nabla_3(r \rho v_\theta)}{r}
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

    \overline{\overline{\tau}} = \mu(\nabla \vec{v} + \nabla \vec{v}^T) +  \lambda \overline{\overline{I}} \nabla \cdot \vec{v}

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

    \tau_{11} &= 2 \mu \bigg(\nabla_1 v_x \bigg)  + \lambda(\nabla \cdot \vec{v}) \\
    \tau_{22} &= 2 \mu \bigg(\nabla_2 v_y \bigg)  + \lambda(\nabla \cdot \vec{v}) \\
    \tau_{33} &= 2 \mu \bigg(\nabla_3 v_z \bigg)  + \lambda(\nabla \cdot \vec{v}) \\
 \\
    \tau_{12} &= \mu \bigg( \nabla_2 v_x + \nabla_1 v_y \bigg) \\
    \tau_{13} &= \mu \bigg( \nabla_3 v_x + \nabla_1 v_z \bigg) \\
    \tau_{23} &= \mu \bigg( \nabla_2 v_z + \nabla_3 v_y \bigg)

.. math::

    \nabla \cdot \vec{v} = \bigg( \frac{\partial v_x}{\partial x} + \frac{\partial v_y}{\partial y} + \frac{\partial v_z}{\partial z} \bigg) = \bigg[ \nabla_1 v_x + \nabla_2 v_y + \nabla_3 v_z \bigg]

Cylindrical
~~~~~~~~~~~

.. math::

    \tau_{11} &= 2 \mu \bigg(\nabla_1 v_r \quad\quad \bigg)       + \lambda(\nabla \cdot \vec{v}) \\
    \tau_{22} &= 2 \mu \bigg(\nabla_2 v_\theta + \frac{v_r}{r} \bigg)  + \lambda(\nabla \cdot \vec{v}) \\
    \tau_{33} &= 2 \mu \bigg(\nabla_3 v_z \quad\quad \bigg)       + \lambda(\nabla \cdot \vec{v}) \\
 \\
    \tau_{12} &= \mu \bigg( \nabla_2 v_r + \nabla_1 v_\theta - \frac{v_\theta}{r} \bigg) \\
    \tau_{13} &= \mu \bigg( \nabla_3 v_r + \nabla_1 v_z \quad\quad \bigg) \\
    \tau_{23} &= \mu \bigg( \nabla_2 v_z + \nabla_3 v_\theta \quad\quad \bigg)

.. math::

    \nabla \cdot \vec{v} = \bigg( \frac{1}{r}\frac{\partial r v_r}{\partial r} + \frac{1}{r}\frac{\partial v_\theta}{\partial \theta} + \frac{\partial v_z}{\partial z}\bigg)  =  \bigg( \frac{\partial v_r}{\partial r} + \frac{1}{r}\frac{\partial v_\theta}{\partial \theta} + \frac{\partial v_z}{\partial z}  +  \frac{v_r}{r} \bigg) = \bigg[\nabla_1 v_r + \nabla_2 v_\theta + \nabla_3 v_z + \frac{v_r}{r}\bigg]










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

    \nabla T(\rho, \rho v_1, \rho v_2, \rho v_3, \rho E) = 
    \frac{\partial T}{\partial \rho}     \nabla \rho + 
    \frac{\partial T}{\partial \rho v_1} \nabla \rho v_1 + 
    \frac{\partial T}{\partial \rho v_2} \nabla \rho v_2 + 
    \frac{\partial T}{\partial \rho v_3} \nabla \rho v_3 + 
    \frac{\partial T}{\partial \rho E}   \nabla \rho E


.. note::

    In ``Cylindrical`` coordinates, we have :math:`\nabla(r \rho v_\theta)` instead 
    of :math:`\nabla(\rho v_\theta)`. The gradient of tangential momentum is 
    computed from the angular momentum gradient as

    .. math::

        \nabla(\rho v_\theta) = 
        \begin{bmatrix}
            \frac{\nabla_1(r \rho v_\theta)}{r} - \frac{v_\theta}{r}, &
            \frac{\nabla_2(r \rho v_\theta)}{r}, &
            \frac{\nabla_3(r \rho v_\theta)}{r}
        \end{bmatrix}





Vorticity
---------

+-------------------+-------------------------------------------------------------------------------------------+
| **Model name:**   | ``Vorticity``                                                                             |
+-------------------+-------------------------------------------------------------------------------------------+
| **Model fields:** | ``Vorticity-1``  ``Vorticity-2``  ``Vorticity-3``                                         |
+-------------------+-------------------------------------------------------------------------------------------+


Vorticity, used in the Spalart-Allmaras turbulence model, is defined as the Curl of 
velocity as

.. math::

    \vec{\omega} = \nabla \times \vec{v}



In Cartesian coordinates, this is computed as

.. math::

    \vec{\omega} = \nabla \times \vec{v} &= \bigg(\frac{\partial v_z}{\partial x} - \frac{\partial v_y}{\partial z}\bigg) \hat{x}  +  \bigg(\frac{\partial v_x}{\partial z} - \frac{\partial v_z}{\partial x}\bigg) \hat{y}  +  \bigg( \frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}\bigg) \hat{z} \\
    &= \bigg[ \bigg(\nabla_2 v_3  -  \nabla_3 v_2\bigg), \bigg(\nabla_3 v_1 - \nabla_1 v_3\bigg), \bigg(\nabla_1 v_2 - \nabla_2 v_1\bigg)\bigg]


In Cylindrical coordinates, the Curl of a vector is given as


.. math::

    \vec{\omega} = \nabla \times \vec{v} &= \bigg(\frac{1}{r}\frac{\partial v_z}{\partial \theta} - \frac{\partial v_\theta}{\partial z}\bigg) \hat{r}  +  \bigg(\frac{\partial v_r}{\partial z} - \frac{\partial v_z}{\partial r}\bigg) \hat{\theta}  +  \frac{1}{r}\bigg( \frac{\partial r v_\theta}{\partial r} - \frac{\partial v_r}{\partial \theta}\bigg) \hat{z} \\
    &= \bigg(\frac{1}{r}\frac{\partial v_z}{\partial \theta} - \frac{\partial v_\theta}{\partial z}\bigg) \hat{r}  +  \bigg(\frac{\partial v_r}{\partial z} - \frac{\partial v_z}{\partial r}\bigg) \hat{\theta}  +  \bigg( \frac{\partial v_\theta}{\partial r} - \frac{1}{r}\frac{\partial v_r}{\partial \theta} + \frac{v_\theta}{r}\bigg) \hat{z} \\
    &= \bigg[ \bigg(\nabla_2 v_3  -  \nabla_3 v_2\bigg), \bigg(\nabla_3 v_1 - \nabla_1 v_3\bigg), \bigg(\nabla_1 v_2 - \nabla_2 v_1 + \frac{v_2}{r}\bigg)\bigg]


In the non-inertial frame for Cylindrical coordinates, the relative 
vorticity is accounted for as

.. math::

    \omega_3 = \omega_3 - 2\omega

Note, that here :math:`\omega_3` is the third component of the vorticity 
vector, while :math:`\omega` is the rate of rotation for the non-inertial frame.










