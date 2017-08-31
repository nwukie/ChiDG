module RANS_bc_diffusion
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,THREE,HALF
    use mod_fluid,              only: omega, gam
    use mod_spalart_allmaras

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: RANS_bc_diffusion_t

!        real(rk)    :: gam = 1.4_rk
!        real(rk)    :: R   = 287.15_rk
!        real(rk)    :: Cp  = 1003.0_rk
!        real(rk)    :: Pr  = 0.72_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type RANS_bc_diffusion_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/20/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(RANS_bc_diffusion_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("RANS BC Diffusion")

        ! Set operator type
        call self%set_operator_type("BC Diffusive Flux")

        ! Set operator equations
        call self%add_primary_field('Density'          )
        call self%add_primary_field('Momentum-1'       )
        call self%add_primary_field('Momentum-2'       )
        call self%add_primary_field('Momentum-3'       )
        call self%add_primary_field('Energy'           )
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(RANS_bc_diffusion_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                    &
            density, mom1, mom2, mom3, energy, density_nutilde,     &
            enthalpy, pressure, nu, nutilde,                        &
            grad1_density, grad2_density, grad3_density,            &
            grad1_mom1,    grad2_mom1,    grad3_mom1,               &
            grad1_mom2,    grad2_mom2,    grad3_mom2,               &
            grad1_mom3,    grad2_mom3,    grad3_mom3,               &
            grad1_energy,  grad2_energy,  grad3_energy,             &
            grad1_nutilde, grad2_nutilde, grad3_nutilde,            &
            grad1_density_nutilde, grad2_density_nutilde, grad3_density_nutilde,    &
            grad1_u,       grad2_u,       grad3_u,                                  &
            grad1_v,       grad2_v,       grad3_v,                                  &
            grad1_w,       grad2_w,       grad3_w,                                  &
            grad1_T,       grad2_T,       grad3_T,                                  &
            mu,            mu_l,          mu_t,                                     &
            mu2,           mu2_l,         mu2_t,                                    &
            k,             k_l,           k_t,                                      &
            shear_11,      shear_22,      shear_33,                                 &
            shear_12,      shear_13,      shear_23,                                 &
            dke_ddensity,  dke_dmom1, dke_dmom2, dke_dmom3, dke_denergy,            &
            du_ddensity,   dv_ddensity,   dw_ddensity,                              &
            du_dmom1,      dv_dmom2,      dw_dmom3,                                 &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,                  &
            dT_ddensity, dT_dmom1, dT_dmom2, dT_dmom3, dT_denergy,                  &
            dnutilde_ddensity,                                                      &
            dnutilde_ddensitynutilde,                                               &
            invdensity, div_velocity, u, v, w, chi, diffusion, f_n1, f_v1,          &
            vorticity_1, vorticity_2, vorticity_3,                                  &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r
        real(rk)                                :: const

        ! Sutherlands Law constants
        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]




        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_field('Density'   ,        'value','boundary')
        mom1            = worker%get_field('Momentum-1',        'value','boundary')
        mom2            = worker%get_field('Momentum-2',        'value','boundary')
        mom3            = worker%get_field('Momentum-3',        'value','boundary')
        energy          = worker%get_field('Energy'    ,        'value','boundary')
        density_nutilde = worker%get_field('Density * NuTilde', 'value','boundary')
        invdensity = ONE/density


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density         = worker%get_field('Density'   ,       'grad1','boundary')
        grad2_density         = worker%get_field('Density'   ,       'grad2','boundary')
        grad3_density         = worker%get_field('Density'   ,       'grad3','boundary')

        grad1_mom1            = worker%get_field('Momentum-1',       'grad1','boundary')
        grad2_mom1            = worker%get_field('Momentum-1',       'grad2','boundary')
        grad3_mom1            = worker%get_field('Momentum-1',       'grad3','boundary')

        grad1_mom2            = worker%get_field('Momentum-2',       'grad1','boundary')
        grad2_mom2            = worker%get_field('Momentum-2',       'grad2','boundary')
        grad3_mom2            = worker%get_field('Momentum-2',       'grad3','boundary')

        grad1_mom3            = worker%get_field('Momentum-3',       'grad1','boundary')
        grad2_mom3            = worker%get_field('Momentum-3',       'grad2','boundary')
        grad3_mom3            = worker%get_field('Momentum-3',       'grad3','boundary')

        grad1_energy          = worker%get_field('Energy'    ,       'grad1','boundary')
        grad2_energy          = worker%get_field('Energy'    ,       'grad2','boundary')
        grad3_energy          = worker%get_field('Energy'    ,       'grad3','boundary')

        grad1_density_nutilde = worker%get_field('Density * NuTilde','grad1','boundary')
        grad2_density_nutilde = worker%get_field('Density * NuTilde','grad2','boundary')
        grad3_density_nutilde = worker%get_field('Density * NuTilde','grad3','boundary')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)

        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if





        !
        ! Compute Pressure, Enthalpy
        !
        !pressure    = (self%gam-ONE)*(energy - HALF*( (mom1*mom1) + (mom2*mom2) + (mom3*mom3) )*invdensity )
        !temperature = pressure/(density*self%R)
        pressure = worker%get_field('Pressure', 'value', 'boundary')


        ! Compute viscosity using Sutherland's Law
        !mu_l  = mu0*((temperature/T0)**(THREE/TWO))*(T0+S)/(temperature+S)
        mu_l  = worker%get_field('Laminar Viscosity', 'value', 'boundary')
        mu2_l = -(TWO/THREE)*mu_l
        k_l   = self%Cp * mu_l / self%Pr



        !
        ! Compute turbulence fields
        !
        nu      = mu_l*invdensity
        nutilde = density_nutilde*invdensity



        !
        ! Compute chi, f_v1
        !
        chi = nutilde/nu
        f_v1 = chi*chi*chi/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)


        !
        ! Compute f_n1
        !
        f_n1 = density
        where (nutilde >= ZERO)
            f_n1 = ONE
        else where
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where




        !
        ! Initialize derivatives, compute mu_t
        !
        mu_t = density
        where (nutilde >= 0)
            mu_t = density * nutilde * f_v1
        else where
            mu_t = ZERO
        end where


        !
        ! Compute: 
        !   - Second Coefficient of Turbulent Viscosity, Stokes' Hypothesis.
        !   - Turbulent Thermal Conductivity, Reynolds' analogy.
        !
        mu2_t = (-TWO/THREE)*mu_t
        k_t   = self%Cp*mu_t/SA_Pr_t




        !
        ! Compute total coefficients
        !
        mu  = mu_l  + mu_t
        mu2 = mu2_l + mu2_t
        k   = k_l   + k_t




        !
        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        !
        dnutilde_ddensity        = -density_nutilde*invdensity*invdensity
        dnutilde_ddensitynutilde =  invdensity


        grad1_nutilde = dnutilde_ddensity * grad1_density  +  dnutilde_ddensitynutilde * grad1_density_nutilde
        grad2_nutilde = dnutilde_ddensity * grad2_density  +  dnutilde_ddensitynutilde * grad2_density_nutilde
        grad3_nutilde = dnutilde_ddensity * grad3_density  +  dnutilde_ddensitynutilde * grad3_density_nutilde




        !############### Computing Temperature Gradient ###################

        !
        ! Compute temperature gradient
        !

        ! Compute velocities
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity

        ! Compute Kinetic Energy Jacobians
        dke_ddensity = -HALF*(u*u + v*v + w*w)
        dke_dmom1    = u
        dke_dmom2    = v
        dke_dmom3    = w

        ! Compute Pressure Jacobians
        dp_ddensity = -(self%gam-ONE)*dke_ddensity
        dp_dmom1    = -(self%gam-ONE)*dke_dmom1
        dp_dmom2    = -(self%gam-ONE)*dke_dmom2
        dp_dmom3    = -(self%gam-ONE)*dke_dmom3
        dp_denergy  =  dp_dmom3    ! Initialize derivatives
        dp_denergy  =  (self%gam-ONE)   ! No negative sign

        ! Compute Temperature Jacobians
        const = ONE/self%R
        dT_ddensity = const*invdensity*dp_ddensity  -  const*invdensity*invdensity*pressure
        dT_dmom1    = const*invdensity*dp_dmom1
        dT_dmom2    = const*invdensity*dp_dmom2
        dT_dmom3    = const*invdensity*dp_dmom3
        dT_denergy  = const*invdensity*dp_denergy

        ! Temperature gradient
        grad1_T = dT_ddensity*grad1_density + dT_dmom1*grad1_mom1 + dT_dmom2*grad1_mom2 + dT_dmom3*grad1_mom3 + dT_denergy*grad1_energy
        grad2_T = dT_ddensity*grad2_density + dT_dmom1*grad2_mom1 + dT_dmom2*grad2_mom2 + dT_dmom3*grad2_mom3 + dT_denergy*grad2_energy
        grad3_T = dT_ddensity*grad3_density + dT_dmom1*grad3_mom1 + dT_dmom2*grad3_mom2 + dT_dmom3*grad3_mom3 + dT_denergy*grad3_energy





        !############### Computing Vorticity ######################

        !
        ! compute velocity jacobians
        !
        du_ddensity = -invdensity*invdensity*mom1
        dv_ddensity = -invdensity*invdensity*mom2
        dw_ddensity = -invdensity*invdensity*mom3

        du_dmom1 = invdensity
        dv_dmom2 = invdensity
        dw_dmom3 = invdensity



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
        grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1
        grad3_u = du_ddensity*grad3_density  +  du_dmom1*grad3_mom1

        grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
        grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2
        grad3_v = dv_ddensity*grad3_density  +  dv_dmom2*grad3_mom2

        grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
        grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
        grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3




        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute vorticity:
            !
            !   vorticity = Curl(U)
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   curl(U) = (dwdx - dvdz)i  +  (dudz - dwdx)j  +  (dvdx - dudy)k
            !
            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u) 



        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   curl(U) = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  (1/r)(d(rv)dr - dudtheta)k
            !           = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  ( dvdr - (1/r)dudtheta  + (v/r) )k
            !           = (grad2_w - grad3_v)i  +  (grad3_u - grad1_w)j  +  ( grad1_v - grad2_u  + (v/r) )k
            !
            !   Note:
            !       grad1_ = d/dr
            !       grad2_ = (1/r)d/dtheta
            !       grad3_ = d/dz
            !
            v = mom2*invdensity

            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u + (v/r)) 


            !
            ! Account for rotation, convert to relative vorticity
            !
            vorticity_3 = vorticity_3 - TWO*omega

        end if










        !################## Computing Shear Stresses #####################

        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   div(V) = dudx + dvdy + dwdz
            !
            div_velocity = grad1_u + grad2_v + grad3_w

            !
            ! Compute shear stress components
            !
            shear_11 = TWO*mu*grad1_u  +  mu2*(div_velocity)
            shear_22 = TWO*mu*grad2_v  +  mu2*(div_velocity)
            shear_33 = TWO*mu*grad3_w  +  mu2*(div_velocity)

            shear_12 = mu*(grad2_u + grad1_v)
            shear_13 = mu*(grad3_u + grad1_w)
            shear_23 = mu*(grad2_w + grad3_v)




        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   div(U) = (1/r)d(r*u)dr   + (1/r)dv/dtheta + dw/dz
            !          = (du/dr + u/r) + (1/r)dv/dtheta + dw/dz
            !
            !   Note:
            !       grad1_u = du/dr
            !       grad2_v = (1/r)dv/dtheta
            !       grad3_w = dw/dz
            !
            u = mom1*invdensity
            v = mom2*invdensity
            div_velocity = (grad1_u + u/r) + grad2_v + grad3_w

            !
            ! Compute shear stress components
            !
            shear_11 = TWO*mu*(grad1_u        )  +  mu2*(div_velocity)
            shear_22 = TWO*mu*(grad2_v + (u/r))  +  mu2*(div_velocity)
            shear_33 = TWO*mu*(grad3_w        )  +  mu2*(div_velocity)

            shear_12 = mu*(grad2_u + grad1_v  - (v/r))
            shear_13 = mu*(grad3_u + grad1_w         )
            shear_23 = mu*(grad2_w + grad3_v         )

        end if





        !=================================================
        ! mass flux
        !=================================================


        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = -shear_11
        flux_2 = -shear_12
        flux_3 = -shear_13

        call worker%integrate_boundary_condition('Momentum-1','Diffusion',flux_1,flux_2,flux_3)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = -shear_12
        flux_2 = -shear_22
        flux_3 = -shear_23

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        call worker%integrate_boundary_condition('Momentum-2','Diffusion',flux_1,flux_2,flux_3)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = -shear_13
        flux_2 = -shear_23
        flux_3 = -shear_33

        call worker%integrate_boundary_condition('Momentum-3','Diffusion',flux_1,flux_2,flux_3)
                             
        !=================================================
        ! energy flux
        !=================================================
        flux_1 = -k*grad1_T  -  (u*shear_11 + v*shear_12 + w*shear_13)
        flux_2 = -k*grad2_T  -  (u*shear_12 + v*shear_22 + w*shear_23)
        flux_3 = -k*grad3_T  -  (u*shear_13 + v*shear_23 + w*shear_33)

        call worker%integrate_boundary_condition('Energy','Diffusion',flux_1,flux_2,flux_3)


        !================================
        !       TURBULENCE FLUX
        !================================
        diffusion = -(ONE/SA_sigma)*(mu_l + f_n1*density_nutilde)

        flux_1 = (diffusion*grad1_nutilde)
        flux_2 = (diffusion*grad2_nutilde)
        flux_3 = (diffusion*grad3_nutilde)

        call worker%integrate_boundary_condition('Density * NuTilde','Diffusion',flux_1,flux_2,flux_3)

    

    end subroutine compute
    !*******************************************************************************




end module RANS_bc_diffusion
