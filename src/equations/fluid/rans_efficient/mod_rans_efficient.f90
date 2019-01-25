module mod_rans_efficient
#include <messenger.h>
    use mod_constants,          only: ZERO, ONE, TWO, THREE, HALF
    use mod_fluid,              only: Rgas, Pr, cp, gam
    use mod_spalart_allmaras,   only: SA_c_v1, SA_Pr_t
    use DNAD_D
    implicit none

    character(1000) :: turbulence_model = 'Spalart Allmaras'
    character(1000) :: viscosity_model  = 'Sutherlands Law'
    character(1000) :: state_equation   = 'Ideal Gas'

    namelist /fluid/ turbulence_model, viscosity_model, state_equation


contains


    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------
    subroutine compute_pressure_temperature(density,mom1,mom2,mom3,energy,pressure,temperature)
        type(AD_D),                 intent(in)      :: density(:)
        type(AD_D),                 intent(in)      :: mom1(:)
        type(AD_D),                 intent(in)      :: mom2(:)
        type(AD_D),                 intent(in)      :: mom3(:)
        type(AD_D),                 intent(in)      :: energy(:)
        type(AD_D), allocatable,    intent(inout)   :: pressure(:)
        type(AD_D), allocatable,    intent(inout)   :: temperature(:)

        pressure = (gam-ONE)*(energy - HALF*( (mom1*mom1) + (mom2*mom2) + (mom3*mom3) )/density )
        temperature = pressure/(density*Rgas)

    end subroutine compute_pressure_temperature
    !*****************************************************************










    !>
    !!
    !!
    !!
    !------------------------------------------------------------------
    subroutine compute_viscosity(density,temperature,density_nutilde,mu_l,mu_t,lambda_l,lambda_t)
        type(AD_D),                 intent(in)      :: density(:)
        type(AD_D),                 intent(in)      :: temperature(:)
        type(AD_D),                 intent(in)      :: density_nutilde(:)
        type(AD_D), allocatable,    intent(inout)   :: mu_l(:)
        type(AD_D), allocatable,    intent(inout)   :: mu_t(:)
        type(AD_D), allocatable,    intent(inout)   :: lambda_l(:)
        type(AD_D), allocatable,    intent(inout)   :: lambda_t(:)

        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]

        type(AD_D), allocatable :: chi(:), f_v1(:), nutilde(:), nu_l(:), chi3(:)


        !
        ! Sutherlands Law for Laminar Viscosity
        !
        mu_l = mu0*((temperature/T0)**(THREE/TWO))*(T0+S)/(temperature+S)
        lambda_l = -(TWO/THREE)*mu_l





        !
        ! Get viscosity: compute nu, nutilde
        !
        nu_l   = mu_l/density
        nutilde = density_nutilde/density


        !
        ! Compute chi, f_v1
        !
        chi = nutilde/nu_l
        chi3 = chi*chi*chi
        f_v1 = chi3/(chi3 + SA_c_v1*SA_c_v1*SA_c_v1)


        !
        ! Initialize derivatives, compute mu_t
        !
        mu_t = density
        where (nutilde >= 0)
            mu_t = density * nutilde * f_v1
        else where
            mu_t = ZERO
        end where
        lambda_t = (-TWO/THREE)*mu_t


    end subroutine compute_viscosity
    !*********************************************************************







    !>
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------
    subroutine compute_thermal_conductivity(mu_l,mu_t,k_l,k_t)
        type(AD_D),                 intent(in)      :: mu_l(:)
        type(AD_D),                 intent(in)      :: mu_t(:)
        type(AD_D), allocatable,    intent(inout)   :: k_l(:)
        type(AD_D), allocatable,    intent(inout)   :: k_t(:)

        k_l = cp * mu_l / Pr
        k_t = cp * mu_t / SA_Pr_t

    end subroutine compute_thermal_conductivity
    !**********************************************************************





    

    !>
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------
    subroutine compute_temperature_gradient(density,u,v,w,p,                            &
                                            grad1_density,grad2_density,grad3_density,  &
                                            grad1_mom1,   grad2_mom1,   grad3_mom1,     &
                                            grad1_mom2,   grad2_mom2,   grad3_mom2,     &
                                            grad1_mom3,   grad2_mom3,   grad3_mom3,     &
                                            grad1_energy, grad2_energy, grad3_energy,   &
                                            grad1_T, grad2_T, grad3_T)
        type(AD_D),                 intent(in)  :: density(:)
        type(AD_D),                 intent(in)  :: u(:)
        type(AD_D),                 intent(in)  :: v(:)
        type(AD_D),                 intent(in)  :: w(:)
        type(AD_D),                 intent(in)  :: p(:)
        type(AD_D),                 intent(in)  :: grad1_density(:)
        type(AD_D),                 intent(in)  :: grad2_density(:)
        type(AD_D),                 intent(in)  :: grad3_density(:)
        type(AD_D),                 intent(in)  :: grad1_mom1(:)
        type(AD_D),                 intent(in)  :: grad2_mom1(:)
        type(AD_D),                 intent(in)  :: grad3_mom1(:)
        type(AD_D),                 intent(in)  :: grad1_mom2(:)
        type(AD_D),                 intent(in)  :: grad2_mom2(:)
        type(AD_D),                 intent(in)  :: grad3_mom2(:)
        type(AD_D),                 intent(in)  :: grad1_mom3(:)
        type(AD_D),                 intent(in)  :: grad2_mom3(:)
        type(AD_D),                 intent(in)  :: grad3_mom3(:)
        type(AD_D),                 intent(in)  :: grad1_energy(:)
        type(AD_D),                 intent(in)  :: grad2_energy(:)
        type(AD_D),                 intent(in)  :: grad3_energy(:)
        type(AD_D), allocatable,    intent(inout)  :: grad1_T(:)
        type(AD_D), allocatable,    intent(inout)  :: grad2_T(:)
        type(AD_D), allocatable,    intent(inout)  :: grad3_T(:)

        type(AD_D), allocatable, dimension(:)   ::                  &
            invdensity,                                             &
            dke_ddensity, dke_dmom1, dke_dmom2, dke_dmom3,          &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,  &
            dT_ddensity, dT_dmom1, dT_dmom2, dT_dmom3, dT_denergy

        real(rk)    :: const


        ! Precompute density inverse
        invdensity = ONE/density


        ! Compute Kinetic Energy Jacobians
        dke_ddensity = -HALF*(u*u + v*v + w*w)
        dke_dmom1    = u
        dke_dmom2    = v
        dke_dmom3    = w


        ! Compute Pressure Jacobians
        dp_ddensity = -(gam-ONE)*dke_ddensity
        dp_dmom1    = -(gam-ONE)*dke_dmom1
        dp_dmom2    = -(gam-ONE)*dke_dmom2
        dp_dmom3    = -(gam-ONE)*dke_dmom3
        dp_denergy  =  dp_dmom3    ! Initialize derivatives
        dp_denergy  =  (gam-ONE)   ! No negative sign


        ! Compute Temperature Jacobians
        const = ONE/Rgas
        dT_ddensity = const*invdensity*dp_ddensity  -  const*invdensity*invdensity*p
        dT_dmom1    = const*invdensity*dp_dmom1
        dT_dmom2    = const*invdensity*dp_dmom2
        dT_dmom3    = const*invdensity*dp_dmom3
        dT_denergy  = const*invdensity*dp_denergy


        ! Compute temperature gradient
        grad1_T = dT_ddensity*grad1_density + dT_dmom1*grad1_mom1 + dT_dmom2*grad1_mom2 + dT_dmom3*grad1_mom3 + dT_denergy*grad1_energy
        grad2_T = dT_ddensity*grad2_density + dT_dmom1*grad2_mom1 + dT_dmom2*grad2_mom2 + dT_dmom3*grad2_mom3 + dT_denergy*grad2_energy
        grad3_T = dT_ddensity*grad3_density + dT_dmom1*grad3_mom1 + dT_dmom2*grad3_mom2 + dT_dmom3*grad3_mom3 + dT_denergy*grad3_energy


    end subroutine compute_temperature_gradient
    !***********************************************************************









    !>
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------
    subroutine compute_shear_stress(coordinate_system,r,                                &
                                    density,grad1_density,grad2_density,grad3_density,  &
                                    mom1,   grad1_mom1,   grad2_mom1,   grad3_mom1,     &
                                    mom2,   grad1_mom2,   grad2_mom2,   grad3_mom2,     &
                                    mom3,   grad1_mom3,   grad2_mom3,   grad3_mom3,     &
                                    energy, grad1_energy, grad2_energy, grad3_energy,   &
                                    mu, lambda,                                         &
                                    shear_11,shear_22,shear_33,shear_12,shear_13,shear_23)
        character(*),               intent(in)      :: coordinate_system
        real(rk),                   intent(in)      :: r(:)
        type(AD_D),                 intent(in)      :: density(:)
        type(AD_D),                 intent(in)      :: grad1_density(:)
        type(AD_D),                 intent(in)      :: grad2_density(:)
        type(AD_D),                 intent(in)      :: grad3_density(:)
        type(AD_D),                 intent(in)      :: mom1(:)
        type(AD_D),                 intent(in)      :: grad1_mom1(:)
        type(AD_D),                 intent(in)      :: grad2_mom1(:)
        type(AD_D),                 intent(in)      :: grad3_mom1(:)
        type(AD_D),                 intent(in)      :: mom2(:)
        type(AD_D),                 intent(in)      :: grad1_mom2(:)
        type(AD_D),                 intent(in)      :: grad2_mom2(:)
        type(AD_D),                 intent(in)      :: grad3_mom2(:)
        type(AD_D),                 intent(in)      :: mom3(:)
        type(AD_D),                 intent(in)      :: grad1_mom3(:)
        type(AD_D),                 intent(in)      :: grad2_mom3(:)
        type(AD_D),                 intent(in)      :: grad3_mom3(:)
        type(AD_D),                 intent(in)      :: energy(:)
        type(AD_D),                 intent(in)      :: grad1_energy(:)
        type(AD_D),                 intent(in)      :: grad2_energy(:)
        type(AD_D),                 intent(in)      :: grad3_energy(:)
        type(AD_D),                 intent(in)      :: mu(:)
        type(AD_D),                 intent(in)      :: lambda(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_11(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_22(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_33(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_12(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_13(:)
        type(AD_D), allocatable,    intent(inout)   :: shear_23(:)

        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u, grad2_u, grad3_u,              &
            grad1_v, grad2_v, grad3_v,              &
            grad1_w, grad2_w, grad3_w,              &
            invdensity, invdensity2, div_velocity, u, v


!        !
!        ! compute velocity jacobians
!        !
!        invdensity = ONE/density
!        du_ddensity  = -invdensity*invdensity*mom1
!        dv_ddensity  = -invdensity*invdensity*mom2
!        dw_ddensity  = -invdensity*invdensity*mom3
!
!        du_dmom1 = invdensity
!        dv_dmom2 = invdensity
!        dw_dmom3 = invdensity
!
!
!
!        !
!        ! compute velocity gradients via chain rule:
!        !
!        !   u = f(rho,rhou)
!        !
!        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
!        !
!        grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
!        grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1
!        grad3_u = du_ddensity*grad3_density  +  du_dmom1*grad3_mom1
!
!        grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
!        grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2
!        grad3_v = dv_ddensity*grad3_density  +  dv_dmom2*grad3_mom2
!
!        grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
!        grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
!        grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3



        !
        ! compute velocity jacobians
        !
        invdensity = ONE/density
        invdensity2 = -invdensity*invdensity
!        du_ddensity  = -invdensity*invdensity*mom1
!        dv_ddensity  = -invdensity*invdensity*mom2
!        dw_ddensity  = -invdensity*invdensity*mom3
!
!        du_dmom1 = invdensity
!        dv_dmom2 = invdensity
!        dw_dmom3 = invdensity



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = invdensity2*mom1*grad1_density  +  invdensity*grad1_mom1
        grad2_u = invdensity2*mom1*grad2_density  +  invdensity*grad2_mom1
        grad3_u = invdensity2*mom1*grad3_density  +  invdensity*grad3_mom1

        grad1_v = invdensity2*mom2*grad1_density  +  invdensity*grad1_mom2
        grad2_v = invdensity2*mom2*grad2_density  +  invdensity*grad2_mom2
        grad3_v = invdensity2*mom2*grad3_density  +  invdensity*grad3_mom2

        grad1_w = invdensity2*mom3*grad1_density  +  invdensity*grad1_mom3
        grad2_w = invdensity2*mom3*grad2_density  +  invdensity*grad2_mom3
        grad3_w = invdensity2*mom3*grad3_density  +  invdensity*grad3_mom3













        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        !if (worker%coordinate_system() == 'Cartesian') then
        if (coordinate_system == 'Cartesian') then

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
            shear_11 = TWO*mu*grad1_u  +  lambda*(div_velocity)
            shear_22 = TWO*mu*grad2_v  +  lambda*(div_velocity)
            shear_33 = TWO*mu*grad3_w  +  lambda*(div_velocity)

            shear_12 = mu*(grad2_u + grad1_v)
            shear_13 = mu*(grad3_u + grad1_w)
            shear_23 = mu*(grad2_w + grad3_v)




        else if (coordinate_system == 'Cylindrical') then

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
            u = mom1/density
            v = mom2/density
            div_velocity = (grad1_u + u/r) + grad2_v + grad3_w

            !
            ! Compute shear stress components
            !
            shear_11 = TWO*mu*(grad1_u        )  +  lambda*(div_velocity)
            shear_22 = TWO*mu*(grad2_v + (u/r))  +  lambda*(div_velocity)
            shear_33 = TWO*mu*(grad3_w        )  +  lambda*(div_velocity)

            shear_12 = mu*(grad2_u + grad1_v  - (v/r))
            shear_13 = mu*(grad3_u + grad1_w         )
            shear_23 = mu*(grad2_w + grad3_v         )

        end if


    end subroutine compute_shear_stress
    !*************************************************************************









    !>
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------
    subroutine compute_vorticity(coordinate_system,r,                                &
                                 density,grad1_density,grad2_density,grad3_density,  &
                                 mom1,   grad1_mom1,   grad2_mom1,   grad3_mom1,     &
                                 mom2,   grad1_mom2,   grad2_mom2,   grad3_mom2,     &
                                 mom3,   grad1_mom3,   grad2_mom3,   grad3_mom3,     &
                                 vorticity_1, vorticity_2, vorticity_3)
        character(*),               intent(in)      :: coordinate_system
        real(rk),                   intent(in)      :: r(:)
        type(AD_D),                 intent(in)      :: density(:)
        type(AD_D),                 intent(in)      :: grad1_density(:)
        type(AD_D),                 intent(in)      :: grad2_density(:)
        type(AD_D),                 intent(in)      :: grad3_density(:)
        type(AD_D),                 intent(in)      :: mom1(:)
        type(AD_D),                 intent(in)      :: grad1_mom1(:)
        type(AD_D),                 intent(in)      :: grad2_mom1(:)
        type(AD_D),                 intent(in)      :: grad3_mom1(:)
        type(AD_D),                 intent(in)      :: mom2(:)
        type(AD_D),                 intent(in)      :: grad1_mom2(:)
        type(AD_D),                 intent(in)      :: grad2_mom2(:)
        type(AD_D),                 intent(in)      :: grad3_mom2(:)
        type(AD_D),                 intent(in)      :: mom3(:)
        type(AD_D),                 intent(in)      :: grad1_mom3(:)
        type(AD_D),                 intent(in)      :: grad2_mom3(:)
        type(AD_D),                 intent(in)      :: grad3_mom3(:)
        type(AD_D), allocatable,    intent(inout)   :: vorticity_1(:)
        type(AD_D), allocatable,    intent(inout)   :: vorticity_2(:)
        type(AD_D), allocatable,    intent(inout)   :: vorticity_3(:)

        type(AD_D), allocatable, dimension(:)   ::  &
            invdensity, du_ddensity, dv_ddensity, dw_ddensity,  &
            du_dmom1, dv_dmom2, dw_dmom3,                       &
            grad1_u, grad2_u, grad3_u,                          &
            grad1_v, grad2_v, grad3_v,                          &
            grad1_w, grad2_w, grad3_w, v


        !
        ! compute velocity jacobians
        !
        invdensity  = ONE/density
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
        if (trim(coordinate_system) == 'Cartesian') then
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



        else if (trim(coordinate_system) == 'Cylindrical') then
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
            v = mom2/density

            vorticity_1 =  (grad2_w - grad3_v)
            vorticity_2 =  (grad3_u - grad1_w)
            vorticity_3 =  (grad1_v - grad2_u + (v/r)) 


        end if


    end subroutine compute_vorticity











































end module mod_rans_efficient
