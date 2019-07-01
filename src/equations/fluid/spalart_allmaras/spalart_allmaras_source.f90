module spalart_allmaras_source
#include<messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use mod_chidg_mpi,          only: GLOBAL_MASTER
    use mod_spalart_allmaras,   only: SA_c_b1, SA_c_b2, SA_kappa, SA_sigma,         &
                                      SA_c_w1, SA_c_w2, SA_c_w3, SA_c_v1, SA_c_v2,  &
                                      SA_c_v3, SA_c_t3, SA_c_t4, SA_c_n1,           &
                                      SA_Pr_t, SA_rlim, SA_b    

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private

    
    !>  Spalart-Allmaras Source Term.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_source_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_source_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_source_operator_t),   intent(inout)      :: self

        ! Set operator name.
        call self%set_name('Spalart-Allmaras Source Operator')

        ! Set operator type.
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations being integrated.
        call self%add_primary_field('Density * NuTilde')

        ! Set auxiliary variables being used.
        call self%add_auxiliary_field('Wall Distance : p-Poisson')

        ! Add Turbulent Egrad2 Viscosity model
        call self%add_model('Spalart Allmaras Turbulent Model Fields')
        call self%add_model('Wall Distance : p-Poisson Normalization')
        call self%add_model('Vorticity')

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_source_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                            &
            rho, rho_nutilde, invrho,                                       &
            mu, nu, lamda,                                                  &
            grad1_rho,          grad2_rho,          grad3_rho,              &
            grad1_rho_nutilde,  grad2_rho_nutilde,  grad3_rho_nutilde,      &
            vorticity2, vorticity, vorticity_bar, vorticity_mod,            &
            f_v1, f_v2, f_t2, f_w, f_n1, production, destruction,           &
            dnutilde_drho, dnutilde_drho_nutilde, chi, r, rbar, g, mu_t,    &
            nutilde, grad1_nutilde, grad2_nutilde, grad3_nutilde,           &
            source, dwall, vorticity_1, vorticity_2, vorticity_3

        real(rk)    :: const, epsilon_vorticity, eps


        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_field('Density'          , 'value', 'element')
        rho_nutilde = worker%get_field('Density * NuTilde', 'value', 'element')



        !
        ! Interpolate solution gradients to quadrature nodes
        !
        grad1_rho         = worker%get_field('Density',           'grad1', 'element')
        grad2_rho         = worker%get_field('Density',           'grad2', 'element')
        grad3_rho         = worker%get_field('Density',           'grad3', 'element')

        grad1_rho_nutilde = worker%get_field('Density * NuTilde', 'grad1', 'element')
        grad2_rho_nutilde = worker%get_field('Density * NuTilde', 'grad2', 'element')
        grad3_rho_nutilde = worker%get_field('Density * NuTilde', 'grad3', 'element')



        !
        ! Interpolate auxiliary field, Wall Distance
        !
        !eps = 1.e-10_rk
        eps = 1.e-11_rk
        !eps = ZERO
        dwall = worker%get_field('Wall Distance', 'value', 'element')
        if (any(ieee_is_nan(dwall(:)%x_ad_))) call write_line('dwall is nan')

        !
        ! Divide by density
        !
        invrho  = ONE/rho
        nutilde = rho_nutilde*invrho



        !
        ! Compute model values
        !
        mu  = worker%get_field('Laminar Viscosity', 'value', 'element')
        nu  = mu*invrho


        !
        ! Compute turbulence viscosity
        !
        chi  = nutilde/nu
        f_v1 = (chi*chi*chi)/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)

        mu_t = rho_nutilde
        mu_t = ZERO
        where(nutilde >= ZERO)
            mu_t = rho_nutilde * f_v1
        end where
        

        !
        ! Compute f_n1
        !
        f_n1 = rho_nutilde
        f_n1 = ONE
        where(nutilde < ZERO)
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where


        !
        ! Compute vorticity and modified vorticity
        !
        vorticity_1 = worker%get_field('Vorticity-1', 'value', 'element')
        vorticity_2 = worker%get_field('Vorticity-2', 'value', 'element')
        vorticity_3 = worker%get_field('Vorticity-3', 'value', 'element')


        vorticity2 =  vorticity_1**TWO  +  vorticity_2**TWO  +  vorticity_3**TWO 
        epsilon_vorticity = 1.e-6_rk
        vorticity = vorticity2
        where(vorticity2 < epsilon_vorticity)
            vorticity = HALF*(epsilon_vorticity + vorticity2/epsilon_vorticity)
        else where
            vorticity = sqrt(vorticity2)
        end where


        f_v2 = ONE - (chi/(ONE+chi*f_v1))
        vorticity_bar = (nutilde/(SA_kappa*SA_kappa*(dwall*dwall + eps)))*f_v2


        vorticity_mod = vorticity
        where (vorticity_bar >= -SA_c_v2*vorticity)
            vorticity_mod = vorticity + vorticity_bar
        else where
            vorticity_mod = vorticity + vorticity*(SA_c_v2*SA_c_v2*vorticity + SA_c_v3*vorticity_bar)/( (SA_c_v3 - TWO*SA_c_v2)*vorticity - vorticity_bar ) 
        end where


        !
        ! Compute f_t2, f_w, g, r
        !
        rbar = nutilde/(vorticity_mod * SA_kappa * SA_kappa * (dwall*dwall + eps) )

        !r = min(SA_rlim,rbar)
        r = SA_rlim - (SA_rlim - rbar)*(atan(SA_b * (SA_rlim - rbar))/PI  + HALF)  + atan(SA_b)/PI  -  HALF

        g = r + SA_c_w2*(r**SIX  -  r)

        f_w = g*(((ONE + SA_c_w3**SIX)/(g**SIX + SA_c_w3**SIX))**(ONE/SIX))

        f_t2 = SA_c_t3*exp(-SA_c_t4*chi*chi)



        !
        ! Compute Production, Destruction
        !
        production = vorticity_mod
        where ( nutilde >= ZERO )
            production = SA_c_b1*(ONE - f_t2)*vorticity_mod*nutilde
        else where
            production = SA_C_b1*(ONE - SA_c_t3)*vorticity*nutilde
        end where


        destruction = vorticity_mod
        where ( nutilde >= ZERO )
            !destruction = (SA_c_w1*f_w - (SA_c_b1/(SA_kappa*SA_kappa))*f_t2) * (nutilde/dwall)**TWO
            destruction = (SA_c_w1*f_w - (SA_c_b1/(SA_kappa*SA_kappa))*f_t2) * (nutilde*nutilde/(dwall*dwall + eps))
        else where
            !destruction = -SA_c_w1 * (nutilde/dwall)**TWO
            destruction = -SA_c_w1 * (nutilde*nutilde/(dwall*dwall + eps))
        end where


        !
        ! Compute jacobian of nutilde
        !
        dnutilde_drho         = -invrho*invrho*rho_nutilde
        dnutilde_drho_nutilde =  invrho


        !
        ! Compute gradient of nutilde
        !
        grad1_nutilde = dnutilde_drho*grad1_rho  +  dnutilde_drho_nutilde*grad1_rho_nutilde
        grad2_nutilde = dnutilde_drho*grad2_rho  +  dnutilde_drho_nutilde*grad2_rho_nutilde
        grad3_nutilde = dnutilde_drho*grad3_rho  +  dnutilde_drho_nutilde*grad3_rho_nutilde


        !========================================================================
        !                       Spalart-Allmaras Source Term
        !========================================================================
        source = -(                                 &
                    -rho*(production-destruction)   &
                    -(SA_c_b2/SA_sigma)*rho*(grad1_nutilde*grad1_nutilde + grad2_nutilde*grad2_nutilde + grad3_nutilde*grad3_nutilde)   &
                    +(ONE/SA_sigma)*(nu + f_n1*nutilde)*(grad1_rho*grad1_nutilde + grad2_rho*grad2_nutilde + grad3_rho*grad3_nutilde)   &
                  )


        call worker%integrate_volume_source('Density * NuTilde',source)


    end subroutine compute
    !*********************************************************************************************************






end module spalart_allmaras_source
