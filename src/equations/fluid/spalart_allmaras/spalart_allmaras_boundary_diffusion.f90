module spalart_allmaras_boundary_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use mod_spalart_allmaras,   only: SA_c_n1, SA_sigma
    implicit none

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_boundary_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_boundary_diffusion_operator_t
    !********************************************************************************


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_boundary_diffusion_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('Spalart-Allmaras Boundary Diffusion Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Diffusive Flux')

        ! Set operator equations
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_boundary_diffusion_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                                   intent(inout)   :: worker
        class(properties_t),                                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            rho_m, rho_nutilde_m, invrho_m,                 &
            rho_p, rho_nutilde_p, invrho_p,                 &
            nutilde_m, chi_m, f_n1_m, nu_l_m, mu_l_m,       &
            nutilde_p, chi_p, f_n1_p, nu_l_p, mu_l_p,       &
            grad1_rho_m, grad1_nutilde_m, grad1_rho_nutilde_m,    &
            grad2_rho_m, grad2_nutilde_m, grad2_rho_nutilde_m,    &
            grad3_rho_m, grad3_nutilde_m, grad3_rho_nutilde_m,    &
            grad1_rho_p, grad1_nutilde_p, grad1_rho_nutilde_p,    &
            grad2_rho_p, grad2_nutilde_p, grad2_rho_nutilde_p,    &
            grad3_rho_p, grad3_nutilde_p, grad3_rho_nutilde_p,    &
            dnutilde_drho_m, dnutilde_drhonutilde_m,        &
            dnutilde_drho_p, dnutilde_drhonutilde_p,        &
            diffusion_m, diffusion_p,                       &
            flux_1_m, flux_2_m, flux_3_m,                   &
            flux_1_p, flux_2_p, flux_3_p,                   &
            flux_1, flux_2, flux_3, integrand


        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m         = worker%get_primary_field_face('Density',          'value', 'face interior')
        rho_p         = worker%get_primary_field_face('Density',          'value', 'face exterior')

        rho_nutilde_m = worker%get_primary_field_face('Density * NuTilde','value', 'face interior')
        rho_nutilde_p = worker%get_primary_field_face('Density * NuTilde','value', 'face exterior')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_rho_m         = worker%get_primary_field_face('Density',          'grad1+lift', 'face interior')
        grad2_rho_m         = worker%get_primary_field_face('Density',          'grad2+lift', 'face interior')
        grad3_rho_m         = worker%get_primary_field_face('Density',          'grad3+lift', 'face interior')
        grad1_rho_p         = worker%get_primary_field_face('Density',          'grad1+lift', 'face exterior')
        grad2_rho_p         = worker%get_primary_field_face('Density',          'grad2+lift', 'face exterior')
        grad3_rho_p         = worker%get_primary_field_face('Density',          'grad3+lift', 'face exterior')


        grad1_rho_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad1+lift', 'face interior')
        grad2_rho_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad2+lift', 'face interior')
        grad3_rho_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad3+lift', 'face interior')
        grad1_rho_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad1+lift', 'face exterior')
        grad2_rho_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad2+lift', 'face exterior')
        grad3_rho_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad3+lift', 'face exterior')





        invrho_m = ONE/rho_m
        invrho_p = ONE/rho_p


        !
        ! Get normal vector
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)




        !
        ! Get model fields:
        !   Viscosity
        !
        mu_l_m = worker%get_model_field_face('Laminar Viscosity','value', 'face interior')
        mu_l_p = worker%get_model_field_face('Laminar Viscosity','value', 'face exterior')

        nu_l_m = mu_l_m*invrho_m
        nu_l_p = mu_l_p*invrho_p


        !
        ! Compute nutilde
        !
        nutilde_m = rho_nutilde_m*invrho_m
        nutilde_p = rho_nutilde_p*invrho_p


        chi_m = nutilde_m/nu_l_m
        chi_p = nutilde_p/nu_l_p



        ! Initialize derivatives first
        f_n1_m = rho_m
        where (nutilde_m >= ZERO)
            f_n1_m = ONE
        else where
            f_n1_m = (SA_c_n1 + chi_m*chi_m*chi_m)/(SA_c_n1 - chi_m*chi_m*chi_m)
        end where

        f_n1_p = rho_p
        where (nutilde_p >= ZERO)
            f_n1_p = ONE
        else where
            f_n1_p = (SA_c_n1 + chi_p*chi_p*chi_p)/(SA_c_n1 - chi_p*chi_p*chi_p)
        end where




        !
        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        !
        dnutilde_drho_m        = -rho_nutilde_m*invrho_m*invrho_m
        dnutilde_drhonutilde_m =  invrho_m

        dnutilde_drho_p        = -rho_nutilde_p*invrho_p*invrho_p
        dnutilde_drhonutilde_p =  invrho_p

        grad1_nutilde_m = dnutilde_drho_m * grad1_rho_m  +  dnutilde_drhonutilde_m * grad1_rho_nutilde_m
        grad2_nutilde_m = dnutilde_drho_m * grad2_rho_m  +  dnutilde_drhonutilde_m * grad2_rho_nutilde_m
        grad3_nutilde_m = dnutilde_drho_m * grad3_rho_m  +  dnutilde_drhonutilde_m * grad3_rho_nutilde_m

        grad1_nutilde_p = dnutilde_drho_p * grad1_rho_p  +  dnutilde_drhonutilde_p * grad1_rho_nutilde_p
        grad2_nutilde_p = dnutilde_drho_p * grad2_rho_p  +  dnutilde_drhonutilde_p * grad2_rho_nutilde_p
        grad3_nutilde_p = dnutilde_drho_p * grad3_rho_p  +  dnutilde_drhonutilde_p * grad3_rho_nutilde_p




        !-----------------------------------------
        !            TURBULENCE FLUX
        !-----------------------------------------
        diffusion_m = -(ONE/SA_sigma)*(mu_l_m + f_n1_m*rho_nutilde_m)
        diffusion_p = -(ONE/SA_sigma)*(mu_l_p + f_n1_p*rho_nutilde_p)

        flux_1_m = diffusion_m*grad1_nutilde_m
        flux_2_m = diffusion_m*grad2_nutilde_m
        flux_3_m = diffusion_m*grad3_nutilde_m

        flux_1_p = diffusion_p*grad1_nutilde_p
        flux_2_p = diffusion_p*grad2_nutilde_p
        flux_3_p = diffusion_p*grad3_nutilde_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !********************************************************************************************












end module spalart_allmaras_boundary_diffusion
