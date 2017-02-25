module spalart_allmaras_bc_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use mod_spalart_allmaras,   only: SA_c_n1, SA_sigma
    implicit none

    private



    !>  Volume diffusion integral for Spalart-Allmaras flux terms.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-----------------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_bc_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_bc_diffusion_operator_t
    !*****************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_bc_diffusion_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Spalart-Allmaras BC Diffusion Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Diffusive Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density * NuTilde")

    end subroutine init
    !*****************************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_bc_diffusion_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            rho, rho_nutilde, invrho,                       &
            nutilde, chi, f_n1, nu_l, mu_l,                 &
            grad1_rho, grad1_nutilde, grad1_rho_nutilde,    &
            grad2_rho, grad2_nutilde, grad2_rho_nutilde,    &
            grad3_rho, grad3_nutilde, grad3_rho_nutilde,    &
            dnutilde_drho, dnutilde_drhonutilde,            &
            flux_1, flux_2, flux_3, diffusion, integrand


        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3



        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_primary_field_face('Density',          'value','boundary')
        rho_nutilde = worker%get_primary_field_face('Density * NuTilde','value','boundary')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_rho         = worker%get_primary_field_face('Density',          'grad1+lift','boundary')
        grad2_rho         = worker%get_primary_field_face('Density',          'grad2+lift','boundary')
        grad3_rho         = worker%get_primary_field_face('Density',          'grad3+lift','boundary')

        grad1_rho_nutilde = worker%get_primary_field_face('Density * NuTilde','grad1+lift','boundary')
        grad2_rho_nutilde = worker%get_primary_field_face('Density * NuTilde','grad2+lift','boundary')
        grad3_rho_nutilde = worker%get_primary_field_face('Density * NuTilde','grad3+lift','boundary')





        invrho = ONE/rho


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
        mu_l = worker%get_model_field_face('Laminar Viscosity','value','boundary')
        nu_l = mu_l*invrho


        !
        ! Compute nutilde, chi
        !
        nutilde = rho_nutilde*invrho
        chi     = nutilde/nu_l



        ! Initialize derivatives first
        f_n1 = rho
        where (nutilde >= ZERO)
            f_n1 = ONE
        else where
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where



        !
        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        !
        dnutilde_drho        = -rho_nutilde*invrho*invrho
        dnutilde_drhonutilde =  invrho


        grad1_nutilde = dnutilde_drho * grad1_rho  +  dnutilde_drhonutilde * grad1_rho_nutilde
        grad2_nutilde = dnutilde_drho * grad2_rho  +  dnutilde_drhonutilde * grad2_rho_nutilde
        grad3_nutilde = dnutilde_drho * grad3_rho  +  dnutilde_drhonutilde * grad3_rho_nutilde


        !-------------------------------------
        !           TURBULENCE FLUX
        !-------------------------------------
        diffusion = -(ONE/SA_sigma)*(mu_l + f_n1*rho_nutilde)


        flux_1 = diffusion*grad1_nutilde
        flux_2 = diffusion*grad2_nutilde
        flux_3 = diffusion*grad3_nutilde

        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3

        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !******************************************************************************************












end module spalart_allmaras_bc_diffusion
