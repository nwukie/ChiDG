module spalart_allmaras_volume_diffusion
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
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_volume_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_volume_diffusion_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_volume_diffusion_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Spalart-Allmaras Volume Diffusion Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Volume Diffusive Flux')

        !
        ! Set operator equations
        !
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
        class(spalart_allmaras_volume_diffusion_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                                   intent(inout)   :: worker
        class(properties_t),                                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                    &
            density, density_nutilde, invdensity,                   &
            nutilde, chi, f_n1, nu_l, mu_l,                         &
            grad1_density, grad1_nutilde, grad1_density_nutilde,    &
            grad2_density, grad2_nutilde, grad2_density_nutilde,    &
            grad3_density, grad3_nutilde, grad3_density_nutilde,    &
            dnutilde_ddensity, dnutilde_ddensitynutilde,            &
            flux_1, flux_2, flux_3, diffusion


        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3



        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_field('Density',           'value', 'element')
        density_nutilde = worker%get_field('Density * NuTilde', 'value', 'element')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density         = worker%get_field('Density',           'grad1')
        grad2_density         = worker%get_field('Density',           'grad2')
        grad3_density         = worker%get_field('Density',           'grad3')

        grad1_density_nutilde = worker%get_field('Density * NuTilde', 'grad1')
        grad2_density_nutilde = worker%get_field('Density * NuTilde', 'grad2')
        grad3_density_nutilde = worker%get_field('Density * NuTilde', 'grad3')





        invdensity = ONE/density


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
        mu_l = worker%get_field('Laminar Viscosity', 'value', 'element')
        nu_l = mu_l*invdensity


        !
        ! Compute nutilde, chi
        !
        nutilde = density_nutilde*invdensity
        chi     = nutilde/nu_l



        ! Initialize derivatives first
        f_n1 = density
        where (nutilde >= ZERO)
            f_n1 = ONE
        else where
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where



        !
        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        !
        dnutilde_ddensity        = -density_nutilde*invdensity*invdensity
        dnutilde_ddensitynutilde =  invdensity


        grad1_nutilde = dnutilde_ddensity * grad1_density  +  dnutilde_ddensitynutilde * grad1_density_nutilde
        grad2_nutilde = dnutilde_ddensity * grad2_density  +  dnutilde_ddensitynutilde * grad2_density_nutilde
        grad3_nutilde = dnutilde_ddensity * grad3_density  +  dnutilde_ddensitynutilde * grad3_density_nutilde




        !================================
        !       TURBULENCE FLUX
        !================================
        diffusion = -(ONE/SA_sigma)*(mu_l + f_n1*density_nutilde)

        flux_1 = diffusion*grad1_nutilde
        flux_2 = diffusion*grad2_nutilde
        flux_3 = diffusion*grad3_nutilde


        call worker%integrate_volume_flux('Density * NuTilde','Diffusion',flux_1,flux_2,flux_3)


    end subroutine compute
    !************************************************************************************************












end module spalart_allmaras_volume_diffusion
