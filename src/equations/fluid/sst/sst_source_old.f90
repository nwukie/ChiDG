module sst_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_source_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_source_operator_t
    !******************************************************************************










contains

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_source_operator_t),   intent(inout)      :: self

        ! Set operator name.
        call self%set_name('SST Source Operator')

        ! Set operator type.
        call self%set_operator_type('Volume Diffusive Operator')
        call self%add_auxiliary_field('Wall Distance : p-Poisson')
        ! Set operator equations being integrated.
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')

    end subroutine init
    !********************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_source_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                            &
            density, density_omega, src_term, omega, k,           &
            omega_source_term,omega_grad_sq,                                     &
            density_k, &
            beta_k, prod_k, dest_k, prod_k_mod,  &
            prod_om, alpha_w, beta_w, sigma_w, sigma_d, &
            source, mu_l, mu_t

        real(rk)    :: const, epsilon_vorticity, eps
        

        !
        ! Interpolate solution to quadrature nodes
        !

        mu_l       = worker%get_field('Laminar Viscosity',           'value', 'element')
        mu_t       = worker%get_field('Turbulent Viscosity',           'value', 'element')

        density         = worker%get_field('Density'          , 'value', 'element')

        ! k source term
        beta_k      = worker%get_field('SST beta_k', 'value', 'element')
        prod_k      = worker%get_field('SST k Production Term', 'value', 'element')
        dest_k      = worker%get_field('SST k Destruction Term', 'value', 'element')
        
        dest_k      = beta_k*dest_k
        prod_k_mod  = prod_k
        ! Limit production
        prod_k_mod  = min(prod_k, 10.0_rk*dest_k)

        ! Omega source term
        prod_om     = worker%get_field('SST k Production Term', 'value', 'element')
        alpha_w     = worker%get_field('SST alpha_w', 'value', 'element')
        beta_w      = worker%get_field('SST beta_w', 'value', 'element')
        sigma_w     = worker%get_field('SST sigma_w', 'value', 'element')
        sigma_d     = worker%get_field('SST sigma_d', 'value', 'element')

        omega       = worker%get_field('Omega', 'value', 'element')

        omega_source_term   = worker%get_field('SST Omega Source Term', 'value', 'element')
        omega_grad_sq       = worker%get_field('Omega Gradient Squared', 'value', 'element')

        


        !========================================================================
        !                       Omega Source Term
        !========================================================================
        source  = alpha_w*density/(mu_t*exp(omega))*prod_om - &
                    beta_w*density*exp(omega) + sigma_d*omega_source_term + &
                    (mu_l + sigma_w*mu_t)*omega_grad_sq

        call worker%integrate_volume_source('Density * Omega',source)
        
        !========================================================================
        !                       k Source Term
        !========================================================================
        source = prod_k_mod - dest_k

        call worker%integrate_volume_source('Density * k',source)



    end subroutine compute
    !*********************************************************************************************************






end module sst_source
