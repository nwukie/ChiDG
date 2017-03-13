module bc_state_spalart_allmaras_inlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_point,             only: point_t
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none
    


    !>  Inflow boundary condition for Spalart-Allmaras turbulent working variable.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: spalart_allmaras_inlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type spalart_allmaras_inlet_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_inlet_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name('Spalart Allmaras Inlet')
        call self%set_family('Inlet')


        !
        ! Add turbulent inlet parameter and default value.
        !
        call self%bcproperties%add('Turbulent Viscosity Ratio',   'Required')
        call self%set_fcn_option('Turbulent Viscosity Ratio', 'val', 3._rk)


    end subroutine init
    !********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(spalart_allmaras_inlet_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D),     allocatable, dimension(:)   ::  &
            density_nutilde_m, density_nutilde_bc,      &
            grad1_density_nutilde_m,                    &
            grad2_density_nutilde_m,                    &
            grad3_density_nutilde_m,                    &
            density_m, mu_m, nu_m

        real(rk),       allocatable, dimension(:)   :: nutilde_nu


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m               = worker%get_primary_field_face('Density',           'value', 'face interior')
        density_nutilde_m       = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')

        grad1_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad3', 'face interior')


        
        !
        ! Get User boundary condition viscosity ratio
        !
        nutilde_nu = self%bcproperties%compute('Turbulent Viscosity Ratio',worker%time(), worker%coords())


        !
        ! Get viscosity from interior
        !
        mu_m = worker%get_model_field_face('Laminar Viscosity', 'value', 'face interior')
        nu_m = mu_m/density_m


        !
        ! Compute boundary condition state
        !
        density_nutilde_bc = density_m * (nutilde_nu * nu_m)



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * NuTilde', density_nutilde_bc,'value')


        !
        ! Store boundary condition gradient - Extrapolate
        !
        call worker%store_bc_state('Density * NuTilde', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * NuTilde', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * NuTilde', grad3_density_nutilde_m, 'grad3')
                                                
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_spalart_allmaras_inlet
