module bc_state_spalart_allmaras_inlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
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
    type, public, extends(bc_state_t) :: bc_state_spalart_allmaras_inlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type bc_state_spalart_allmaras_inlet_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(bc_state_spalart_allmaras_inlet_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name("Spalart Allmaras Inlet")
        call self%set_family("Inlet")


        !
        ! Add turbulent inlet parameter and default value.
        !
        call self%bcproperties%add('Turbulent Viscosity Ratio',   'Required')
        call self%set_fcn_option('Turbulent Viscosity Ratio', 'val', 1._rk)


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
        class(bc_state_spalart_allmaras_inlet_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                  &
            rhoNutilde_m, rhoNutilde_bc,                            &
            drhoNutilde_dx_m, drhoNutilde_dy_m, drhoNutilde_dz_m

        real(rk),   allocatable, dimension(:)   :: nutilde_nu

        !
        ! Interpolate interior solution to quadrature nodes
        !
        rhoNutilde_m     = worker%get_primary_field_face("Density * NuTilde", 'value', 'face interior')
        drhoNutilde_dx_m = worker%get_primary_field_face("Density * NuTilde", 'ddx',   'face interior')
        drhoNutilde_dy_m = worker%get_primary_field_face("Density * NuTilde", 'ddy',   'face interior')
        drhoNutilde_dz_m = worker%get_primary_field_face("Density * NuTilde", 'ddz',   'face interior')



        !
        ! Get User boundary condition viscosity ratio
        !
        coords = worker%coords()
        time   = worker%time()
        nutilde_nu = self%bcproperties%compute("Turbulent Viscosity Ratio",time, coords)


        !
        ! Compute boundary condition state
        !
        rhoNutilde_bc = nutilde_nu * rhoNutilde_m



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state("Density * NuTilde", rhoNutilde_bc,'value')


        !
        ! Store boundary condition gradient - Extrapolate
        !
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dx_m, 'ddx')
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dy_m, 'ddy')
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dz_m, 'ddz')
                                                
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_spalart_allmaras_inlet
