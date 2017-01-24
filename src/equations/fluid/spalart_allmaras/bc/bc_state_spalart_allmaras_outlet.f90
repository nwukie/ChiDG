module bc_state_spalart_allmaras_outlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none
    


    !>  Outlet boundary condition for Spalart-Allmaras turbulent working variable.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: spalart_allmaras_outlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type spalart_allmaras_outlet_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_outlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Spalart Allmaras Outlet")
        call self%set_family("Outlet")

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
        class(spalart_allmaras_outlet_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            rhoNutilde_m, drhoNutilde_dx_m, drhoNutilde_dy_m, drhoNutilde_dz_m


        !
        ! Interpolate interior solution to quadrature nodes
        !
        rhoNutilde_m     = worker%get_primary_field_face("Density * NuTilde", 'value', 'face interior')
        drhoNutilde_dx_m = worker%get_primary_field_face("Density * NuTilde", 'ddx',   'face interior')
        drhoNutilde_dy_m = worker%get_primary_field_face("Density * NuTilde", 'ddy',   'face interior')
        drhoNutilde_dz_m = worker%get_primary_field_face("Density * NuTilde", 'ddz',   'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state("Density * NuTilde", rhoNutilde_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        drhoNutilde_dx_m = ZERO
        drhoNutilde_dy_m = ZERO
        drhoNutilde_dz_m = ZERO
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dx_m, 'ddx')
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dy_m, 'ddy')
        call worker%store_bc_state("Density * NuTilde", drhoNutilde_dz_m, 'ddz')
                                                
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_spalart_allmaras_outlet
