module bc_state_spalart_allmaras_outlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
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
        call self%set_name('Spalart Allmaras Outlet')
        call self%set_family('Outlet')

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
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(spalart_allmaras_outlet_t),   intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_nutilde_m, grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_nutilde_m       = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        grad1_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * NuTilde', density_nutilde_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_nutilde_m = ZERO
        grad2_density_nutilde_m = ZERO
        grad3_density_nutilde_m = ZERO
        call worker%store_bc_state('Density * NuTilde', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * NuTilde', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * NuTilde', grad3_density_nutilde_m, 'grad3')
                                                
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_spalart_allmaras_outlet
