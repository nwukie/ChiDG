module bc_state_sst_outlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none
    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: sst_outlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type sst_outlet_t
    !*******************************************************************************************




contains


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_outlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('SST Outlet')
        call self%set_family('Outlet')

    end subroutine init
    !********************************************************************************







    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(sst_outlet_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            distance_elem, density_m, mu_l, &
            density_k_m, grad1_density_k_m, grad2_density_k_m, grad3_density_k_m, & 
            density_omega_m, grad1_density_omega_m, grad2_density_omega_m, grad3_density_omega_m

        real(rk) :: distance

        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_omega_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_omega_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_omega_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_omega_m = worker%get_field('Density * Omega', 'grad3', 'face interior')




        !
        ! Store boundary condition gradient - Extrapolate
        !
        call worker%store_bc_state('Density * Omega', density_omega_m, 'value')
        call worker%store_bc_state('Density * Omega', grad1_density_omega_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_omega_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_omega_m, 'grad3')
                                                

        ! R_11
        density_k_m       = worker%get_field('Density * k', 'value', 'face interior')
        grad1_density_k_m = worker%get_field('Density * k', 'grad1', 'face interior')
        grad2_density_k_m = worker%get_field('Density * k', 'grad2', 'face interior')
        grad3_density_k_m = worker%get_field('Density * k', 'grad3', 'face interior')

        call worker%store_bc_state('Density * k', density_k_m, 'value')
        call worker%store_bc_state('Density * k', grad1_density_k_m, 'grad1')
        call worker%store_bc_state('Density * k', grad2_density_k_m, 'grad2')
        call worker%store_bc_state('Density * k', grad3_density_k_m, 'grad3')

        
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_sst_outlet
