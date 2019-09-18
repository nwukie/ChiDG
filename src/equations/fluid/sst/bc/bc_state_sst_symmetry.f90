module bc_state_sst_symmetry
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none
    


    !>  Symmetry boundary condition for Spalart-Allmaras turbulent working variable.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: sst_symmetry_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type sst_symmetry_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_symmetry_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('SST Symmetry')
        call self%set_family('Symmetry')

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
        class(sst_symmetry_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_m, grad1_density_m, grad2_density_m, grad3_density_m, normal_grad
            
        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3

        !
        ! Get unit normal vector
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m       = worker%get_field('Density * k', 'value', 'face interior')
        grad1_density_m = worker%get_field('Density * k', 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density * k', 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density * k', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * k', density_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        normal_grad = grad1_density_m*unorm_1 + grad2_density_m*unorm_2 + grad3_density_m*unorm_3
        grad1_density_m = grad1_density_m-normal_grad*unorm_1
        grad2_density_m = grad2_density_m-normal_grad*unorm_2
        grad3_density_m = grad3_density_m-normal_grad*unorm_3
        call worker%store_bc_state('Density * k', grad1_density_m, 'grad1')
        call worker%store_bc_state('Density * k', grad2_density_m, 'grad2')
        call worker%store_bc_state('Density * k', grad3_density_m, 'grad3')
                                                
        ! Omega;
        density_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density * Omega', 'grad3', 'face interior')


        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Omega', density_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        normal_grad = grad1_density_m*unorm_1 + grad2_density_m*unorm_2 + grad3_density_m*unorm_3
        grad1_density_m = grad1_density_m-normal_grad*unorm_1
        grad2_density_m = grad2_density_m-normal_grad*unorm_2
        grad3_density_m = grad3_density_m-normal_grad*unorm_3
        call worker%store_bc_state('Density * Omega', grad1_density_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_m, 'grad3')
                                                



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_sst_symmetry
