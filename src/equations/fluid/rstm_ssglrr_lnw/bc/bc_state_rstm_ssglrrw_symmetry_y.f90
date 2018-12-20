module bc_state_rstm_ssglrrw_symmetry_y
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
    type, public, extends(bc_state_t) :: rstm_ssglrrw_symmetry_y_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type rstm_ssglrrw_symmetry_y_t
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
        class(rstm_ssglrrw_symmetry_y_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Symmetry-y')
        call self%set_family('Symmetry')

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
        class(rstm_ssglrrw_symmetry_y_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_reynolds_m, grad1_density_reynolds_m, grad2_density_reynolds_m, grad3_density_reynolds_m, & 
            density_nutilde_m, grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_nutilde_m       = worker%get_field('Density * Omega', 'value', 'face interior')
        grad1_density_nutilde_m = worker%get_field('Density * Omega', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_field('Density * Omega', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_field('Density * Omega', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Omega', density_nutilde_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_nutilde_m = ZERO
        grad2_density_nutilde_m = ZERO
        grad3_density_nutilde_m = ZERO
        call worker%store_bc_state('Density * Omega', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * Omega', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * Omega', grad3_density_nutilde_m, 'grad3')
                                                


        ! R_11
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-11', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-11', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Reynolds-11', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-11', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-11', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-11', grad3_density_reynolds_m, 'grad3')
        
        ! R_22
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-22', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-22', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Reynolds-22', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-22', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-22', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-22', grad3_density_reynolds_m, 'grad3')
 
        ! R_33
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-33', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-33', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Reynolds-33', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-33', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-33', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-33', grad3_density_reynolds_m, 'grad3')
 
        ! R_12
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-12', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-12', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-12', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        call worker%store_bc_state('Density * Reynolds-12', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-12', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-12', grad3_density_reynolds_m, 'grad3')
 

        ! R_13
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-13', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-13', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        call worker%store_bc_state('Density * Reynolds-13', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_reynolds_m = ZERO
        grad2_density_reynolds_m = ZERO
        grad3_density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-13', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-13', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-13', grad3_density_reynolds_m, 'grad3')
 

        ! R_23
        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_reynolds_m       = worker%get_field('Density * Reynolds-23', 'value', 'face interior')
        grad1_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad1', 'face interior')
        grad2_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad2', 'face interior')
        grad3_density_reynolds_m = worker%get_field('Density * Reynolds-23', 'grad3', 'face interior')



        !
        ! Store boundary condition state - Extrapolate
        !
        density_reynolds_m = ZERO
        call worker%store_bc_state('Density * Reynolds-23', density_reynolds_m,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        call worker%store_bc_state('Density * Reynolds-23', grad1_density_reynolds_m, 'grad1')
        call worker%store_bc_state('Density * Reynolds-23', grad2_density_reynolds_m, 'grad2')
        call worker%store_bc_state('Density * Reynolds-23', grad3_density_reynolds_m, 'grad3')
 
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_rstm_ssglrrw_symmetry_y
