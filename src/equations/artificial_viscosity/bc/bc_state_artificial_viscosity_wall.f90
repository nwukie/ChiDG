module bc_state_artificial_viscosity_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: artificial_viscosity_wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type artificial_viscosity_wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_wall_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name('Artificial Viscosity Wall')
        call self%set_family('Wall')


    end subroutine init
    !********************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(artificial_viscosity_wall_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            eps_m,      eps_bc,                     &
            deps_dx_m,  deps_dy_m,  deps_dz_m




        !
        ! Interpolate interior solution to quadrature nodes
        !
        eps_m     = worker%get_field('Artificial Viscosity' ,'value', 'face interior')
        deps_dx_m = worker%get_field('Artificial Viscosity' ,'grad1', 'face interior')
        deps_dy_m = worker%get_field('Artificial Viscosity' ,'grad2', 'face interior')
        deps_dz_m = worker%get_field('Artificial Viscosity' ,'grad3', 'face interior')


        !
        ! Initialize arrays
        !
        eps_bc = eps_m


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Artificial Viscosity' ,eps_bc, 'value')


        deps_dx_m = ZERO
        deps_dy_m = ZERO
        deps_dz_m = ZERO
        call worker%store_bc_state('Artificial Viscosity' ,deps_dx_m, 'grad1')
        call worker%store_bc_state('Artificial Viscosity' ,deps_dy_m, 'grad2')
        call worker%store_bc_state('Artificial Viscosity' ,deps_dz_m, 'grad3')
                                                


    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_artificial_viscosity_wall
