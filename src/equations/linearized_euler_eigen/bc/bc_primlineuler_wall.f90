module bc_primlineuler_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO

    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    implicit none


    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: primlineuler_wall_t

    contains
    
        procedure   :: init
        procedure   :: compute_bc_state

    end type primlineuler_wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self)
        class(primlineuler_wall_t),    intent(inout)   :: self

        !
        ! Set name
        ! 
        call self%set_name('primlineuler_wall')
        call self%set_family('Wall')


    end subroutine init
    !*******************************************************************************************









    !>  Wall for linearized euler.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(primlineuler_wall_t),     intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            rho_r,      u_r,     v_r,     w_r,     p_r,                 &
            rho_bc,     u_bc,    v_bc,    w_bc,    p_bc,                &
            grad1_rho_m, grad1_u_m, grad1_v_m, grad1_w_m, grad1_p_m,    &
            grad2_rho_m, grad2_u_m, grad2_v_m, grad2_w_m, grad2_p_m,    &
            grad3_rho_m, grad3_u_m, grad3_v_m, grad3_w_m, grad3_p_m



        !
        ! Interpolate interior solution to quadrature nodes
        !
        rho_r = worker%get_field('Density',    'value', 'face interior')
        u_r   = worker%get_field('Velocity-1', 'value', 'face interior')
        v_r   = worker%get_field('Velocity-2', 'value', 'face interior')
        w_r   = worker%get_field('Velocity-3', 'value', 'face interior')
        p_r   = worker%get_field('Pressure',   'value', 'face interior')



        !
        ! Interpolate gradients
        !
        grad1_rho_m = worker%get_field('Density',  'grad1', 'face interior')
        grad2_rho_m = worker%get_field('Density',  'grad2', 'face interior')
        grad3_rho_m = worker%get_field('Density',  'grad3', 'face interior')

        grad1_u_m = worker%get_field('Velocity-1', 'grad1', 'face interior')
        grad2_u_m = worker%get_field('Velocity-1', 'grad2', 'face interior')
        grad3_u_m = worker%get_field('Velocity-1', 'grad3', 'face interior')

        grad1_v_m = worker%get_field('Velocity-2', 'grad1', 'face interior')
        grad2_v_m = worker%get_field('Velocity-2', 'grad2', 'face interior')
        grad3_v_m = worker%get_field('Velocity-2', 'grad3', 'face interior')

        grad1_w_m = worker%get_field('Velocity-3', 'grad1', 'face interior')
        grad2_w_m = worker%get_field('Velocity-3', 'grad2', 'face interior')
        grad3_w_m = worker%get_field('Velocity-3', 'grad3', 'face interior')

        grad1_p_m = worker%get_field('Pressure',   'grad1', 'face interior')
        grad2_p_m = worker%get_field('Pressure',   'grad2', 'face interior')
        grad3_p_m = worker%get_field('Pressure',   'grad3', 'face interior')



        !
        ! Initialize storage
        !
        rho_bc = rho_r
        u_bc   = rho_r
        v_bc   = rho_r
        w_bc   = rho_r
        p_bc   = rho_r


        !
        ! Set boundary state
        !
        rho_bc = rho_r
        u_bc   = ZERO
        v_bc   = ZERO
        w_bc   = ZERO
        p_bc   = p_r



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , rho_bc, 'value')
        call worker%store_bc_state('Velocity-1', u_bc,   'value')
        call worker%store_bc_state('Velocity-2', v_bc,   'value')
        call worker%store_bc_state('Velocity-3', w_bc,   'value')
        call worker%store_bc_state('Pressure'  , p_bc,   'value')





        grad1_rho_m = ZERO
        grad2_rho_m = ZERO
        grad3_rho_m = ZERO
        call worker%store_bc_state('Density'   , grad1_rho_m, 'grad1')
        call worker%store_bc_state('Density'   , grad2_rho_m, 'grad2')
        call worker%store_bc_state('Density'   , grad3_rho_m, 'grad3')
                                                
        call worker%store_bc_state('Velocity-1', grad1_u_m,   'grad1')
        call worker%store_bc_state('Velocity-1', grad2_u_m,   'grad2')
        call worker%store_bc_state('Velocity-1', grad3_u_m,   'grad3')
                                                
        call worker%store_bc_state('Velocity-2', grad1_v_m,   'grad1')
        call worker%store_bc_state('Velocity-2', grad2_v_m,   'grad2')
        call worker%store_bc_state('Velocity-2', grad3_v_m,   'grad3')
                                                
        call worker%store_bc_state('Velocity-3', grad1_w_m,   'grad1')
        call worker%store_bc_state('Velocity-3', grad2_w_m,   'grad2')
        call worker%store_bc_state('Velocity-3', grad3_w_m,   'grad3')

        grad1_p_m = ZERO
        grad2_p_m = ZERO
        grad3_p_m = ZERO
        call worker%store_bc_state('Pressure',   grad1_p_m,   'grad1')
        call worker%store_bc_state('Pressure',   grad2_p_m,   'grad2')
        call worker%store_bc_state('Pressure',   grad3_p_m,   'grad3')






    end subroutine compute_bc_state
    !*********************************************************************************************************






end module bc_primlineuler_wall
