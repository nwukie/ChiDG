module bc_state_HP_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    implicit none


    !>
    !!  
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   11/09/2018
    !!
    !-----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: HP_wall_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type HP_wall_t
    !******************************************************************************************




contains



    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(HP_wall_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('HP Wall')
        call self%set_family('Scalar')


    end subroutine init
    !******************************************************************************************












    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      eqnset  Equation Set type governing the current domain
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !---------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(HP_wall_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM

        type(AD_D), allocatable, dimension(:)   ::  &
            u_bc, p_bc, q_bc, r_bc,                 &
            grad1_u_bc, grad2_u_bc, grad3_u_bc,     &
            grad1_p_bc, grad2_p_bc, grad3_p_bc,     &
            grad1_q_bc, grad2_q_bc, grad3_q_bc,     &
            grad1_r_bc, grad2_r_bc, grad3_r_bc
                    

        ! Get u_m from face interior to initialize derivatives
        u_bc       = worker%get_field('u','value', 'face interior')
        grad1_u_bc = worker%get_field('u','grad1', 'face interior')
        grad2_u_bc = worker%get_field('u','grad2', 'face interior')
        grad3_u_bc = worker%get_field('u','grad3', 'face interior')

        p_bc       = worker%get_field('p','value', 'face interior')
        grad1_p_bc = worker%get_field('p','grad1', 'face interior')
        grad2_p_bc = worker%get_field('p','grad2', 'face interior')
        grad3_p_bc = worker%get_field('p','grad3', 'face interior')

        q_bc       = worker%get_field('q','value', 'face interior')
        grad1_q_bc = worker%get_field('q','grad1', 'face interior')
        grad2_q_bc = worker%get_field('q','grad2', 'face interior')
        grad3_q_bc = worker%get_field('q','grad3', 'face interior')

        r_bc       = worker%get_field('r','value', 'face interior')
        grad1_r_bc = worker%get_field('r','grad1', 'face interior')
        grad2_r_bc = worker%get_field('r','grad2', 'face interior')
        grad3_r_bc = worker%get_field('r','grad3', 'face interior')


        ! ZERO dirichlet condition on walls
        u_bc = ZERO


        ! Store boundary condition state, Value
        call worker%store_bc_state('u', u_bc, 'value')
        call worker%store_bc_state('p', p_bc, 'value')
        call worker%store_bc_state('q', q_bc, 'value')
        call worker%store_bc_state('r', r_bc, 'value')


        ! Store boundary condition state, gradient
        call worker%store_bc_state('u', grad1_u_bc, 'grad1')
        call worker%store_bc_state('u', grad2_u_bc, 'grad2')
        call worker%store_bc_state('u', grad3_u_bc, 'grad3')

        call worker%store_bc_state('p', grad1_p_bc, 'grad1')
        call worker%store_bc_state('p', grad2_p_bc, 'grad2')
        call worker%store_bc_state('p', grad3_p_bc, 'grad3')

        call worker%store_bc_state('q', grad1_q_bc, 'grad1')
        call worker%store_bc_state('q', grad2_q_bc, 'grad2')
        call worker%store_bc_state('q', grad3_q_bc, 'grad3')

        call worker%store_bc_state('r', grad1_r_bc, 'grad1')
        call worker%store_bc_state('r', grad2_r_bc, 'grad2')
        call worker%store_bc_state('r', grad3_r_bc, 'grad3')



    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_HP_wall
