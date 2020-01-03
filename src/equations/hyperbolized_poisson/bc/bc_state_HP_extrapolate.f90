module bc_state_HP_extrapolate
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
    type, public, extends(bc_state_t) :: HP_extrapolate_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type HP_extrapolate_t
    !******************************************************************************************




contains



    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(HP_extrapolate_t),  intent(inout)   :: self

        call self%set_name('HP Extrapolate')
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
        class(HP_extrapolate_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM

        type(AD_D), allocatable, dimension(:)   ::  &
            u_bc, p_bc, q_bc, r_bc,                 &
            grad1_u_bc, grad2_u_bc, grad3_u_bc,     &
            grad1_p_bc, grad2_p_bc, grad3_p_bc,     &
            grad1_q_bc, grad2_q_bc, grad3_q_bc,     &
            grad1_r_bc, grad2_r_bc, grad3_r_bc, gradu_normal, gradu_normal_old, gradu_normal_new, gradu_tang, &
            pt_bc, qt_bc, rt_bc, mag, unorm_1, unorm_2, unorm_3
                    


        !
        ! Get u_m from face interior to initialize derivatives
        !
        u_bc       = worker%get_field('u','value', 'face interior')
        p_bc       = worker%get_field('p','value', 'face interior')
        q_bc       = worker%get_field('q','value', 'face interior')
        r_bc       = worker%get_field('r','value', 'face interior')

        grad1_u_bc = worker%get_field('u','grad1', 'face interior')
        grad2_u_bc = worker%get_field('u','grad2', 'face interior')
        grad3_u_bc = worker%get_field('u','grad3', 'face interior')

        grad1_p_bc = worker%get_field('p','grad1', 'face interior')
        grad2_p_bc = worker%get_field('p','grad2', 'face interior')
        grad3_p_bc = worker%get_field('p','grad3', 'face interior')

        grad1_q_bc = worker%get_field('q','grad1', 'face interior')
        grad2_q_bc = worker%get_field('q','grad2', 'face interior')
        grad3_q_bc = worker%get_field('q','grad3', 'face interior')

        grad1_r_bc = worker%get_field('r','grad1', 'face interior')
        grad2_r_bc = worker%get_field('r','grad2', 'face interior')
        grad3_r_bc = worker%get_field('r','grad3', 'face interior')


! Homogeneous Neumann
        ! Get normal component of gradient
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)
        gradu_normal = p_bc*unorm_1 + q_bc*unorm_2 + r_bc*unorm_3

        ! Subtract normal component
        p_bc = p_bc - gradu_normal*unorm_1
        q_bc = q_bc - gradu_normal*unorm_2
        r_bc = r_bc - gradu_normal*unorm_3


!!! Extrapolate tangential gradient
!        ! Get normal component of gradient
!        unorm_1 = worker%unit_normal(1)
!        unorm_2 = worker%unit_normal(2)
!        unorm_3 = worker%unit_normal(3)
!        gradu_normal_old = p_bc*unorm_1 + q_bc*unorm_2 + r_bc*unorm_3
!
!
!        ! Subtract from total to get tangential
!        pt_bc = p_bc - gradu_normal_old*unorm_1
!        qt_bc = q_bc - gradu_normal_old*unorm_2
!        rt_bc = r_bc - gradu_normal_old*unorm_3
!
!
!        ! Magnitude of tangential
!        gradu_tang = sqrt(pt_bc*pt_bc + qt_bc*qt_bc + rt_bc*rt_bc + 1.e-6_rk)
!
!
!        ! If we want the magnitude of the gradient to be 1, and we are extrapolating the tangential part, scale the normal part
!        ! so that the total equals 1
!        gradu_normal_new = ONE - gradu_tang
!
!
!
!        ! Add scaled normal component
!!        p_bc = p_bc*(gradu_normal_new/gradu_normal_old)
!!        q_bc = q_bc*(gradu_normal_new/gradu_normal_old)
!!        r_bc = r_bc*(gradu_normal_new/gradu_normal_old)
!
!        p_bc = pt_bc + gradu_normal_new*worker%unit_normal(1)
!        q_bc = qt_bc + gradu_normal_new*worker%unit_normal(2)
!        r_bc = rt_bc + gradu_normal_new*worker%unit_normal(3)
!
!
!
!
!! Scaled magnitude
!        ! Get normal component of gradient
!        unorm_1 = worker%unit_normal(1)
!        unorm_2 = worker%unit_normal(2)
!        unorm_3 = worker%unit_normal(3)
!        gradu_normal_old = p_bc*unorm_1 + q_bc*unorm_2 + r_bc*unorm_3
!
!        mag = sqrt(p_bc*p_bc + q_bc*q_bc + r_bc*r_bc + 1.e-6_rk)
!
!        ! Subtract from total to get tangential
!        p_bc = p_bc/mag
!        q_bc = q_bc/mag
!        r_bc = r_bc/mag




        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state('u', u_bc, 'value')
        call worker%store_bc_state('p', p_bc, 'value')
        call worker%store_bc_state('q', q_bc, 'value')
        call worker%store_bc_state('r', r_bc, 'value')


        !
        ! Store boundary condition state, gradient
        !
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






end module bc_state_HP_extrapolate
