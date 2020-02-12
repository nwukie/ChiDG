module bc_state_scalar_derivative
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    use ieee_arithmetic
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: scalar_derivative_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type scalar_derivative_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(scalar_derivative_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Scalar Derivative')
        call self%set_family('Scalar')


        !
        ! Add functions
        !
        call self%bcproperties%add('Normal Gradient','Required')         ! add StaticPressure


        !
        ! Add parameters
        !


    end subroutine init
    !******************************************************************************************












    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/15/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(scalar_derivative_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM


        type(AD_D), allocatable, dimension(:)   :: u_bc, grad1u_bc, grad2u_bc, grad3u_bc, grad1u_t, grad2u_t, grad3u_t, normal_gradient

        ! Get 'u' value from face interior to extrapolate
        u_bc      = worker%get_field('u', 'value', 'face interior')
        grad1u_bc = worker%get_field('u', 'grad1', 'face interior')
        grad2u_bc = worker%get_field('u', 'grad2', 'face interior')
        grad3u_bc = worker%get_field('u', 'grad3', 'face interior')


        ! Initialize derivative arrays

        ! Normal gradient by taking advantage of the fact that only the normal component
        ! is felt in the flux anyways since it is dotted with the normal vector. So, just
        ! compose the boundary gradient from the normal component of the gradient.
        normal_gradient = self%bcproperties%compute("Normal Gradient",worker%time(),worker%coords())
        grad1u_bc = normal_gradient*worker%unit_normal(1)
        grad2u_bc = normal_gradient*worker%unit_normal(2)
        grad3u_bc = normal_gradient*worker%unit_normal(3)



!        ! Normal gradient by extrapolating tangential and adding back normal
!        normal_gradient = self%bcproperties%compute("Normal Gradient",worker%time(),worker%coords())
!
!
!        ! Subtract the normal component of the gradient to leave the tangential part
!        grad1u_t = grad1u_bc - (grad1u_bc*worker%unit_normal(1))*worker%unit_normal(1)
!        grad2u_t = grad2u_bc - (grad2u_bc*worker%unit_normal(2))*worker%unit_normal(2)
!        grad3u_t = grad3u_bc - (grad3u_bc*worker%unit_normal(3))*worker%unit_normal(3)
!
!
!        ! Add back the prescribed normal component of the gradient to the tangential part to get the total
!        grad1u_bc = grad1u_t + normal_gradient*worker%unit_normal(1)
!        grad2u_bc = grad2u_t + normal_gradient*worker%unit_normal(2)
!        grad3u_bc = grad3u_t + normal_gradient*worker%unit_normal(3)


        ! Store
        call worker%store_bc_state('u', u_bc,      'value')
        call worker%store_bc_state('u', grad1u_bc, 'grad1')
        call worker%store_bc_state('u', grad2u_bc, 'grad2')
        call worker%store_bc_state('u', grad3u_bc, 'grad3')


    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_scalar_derivative
