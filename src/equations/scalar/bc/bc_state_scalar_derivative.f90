module bc_state_scalar_derivative
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
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
        call self%bcproperties%add('Derivative','Required')         ! add StaticPressure


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
    subroutine compute_bc_state(self,worker,prop)
        class(scalar_derivative_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: iu

        type(AD_D), allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc


        !
        ! Get equation index
        !
        iu = prop%get_primary_field_index("u")



        !
        ! Get 'u' value from face interior to extrapolate
        !
        u_bc = worker%get_primary_field_face('u', 'value', 'face interior')



        !
        ! Initialize derivative arrays
        !
        dudx_bc = ZERO*worker%get_primary_field_face('u', 'grad1','face interior')
        dudy_bc = ZERO*dudx_bc
        dudz_bc = ZERO*dudx_bc







        !
        ! Get derivative value
        !
        dudx_bc = self%bcproperties%compute("Derivative",worker%time(),worker%coords())








        call worker%store_bc_state('u', u_bc,    'value')
        call worker%store_bc_state('u', dudx_bc, 'grad1')
        call worker%store_bc_state('u', dudy_bc, 'grad2')
        call worker%store_bc_state('u', dudz_bc, 'grad3')






    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_scalar_derivative
