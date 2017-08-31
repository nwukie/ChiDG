module bc_state_mesh_motion_derivative
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    use mpi_f08,            only: mpi_comm
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: mesh_motion_derivative_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type mesh_motion_derivative_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(mesh_motion_derivative_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Mesh Motion Derivative')
        call self%set_family('Mesh Motion')


        !
        ! Add functions
        !
        call self%bcproperties%add('Derivative1','Required')         ! add StaticPressure
        call self%bcproperties%add('Derivative2','Required')         ! add StaticPressure
        call self%bcproperties%add('Derivative3','Required')         ! add StaticPressure


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
    subroutine compute_bc_state(self,worker,prop, bc_COMM)
        class(mesh_motion_derivative_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM

        ! Equation indices
        integer(ik)     :: iu

        type(AD_D), allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc



        !
        ! Get 'u' value from face interior to extrapolate
        !
        u_bc = worker%get_field('grid_displacement1', 'value', 'face interior')



        !
        ! Initialize derivative arrays
        !
        dudx_bc = ZERO*worker%get_field('grid_displacement1', 'grad1', 'face interior')
        dudy_bc = ZERO*dudx_bc
        dudz_bc = ZERO*dudx_bc







        !
        ! Get derivative value
        !
        dudx_bc = self%bcproperties%compute("Derivative1",worker%time(),worker%coords())








        call worker%store_bc_state('grid_displacement1', u_bc,    'value')
        call worker%store_bc_state('grid_displacement1', dudx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement1', dudy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement1', dudz_bc, 'grad3')




        !
        ! Get 'u' value from face interior to extrapolate
        !
        u_bc = worker%get_field('grid_displacement2', 'value', 'face interior')



        !
        ! Initialize derivative arrays
        !
        dudx_bc = ZERO*worker%get_field('grid_displacement2', 'grad1','face interior')
        dudy_bc = ZERO*dudx_bc
        dudz_bc = ZERO*dudx_bc







        !
        ! Get derivative value
        !
        dudx_bc = self%bcproperties%compute("Derivative2",worker%time(),worker%coords())








        call worker%store_bc_state('grid_displacement2', u_bc,    'value')
        call worker%store_bc_state('grid_displacement2', dudx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement2', dudy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement2', dudz_bc, 'grad3')




        !
        ! Get 'u' value from face interior to extrapolate
        !
        u_bc = worker%get_field('grid_displacement3', 'value', 'face interior')



        !
        ! Initialize derivative arrays
        !
        dudx_bc = ZERO*worker%get_field('grid_displacement3', 'grad1','face interior')
        dudy_bc = ZERO*dudx_bc
        dudz_bc = ZERO*dudx_bc







        !
        ! Get derivative value
        !
        dudx_bc = self%bcproperties%compute("Derivative3",worker%time(),worker%coords())








        call worker%store_bc_state('grid_displacement3', u_bc,    'value')
        call worker%store_bc_state('grid_displacement3', dudx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement3', dudy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement3', dudz_bc, 'grad3')




    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_mesh_motion_derivative
