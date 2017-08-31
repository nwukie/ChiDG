mdule bc_state_mesh_motion_traction
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_point,         only: point_t
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
    type, public, extends(bc_state_t) :: mesh_motion_traction_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type mesh_motion_traction_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(mesh_motion_traction_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Mesh Motion Traction')
        call self%set_family('Mesh Motion')


        !
        ! Add functions
        !
        call self%bcproperties%add('Strain11','Required')
        call self%bcproperties%add('Strain22','Required')
        call self%bcproperties%add('Strain33','Required')
        call self%bcproperties%add('Strain12','Required')
        call self%bcproperties%add('Strain23','Required')
        call self%bcproperties%add('Strain31','Required')


        call self%set_fcn_option('Strain11','val',ZERO)
        call self%set_fcn_option('Strain22','val',ZERO)
        call self%set_fcn_option('Strain33','val',ZERO)
        call self%set_fcn_option('Strain12','val',ZERO)
        call self%set_fcn_option('Strain23','val',ZERO)
        call self%set_fcn_option('Strain32','val',ZERO)
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
    subroutine compute_bc_state(self,worker,prop)
        class(mesh_motion_traction_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: iu1, iu2, iu3

        integer(ik)     :: igq
        real(rk)        :: time

        type(AD_D),     allocatable, dimension(:)   :: u_bc, u1_bc, du1dx_bc, du1dy_bc, du1dz_bc,  &
                                                        u2_bc, du2dx_bc, du2dy_bc, du2dz_bc,  &
                                                        u3_bc, du3dx_bc, du3dy_bc, du3dz_bc

        real(rk),     allocatable, dimension(:)   :: u_input
        type(point_t),  allocatable, dimension(:)   :: coords



        !
        ! Get equation index
        !
        !iu1 = prop%get_primary_field_index("grid_displacement1")
        !iu2 = prop%get_primary_field_index("grid_displacement2")
        !iu3 = prop%get_primary_field_index("grid_displacement3")



        !
        ! Get u_m from face interior to initialize derivatives
        !
        u1_bc    = worker%get_field('grid_displacement1','value', 'face interior')
        du1dx_bc = worker%get_field('grid_displacement1','grad1', 'face interior')
        du1dy_bc = worker%get_field('grid_displacement1','grad2', 'face interior')
        du1dz_bc = worker%get_field('grid_displacement1','grad3', 'face interior')

        u2_bc    = worker%get_field('grid_displacement2','value', 'face interior')
        du2dx_bc = worker%get_field('grid_displacement2','grad1', 'face interior')
        du2dy_bc = worker%get_field('grid_displacement2','grad2', 'face interior')
        du2dz_bc = worker%get_field('grid_displacement2','grad3', 'face interior')

        u3_bc    = worker%get_field('grid_displacement3','value', 'face interior')
        du3dx_bc = worker%get_field('grid_displacement3','grad1', 'face interior')
        du3dy_bc = worker%get_field('grid_displacement3','grad2', 'face interior')
        du3dz_bc = worker%get_field('grid_displacement3','grad3', 'face interior')


        !
        ! Get derivative value from boundary condition parameter
        !
        coords = worker%coords()
        time   = worker%time()

        u_bc = u1_bc
        u_bc = ZERO
        !GD1
        u_input   = self%bcproperties%compute("Value1",time,coords)
        u_bc   = u_input 

        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state('grid_displacement1', u_bc, 'value')



        !
        ! Store boundary condition state, gradient
        !
        call worker%store_bc_state('grid_displacement1', du1dx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement1', du1dy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement1', du1dz_bc, 'grad3')


        !GD2

        u_input   = self%bcproperties%compute("Value2",time,coords)
        u_bc   = u_input 




        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state('grid_displacement2', u_bc, 'value')



        !
        ! Store boundary condition state, gradient
        !
        call worker%store_bc_state('grid_displacement2', du2dx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement2', du2dy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement2', du2dz_bc, 'grad3')

        !GD3


        u_input   = self%bcproperties%compute("Value3",time,coords)
        u_bc   = u_input 



        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state('grid_displacement3', u_bc, 'value')



        !
        ! Store boundary condition state, gradient
        !
        call worker%store_bc_state('grid_displacement3', du3dx_bc, 'grad1')
        call worker%store_bc_state('grid_displacement3', du3dy_bc, 'grad2')
        call worker%store_bc_state('grid_displacement3', du3dz_bc, 'grad3')


    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_mesh_motion_traction
