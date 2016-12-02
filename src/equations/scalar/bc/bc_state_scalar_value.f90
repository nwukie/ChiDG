module bc_state_scalar_value
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_point,         only: point_t
    use DNAD_D
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: scalar_value_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type scalar_value_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(scalar_value_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Scalar Value')
        call self%set_family('Scalar')


        !
        ! Add functions
        !
        call self%bcproperties%add('Value','Required')


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
        class(scalar_value_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: iu

        integer(ik)     :: igq
        real(rk)        :: time

        type(AD_D),     allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc
        type(point_t),  allocatable, dimension(:)   :: coords



        !
        ! Get equation index
        !
        iu = prop%get_primary_field_index("u")



        !
        ! Get u_m from face interior to initialize derivatives
        !
        u_bc    = worker%get_primary_field_face(iu, 'value', 'face interior')
        dudx_bc = worker%get_primary_field_face(iu, 'ddx',   'face interior')
        dudy_bc = worker%get_primary_field_face(iu, 'ddy',   'face interior')
        dudz_bc = worker%get_primary_field_face(iu, 'ddz',   'face interior')


        !
        ! Get derivative value from boundary condition parameter
        !
        coords = worker%coords()
        time   = worker%time()
        u_bc   = self%bcproperties%compute("Value",time,coords)



        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state(iu, u_bc, 'value')



        !
        ! Store boundary condition state, gradient
        !
        call worker%store_bc_state(iu, dudx_bc, 'ddx')
        call worker%store_bc_state(iu, dudy_bc, 'ddy')
        call worker%store_bc_state(iu, dudz_bc, 'ddz')


    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_scalar_value
