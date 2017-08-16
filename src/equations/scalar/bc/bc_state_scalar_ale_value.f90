module bc_state_scalar_ale_value
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_point,         only: point_t
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
    type, public, extends(bc_state_t) :: scalar_ale_value_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type scalar_ale_value_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(scalar_ale_value_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('Scalar ALE Value')
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
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(scalar_ale_value_t),      intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM

        integer(ik)     :: igq
        real(rk)        :: time

        type(AD_D),     allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc
        type(point_t),  allocatable, dimension(:)   :: coords


        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid

        logical, allocatable, dimension(:)  :: inflow, outflow


        type(AD_D), allocatable, dimension(:,:)   ::  &
            gradu


        !
        ! Get u_m from face interior to initialize derivatives
        !
        u_bc  = worker%get_primary_field_value_ale_face('u', 'face interior')
        gradu = worker%get_primary_field_grad_ale_face('u', 'gradient', 'face interior')
        
        dudx_bc = gradu(:,1) 
        dudy_bc = gradu(:,2)
        dudz_bc = gradu(:,3)


        !
        ! Get derivative value from boundary condition parameter
        !
        coords = worker%coords()
        time   = worker%time()
        u_bc   = self%bcproperties%compute("Value",time,coords)




        !
        ! Store boundary condition state, Value
        !
        call worker%store_bc_state('u', u_bc, 'value')



        !
        ! Store boundary condition state, gradient
        !
        call worker%store_bc_state('u', dudx_bc, 'grad1')
        call worker%store_bc_state('u', dudy_bc, 'grad2')
        call worker%store_bc_state('u', dudz_bc, 'grad3')


    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_scalar_ale_value
