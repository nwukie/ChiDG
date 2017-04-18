module bc_empty
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use mpi_f08,            only: mpi_comm
    implicit none
    


    !> No boundary condition 
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: empty_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type empty_t
    !-------------------------------------------------------------------------------------------




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(empty_t), intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("empty")


    end subroutine init
    !********************************************************************************











    !> Specialized compute routine for Euler Slip-Wall Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(empty_t),                 intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


    end subroutine compute_bc_state
    !********************************************************************************************






end module bc_empty
