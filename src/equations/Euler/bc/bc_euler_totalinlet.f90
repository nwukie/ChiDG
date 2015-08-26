module bc_euler_totalinlet
    use mod_kinds,          only: rk,ik
    use atype_bc,           only: bc_t
    use atype_solverdata,   only: solverdata_t
    use type_mesh,          only: mesh_t
    use atype_equationset,  only: equationset_t


    !> Total condition inlet boundary condition
    !!      - Specify total pressure(pt) and total temperature(tt) 
    !!        at the inlet boundary.
    !!  
    !!  @author Nathan A. Wukie
    !--------------------------------------------------
    type, public, extends(bc_t) :: euler_totalinlet_t

        real(rk)    :: pt   !> Inlet total pressure
        real(rk)    :: tt   !> Inlet total temperature




    contains
        procedure :: compute    !> bc implementation
    end type euler_totalinlet_t




contains

    !> Specialized compute routine for Inlet Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      eqnset  Equation Set type governing the current domain
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !----------------------------------------------------------
    subroutine compute(self,eqnset,mesh,sdata,ielem,iface,iblk)
        class(euler_totalinlet_t),    intent(inout)   :: self
        class(equationset_t),   intent(in)      :: eqnset
        type(mesh_t),           intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: iblk



















    end subroutine






end module bc_euler_totalinlet
