module bc_linearadvection_extrapolate
    use mod_kinds,          only: rk,ik
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !--------------------------------------------------
    type, public, extends(bc_t) :: linearadvection_extrapolate_t



    contains
        procedure :: compute    !> bc implementation
    end type linearadvection_extrapolate_t




contains

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
    !----------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk)
        class(linearadvection_extrapolate_t),   intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        integer(ik),                            intent(in)      :: idom
        integer(ik),                            intent(in)      :: ielem
        integer(ik),                            intent(in)      :: iface
        integer(ik),                            intent(in)      :: iblk













    end subroutine






end module bc_linearadvection_extrapolate
