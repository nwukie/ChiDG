module bc_state_linearadvection_extrapolate
    use mod_kinds,          only: rk,ik

    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !-----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: linearadvection_extrapolate_t



    contains

        procedure   :: init
        procedure   :: compute_bc_state    !> bc implementation

    end type linearadvection_extrapolate_t
    !******************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)    
        class(linearadvection_extrapolate_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('linearadvection_extrapolate')


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
        class(linearadvection_extrapolate_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop






    end subroutine compute_bc_state
    !*********************************************************************************************






end module bc_state_linearadvection_extrapolate
