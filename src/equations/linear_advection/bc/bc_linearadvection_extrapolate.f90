module bc_linearadvection_extrapolate
    use mod_kinds,          only: rk,ik

    use type_bc_operator,   only: bc_operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !-----------------------------------------------------------------------------------------
    type, public, extends(bc_operator_t) :: linearadvection_extrapolate_t



    contains

        procedure   :: init
        procedure   :: compute    !> bc implementation

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


        !
        ! Set operator equations
        !
        call self%set_equation("u")


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
    subroutine compute(self,worker,prop)
        class(linearadvection_extrapolate_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop






    end subroutine compute
    !*********************************************************************************************






end module bc_linearadvection_extrapolate
