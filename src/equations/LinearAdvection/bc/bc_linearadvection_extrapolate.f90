module bc_linearadvection_extrapolate
    use mod_kinds,          only: rk,ik
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t

    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !--------------------------------------------------
    type, public, extends(bc_t) :: linearadvection_extrapolate_t



    contains

        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type linearadvection_extrapolate_t




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(linearadvection_extrapolate_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('linearadvection_extrapolate')


        !
        ! Add functions
        !


        !
        ! Add parameters
        !


    end subroutine add_options
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
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(linearadvection_extrapolate_t),   intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        type(face_info_t),                      intent(in)      :: face
        type(function_info_t),                  intent(in)      :: fcn






    end subroutine compute
    !*********************************************************************************************






end module bc_linearadvection_extrapolate
