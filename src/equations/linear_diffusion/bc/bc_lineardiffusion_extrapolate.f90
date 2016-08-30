module bc_lineardiffusion_extrapolate
    use mod_kinds,          only: rk,ik
    use type_bc,            only: bc_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none



    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: lineardiffusion_extrapolate_t



    contains

        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type lineardiffusion_extrapolate_t
    !****************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/16/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(lineardiffusion_extrapolate_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('lineardiffusion_extrapolate')


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
    subroutine compute(self,worker,prop)
        class(lineardiffusion_extrapolate_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop






    end subroutine compute
    !*********************************************************************************************






end module bc_lineardiffusion_extrapolate
