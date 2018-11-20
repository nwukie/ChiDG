module precon_identity
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_matrix,       only: chidg_matrix_t
    use type_chidg_vector,       only: chidg_vector_t

    use mod_inv,    only: inv


    !> Identity preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_identity_t


    contains
        procedure   :: update
        procedure   :: apply

    end type precon_identity_t




contains


    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_identity_t), intent(inout)   :: self
        !type(blockmatrix_t),    intent(in)      :: A
        !type(blockvector_t),    intent(in)      :: b
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: b



        
    end subroutine update








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    function apply(self,A,v,z_old) result(z)
        class(precon_identity_t), intent(inout)   :: self
        type(chidg_matrix_t),      intent(in)      :: A
        type(chidg_vector_t),      intent(in)      :: v
        type(chidg_vector_t),      intent(in), optional :: z_old

        !type(blockvector_t) :: z
        type(chidg_vector_t) :: z


        !
        ! Identity preconditioner. Do nothing, copy incoming vector to outgoing vector without modification
        !
        z = v



    end function apply













end module precon_identity
