module type_preconditioner
    use mod_kinds,          only: rk, ik
    use type_chidg_matrix,   only: chidg_matrix_t
    use type_chidg_vector,   only: chidg_vector_t
    use type_chidg_data,    only: chidg_data_t
    use type_timer,         only: timer_t

    implicit none




    !> An abstract preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------
    type, public, abstract :: preconditioner_t

        type(timer_t)       :: timer
        logical             :: initialized = .false.

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply


    end type preconditioner_t


contains




    !> Default empty initialization routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(preconditioner_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data




    end subroutine init
    !**************************************************************************************









    !> This routine computes the preconditioner so that it can be applied
    !! in the form of a matrix-vector product.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(preconditioner_t),    intent(inout)   :: self
        type(chidg_matrix_t),        intent(in)      :: A
        type(chidg_vector_t),        intent(in)      :: b

    end subroutine update
    !************************************************************************************







    !> This routine applies the preconditioner in the form of the matrix-vector product
    !! z = Minv*v
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(preconditioner_t),    intent(inout)   :: self
        type(chidg_matrix_t),        intent(in)      :: A
        type(chidg_vector_t),        intent(in)      :: v

        type(chidg_vector_t)     :: z


        !
        ! Default - no preconditioning
        !
        z = v

    end function apply
    !************************************************************************************













end module type_preconditioner
