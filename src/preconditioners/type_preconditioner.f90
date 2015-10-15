module type_preconditioner
    use mod_kinds,          only: rk, ik
    use type_domain,        only: domain_t
    use type_blockmatrix,   only: blockmatrix_t
    use type_blockvector,   only: blockvector_t
    use type_chidgData,     only: chidgData_t

    implicit none




    !> An abstract preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------
    type, public, abstract :: preconditioner_t


        logical  :: initialized = .false.

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
        type(chidgData_t),          intent(inout)   :: data




    end subroutine init









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
        type(blockmatrix_t),        intent(in)      :: A
        type(blockvector_t),        intent(in)      :: b


    end subroutine update







    !> This routine applies the preconditioner in the form of the matrix-vector product
    !! z = Minv*v
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(preconditioner_t),    intent(inout)   :: self
        type(blockmatrix_t),        intent(in)      :: A
        type(blockvector_t),        intent(in)      :: v

        type(blockvector_t)     :: z


        !
        ! Default - no preconditioning
        !
        z = v

    end function













end module type_preconditioner
