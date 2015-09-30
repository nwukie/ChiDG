module atype_preconditioner
    use mod_kinds,  only: rk, ik

    implicit none




    !> An abstract preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------
    type, public, abstract :: preconditioner_t

        type(blockmatrix_t)     :: P
        type(blockmatrix_t)     :: Pinv


    contains
        procedure   :: apply


    end type preconditioner_t


contains






    subroutine apply(self,A,x)













end module atype_preconditioner
