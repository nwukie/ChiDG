module type_blockjacobi
    use mod_kinds,          only: rk, ik


    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------
    type, public :: blockjacobi_t



    contains


    end type blockjacobi_t







contains


    subroutine solver(self,A,x,b)
       class(blockjacobi_t),    intent(in)      :: self
       type(blockmatrix_t),     intent(in)      :: A
       type(blockvector_t),     intent(inout)   :: x
       type(blockvector_t),     intent(in)      :: b 





    end subroutine solver



end module
