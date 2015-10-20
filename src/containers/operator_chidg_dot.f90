module operator_chidg_dot
    use mod_kinds,          only: rk, ik
    
    use type_chidgVector,   only: chidgVector_t
    use operator_block_dot, only: block_dot
    implicit none


contains



    !> Compute vector-vector dot product from two chidgVector_t types.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    function dot(a,b) result(res)
        type(chidgVector_t),    intent(in)  :: a
        type(chidgVector_t),    intent(in)  :: b

        real(rk)    :: res
        integer(ik) :: idom, ielem

        res = ZERO

        !
        ! Compute vector dot-product
        !
        do idom = 1,size(self%dom)
            
            res = res + block_dot(a%dom(idom),b%dom(idom))

        end do




    end function



end module operator_chidg_dot
