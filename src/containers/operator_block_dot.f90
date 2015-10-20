module operator_block_dot
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    use type_blockvector,   only: blockvector_t

    implicit none

    


contains





    !> Compute vector-vector dot product from two blockvector_t types.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    function block_dot(a,b) result(res)
        type(blockvector_t),    intent(in)  :: a
        type(blockvector_t),    intent(in)  :: b

        real(rk)    :: res
        integer(ik) :: ielem

        res = ZERO

        !
        ! Compute vector dot-product
        !
        do ielem = 1,size(a%lvecs)

            res = res + dot_product(a%lvecs(ielem)%vec, b%lvecs(ielem)%vec)

        end do  ! ielem
        


    end function










end module operator_block_dot

