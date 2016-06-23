module operator_block_dot
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    use type_blockvector,   only: blockvector_t

    implicit none

    


contains





    !> Compute vector-vector dot product from two blockvector_t types.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
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
        do ielem = 1,size(a%vecs)

            res = res + dot_product(a%vecs(ielem)%vec, b%vecs(ielem)%vec)

        end do  ! ielem
        


    end function block_dot
    !*****************************************************************************










end module operator_block_dot

