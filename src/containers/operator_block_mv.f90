module operator_block_mv
    use mod_kinds,          only: rk, ik
    use type_blockmatrix,   only: blockmatrix_t
    use type_blockvector,   only: blockvector_t

    implicit none



    public operator(*)
    interface operator(*)
        module procedure MULT_blockmatrix_blockvector
    end interface


contains




    !> This function implements the extremely important matrix-vector multiplication
    !! operation Ax
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------
    function MULT_blockmatrix_blockvector(A,x) result(res)
        type(blockmatrix_t), intent(in) :: A
        type(blockvector_t), intent(in) :: x 

        type(blockvector_t) :: res
        integer(ik)         :: ielem, iblk, iparent

        res = x
        call res%clear


        !
        ! Compute Ax
        !
        do ielem = 1,size(A%lblks,1)
            do iblk = 1,size(A%lblks,2)

                if (allocated(A%lblks(ielem,iblk)%mat)) then
                    iparent = A%lblks(ielem,iblk)%eparent()
                    res%lvecs(ielem)%vec = res%lvecs(ielem)%vec + matmul(A%lblks(ielem,iblk)%mat,x%lvecs(iparent)%vec)
                end if

            end do  ! iblk
        end do  ! ielem




    end function MULT_blockmatrix_blockvector








end module operator_block_mv
