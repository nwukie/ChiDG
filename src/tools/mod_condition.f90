module mod_condition
    use mod_types, only: rk,ik

    use mod_inv

    implicit none

contains

    subroutine condition_number(dim1,dim2,mat,rcond)
        integer(kind=ik), intent(in)                        :: dim1,dim2
        real(kind=rk),    intent(in), dimension(dim1,dim2)  :: mat


        real(kind=rk)                                   :: mat_norm,invmat_norm,rcond
        integer(kind=ik)                                :: i, info, j, n, idim1,idim2
        real(kind=rk),      dimension(dim1,dim2)        :: invmat
        real(kind=rk),      dimension(dim2)             :: vector_sum_mat, vector_sum_invmat

        invmat = inv(mat)

        ! compute the 1-norm of each matrix
        do idim2 = 1,dim2
            vector_sum_mat(idim2)    = sum(abs(mat(:,idim2)))
            vector_sum_invmat(idim2) = sum(abs(invmat(:,idim2)))
        end do

        mat_norm    = maxval(vector_sum_mat)
        invmat_norm = maxval(vector_sum_invmat)

        ! compute the matrix condition number
        rcond = mat_norm*invmat_norm

    end subroutine condition_number

end module mod_condition
