module mod_dgemv
    use mod_kinds,  only: rk, ik
    implicit none





    ! Import BLAS dot product
    external DAXPY
!    external DDOT
!    real(rdouble) :: DDOT





    !integer(ik),    parameter   :: block_size = 40
    !integer(ik),    parameter   :: block_size = 135
    integer(ik),    parameter   :: block_size = 320

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/26/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine chidg_matmul(A,x,y)
        real(rk),   intent(in)      :: A(:,:)
        real(rk),   intent(in)      :: x(:)
        real(rk),   intent(inout)   :: y(:)

        integer(ik) :: icol, nrows, ncols, i,j, ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8

        real(rk), dimension(block_size)    :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, restmp, &
                                             Atmp1, Atmp2, Atmp3, Atmp4, Atmp5, Atmp6, Atmp7, Atmp8

        nrows = size(A,1)
        ncols = size(A,2)


!        do icol = 1,ncols
!            call DAXPY(nrows, x(icol), A(:,icol),1,y,1)
!        end do


!        !do j = 1,ncols
!        !    do i = 1,nrows
!        do j = 1,block_size
!            do i = 1,block_size
!                y(i) = y(i) + x(j)*A(i,j)
!            end do
!        end do

        restmp = 0._rk

        do icol = 1,int(real(block_size)/real(8))
            ind1 = 1 + 8*(icol-1)
            ind2 = 2 + 8*(icol-1)
            ind3 = 3 + 8*(icol-1)
            ind4 = 4 + 8*(icol-1)
            ind5 = 5 + 8*(icol-1)
            ind6 = 6 + 8*(icol-1)
            ind7 = 7 + 8*(icol-1)
            ind8 = 8 + 8*(icol-1)

            Atmp1 = A(:,ind1)
            Atmp2 = A(:,ind2)
            Atmp3 = A(:,ind3)
            Atmp4 = A(:,ind4)
            Atmp5 = A(:,ind5)
            Atmp6 = A(:,ind6)
            Atmp7 = A(:,ind7)
            Atmp8 = A(:,ind8)

            tmp1 = x(ind1)*Atmp1
            tmp2 = x(ind2)*Atmp2
            tmp3 = x(ind3)*Atmp3
            tmp4 = x(ind4)*Atmp4
            tmp5 = x(ind5)*Atmp5
            tmp6 = x(ind6)*Atmp6
            tmp7 = x(ind7)*Atmp7
            tmp8 = x(ind8)*Atmp8

            restmp = restmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8
        end do

        y = y + restmp

    end subroutine chidg_matmul
    !*****************************************************************************************




end module mod_dgemv
