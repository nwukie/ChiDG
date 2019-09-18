!>  Module contains procedures for analyzing the eigensystem for a
!!  matrix.
!!
!!
!!  Procedures          Description
!!  --------------------------------------
!!  eigen_solver        Compute eigenvalues/eigenvectors of a matrix.
!!  eigenvalues         Compute the eigenvalues of a matrix.
!!
!!
!!  @author Nathan A. Wukie
!!  @date   1/17/2018
!!
!!
!---------------------------------------------------------------------
module mod_eigen
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use DNAD_D
    implicit none

    ! External procedures defined in LAPACK
    external DGETRI
    external DGETRF
    external DGECON
    external DLANGE
    external DGEEV
    real(rk) :: DLANGE



contains




    !>  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2018
    !!
    !!
    !-----------------------------------------------------------------
    subroutine eigen_solver(A)
        type(AD_D), intent(in), dimension(:,:)  :: A

        real(rk),   allocatable, dimension(:,:) :: Areal
        integer(ik)                             :: ierr

        !
        ! Allocate real storage and copy real entries from A
        !
        allocate(Areal(size(A,1),size(A,2)), stat=ierr)
        if (ierr /= 0) call AllocationError
        Areal = A%x_ad_

        
!        !
!        ! Now, we want to compute the derivatives of the eigensolver
!        ! with respect to the solution components being differentiated.
!        ! So, we need to perturb
!        !
!        do iperturb = 1,size(A(1,1)%x_ad_,1)
!
!
!        end do


        !
        !
        !
        ! call LAPACK eigensolver 


    end subroutine eigen_solver
    !*****************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!  TODO: Try to find URL for citation.
    !!
    !---------------------------------------------------------------------------
    subroutine eigenvalues(A,wr,wi)
        real(rk), dimension(:,:), intent(in)    :: A
        real(rk), dimension(:),   intent(inout) :: wr, wi
        real(rk), dimension(size(A,1),size(A,2)) :: Acopy

        real(rk), dimension(size(A,1)*4)   :: work    ! work array for LAPACK
        real(rk), dimension(size(A,1),size(A,2)) :: junk
        integer :: n, info, nwork



        !
        ! Store A in Acopy to prevent it from being overwritten by LAPACK
        !
        Acopy = A
        n = size(A,1)
        nwork = size(work,1)


        !
        ! Estimate eigenvalues
        !
        call DGEEV('N','N',n, Acopy, n, wr, wi, junk, n, junk, n, work, nwork, info)



    end subroutine eigenvalues
    !******************************************************************************






end module mod_eigen
