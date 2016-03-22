module mod_eigenvalues
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ONE

    use mod_inv,    only: inv

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
        print*, 'eigenvalues DGEEV info ', info



    end subroutine eigenvalues
    !******************************************************************************












end module mod_eigenvalues
