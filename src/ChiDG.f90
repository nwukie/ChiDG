!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie


program chidg
    use mod_kinds,              only: ik,rk
    use cgns
    use hdf5_dgtools
    use messenger

    !-----------------------------------------------------------
    !-----------------------------------------------------------
    ! Variable declarations
    implicit none

    !======================================================
    !
    !   Initialization
    !
    !======================================================
    call plot3d_to_hdf5()



    !======================================================
    !
    !   Solve
    !
    !======================================================







    !=======================================================
    !
    !   Post-execution tasks
    !
    !=======================================================

!    print*, "Completed in: ", elapsed, " Seconds"
!    print*, "======================="
!    print*, "Done!"
!    print*, "======================="

end program chidg
