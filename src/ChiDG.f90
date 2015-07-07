!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie


program chidg
    use mod_kinds,          only: ik,rk
    use mod_hdfio,              only: read_grid_hdf5
    use type_domain,            only: domain_t
    use messenger

    !-----------------------------------------------------------
    !-----------------------------------------------------------
    ! Variable declarations
    implicit none
    type(domain_t), allocatable :: domains(:)

    !======================================================
    !
    !   Initialization
    !
    !======================================================




    call read_grid_hdf5('test',domains)



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
