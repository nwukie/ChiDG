!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie
!!
!!



program dummy_driver
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, ONE, TWO, THREE, FOUR, FIVE, ZERO
    
    !
    ! Variable declarations
    !
    implicit none
    integer(ik)                         :: nterms_c, idom

end program dummy_driver
