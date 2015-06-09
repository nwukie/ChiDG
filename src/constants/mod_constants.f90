!> Module containing general constants that can be used consistently
!! throughout the program.
!!
!! These include:
!!      - Mathematical constants (pi, e, etc.)
!!      - Floating-point values for numbers and fractions
!!          (to ensure the correct floating point precision is used consistently)
!!      - Integers for consistent indexing

module mod_constants
    use mod_kinds, only: rk,ik,ishort
    implicit none

    !> Mathematical constants
    real(rk), parameter     :: PI       = 3.14159265358979323846_rk

    !> Floating-point numbers
    real(rk), parameter     :: ZERO     = 0._rk
    real(rk), parameter     :: ONE      = 1._rk
    real(rk), parameter     :: TWO      = 2._rk
    real(rk), parameter     :: THREE    = 3._rk
    real(rk), parameter     :: FOUR     = 4._rk
    real(rk), parameter     :: FIVE     = 5._rk
    real(rk), parameter     :: SIX      = 6._rk
    real(rk), parameter     :: SEVEN    = 7._rk
    real(rk), parameter     :: EIGHT    = 8._rk
    real(rk), parameter     :: NINE     = 9._rk


    !> Floating-point fractions
    real(rk), parameter     :: HALF     = 1._rk/2._rk
    real(rk), parameter     :: THIRD    = 1._rk/3._rk
    real(rk), parameter     :: FOURTH   = 1._rk/4._rk
    real(rk), parameter     :: FIFTH    = 1._rk/5._rk
    real(rk), parameter     :: SIXTH    = 1._rk/6._rk
    real(rk), parameter     :: SEVENTH  = 1._rk/7._rk
    real(rk), parameter     :: EIGHTH   = 1._rk/8._rk
    real(rk), parameter     :: NINTH    = 1._rk/9._rk

    real(rk), parameter     :: TWO_THR  = 2._rk/3._rk
    real(rk), parameter     :: THR_FIF  = 3._rk/5._rk

    !> Indexing for faces directions and jacobian blocks
    integer(ishort), parameter :: XI_MIN    = 1
    integer(ishort), parameter :: XI_MAX    = 2
    integer(ishort), parameter :: ETA_MIN   = 3
    integer(ishort), parameter :: ETA_MAX   = 4
    integer(ishort), parameter :: ZETA_MIN  = 5
    integer(ishort), parameter :: ZETA_MAX  = 6
    integer(ishort), parameter :: DIAG      = 7

    !> Constants for spatial definition
    integer(ishort), parameter :: SPACEDIM  = 3
    integer(ishort), parameter :: NFACES    = 6

    !> Spatial direction
    integer(ishort), parameter :: XI_DIR    = 1
    integer(ishort), parameter :: ETA_DIR   = 2
    integer(ishort), parameter :: ZETA_DIR  = 3


    ! QUADRATURE CONSTANTS
    ! specify number of quadrature orders to initialize
    integer(ishort), parameter :: NGQ = 1

    ! INPUT/OUTPUT CONSTANTS
    integer(ishort), parameter :: RES_MESH_OUT = 5
    integer(ishort), parameter :: RES_SOLN_OUT = 5
    integer(ishort), parameter :: MAXBLOCKS    = 200







end module mod_constants
