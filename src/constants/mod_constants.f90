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


    !
    ! Floating-point mathematical constants
    !
    real(rk), parameter     :: PI       = 3.14159265358979323846_rk     !< Mathematical constant pi


    !
    ! Floating-point numbers
    !
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
    real(rk), parameter     :: TEN      = 10._rk



    !
    ! Floating-point fractions
    !
    real(rk), parameter     :: HALF     = 0.5_rk
    real(rk), parameter     :: THIRD    = 1._rk/3._rk
    real(rk), parameter     :: FOURTH   = 0.25_rk
    real(rk), parameter     :: FIFTH    = 1._rk/5._rk
    real(rk), parameter     :: SIXTH    = 0.2_rk
    real(rk), parameter     :: SEVENTH  = 1._rk/7._rk
    real(rk), parameter     :: EIGHTH   = 1._rk/8._rk
    real(rk), parameter     :: NINTH    = 1._rk/9._rk
    real(rk), parameter     :: TENTH    = 0.1_rk

    real(rk), parameter     :: TWO_THR  = 2._rk/3._rk
    real(rk), parameter     :: THR_FIV  = 3._rk/5._rk



    !
    ! Constants for spatial definition
    !
    integer(ik), parameter :: SPACEDIM  = 3
    integer(ik), parameter :: NFACES    = 6



    !
    ! Spatial direction
    !
    integer(ik), parameter :: X_DIR     = 1
    integer(ik), parameter :: y_DIR     = 2
    integer(ik), parameter :: Z_DIR     = 3

    integer(ik), parameter :: XI_DIR    = 1
    integer(ik), parameter :: ETA_DIR   = 2
    integer(ik), parameter :: ZETA_DIR  = 3




    !
    ! Indexing for faces directions and jacobian blocks
    !
    integer(ik), parameter :: XI_MIN    = 1
    integer(ik), parameter :: XI_MAX    = 2
    integer(ik), parameter :: ETA_MIN   = 3
    integer(ik), parameter :: ETA_MAX   = 4
    integer(ik), parameter :: ZETA_MIN  = 5
    integer(ik), parameter :: ZETA_MAX  = 6
    integer(ik), parameter :: DIAG      = 7




    !
    ! Face types. These should be distinct from the above 'face directions'
    !
    integer(ik), parameter :: INTERIOR  =  0     ! interior face
    integer(ik), parameter :: BOUNDARY  = -1     ! boundary condition type
    integer(ik), parameter :: CHIMERA   = -2     ! Chimera face
    integer(ik), parameter :: ORPHAN    = -3     ! orphan face - has no identity. Every face needs an identity.




    !
    ! Interpolation type
    !
    integer(ik), parameter :: LOCAL     = 0         ! Interpolate from current element
    integer(ik), parameter :: NEIGHBOR  = 1         ! Interpolate from neighbor element





    !
    ! QUADRATURE CONSTANTS: specify number of quadrature orders to initialize
    !
    integer(ik), parameter :: NGQ = 1





    !
    ! IO constants
    !
    integer(ik),   parameter :: MAXBLOCKS      = 200
    character(4),  parameter :: IO_DESTINATION = 'both'







end module mod_constants
