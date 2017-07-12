module mod_fluid
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE
    implicit none


    ! Fluid equation constants
    real(rk), parameter :: Rgas  = 287.15_rk            ! Specific gas constant
    real(rk), parameter :: gam   = 1.4_rk               ! Ratio of specific heats
    real(rk), parameter :: cp    = Rgas*gam/(gam-ONE)   ! Specific heat capacity, constant pressure
    real(rk)            :: omega = 0._rk                ! Rotation rate (rad/s)


end module mod_fluid
