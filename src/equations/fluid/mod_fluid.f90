module mod_fluid
#include <messenger.h>
    use mod_kinds,      only: rk
    use mod_constants,  only: ONE
    implicit none

    ! Gas Constants
    real(rk), parameter :: Rgas = 287.15_rk
    real(rk), parameter :: gam  = 1.4_rk
    real(rk), parameter :: cp   = Rgas*gam/(gam-ONE)


    !real(rk)    :: omega = 366.5191_rk
    real(rk)    :: omega = 0._rk




end module mod_fluid
