module type_r3
    use mod_kinds,      only: rk, ik
    implicit none


    type, public :: r3_t

        integer(ik)         :: ID
        real(rk)            :: data(3)

    end type r3_t


end module type_r3
