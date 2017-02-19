module mod_time_HB
    use mod_kinds,  only: rk,ik
    implicit none

    real(rk), allocatable        :: D(:,:)
    logical                      :: harmonic_balance = .false.


contains



    !>  Called in type_HB, this subroutine copies D computed there and makes it
    !!  available to the code
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   2/17/2017 (Matteo's nonna's birhtday: 89! Woot! Woot!)
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_pseudo_spectral_operator(D_orig)
        real(rk),   intent(in)          :: D_orig(:,:)


        D = D_orig

        harmonic_balance = .true.


    end subroutine get_pseudo_spectral_operator
    !***********************************************************************************




















end module mod_time_HB
