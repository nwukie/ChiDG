module type_equation
    use mod_kinds,      only: rk,ik
    implicit none
    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------
    type, public :: equation_t
        character(len=20)  :: name
        integer(ik)        :: ind
    end type equation_t
    !***********************************************************************


end module type_equation
