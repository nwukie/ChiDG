module type_equation
    use mod_kinds,      only: rk,ik
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    type, public :: equation_t
        character(len=:), allocatable   :: name
        integer(ik)                     :: ind

    contains
        procedure   :: set_name
        procedure   :: set_index
    end type equation_t
    !***********************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(equation_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        self%name = string

    end subroutine set_name
    !***********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_index(self,ind)
        class(equation_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        self%ind = ind

    end subroutine set_index
    !***********************************************************************************






end module type_equation
