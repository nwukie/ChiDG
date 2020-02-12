module type_field
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    type, public :: field_t
        character(:),   allocatable :: name
        integer(ik)                 :: cache_index
        integer(ik)                 :: functional_ID ! Functional ID associate with the adjoint field
    contains
        procedure   :: set_name
        procedure   :: get_name
        procedure   :: set_cache_index
        procedure   :: get_cache_index
        procedure   :: set_functional_ID
        procedure   :: get_functional_ID
    end type field_t
    !***********************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(field_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        ! Check for blank string
        if (trim(string) == '') call chidg_signal(FATAL,"field%set_name: empty string.")

        self%name = string

    end subroutine set_name
    !***********************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !-----------------------------------------------------------------------------------
    function get_name(self) result(name_)
        class(field_t), intent(in)  :: self

        character(:),   allocatable :: name_

        if (allocated(self%name)) then
            name_ = self%name
        else
            call chidg_signal(FATAL,"field%get_name: The field has not been set with a name.")
        end if


    end function get_name
    !***********************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_cache_index(self,cache_index)
        class(field_t), intent(inout)   :: self
        integer(ik),    intent(in)      :: cache_index

        self%cache_index = cache_index

    end subroutine set_cache_index
    !***********************************************************************************


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !-----------------------------------------------------------------------------------
    function get_cache_index(self) result(cache_index)
        class(field_t), intent(in)  :: self

        integer(ik) :: cache_index

        cache_index = self%cache_index

    end function get_cache_index
    !***********************************************************************************



    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_functional_ID(self,func_id)
        class(field_t), intent(inout)   :: self
        integer(ik),    intent(in)      :: func_id

        self%functional_ID = func_id

    end subroutine set_functional_ID
    !***********************************************************************************



    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !!
    !-----------------------------------------------------------------------------------
    function get_functional_ID(self) result(func_id)
        class(field_t), intent(in)  :: self

        integer(ik) :: func_id

        func_id = self%functional_ID

    end function get_functional_ID
    !***********************************************************************************



end module type_field
