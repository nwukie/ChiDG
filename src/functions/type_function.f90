module type_function
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use type_dict,  only: dict_t
    implicit none
    private


    !>  Function class.
    !!
    !!  \f$     f = f(t,\vec{x})    \f$
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, abstract :: function_t

        character(len=:),   allocatable             :: name
        type(dict_t)                                :: dict

    contains

        procedure                                   :: set          !< Set function value
        procedure(init_interface),      deferred    :: init         !< Initialize function and register options
        procedure(compute_interface),   deferred    :: compute      !< Elemental function definition

        !procedure                                   :: list_options !< Return self%dict contents
    end type function_t
    !********************************************************************************************





    abstract interface
        impure elemental function compute_interface(self,time,coord)
            use type_point, only: point_t
            use mod_kinds,  only: rk
            import function_t

            class(function_t),  intent(inout)   :: self
            real(rk),           intent(in)      :: time
            type(point_t),      intent(in)      :: coord
            real(rk)                            :: compute_interface
        end function
    end interface




    abstract interface
        subroutine init_interface(self)
            import function_t

            class(function_t),  intent(inout)  :: self

        end subroutine
    end interface




contains




    !> Procedure for setting function parameters
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  key     String indicating a dictionary entry to set
    !!  @parma[in]  val     Real value to be associated with the key
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set(self,key,val)
        class(function_t),  intent(inout)    :: self
        character(*),       intent(in)       :: key
        real(rk),           intent(in)       :: val
        
        logical     :: key_exists
        real(rk)    :: test

        !
        ! Check that key exists in the dictionary as a valid option
        !
        key_exists = self%dict%contains(key)



        !
        ! Set dictionary key/value pair
        !
        if (key_exists) then

            call self%dict%set(key,val)

        else


            call write_line("Invalid function option '"//key//"' in "//self%name//"%set('"//key//"',",val,')', delimiter='')
            call write_line(' ')
            call write_line('Valid function options are:')
            call write_line(' ')
            call write_line('key','current value')
            call write_line('----------------------------------------------')
            call self%dict%print()

            call chidg_signal_one(FATAL,"function%set: key is not a valid option for the current function",key)

        end if

    end subroutine set
    !*********************************************************************************************














end module type_function
