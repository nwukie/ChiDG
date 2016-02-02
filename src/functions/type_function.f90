module type_function
    use mod_kinds,  only: rk, ik
    use type_dict,  only: dict_t
    implicit none
    private


    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, abstract :: function_t


    contains

        procedure(order_interface), deferred :: order       !< Returns the order of the function
        procedure(calc_interface),  deferred :: calc        !< Elemental function definition
        procedure                            :: set         !< Set function value

    end type function_t
    !********************************************************************************************




    abstract interface
        function order_interface(self)
            use mod_kinds,  only: ik
            import function_t

            class(function_t),  intent(in)  :: self
            integer(ik)                     :: order_interface
        end function
    end interface



    abstract interface
        elemental function calc_interface(self,pt)
            use type_point, only: point_t
            use mod_kinds,  only: rk
            import function_t

            class(function_t),  intent(in)  :: self
            type(point_t),      intent(in)  :: pt
            real(rk)                        :: calc_interface
        end function
    end interface








contains





!    !>
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------------------------
!    subroutine function_options



    !> Procedure for setting function parameters
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set(self,valstring,val)
        class(function_t),  intent(inout)    :: self
        character(*),       intent(in)       :: valstring
        real(rk),           intent(in)       :: val
        

    end subroutine set
    !*********************************************************************************************














end module type_function
