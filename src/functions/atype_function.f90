module atype_function
    use mod_kinds,  only: rk,ik
    implicit none
    private

    type, public, abstract :: function_t

    contains
        procedure(order_interface), deferred :: order       !> Returns the order of the implemented function so the proper integration rule can be used later on
        procedure(calc_interface),  deferred :: calc        !> Elemental function definition
    end type function_t




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


end module atype_function
