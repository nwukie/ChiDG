module type_equation_builder
    use mod_string,         only: string_t
    use type_equation_set,  only: equation_set_t
    implicit none



    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    type, public, abstract :: equation_builder_t

        type(string_t)  :: name

    contains

        procedure(init_interface),  deferred    :: init
        procedure(build_interface), deferred    :: build

        procedure                               :: set_name
        procedure                               :: get_name

    end type equation_builder_t
    !*********************************************************************************************
    abstract interface
        subroutine init_interface(self)
            import equation_builder_t

            class(equation_builder_t),  intent(inout)   :: self
        end subroutine
    end interface

    abstract interface
        function build_interface(self,blueprint) result(eqnset)
            import equation_builder_t
            import equation_set_t

            class(equation_builder_t),  intent(in)  :: self
            character(len=*),           intent(in)  :: blueprint
            type(equation_set_t)    :: eqnset
        end function
    end interface
    !*********************************************************************************************






contains








    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(equation_builder_t),  intent(inout)   :: self
        character(len=*),           intent(in)      :: string

        call self%name%set(string)

    end subroutine set_name
    !*********************************************************************************************
    

    

    !>
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function get_name(self) result(builder_name)
        class(equation_builder_t),  intent(in)  :: self

        character(len=:), allocatable   :: builder_name

        builder_name = self%name%get()

    end function get_name
    !*********************************************************************************************










end module type_equation_builder
