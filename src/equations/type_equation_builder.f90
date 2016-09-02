module type_equation_builder
    use mod_string,         only: string_t
    use type_equation_set,  only: equation_set_t
    implicit none



    !>  This is a type that constructs an equation set from pre-specified blueprint.
    !!
    !!  init
    !!      - Shall be used to set the name of the builder when it is extended. That way 
    !!        it can be found in a list by it's string.
    !!
    !!  build
    !!      - Shall be used to allocate an equation_set_t and call its initialization routines,
    !!        add operators etc. for some standard equation set. For example, the euler equations.
    !!        
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
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




    !---------------------------------------------------------------------------------------------
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








    !>  Set the name of the equation builder
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
    !!
    !---------------------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(equation_builder_t),  intent(inout)   :: self
        character(len=*),           intent(in)      :: string

        call self%name%set(string)

    end subroutine set_name
    !*********************************************************************************************
    

    

    !>  Return the name of the equation builder
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
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
