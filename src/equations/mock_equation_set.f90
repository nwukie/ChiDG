module mock_equation_set
    use type_equation_set,  only: equation_set_t
    
    type, extends(equation_set_t) :: mock_equation_set_t

    contains
        procedure :: init
    end type

contains

    subroutine init(self)
        class(mock_equation_set_t),  intent(inout)   :: self
    end subroutine

end module mock_equation_set
