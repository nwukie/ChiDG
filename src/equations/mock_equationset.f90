module mock_equationset
    use type_equationset,  only: equationset_t
    
    type, extends(equationset_t) :: mock_equationset_t

    contains
        procedure :: init
    end type

contains

    subroutine init(self)
        class(mock_equationset_t),  intent(inout)   :: self
    end subroutine

end module
