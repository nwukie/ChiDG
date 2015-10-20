module type_equationset_wrapper
    use type_equationset,   only: equationset_t


    type, public :: equationset_wrapper_t

        class(equationset_t), allocatable   :: item

    end type equationset_wrapper_t


end module type_equationset_wrapper
