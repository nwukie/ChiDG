module type_element_location
    use mod_kinds,  only: ik



    type, public :: element_location_t

        integer(ik) :: idomain
        integer(ik) :: ielement

    end type element_location_t




end module type_element_location
