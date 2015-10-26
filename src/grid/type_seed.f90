module type_seed
    use mod_kinds,  only: ik



    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------
    type, public :: seed_t

        integer(ik) :: idom
        integer(ik) :: ielem

    end type seed_t








end module type_seed
