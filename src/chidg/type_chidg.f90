module type_chidg
    use mod_equations,  only: initialize_equations
    use mod_grid,       only: initialize_grid


    implicit none


    type, public    :: chidg_t


    contains
        procedure   :: init
    end type



contains


    subroutine init(self)
        class(chidg_t),  intent(inout)   :: self

        call initialize_equations()
        call initialize_grid()

    end subroutine






end module type_chidg
