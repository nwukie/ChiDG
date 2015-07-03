module atype_fluid
    implicit none
    private

    type, public, abstract :: fluid_t
        private


    contains
        procedure, deferred :: compute_pressure
        procedure, deferred :: compute_gamma
        procedure, deferred :: compute_cp
        procedure, deferred :: compute_rgas


        final :: destructor
    end type fluid_t

contains
    
    subroutine destructor(self)
        type(atype_fluid), intent(in) :: self
    end subroutine

end module atype_fluid
