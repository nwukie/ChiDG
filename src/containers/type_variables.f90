module type_variables
    use mod_types,   only: rk,ik
    implicit none
    private

    type, public :: variables_t
        real(kind=rk), dimension(:,:), allocatable  :: vals
    contains
        procedure :: init
        final :: destructor
    end type variables_t

contains

    subroutine init(self,nterms_sol,neqns)
        class(variables_t),  intent(inout)  :: self
        integer(kind=ik),   intent(in)      :: nterms_sol,neqns

        allocate(self%vals(nterms_sol,neqns))

        ! Initialize to 0
        self%vals = 0._rk
    end subroutine
    
    subroutine destructor(self)
        type(variables_t), intent(inout) :: self
        deallocate(self%vals)
    end subroutine

end module type_variables
