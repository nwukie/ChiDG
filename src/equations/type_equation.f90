module type_equation
    use mod_kinds,      only: rk,ik

    implicit none
    private

    type, public :: equation_t
        character(len=20)  :: name
        integer(ik)        :: ind
!    contains
!        final :: destructor
    end type equation_t

contains
    
!    elemental subroutine destructor(self)
!        type(equation_t), intent(in) :: self
!    end subroutine

end module type_equation
