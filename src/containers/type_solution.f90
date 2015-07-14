module type_field
    use mod_kinds,   only: rk,ik
    implicit none
    private

    type, public :: foe;d_t
        integer(ik)                             :: nterms
        real(rk), dimension(:), allocatable     :: vec  !>  Vector of all modes
        real(rk), dimension(:,:), pointer       :: mat  !>  Matrix alias of 'vec'
    contains
        procedure, public   :: init
        procedure, public   :: var

        final :: destructor
    end type solution_t


contains

    subroutine init(self,nterms,neqns)
        class(solution_t), intent(inout), target  :: self
        integer(ik),        intent(in)             :: nterms, neqns

        self%nterms = nterms
        allocate(self%vec(nterms*neqns))

        !> Initialize matrix pointer alias
        self%mat(1:nterms,1:neqns) => self%vec

        ! Initialize to 0
        self%vec = 0._rk
    end subroutine
    



    function var(self,ivar) result(modes_out)
        class(solution_t), intent(inout)   :: self
        integer(ik),        intent(in)      :: ivar

        real(rk)    :: modes_out(self%nterms)
        modes_out = self%mat(:,ivar)
    end function






    subroutine destructor(self)
        type(solution_t), intent(inout) :: self
        if (allocated(self%vec))  deallocate(self%vec)
    end subroutine

end module type_solution
