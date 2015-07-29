module type_expansion
    use mod_kinds,   only: rk,ik
    implicit none
    private

    type, public :: expansion_t
        integer(ik)                             :: nterms
        real(rk), dimension(:), allocatable     :: vec  !>  Vector of all modes
        real(rk), dimension(:,:), pointer       :: mat  !>  Matrix alias of 'vec'
    contains
        procedure, public   :: init
        procedure, public   :: var

        final :: destructor
    end type expansion_t


contains

    subroutine init(self,nterms,neqns)
        class(expansion_t), intent(inout), target  :: self
        integer(ik),        intent(in)             :: nterms, neqns

        self%nterms = nterms
        allocate(self%vec(nterms*neqns))

        !> Initialize matrix pointer alias
        self%mat(1:nterms,1:neqns) => self%vec

        !> Initialize to 0
        self%vec = 0._rk
    end subroutine
    



    function var(self,ivar) result(modes_out)
        class(expansion_t), intent(inout)   :: self
        integer(ik),        intent(in)      :: ivar

        real(rk)    :: modes_out(self%nterms)
        modes_out = self%mat(:,ivar)
    end function






    subroutine destructor(self)
        type(expansion_t), intent(inout) :: self
        if (allocated(self%vec))  deallocate(self%vec)
    end subroutine

end module type_expansion
