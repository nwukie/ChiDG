module type_expansion
    use mod_kinds,   only: rk,ik
    implicit none
    private

    type, public :: expansion_t
        integer(ik)                             :: neqns
        integer(ik)                             :: nterms
        real(rk), dimension(:), allocatable     :: vec              !>  Vector of all modes
        real(rk), dimension(:,:), pointer       :: mat => null()    !>  Matrix alias of 'vec'
    contains
        procedure, public   :: init
        procedure, public   :: var

        final :: destructor
    end type expansion_t






    !-------------------    OPERATORS   ---------------------
    public operator (*)
    interface operator (*)
        module procedure mult_real_exp    ! real * expansion, ELEMENTAL
        module procedure mult_exp_real    ! expansion * real, ELEMENTAL
    end interface


    public operator (/)
    interface operator (/)
        module procedure div_real_exp    ! real / expansion, ELEMENTAL
        module procedure div_exp_real    ! expansion / real, ELEMENTAL
    end interface





    public operator (+)
    interface operator (+)
        module procedure add_exp_exp    ! expansion + expansion, ELEMENTAL
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_exp_exp    ! expansion - expansion, ELEMENTAL
    end interface








contains

    subroutine init(self,nterms,neqns)
        class(expansion_t), intent(inout), target  :: self
        integer(ik),        intent(in)             :: nterms, neqns

        self%neqns  = neqns
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















    !---------------------------    OPERATORS    ----------------------------
    elemental function mult_real_exp(left,right) result(res)
        real(rk),           intent(in)  :: left
        type(expansion_t),  intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns   = right%neqns
        res%nterms  = right%nterms

        res%vec = left * right%vec
        res%mat(1:res%nterms,1:res%neqns) => res%vec

    end function

    elemental function mult_exp_real(left,right) result(res)
        type(expansion_t),  intent(in)  :: left
        real(rk),           intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns   = left%neqns
        res%nterms  = left%nterms

        res%vec = left%vec * right
        res%mat(1:res%nterms,1:res%neqns) => res%vec

    end function






    elemental function div_real_exp(left,right) result(res)
        real(rk),           intent(in)  :: left
        type(expansion_t),  intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns   = right%neqns
        res%nterms  = right%nterms

        res%vec = left / right%vec
        res%mat(1:res%nterms,1:res%neqns) => res%vec

    end function

    elemental function div_exp_real(left,right) result(res)
        type(expansion_t),  intent(in)  :: left
        real(rk),           intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns   = left%neqns
        res%nterms  = left%nterms

        res%vec = left%vec / right
        res%mat(1:res%nterms,1:res%neqns) => res%vec

    end function












    elemental function add_exp_exp(left,right) result(res)
        type(expansion_t),  intent(in)  :: left
        type(expansion_t),  intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns = left%neqns
        res%nterms = left%nterms

        res%vec = left%vec + right%vec
        res%mat(1:res%nterms,1:res%neqns) => res%vec
    end function





    elemental function sub_exp_exp(left,right) result(res)
        type(expansion_t),  intent(in)  :: left
        type(expansion_t),  intent(in)  :: right
        type(expansion_t), target   :: res

        res%neqns = left%neqns
        res%nterms = left%nterms

        res%vec = left%vec - right%vec
        res%mat(1:res%nterms,1:res%neqns) => res%vec
    end function







    subroutine destructor(self)
        type(expansion_t), intent(inout) :: self
        if (allocated(self%vec))  deallocate(self%vec)
    end subroutine

end module type_expansion
