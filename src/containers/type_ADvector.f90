module type_ADvector
    use mod_kinds,   only: rk,ik,rdouble,ilong
    use DNAD_D
    implicit none
    private

    type, public :: ADvector_t
        type(AD_D), dimension(:), allocatable  :: vals
    contains
        procedure :: init
        final :: destructor
    end type ADvector_t


    public assignment (=)
    interface assignment (=)
        module procedure assign_real_ADvec  ! ADvec = real(:)
        module procedure assign_long_ADvec  ! ADvec = integer(:)
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_ADvec_ADvec   ! ADvec + ADvec
        module procedure add_ADvec_REAL    ! ADvec + REAL
        module procedure add_ADvec_INTEGER ! ADvec + INTEGER
        module procedure add_REAL_ADvec    ! REAL + ADvec
        module procedure add_INTEGER_ADvec ! INTEGER + ADvec
    end interface




contains

    !> Methods
    !------------------------------------------------------------
    subroutine init(self,nreal,nderiv)
        class(ADvector_t),   intent(inout)  :: self
        integer(ik),         intent(in)     :: nreal, nderiv

        integer(ik)     :: ival

        ! Allocate number of AD values
        allocate(self%vals(nreal))


        ! Loop through and allocate space for derivatives
        do ival = 1,nreal
            self%vals(ival) = AD_D(nderiv)
        end do


        ! Initialize to 0
        self%vals = 0._rk
    end subroutine

    



    !-----------------------------------------------------------
    !> Operators
    !-----------------------------------------------------------
    !> ASSIGNMENT
    !
    !> ADvec(:) = real(:)
    pure subroutine assign_real_ADvec(u,r)
        type(ADvector_t), intent(inout) :: u
        real(rk),         intent(in)    :: r(:)

        u%vals = r
    end subroutine assign_real_ADvec

    !> ADvec(:) = long_int(:)
    pure subroutine assign_long_ADvec(u,i)
        type(ADvector_t), intent(inout) :: u
        integer(ilong),   intent(in)    :: i(:)

        u%vals = i
    end subroutine assign_long_ADvec

    !> ADDITION
    !
    !> ADvec = ADvec + ADvec
    pure function add_ADvec_ADvec(u,v) result(res)
        type(ADvector_t), intent(in) :: u,v
        type(ADvector_t) :: res

        res%vals = u%vals + v%vals
    end function

    !> ADvec = ADvec + REAL(:)
    pure function add_ADvec_REAL(u,v) result(res)
        type(ADvector_t), intent(in) :: u
        real(rdouble),    intent(in) :: v(:)
        type(ADvector_t) :: res

        res%vals = u%vals + v
    end function

    !> ADvec = REAL(:) + ADvec
    pure function add_REAL_ADvec(v,u) result(res)
        real(rdouble),    intent(in) :: v(:)
        type(ADvector_t), intent(in) :: u
        type(ADvector_t) :: res

        res%vals = u%vals + v
    end function

    !> ADvec = ADvec + INTEGER(:)
    pure function add_ADvec_INTEGER(u,v) result(res)
        type(ADvector_t), intent(in) :: u
        integer(ilong),   intent(in) :: v(:)
        type(ADvector_t) :: res

        res%vals = u%vals + v
    end function

    !> ADvec = INTEGER(:) + ADvec
    pure function add_INTEGER_ADvec(v,u) result(res)
        integer(ilong),   intent(in) :: v(:)
        type(ADvector_t), intent(in) :: u
        type(ADvector_t) :: res

        res%vals = u%vals + v
    end function









    !> Destructor
    !-----------------------------------------------------------
    subroutine destructor(self)
        type(ADvector_t), intent(inout) :: self

        ! Could this cause memory leak on allocated derivatives?
        deallocate(self%vals)
    end subroutine

end module type_ADvector
