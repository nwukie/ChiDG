module type_test
    use mod_kinds,      only: ik,rk
    use mod_constants,  only: ZERO

    implicit none
    private

    !> Point data type, containing three spatial coordinates
    !! @author Nathan A. Wukie
    type, public :: test_t

        real(rk) :: c1_ = ZERO, &
                    c2_ = ZERO, &
                    c3_ = ZERO

        real(rk), allocatable   :: comp(:)
    contains

        final :: destructor
    end type test_t

contains

    !> Default destructor for point_t
    subroutine destructor(self)
        type(test_t), intent(inout) :: self

        if (allocated(self%comp)) deallocate(self%comp)
    end subroutine

end module type_test
