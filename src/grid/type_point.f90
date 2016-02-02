module type_point
    use mod_kinds,      only: ik,rk
    use mod_constants,  only: ZERO
    implicit none

    private




    !>  Point data type, containing three spatial coordinates
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------------
    type, public :: point_t

        real(rk) :: c1_ = ZERO, &
                    c2_ = ZERO, &
                    c3_ = ZERO

    contains

        procedure :: set
        procedure :: x
        procedure :: y
        procedure :: z
        procedure :: xi
        procedure :: eta
        procedure :: zeta

    end type point_t
    !*************************************************************************************************************






contains




    !> Set coordinates for point
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] x First coordinate value
    !!  @param[in] y Second coordinate value
    !!  @param[in] z Third coordinate value
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine set(self, c1, c2, c3)
        class(point_t), intent(inout) :: self
        real(rk) :: c1,c2,c3

        self%c1_ = c1
        self%c2_ = c2
        self%c3_ = c3

    end subroutine set
    !*************************************************************************************************************
    





    !> Set first coordinate.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] x_in     First coordinate value
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine x(self,x_in)
        class(point_t), intent(inout) :: self
        real(rk) :: x_in

        self%c1_ = x_in

    end subroutine x
    !*************************************************************************************************************







    !> Set second coordinate.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] y_in     Second coordinate value
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine y(self,y_in)
        class(point_t), intent(inout) :: self
        real(rk) :: y_in

        self%c2_ = y_in

    end subroutine y
    !*************************************************************************************************************







    !> Set third coordinate
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] z_in     Third coordinate value
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine z(self,z_in)
        class(point_t), intent(inout) :: self
        real(rk) :: z_in

        self%c3_ = z_in

    end subroutine z
    !*************************************************************************************************************





    !> type-bound procedure for setting first coordinate
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in] xi First coordinate value
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine xi(self,xi_in)
        class(point_t), intent(inout) :: self
        real(rk) :: xi_in

        self%c1_ = xi_in
    end subroutine
    !*********************************************************************************************************





    !> type-bound procedure for setting second computational coordinate
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !! @param[in] eta Second coordinate value
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine eta(self,eta_in)
        class(point_t), intent(inout) :: self
        real(rk) :: eta_in

        self%c2_ = eta_in
    end subroutine
    !*********************************************************************************************************





    !> type-bound procedure for setting third computational coordinate
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !! @param[in] zeta Third coordinate value
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine zeta(self,zeta_in)
        class(point_t), intent(inout) :: self
        real(rk) :: zeta_in

        self%c3_ = zeta_in
    end subroutine
    !*********************************************************************************************************






end module type_point
