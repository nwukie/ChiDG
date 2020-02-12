module type_point_ad
#include <messenger.h>
    use mod_kinds,      only: ik,rk
    use mod_constants,  only: ZERO, VALID_POINT, INVALID_POINT
    use DNAD_D
    implicit none

    private




    !>  Point data type, containing three spatial coordinates
    !!  and derivatives
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: point_ad_t

        integer(ik) :: status = 0   ! Status code. Could be used to indicate if the point is
                                    ! valid. For example, if the chimera routine was searching
                                    ! for a donor point, this could indicate if the point
                                    ! returned was actually included in the element, or if the
                                    ! routine failed. Other than that, non-consequential.

!        real(rk)    :: c1_ = ZERO, &
!                       c2_ = ZERO, &
!                       c3_ = ZERO

        type(AD_D)  :: c1_ , &
                       c2_ , &
                       c3_ 
    contains

        generic     :: set => set_from_rk,set_from_ad
        procedure   :: set_from_rk
        procedure   :: set_from_ad
        procedure   :: x
        procedure   :: y
        procedure   :: z
        procedure   :: xi
        procedure   :: eta
        procedure   :: zeta

        procedure   :: add_x
        procedure   :: add_y
        procedure   :: add_z

        procedure   :: valid

    end type point_ad_t
    !******************************************************************************************




    !
    ! Constructors
    !
    interface point_ad_t
        module procedure construct_point_from_ad
        module procedure construct_points_from_ads
        module procedure construct_point_from_rk
        module procedure construct_points_from_rks
    end interface


contains




    !>  Construct a single point_ad_t from a Rank1 AD_D array of size 3.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!
    !---------------------------------------------------------------------------------
    function construct_point_from_ad(point_ad) result(point_)
        type(AD_D),   intent(in)  :: point_ad(:)

        type(point_ad_t)   :: point_

        point_%c1_ = point_ad(1)
        point_%c2_ = point_ad(2)
        point_%c3_ = point_ad(3)

    end function construct_point_from_ad
    !*********************************************************************************



    !>  Construct an array of point_ad_t's from a Rank2 AD_D array of size (npts,3)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !---------------------------------------------------------------------------------
    function construct_points_from_ads(points_ad) result(points_)
        type(AD_D),   intent(in)  :: points_ad(:,:)

        integer(ik)                 :: ierr, ipoint
        type(point_ad_t),  allocatable :: points_(:)

        allocate(points_(size(points_ad,1)),stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipoint = 1,size(points_ad,1)
            points_(ipoint)%c1_ = points_ad(ipoint,1)
            points_(ipoint)%c2_ = points_ad(ipoint,2)
            points_(ipoint)%c3_ = points_ad(ipoint,3)
        end do

    end function construct_points_from_ads
    !*********************************************************************************





    !>  Construct a single point_ad_t from a Rank1 rk array of size 3.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/11/2018
    !!
    !---------------------------------------------------------------------------------
    function construct_point_from_rk(point_rk) result(point_)
        real(rk),   intent(in)  :: point_rk(:)

        type(point_ad_t)   :: point_

        call point_%set(point_rk(1),point_rk(2),point_rk(3))

    end function construct_point_from_rk
    !*********************************************************************************



    !>  Construct an array of point_ad_t's from a Rank2 rk array of size (npts,3)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/1/2018
    !!
    !---------------------------------------------------------------------------------
    function construct_points_from_rks(points_rk) result(points_)
        real(rk),   intent(in)  :: points_rk(:,:)

        integer(ik)                 :: ierr, ipoint
        type(point_ad_t),  allocatable :: points_(:)

        allocate(points_(size(points_rk,1)),stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipoint = 1,size(points_rk,1)
            call points_(ipoint)%set(points_rk(ipoint,1),points_rk(ipoint,2),points_rk(ipoint,3))
        end do

    end function construct_points_from_rks
    !*********************************************************************************









    !> Set coordinates for point, given AD_D coordiantes
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] x First coordinate value
    !!  @param[in] y Second coordinate value
    !!  @param[in] z Third coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_from_ad(self, c1, c2, c3)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: c1,c2,c3

        self%c1_ = c1
        self%c2_ = c2
        self%c3_ = c3

    end subroutine set_from_ad
    !******************************************************************************************







    !>  Set coordinates for point, given real coordinates
    !!  If so, it means that the user is not interested in derivatives.
    !!  Therefore, initialize derivatives with 0.
    !!
    !!  WARNING: make sure to use only the real part and not access the derivatives
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] x First coordinate value
    !!  @param[in] y Second coordinate value
    !!  @param[in] z Third coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_from_rk(self, c1_rk, c2_rk, c3_rk)
        class(point_ad_t), intent(inout) :: self
        real(rk),          intent(in)    :: c1_rk, c2_rk, c3_rk
        
        type(AD_D)   :: c1_ad, c2_ad, c3_ad
        
        ! 
        ! Initialize derivatives
        !
        c1_ad = AD_D(0)
        c2_ad = AD_D(0)
        c3_ad = AD_D(0)

        !
        ! Assign real parts
        !
        c1_ad = c1_rk 
        c2_ad = c2_rk 
        c3_ad = c3_rk 


        self%c1_ = c1_ad
        self%c2_ = c2_ad
        self%c3_ = c3_ad

    end subroutine set_from_rk
    !******************************************************************************************
    





    !> Set first coordinate.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] x_in     First coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine x(self,x_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: x_in

        self%c1_ = x_in

    end subroutine x
    !******************************************************************************************







    !> Set second coordinate.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] y_in     Second coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine y(self,y_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: y_in

        self%c2_ = y_in

    end subroutine y
    !******************************************************************************************







    !> Set third coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] z_in     Third coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine z(self,z_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: z_in

        self%c3_ = z_in

    end subroutine z
    !******************************************************************************************





    !> type-bound procedure for setting first coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !!  @param[in] xi First coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine xi(self,xi_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: xi_in

        self%c1_ = xi_in
    end subroutine
    !******************************************************************************************





    !> type-bound procedure for setting second computational coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !! @param[in] eta Second coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine eta(self,eta_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: eta_in

        self%c2_ = eta_in
    end subroutine
    !******************************************************************************************





    !> type-bound procedure for setting third computational coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/31/2018
    !!
    !! @param[in] zeta Third coordinate value
    !!
    !------------------------------------------------------------------------------------------
    subroutine zeta(self,zeta_in)
        class(point_ad_t), intent(inout) :: self
        type(AD_D) :: zeta_in

        self%c3_ = zeta_in

    end subroutine
    !******************************************************************************************








    !>  Add to x-coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/9/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_x(self, inc_x)
        class(point_ad_t), intent(inout)   :: self
        type(AD_D),       intent(in)      :: inc_x


        self%c1_ = self%c1_ + inc_x


    end subroutine add_x
    !******************************************************************************************








    !>  Add to y-coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/9/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_y(self, inc_y)
        class(point_ad_t), intent(inout)   :: self
        type(AD_D),       intent(in)      :: inc_y


        self%c2_ = self%c2_ + inc_y


    end subroutine add_y
    !******************************************************************************************








    !>  Add to x-coordinate
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/9/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine add_z(self, inc_z)
        class(point_ad_t), intent(inout)   :: self
        type(AD_D),       intent(in)      :: inc_z


        self%c3_ = self%c3_ + inc_z


    end subroutine add_z
    !*****************************************************************************************








    !>  Return status of the point.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/23/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function valid(self)  result(validity)
        class(point_ad_t), intent(in)  :: self

        logical :: validity

        validity = (self%status == VALID_POINT)

    end function valid
    !*****************************************************************************************





end module type_point_ad
