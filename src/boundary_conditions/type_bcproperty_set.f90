module type_bcproperty_set
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ONE, dBC_DIFF
    use type_bcproperty,    only: bcproperty_t
    use type_point,         only: point_t
    use type_point_ad,      only: point_ad_t
    use type_function,      only: function_t
    use type_function_info, only: function_info_t
    use DNAD_D
    implicit none







    !>  A class containing a set of bcproperty_t instances. This class manages the 
    !!  addition of properties, setting their functions, setting their options.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !--------------------------------------------------------------------------------
    type, public :: bcproperty_set_t

        integer(ik)                         :: nproperties_ = 0 !< Number of functions
        type(bcproperty_t),     allocatable :: bcprop(:)        !< boundary property list

    contains

        ! bcfcn procedures
        procedure   :: add                  ! Procedure for adding bcfunction_t's to the list


        procedure   :: get_nproperties      ! Return the number of properties in the set.
        procedure   :: get_property_index   ! Return the index of a property in the set, given an identifying string.
        procedure   :: get_property_name    ! Return a string of the property, given an index in the set.

        ! bcfunction%fcn procedures
        procedure   :: set_fcn              ! Procedure for setting the particular function associated with bcfunction_t
        procedure   :: set_fcn_option       ! Procedure for setting options of an allocated function; bcfunction_t%fcn

        ! compute bcfunction%fcn
        procedure   :: compute_real         ! Compute values from a bcfunction_t definition.
        procedure   :: compute_ad           ! Compute values from a bcfunction_t definition.

        ! compute generic procedure
        generic     :: compute => compute_real,compute_ad

    end type bcproperty_set_t
    !********************************************************************************



contains


    !>  Add a bcfunction_t to the array self%bcfcn(:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  bcfcn   String specifying the name of the bcfunction_t. Ex: 'Static Pressure'
    !!  @param[in]  ftype   Character string specifying the requirement type. 'Required' or 'Optional'
    !!
    !--------------------------------------------------------------------------------
    subroutine add(self,bcprop,ftype)
        class(bcproperty_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: bcprop
        character(*),               intent(in)      :: ftype

        type(bcproperty_t), allocatable :: temp_bcprop(:)
        integer(ik)                     :: ierr, ifcn

        ! Increment nfunctions
        self%nproperties_ = self%nproperties_ + 1
        ifcn = self%nproperties_

        ! Resize array storage
        allocate(temp_bcprop(self%nproperties_), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Copy previously initialized instances to new array.
        if (self%nproperties_ > 1) temp_bcprop(1:size(self%bcprop)) = self%bcprop(1:size(self%bcprop))

        ! Initialize new bcproperty
        temp_bcprop(ifcn)%name_  = bcprop
        temp_bcprop(ifcn)%type_  = ftype

        ! Set default function to f(t,x,y,z) = constant
        call temp_bcprop(ifcn)%set('function','constant')

        ! Move resized temp allocation back to self%bcprop(:)
        call move_alloc(temp_bcprop,self%bcprop)

    end subroutine add
    !********************************************************************************







    !>  Procedure for setting a concrete function in a particular bcfunction_t object.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  bcfcn   String containing the name of the bcfunction_t to be modified.
    !!  @param[in]  fcn     String indicating the type of function to be set in the bcfunction_t
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_fcn(self,bcprop,fcn)
        class(bcproperty_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: bcprop
        character(*),               intent(in)      :: fcn

        integer(ik)     :: ind
        
        ! Get index of bcfcn in self%bcfcn(:)
        ind = self%get_property_index(bcprop)

        ! Set concrete function
        call self%bcprop(ind)%set('function',fcn)

    end subroutine set_fcn
    !************************************************************************************






    !>  Procedure for setting options of a concrete function in a particular bcfunction_t object.
    !!  bcfunction%fcn
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcfcn   String containing the name of the bcfunction_t to be modified.
    !!  @param[in]  option  String indicating the option to be set.
    !!  @param[in]  val     Real value to be set for the specified option.
    !!
    !-------------------------------------------------------------------------------------
    subroutine set_fcn_option(self,bcprop,key,val)
        class(bcproperty_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: bcprop
        character(*),               intent(in)      :: key
        real(rk),                   intent(in)      :: val

        integer(ik) :: ifcn

        ! Get index of bcfcn in self%bcfcn(:)
        ifcn = self%get_property_index(bcprop)

        ! Set function option
        call self%bcprop(ifcn)%fcn%set_option(key,val)

    end subroutine set_fcn_option
    !*************************************************************************************






    !>  Procedure for returning the index of a bcfunction_t in self%bcfcn(:) specified
    !!  by an identifying character string. Ex: 'Static Pressure'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  bcfcn   String specifying the name of a particular bcfunction_t
    !!  @result     ind     Index of the given bcfunction_t in self%bcfcn(:)
    !!
    !------------------------------------------------------------------------------------
    function get_property_index(self,bcprop) result(ind)
        class(bcproperty_set_t),    intent(in)  :: self
        character(*),               intent(in)  :: bcprop

        integer(ik) :: ind, ifcn
        logical     :: name_matches

        
        ! Set known error value
        ind = 0

        
        ! Search through functions
        do ifcn = 1,size(self%bcprop)

            ! Check for matching name
            name_matches = ( trim(bcprop) == trim(self%bcprop(ifcn)%get_name()) )

            if ( name_matches ) then
                ind = ifcn
                exit
            end if

        end do !ifcn

        ! Check that 'ind' has been set.
        if ( ind == 0 ) call chidg_signal_one(FATAL,"bcfunction_set%get_fcn_index: function key not found",bcprop)

    end function get_property_index
    !*************************************************************************************






    !>  Given the index of a property in the list, return the property name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!  @param[in]  iprop   Index of a property in the set to be queried.
    !!  @result     pname   String containing the name of the specified property.
    !!
    !-------------------------------------------------------------------------------------
    function get_property_name(self,iprop) result(pname)
        class(bcproperty_set_t),    intent(in)  :: self
        integer(ik),                intent(in)  :: iprop

        character(len=:),   allocatable :: pname

        pname = self%bcprop(iprop)%get_name()

    end function get_property_name
    !*************************************************************************************






    !>  Return the number of properties in the set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_nproperties(self) result(nprop)
        class(bcproperty_set_t),    intent(in)  :: self

        integer(ik) :: nprop

        nprop = self%nproperties_

    end function get_nproperties
    !****************************************************************************************





!    !>  Procedure to compute values from a bcfunction_t definition.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/3/2016
!    !!
!    !!  @param[in]  bcprop  String of the property to be computed.
!    !!  @param[in]  time    Real value of the global time.
!    !!  @param[in]  coord   point_t instance containing the coordinates.
!    !!  @result     val     Real value of the function. f = f(t,coord)
!    !!
!    !---------------------------------------------------------------------------------------
!    impure elemental function compute(self,bcprop,time,coord) result(val)
!        class(bcproperty_set_t),    intent(inout)   :: self
!        character(*),               intent(in)      :: bcprop
!        real(rk),                   intent(in)      :: time
!        type(point_t),              intent(in)      :: coord
!
!        integer(ik) :: ifcn
!        real(rk)    :: val
!        logical     :: fcn_set
!
!        ! Get index of bcfunction_t specified by bcfcn in self%bcfcn(:)
!        ifcn = self%get_property_index(bcprop)
!
!        ! Check function status
!        fcn_set = self%bcprop(ifcn)%status()
!        
!        ! Compute function
!        if (fcn_set) then
!            val = self%bcprop(ifcn)%fcn%compute(time,coord)
!        else
!            call chidg_signal(FATAL,"bcfunction_set%compute: bcfcn%fcn not allocated")
!        end if
!
!    end function compute
!    !***************************************************************************************






    !>  Procedure to compute values from a bcfunction_t definition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  bcprop  String of the property to be computed.
    !!  @param[in]  time    Real value of the global time.
    !!  @param[in]  coord   point_t instance containing the coordinates.
    !!  @result     val     Real value of the function. f = f(t,coord)
    !!
    !---------------------------------------------------------------------------------------
    impure elemental function compute_real(self,bcprop,time,coord_rk) result(val_rk)
        class(bcproperty_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: bcprop
        real(rk),                   intent(in)      :: time
        type(point_t),              intent(in)      :: coord_rk

        real(rk)    :: val_rk
        
        integer(ik)         :: ifcn
        logical             :: fcn_set
        type(AD_D)          :: val, x_ad, y_ad, z_ad
        type(point_ad_t)    :: coord_ad

        ! Get index of bcfunction_t specified by bcfcn in self%bcfcn(:)
        ifcn = self%get_property_index(bcprop)

        ! Check function status
        fcn_set = self%bcprop(ifcn)%status()
        
        ! Initialize dummy AD_D type for coordinates
        ! This is because fcn%compute accepts coordinates in AD_D form
        x_ad = AD_D(0)
        y_ad = AD_D(0)
        z_ad = AD_D(0)

        x_ad = coord_rk%c1_
        y_ad = coord_rk%c2_
        z_ad = coord_rk%c3_

        call coord_ad%set(x_ad,y_ad,z_ad)

        ! Compute function
        if (fcn_set) then
            val    = self%bcprop(ifcn)%fcn%compute(time,coord_ad)
            val_rk = val%x_ad_
        else
            call chidg_signal(FATAL,"bcfunction_set%compute: bcfcn%fcn not allocated")
        end if

    end function compute_real
    !***************************************************************************************








    !>  Procedure to compute values from a bcfunction_t definition.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/31/2018
    !!
    !!  @param[in]  bcprop  String of the property to be computed.
    !!  @param[in]  time    Real value of the global time.
    !!  @param[in]  coord   point_add_t instance containing the coordinates.
    !!  @result     val     Real value of the function. f = f(t,coord)
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !!  Added initialization of derivatives for BC parameter sensitivities.
    !!
    !---------------------------------------------------------------------------------------
    impure elemental function compute_ad(self,bcprop,time,coord,fcn_info) result(val)
        class(bcproperty_set_t),    intent(inout)           :: self
        character(*),               intent(in)              :: bcprop
        real(rk),                   intent(in)              :: time
        type(point_ad_t),           intent(in)              :: coord
        type(function_info_t),      intent(in),    optional :: fcn_info

        integer(ik) :: ifcn
        type(AD_D)  :: val
        logical     :: fcn_set, diff_me

        ! Get index of bcfunction_t specified by bcfcn in self%bcfcn(:)
        ifcn = self%get_property_index(bcprop)

        ! Check function status
        fcn_set = self%bcprop(ifcn)%status()
        
        ! Compute function
        if (fcn_set) then
            val = self%bcprop(ifcn)%fcn%compute(time,coord)
        else
            call chidg_signal(FATAL,"bcfunction_set%compute: bcfcn%fcn not allocated")
        end if

        ! Initialize derivatives if BC linearization
        !diff_me = (present(fcn_info)          .and. &
        !           fcn_info%dtype == dBC_DIFF .and. &
        !           fcn_info%bc_group_match    .and. &
        !           fcn_info%bc_param == bcprop     )
        diff_me = .false.
        if (present(fcn_info)) then
            diff_me = (fcn_info%dtype == dBC_DIFF .and. &
                       fcn_info%bc_group_match    .and. &
                       fcn_info%bc_param == bcprop     )
        end if

        if (diff_me) then
            val%xp_ad_ = ONE
        end if

    end function compute_ad
    !***************************************************************************************





end module type_bcproperty_set
