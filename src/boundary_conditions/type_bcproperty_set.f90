module type_bcproperty_set
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_bcproperty,    only: bcproperty_t
    use type_point,         only: point_t
    use type_function,      only: function_t
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
        procedure   :: add              !< Procedure for adding bcfunction_t's to the list


        procedure   :: get_nproperties          !< Return the number of properties in the set.
        procedure   :: get_property_index       !< Return the index of a property in the set, given an identifying string.
        procedure   :: get_property_name        !< Return a string of the property, given an index in the set.

        ! bcfunction%fcn procedures
        procedure   :: set_fcn          !< Procedure for setting the particular function associated with bcfunction_t
        procedure   :: set_fcn_option   !< Procedure for setting options of an allocated function; bcfunction_t%fcn

        ! compute bcfunction%fcn
        procedure   :: compute          !< Compute values from a bcfunction_t definition.

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

        !
        ! Increment nfunctions
        !
        self%nproperties_ = self%nproperties_ + 1
        ifcn = self%nproperties_


        !
        ! Resize array storage
        !
        allocate(temp_bcprop(self%nproperties_), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy previously initialized instances to new array.
        !
        if (self%nproperties_ > 1) then
            temp_bcprop(1:size(self%bcprop)) = self%bcprop(1:size(self%bcprop))
        end if


        !
        ! Initialize new bcproperty
        !
        temp_bcprop(ifcn)%name_  = bcprop
        temp_bcprop(ifcn)%type_  = ftype


        !
        ! Set default function to f(t,x,y,z) = constant
        !
        call temp_bcprop(ifcn)%set('function','constant')



        !
        ! Move resized temp allocation back to self%bcprop(:)
        !
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
        

        !
        ! Get index of bcfcn in self%bcfcn(:)
        !
        ind = self%get_property_index(bcprop)


        !
        ! Set concrete function
        !
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

        !
        ! Get index of bcfcn in self%bcfcn(:)
        !
        ifcn = self%get_property_index(bcprop)


        !
        ! Set function option
        !
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

        
        !
        ! Set known error value
        !
        ind = 0

        
        !
        ! Search through functions
        !
        do ifcn = 1,size(self%bcprop)

        
            !
            ! Check for matching name
            !
            name_matches = ( trim(bcprop) == trim(self%bcprop(ifcn)%get_name()) )


            if ( name_matches ) then
                ind = ifcn
                exit
            end if

        end do !ifcn

        
        !
        ! Check that 'ind' has been set.
        !
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
    impure elemental function compute(self,bcprop,time,coord) result(val)
        class(bcproperty_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: bcprop
        real(rk),                   intent(in)      :: time
        type(point_t),              intent(in)      :: coord

        integer(ik) :: ifcn
        real(rk)    :: val
        logical     :: fcn_set

        !
        ! Get index of bcfunction_t specified by bcfcn in self%bcfcn(:)
        !
        ifcn = self%get_property_index(bcprop)


        !
        ! Check function status
        !
        fcn_set = self%bcprop(ifcn)%status()
        


        !
        ! Compute function
        !
        if (fcn_set) then
            val = self%bcprop(ifcn)%fcn%compute(time,coord)
        else
            call chidg_signal(FATAL,"bcfunction_set%compute: bcfcn%fcn not allocated")
        end if

    end function compute
    !***************************************************************************************

















end module type_bcproperty_set
