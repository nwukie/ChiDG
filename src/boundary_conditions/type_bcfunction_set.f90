module type_bcfunction_set
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_bcfunction,    only: bcfunction_t
    use type_point,         only: point_t
    implicit none







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !--------------------------------------------------------------------------------
    type, public :: bcfunction_set_t

        integer(ik)                         :: nfunctions_ = 0  !< Number of functions
        type(bcfunction_t),     allocatable :: bcfcn(:)         !< Bounday function list

    contains

        ! bcfcn procedures
        procedure   :: add              !< Procedure for adding bcfunction_t's to the list
        procedure   :: get_bcfcn_index  !< Procedure for returning the index of a particular bcfunction_t in bcfun(:)

        ! bcfunction%fcn procedures
        procedure   :: set_fcn          !< Procedure for setting the particular function associated with bcfunction_t
        procedure   :: set_fcn_option   !< Procedure for setting options of an allocated function; bcfunction_t%fcn

        ! compute bcfunction%fcn
        procedure   :: compute          !< Compute values from a bcfunction_t definition.

    end type bcfunction_set_t
    !********************************************************************************



contains


    !>  Add a bcfunction_t to the array self%bcfcn(:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  fname   String specifying the name of the bcfunction_t. Ex: 'Static Pressure'
    !!  @param[in]  ftype   Character string specifying the requirement type. 'Required' or 'Optional'
    !!
    !--------------------------------------------------------------------------------
    subroutine add(self,fname,ftype)
        class(bcfunction_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: fname
        character(*),               intent(in)      :: ftype

        type(bcfunction_t), allocatable :: temp_bcfcn(:)
        integer(ik)                     :: ierr, ifcn

        !
        ! Increment nfunctions
        !
        self%nfunctions_ = self%nfunctions_ + 1
        ifcn = self%nfunctions_


        !
        ! Resize array storage
        !
        allocate(temp_bcfcn(self%nfunctions_), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy previously initialized instances to new array.
        !
        if (self%nfunctions_ > 1) then
            temp_bcfcn(1:size(self%bcfcn)) = self%bcfcn(1:size(self%bcfcn))
        end if


        !
        ! Initialize new bcfunction
        !
        temp_bcfcn(ifcn)%name_  = fname
        temp_bcfcn(ifcn)%type_  = ftype


        !
        ! Move resized temp allocation back to self%bcfcn(:)
        !
        call move_alloc(temp_bcfcn,self%bcfcn)


    end subroutine add
    !********************************************************************************







    !>  Procedure for setting a concrete function in a particular bcfunction_t object.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  fname   String containing the name of the bcfunction_t to be modified.
    !!  @param[in]  fstr    String indicating the type of function to be set in the bcfunction_t
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_fcn(self,fname,fstr)
        class(bcfunction_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: fname
        character(*),               intent(in)      :: fstr

        integer(ik)     :: ind
        

        !
        ! Get index of fname in self%bcfcn(:)
        !
        ind = self%get_bcfcn_index(fname)


        !
        ! Set concrete function
        !
        call self%bcfcn(ind)%set('function',fstr)


    end subroutine set_fcn
    !************************************************************************************






    !>  Procedure for setting options of a concrete function in a particular bcfunction_t object.
    !!  bcfunction%fcn
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine set_fcn_option(self,fname,foption,val)
        class(bcfunction_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: fname
        character(*),               intent(in)      :: foption
        real(rk),                   intent(in)      :: val

        integer(ik) :: ifcn

        !
        ! Get index of fname in self%bcfcn(:)
        !
        ifcn = self%get_bcfcn_index(fname)


        !
        ! Set function option
        !
        call self%bcfcn(ifcn)%fcn%set(foption,val)



    end subroutine set_fcn_option
    !*************************************************************************************













    !>  Procedure for returning the index of a bcfunction_t in self%bcfcn(:) specified
    !!  by an identifying character string. Ex: 'Static Pressure'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  fname   String specifying the name of a particular bcfunction_t
    !!  @result     ind     Index of the given bcfunction_t in self%bcfcn(:)
    !!
    !------------------------------------------------------------------------------------
    function get_bcfcn_index(self,fname) result(ind)
        class(bcfunction_set_t),    intent(in)  :: self
        character(*),               intent(in)  :: fname

        integer(ik) :: ind, ifcn
        logical     :: name_matches

        
        !
        ! Set known error value
        !
        ind = 0

        
        !
        ! Search through functions
        !
        do ifcn = 1,size(self%bcfcn)

        
            !
            ! Check for matching name
            !
            name_matches = ( trim(fname) == trim(self%bcfcn(ifcn)%name_) )


            if ( name_matches ) then
                ind = ifcn
                exit
            end if

        end do !ifcn

        
        !
        ! Check that 'ind' has been set.
        !
        if ( ind == 0 ) call chidg_signal_one(FATAL,"bcfunction_set%get_fcn_index: function key not found",fname)


    end function get_bcfcn_index
    !*************************************************************************************










    !>  Procedure to compute values from a bcfunction_t definition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    impure elemental function compute(self,fname,time,coord) result(val)
        class(bcfunction_set_t),    intent(inout)   :: self
        character(*),               intent(in)      :: fname
        real(rk),                   intent(in)      :: time
        type(point_t),              intent(in)      :: coord

        integer(ik) :: ifcn
        real(rk)    :: val
        logical     :: fcn_set

        !
        ! Get index of bcfunction_t specified by fname in self%bcfcn(:)
        !
        ifcn = self%get_bcfcn_index(fname) 


        !
        ! Check function status
        !
        fcn_set = self%bcfcn(ifcn)%status()


        !
        ! Compute function
        !
        if (fcn_set) then
            val = self%bcfcn(ifcn)%fcn%compute(time,coord)
        else
            call chidg_signal(FATAL,"bcfunction_set%compute: bcfcn%fcn not allocated")
        end if

    end function compute
    !***************************************************************************************









end module type_bcfunction_set
