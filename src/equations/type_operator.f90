module type_operator
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: BOUNDARY_ADVECTIVE_FLUX, VOLUME_ADVECTIVE_FLUX, &
                                  BOUNDARY_DIFFUSIVE_FLUX, VOLUME_DIFFUSIVE_FLUX, BC_FLUX
    use mod_string,         only: string_t, string_to_upper
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @note   Refactoring
    !!
    !!  Valid Operator Types are:
    !!      BOUNDARY_ADVECTIVE_FLUX
    !!      BOUNDARY_DIFFUSIVE_FLUX
    !!      VOLUME_ADVECTIVE_FLUX
    !!      VOLUME_DIFFUSIVE_FLUX
    !!      BC_FLUX
    !!
    !----------------------------------------------------------------------------------
    type, abstract, public :: operator_t

        integer(ik)                     :: operator_type
        character(:),       allocatable :: name

        type(string_t),     allocatable :: primary_fields(:)
        type(string_t),     allocatable :: auxiliary_fields(:)

    contains

        procedure(self_interface),      deferred :: init
        procedure(compute_interface),   deferred :: compute

        procedure   :: set_name
        procedure   :: get_name

        procedure   :: set_operator_type
        procedure   :: get_operator_type

        procedure   :: add_primary_field
        procedure   :: add_auxiliary_field

        procedure   :: nprimary_fields
        procedure   :: nauxiliary_fields

    end type operator_t
    !**********************************************************************************




    abstract interface
        subroutine compute_interface(self,worker,prop)
            import operator_t
            import chidg_worker_t
            import properties_t

            class(operator_t),      intent(inout)   :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
        end subroutine
    end interface




    abstract interface
        subroutine self_interface(self)
            import operator_t

            class(operator_t),      intent(inout)   :: self
        end subroutine
    end interface






contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine set_name(self, string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        self%name = trim(string)

    end subroutine set_name
    !****************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_name(self) result(string)
        class(operator_t),  intent(inout)   :: self

        character(len=:), allocatable :: string

        string = trim(self%name)

    end function get_name
    !****************************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine set_operator_type(self, string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        integer(ik)                     :: ind
        character(len=:),   allocatable :: upper
        logical                         :: volume, boundary, advective, diffusive, bc
        


        ! Convert to upper case
        upper = string_to_upper(string)


        ! Get operator string key words
        ind = index(upper,"VOLUME")
        volume    = (ind /= 0)

        ind = index(upper,"BOUNDARY")
        boundary  = (ind /= 0)

        ind = index(upper,"BC")
        bc        = (ind /= 0)

        ind = index(upper,"ADVECTIVE")
        advective = (ind /= 0)

        ind = index(upper,"DIFFUSIVE")
        diffusive = (ind /= 0)

        
        ! Set operator type
        if (volume .and. advective) then
            self%operator_type = VOLUME_ADVECTIVE_FLUX
        else if (volume .and. diffusive) then
            self%operator_type = VOLUME_DIFFUSIVE_FLUX
        else if (boundary .and. advective) then
            self%operator_type = BOUNDARY_ADVECTIVE_FLUX
        else if (boundary .and. diffusive) then
            self%operator_type = BOUNDARY_DIFFUSIVE_FLUX
        else if ( (bc .and. advective) .or. (bc .and. diffusive) ) then
            self%operator_type = BC_FLUX
        else
            call chidg_signal(FATAL,"operator%set_operator_type: Did not recognize a valid operator type")
        end if


    end subroutine set_operator_type
    !**************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    function get_operator_type(self) result(otype)
        class(operator_t),  intent(in)  :: self

        integer(ik) :: otype

        otype = self%operator_type

    end function get_operator_type
    !***************************************************************************************************










    !>  Add a primary field that the operator is integrating.
    !!
    !!  This is a field that is being solved for. Something like Density, or Momentum.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine add_primary_field(self,string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        integer(ik)     :: ierr, ieq
        type(string_t), allocatable  :: temp(:)


        !
        ! Extend fields if necessary
        !
        if (allocated(self%primary_fields)) then
            
            allocate(temp(size(self%primary_fields) + 1), stat=ierr)
            do ieq = 1,size(self%primary_fields)
                temp(ieq) = self%primary_fields(ieq)
            end do
            
        else
            allocate(temp(1), stat=ierr)
        end if
        if (ierr /= 0) call AllocationError


        !
        ! Set new variable
        !
        call temp(size(temp))%set(string)


        !
        ! Copy temp back to self
        !
        self%primary_fields = temp


    end subroutine add_primary_field
    !***************************************************************************************************









    !>  Add an auxiliary field that the operator is using.
    !!
    !!  This might be something like a wall distance field or a blockage field. These
    !!  aren't being solved for, but they are being used, maybe in a source term or something like that.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine add_auxiliary_field(self,string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        integer(ik)     :: ierr, ieq
        type(string_t), allocatable  :: temp(:)


        !
        ! Extend fields if necessary
        !
        if (allocated(self%auxiliary_fields)) then
            
            allocate(temp(size(self%auxiliary_fields) + 1), stat=ierr)
            do ieq = 1,size(self%auxiliary_fields)
                temp(ieq) = self%auxiliary_fields(ieq)
            end do
            
        else
            allocate(temp(1), stat=ierr)
        end if
        if (ierr /= 0) call AllocationError


        !
        ! Set new variable
        !
        call temp(size(temp))%set(string)


        !
        ! Copy temp back to self
        !
        self%auxiliary_fields = temp


    end subroutine add_auxiliary_field
    !********************************************************************************************








    !>  Return the number of primary fields.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !--------------------------------------------------------------------------------------------
    function nprimary_fields(self) result(nfields)
        class(operator_t),  intent(in)  :: self

        integer(ik) :: nfields


        if (allocated(self%primary_fields)) then
            nfields = size(self%primary_fields)
        else
            nfields = 0
        end if

    end function nprimary_fields
    !********************************************************************************************







    !>  Return the number of auxiliary fields.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !--------------------------------------------------------------------------------------------
    function nauxiliary_fields(self) result(nfields)
        class(operator_t),  intent(in)  :: self

        integer(ik) :: nfields


        if (allocated(self%auxiliary_fields)) then
            nfields = size(self%auxiliary_fields)
        else
            nfields = 0
        end if

    end function nauxiliary_fields
    !********************************************************************************************















end module type_operator
