module type_properties
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_field,     only: field_t
    use type_fluid,     only: fluid_t
    use type_scalar,    only: scalar_t
    implicit none




    !>  Base properties type for storing equations, material definitions, 
    !!  and miscellaneous data pertaining to a particular equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!  @date   11/18/2016  auxiliary fields
    !!
    !---------------------------------------------------------------------------------------------
    type, public :: properties_t
        
        ! Fields
        type(field_t),   allocatable :: primary_fields(:)
        type(field_t),   allocatable :: auxiliary_fields(:)
        type(field_t),   allocatable :: model_fields(:)

        ! Materials
        class(fluid_t),     allocatable :: fluid
        class(scalar_t),    allocatable :: scalar

    contains

        procedure   :: add_fluid
        procedure   :: add_scalar

        procedure   :: add_primary_field
        procedure   :: add_auxiliary_field
        procedure   :: add_model_field

        procedure   :: get_primary_field_name
        procedure   :: get_primary_field_index
        procedure   :: get_auxiliary_field_name
        procedure   :: get_auxiliary_field_index
        procedure   :: get_model_field_name
        procedure   :: get_model_field_index

        procedure   :: nprimary_fields
        procedure   :: nauxiliary_fields
        procedure   :: nmodel_fields

    end type properties_t
    !*********************************************************************************************






contains






    !>  Add an equation to the list.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_primary_field(self,field_string)
        class(properties_t),    intent(inout)   :: self
        character(*),           intent(in)      :: field_string

        type(field_t),   allocatable :: temp_fields(:)
        integer(ik) :: ieq, ierr, ind
        logical     :: already_added


        !
        ! Check if equation was already added by another function
        !
        ind = self%get_primary_field_index(field_string)
        already_added = (ind /= 0)


        !
        ! Add equation if necessary
        !
        if (.not. already_added) then


            !
            ! If there are already equations allocated, reallocate and add new equation
            !
            if (allocated(self%primary_fields)) then

                ! Allocate temp field array with one extra slot for new field.
                allocate(temp_fields(self%nprimary_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ieq = 1,self%nprimary_fields()
                    temp_fields(ieq) = self%primary_fields(ieq)
                end do

            !
            ! If there are no equations allocated, allocate one slot and set data
            !
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if


            ! Add new field to last slot
            call temp_fields(size(temp_fields))%set_name(field_string)

            ! Move temporary allocation to data type
            call move_alloc(from=temp_fields, to=self%primary_fields)

        end if


    end subroutine add_primary_field
    !**********************************************************************************************









    !>  Add an auxiliary field to the list.
    !!
    !!  This would be something like a 'Wall Distance' field or 'Blockage' field. These
    !!  are not being solved for, but are used in some way by the operator_t's. Maybe
    !!  in a source term, for example.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_auxiliary_field(self,field_string)
        class(properties_t),    intent(inout)   :: self
        character(*),           intent(in)      :: field_string

        type(field_t),   allocatable :: temp_fields(:)
        integer(ik) :: ifield, ierr, ind
        logical     :: already_added


        !
        ! Check if equation was already added by another function
        !
        ind = self%get_auxiliary_field_index(field_string)
        already_added = (ind /= 0)


        !
        ! Add equation if necessary
        !
        if (.not. already_added) then


            !
            ! If there are already equations allocated, reallocate and add new equation
            !
            if (allocated(self%auxiliary_fields)) then

                ! Allocate temp field array with one extra slot for new field
                allocate(temp_fields(self%nauxiliary_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ifield = 1,self%nauxiliary_fields()
                    temp_fields(ifield) = self%auxiliary_fields(ifield)
                end do


            !
            ! If there are no equations allocated, allocate one slot and set data
            !
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if


            !
            ! Set new field at end
            !
            call temp_fields(size(temp_fields))%set_name(field_string)

            !
            ! Move temporary allocation
            !
            call move_alloc(from=temp_fields, to=self%auxiliary_fields)

        end if


    end subroutine add_auxiliary_field
    !*******************************************************************************************








    !>  Add a model field to the list.
    !!
    !!  This would be something like 'Pressure' or 'Viscosity' getting computed from a model. 
    !!  These are not being solved for, but are used in some way by the operator_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_model_field(self,field_string)
        class(properties_t),    intent(inout)   :: self
        character(*),           intent(in)      :: field_string

        type(field_t),   allocatable    :: temp_fields(:)
        integer(ik)                     :: ifield, ierr, ind
        logical                         :: already_added


        !
        ! Check if equation was already added by another function
        !
        ind = self%get_model_field_index(field_string)
        already_added = (ind /= 0)


        !
        ! Add equation if necessary
        !
        if (.not. already_added) then


            !
            ! If there are already equations allocated, reallocate and add new equation
            !
            if (allocated(self%model_fields)) then

                ! Allocate temp field array with one extra slot for new field
                allocate(temp_fields(self%nmodel_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ifield = 1,self%nmodel_fields()
                    temp_fields(ifield) = self%model_fields(ifield)
                end do


            !
            ! If there are no equations allocated, allocate one slot and set data
            !
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if


            !
            ! Set new field at end
            !
            call temp_fields(size(temp_fields))%set_name(field_string)

            !
            ! Move temporary allocation
            !
            call move_alloc(from=temp_fields, to=self%model_fields)

        end if


    end subroutine add_model_field
    !*******************************************************************************************









    !>  Given an equation index, return the equation name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function get_primary_field_name(self,field_index) result(field_name)
        class(properties_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: user_msg, field_name


        ! Check bounds
        user_msg = "properties%get_primary_field_name: Incoming index to return a field name is &
                    out of bounds."
        if (field_index > self%nprimary_fields()) call chidg_signal_one(FATAL,user_msg,field_index)


        ! Get name
        field_name = self%primary_fields(field_index)%get_name()


    end function get_primary_field_name
    !*******************************************************************************************








    !>  Given an auxiliary field index, return the auxiliary field name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function get_auxiliary_field_name(self,field_index) result(field_name)
        class(properties_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: user_msg, field_name

        ! Check bounds
        user_msg = "properties%get_auxiliarya_field_name: Incoming index to return an auxiliary &
                    field is out of bounds."
        if (field_index > self%nauxiliary_fields()) call chidg_signal_one(FATAL,user_msg,field_index)


        ! Get name
        field_name = self%auxiliary_fields(field_index)%get_name()


    end function get_auxiliary_field_name
    !*******************************************************************************************







    !>  Given a model field index, return the model field name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function get_model_field_name(self,field_index) result(field_name)
        class(properties_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: user_msg, field_name

        ! Check bounds
        user_msg = "properties%get_model_field_name: Incoming index to return an auxiliary &
                    field is out of bounds."
        if (field_index > self%nmodel_fields()) call chidg_signal_one(FATAL,user_msg,field_index)


        ! Get name
        field_name = self%model_fields(field_index)%get_name()


    end function get_model_field_name
    !*******************************************************************************************









    !>  Search for a equation string in the self%primary_fields list. If found, return equation index.
    !!  A set of equations could be stored in any order. So, when an equation is initialized, it
    !!  is initialized with an index indicating its location in the set. That index is used to 
    !!  access the correct solution data values.
    !!
    !!  Return 0 if not found.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !!  @param[in]  field_string   Character string identifying the desired variable
    !!
    !--------------------------------------------------------------------------------------------
    function get_primary_field_index(self,field_string) result(field_index)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: field_string

        integer(ik) :: field_index, ifield
        logical     :: found = .false.


        field_index = 0


        ! Search for character string in self%primary_fields array. If found set index
        do ifield = 1,self%nprimary_fields()
            if (field_string == self%primary_fields(ifield)%name) then
                field_index = ifield
                found = .true.
                exit
            end if
        end do


    end function get_primary_field_index
    !*******************************************************************************************









    !> Search for a equation string in the self%primary_fields list. If found, return equation index.
    !! A set of equations could be stored in any order. So, when an equation is initialized, it
    !! is initialized with an index indicating its location in the set. That index is used to 
    !! access the correct solution data values.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/20/2016
    !!
    !!  @param[in]  field_string   Character string identifying the desired variable
    !!
    !------------------------------------------------------------------------------------------
    function get_auxiliary_field_index(self,field_string) result(field_index)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: field_string

        integer(ik) :: field_index, ifield
        logical     :: found = .false.


        field_index = 0


        ! Search for character string in self%auxiliary_fields array. If found set index
        do ifield = 1,self%nauxiliary_fields()
            if (field_string == self%auxiliary_fields(ifield)%name) then
                field_index = ifield
                found = .true.
                exit
            end if
        end do


    end function get_auxiliary_field_index
    !******************************************************************************************









    !>  Search for a field string in the self%model_fields list. If found, return index of the field.
    !!  A set of fields could be stored in any order. So, when a field is initialized, it
    !!  is initialized with an index indicating its location in the set. That index is used to 
    !!  access the correct solution data values.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!  @param[in]  field_string   Character string identifying the desired variable
    !!
    !------------------------------------------------------------------------------------------
    function get_model_field_index(self,field_string) result(field_index)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: field_string

        integer(ik) :: field_index, ifield
        logical     :: found = .false.


        field_index = 0


        ! Search for character string in self%auxiliary_fields array. If found set index
        do ifield = 1,self%nmodel_fields()
            if (field_string == self%model_fields(ifield)%name) then
                field_index = ifield
                found = .true.
                exit
            end if
        end do


    end function get_model_field_index
    !******************************************************************************************













    !>  Return number of primary equations that have been added.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/31/2016
    !!
    !------------------------------------------------------------------------------------------
    function nprimary_fields(self) result(nfields)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%primary_fields)) then
            nfields = size(self%primary_fields)
        else
            nfields = 0
        end if

    end function nprimary_fields
    !******************************************************************************************









    !>  Return number of auxiliary fields that have been added.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !------------------------------------------------------------------------------------------
    function nauxiliary_fields(self) result(nfields)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%auxiliary_fields)) then
            nfields = size(self%auxiliary_fields)
        else
            nfields = 0
        end if

    end function nauxiliary_fields
    !******************************************************************************************








    !>  Return number of model fields that have been added.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !------------------------------------------------------------------------------------------
    function nmodel_fields(self) result(nfields)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%model_fields)) then
            nfields = size(self%model_fields)
        else
            nfields = 0
        end if

    end function nmodel_fields
    !******************************************************************************************










    !>  Add a fluid definition to the properties type
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_fluid(self,fluid)
        class(properties_t),    intent(inout)   :: self
        class(fluid_t),         intent(in)      :: fluid

        integer(ik) :: ierr

        ! Allocate new material definition
        if (allocated(self%fluid)) deallocate(self%fluid)
        allocate(self%fluid, source=fluid, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine add_fluid
    !******************************************************************************************








    !>  Add a fluid definition to the properties type
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_scalar(self,scalar)
        class(properties_t),    intent(inout)   :: self
        class(scalar_t),        intent(in)      :: scalar

        integer(ik) :: ierr

        ! Allocate new material definition
        if (allocated(self%scalar)) deallocate(self%scalar)
        allocate(self%scalar, source=scalar, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine add_scalar
    !******************************************************************************************














end module type_properties
