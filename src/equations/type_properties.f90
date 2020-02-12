module type_properties
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_field,     only: field_t
    implicit none




    !>  Base properties type for storing equations, material definitions, 
    !!  and miscellaneous data pertaining to a particular equation set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!  @date   11/18/2016  auxiliary fields
    !!  @date   7/6/2017    io fields
    !!
    !---------------------------------------------------------------------------------------------
    type, public :: properties_t
        
        ! Fields
        type(field_t),  allocatable :: primary_fields(:)
        type(field_t),  allocatable :: auxiliary_fields(:)
        type(field_t),  allocatable :: adjoint_fields(:)
        type(field_t),  allocatable :: model_fields(:)
        type(field_t),  allocatable :: io_fields(:)

    contains

        procedure   :: add_primary_field
        procedure   :: add_auxiliary_field
        procedure   :: add_model_field
        procedure   :: add_adjoint_fields
        procedure   :: add_io_field

        procedure   :: clear_io_fields
        procedure   :: clear_primary_fields
        procedure   :: clear_model_fields
        procedure   :: clear_adjoint_fields


        procedure   :: get_primary_field_name
        procedure   :: get_primary_field_index
        procedure   :: get_auxiliary_field_name
        procedure   :: get_auxiliary_field_index
        procedure   :: get_adjoint_field_name
        procedure   :: get_adjoint_field_index
        procedure   :: get_model_field_name
        procedure   :: get_model_field_index
        procedure   :: get_io_field_name
        procedure   :: get_io_field_index

        procedure   :: nprimary_fields
        procedure   :: nauxiliary_fields
        procedure   :: nadjoint_fields
        procedure   :: nmodel_fields
        procedure   :: nio_fields

        procedure   :: add_adjoint_fields_to_io_fields

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

        ! Check if equation was already added by another function
        ind = self%get_primary_field_index(field_string)
        already_added = (ind /= 0)

        ! Add equation if necessary
        if (.not. already_added) then

            ! If there are already equations allocated, reallocate and add new equation
            if (allocated(self%primary_fields)) then

                ! Allocate temp field array with one extra slot for new field.
                allocate(temp_fields(self%nprimary_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ieq = 1,self%nprimary_fields()
                    temp_fields(ieq) = self%primary_fields(ieq)
                end do

            ! If there are no equations allocated, allocate one slot and set data
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

            ! Add new field to last slot
            call temp_fields(size(temp_fields))%set_name(field_string)
            call temp_fields(size(temp_fields))%set_cache_index(size(temp_fields))

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

        ! Check if equation was already added by another function
        ind = self%get_auxiliary_field_index(field_string)
        already_added = (ind /= 0)

        ! Add equation if necessary
        if (.not. already_added) then

            ! If there are already equations allocated, reallocate and add new equation
            if (allocated(self%auxiliary_fields)) then

                ! Allocate temp field array with one extra slot for new field
                allocate(temp_fields(self%nauxiliary_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ifield = 1,self%nauxiliary_fields()
                    temp_fields(ifield) = self%auxiliary_fields(ifield)
                end do

            ! If there are no equations allocated, allocate one slot and set data
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

            ! Set new field at end
            call temp_fields(size(temp_fields))%set_name(field_string)
            call temp_fields(size(temp_fields))%set_cache_index(size(temp_fields))

            ! Move temporary allocation
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


        ! Check if equation was already added by another function
        ind = self%get_model_field_index(field_string)
        already_added = (ind /= 0)


        ! Add equation if necessary
        if (.not. already_added) then

            ! If there are already equations allocated, reallocate and add new equation
            if (allocated(self%model_fields)) then

                ! Allocate temp field array with one extra slot for new field
                allocate(temp_fields(self%nmodel_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ifield = 1,self%nmodel_fields()
                    temp_fields(ifield) = self%model_fields(ifield)
                end do


            ! If there are no equations allocated, allocate one slot and set data
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

            ! Set new field at end
            call temp_fields(size(temp_fields))%set_name(field_string)

            ! Move temporary allocation
            call move_alloc(from=temp_fields, to=self%model_fields)

        end if


    end subroutine add_model_field
    !*******************************************************************************************





    !>  Add a io field to the list.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/6/2017
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_io_field(self,field_string,ifunc,ivar)
        class(properties_t),        intent(inout)   :: self
        character(*),               intent(in)      :: field_string
        integer(ik),    optional,   intent(in)      :: ifunc
        integer(ik),    optional,   intent(in)      :: ivar


        type(field_t),   allocatable    :: temp_fields(:)
        integer(ik)                     :: ifield, ierr, ind
        logical                         :: already_added


        ! Check if equation was already added by another function
        ind = self%get_io_field_index(field_string)
        already_added = (ind /= 0)

        ! Add equation if necessary
        if (.not. already_added) then

            ! If there are already equations allocated, reallocate and add new equation
            if (allocated(self%io_fields)) then

                ! Allocate temp field array with one extra slot for new field
                allocate(temp_fields(self%nio_fields() + 1), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Copy current fields to first temp slots
                do ifield = 1,self%nio_fields()
                    temp_fields(ifield) = self%io_fields(ifield)
                end do

            ! If there are no equations allocated, allocate one slot and set data
            else
                allocate(temp_fields(1), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

            ! Set new field at end
            call temp_fields(size(temp_fields))%set_name(field_string)

            ! If it is an adjoint field set the correspondent functional ID
            if (present(ifunc)) then
                call temp_fields(size(temp_fields))%set_functional_ID(ifunc)
            end if
            
            ! If it is an adjoint field set the correspondent cache index
            if (present(ivar)) then
                call temp_fields(size(temp_fields))%set_cache_index(ivar)
            end if

            ! Move temporary allocation
            call move_alloc(from=temp_fields, to=self%io_fields)

        end if


    end subroutine add_io_field
    !*******************************************************************************************





    !>  Add adjoint fields to the list. No need to pass in the adjoint_field name since each 
    !!  adjoint variable corresponds to a primary variable.
    !!  This is mainly used to add adjoint fields in 'chidg post'
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/4/2017
    !!
    !!  Add capability for io of multiple objective function adjoint variables
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_adjoint_fields(self,nfuncs)
        class(properties_t),    intent(inout)   :: self
        integer(ik),            intent(in)      :: nfuncs

        type(field_t),   allocatable :: temp_fields(:)
        character(:),    allocatable :: primary_f, adjoint_f,appendix
        integer(ik)                  :: ieq, ierr, ind, iadj,ifunc, nfunc, primary_var_index
        logical                      :: already_added
        character(1)                 :: ifunc_string

        ! Loop through all the functionals and store the adjoint variables for each of them
        do ifunc = 1,nfuncs
            
            ! Define the appendix for the ith functional
            write(ifunc_string, "(I1)") ifunc
            appendix = '_' // ifunc_string

            ! For each primary field, create the correspondent adjoint field
            do ieq = 1,self%nprimary_fields()
                
                ! Construct adjoint fields name
                primary_f = self%get_primary_field_name(ieq)
                adjoint_f = 'adjoint_'//trim(primary_f)//trim(appendix)
                

                ! Retrieve correspondent primary field var_index
                primary_var_index = self%primary_fields(ieq)%get_cache_index()
            

                ! Check if adjoint field was already added by another function
                ind = self%get_adjoint_field_index(adjoint_f)
                already_added = (ind /= 0)


                ! Add field if necessary
                if (.not. already_added) then

                    ! If there are already adjoint fields allocated, reallocate and add new field
                    if (allocated(self%adjoint_fields)) then

                        ! Allocate temp field array with one extra slot for new field.
                        allocate(temp_fields(self%nadjoint_fields() + 1), stat=ierr)
                        if (ierr /= 0) call AllocationError

                        ! Copy current fields to first temp slots
                        do iadj = 1,self%nadjoint_fields()
                            temp_fields(iadj) = self%adjoint_fields(iadj)
                        end do

                    ! If there are no adjoint fields allocated, allocate one slot and set data
                    else
                        allocate(temp_fields(1), stat=ierr)
                        if (ierr /= 0) call AllocationError

                    end if


                    ! Add new field to last slot
                    call temp_fields(size(temp_fields))%set_name(adjoint_f)
                    call temp_fields(size(temp_fields))%set_functional_ID(ifunc)
                    call temp_fields(size(temp_fields))%set_cache_index(primary_var_index)

                    ! Move temporary allocation to data type
                    call move_alloc(from=temp_fields, to=self%adjoint_fields)

                end if
            

            end do !primary fields

        end do !ifunc

    end subroutine add_adjoint_fields
    !**********************************************************************************************




    !>  Add a adjoint fields to io fields
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/9/2018
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_adjoint_fields_to_io_fields(self)
        class(properties_t),        intent(inout)   :: self

        integer(ik)                     :: ifield_a,field_ifunc,field_ivar
        character(:),   allocatable     :: field_name
        
        do ifield_a = 1,self%nadjoint_fields()

            field_name  = self%adjoint_fields(ifield_a)%get_name()
            field_ifunc = self%adjoint_fields(ifield_a)%get_functional_ID()
            field_ivar  = self%adjoint_fields(ifield_a)%get_cache_index()
            call self%add_io_field(trim(field_name),field_ifunc,field_ivar)

        end do


    end subroutine add_adjoint_fields_to_io_fields
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





    !>  Given an adjoint field index, return the adjoint field name.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/4/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    function get_adjoint_field_name(self,field_index) result(field_name)
        class(properties_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: user_msg, field_name

        ! Check bounds
        user_msg = "properties%get_adjoint_field_name: Incoming index to return a field name is &
                    out of bounds."
        if (field_index > self%nadjoint_fields()) call chidg_signal_one(FATAL,user_msg,field_index)

        ! Get name
        field_name = self%adjoint_fields(field_index)%get_name()

    end function get_adjoint_field_name
    !*******************************************************************************************





    !>  Given a io field index, return the io field name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/6/2017
    !!
    !-------------------------------------------------------------------------------------------
    function get_io_field_name(self,field_index) result(field_name)
        class(properties_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: field_index

        character(:),   allocatable :: user_msg, field_name

        ! Check bounds
        user_msg = "properties%get_io_field_name: Incoming index to return an auxiliary &
                    field is out of bounds."
        if (field_index > self%nio_fields()) call chidg_signal_one(FATAL,user_msg,field_index)

        ! Get name
        field_name = self%io_fields(field_index)%get_name()

    end function get_io_field_name
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





    !>  Search for a equation string in the self%adjoint_fields list. If found, return equation index.
    !!  A set of equations could be stored in any order. So, when an equation is initialized, it
    !!  is initialized with an index indicating its location in the set. That index is used to 
    !!  access the correct solution data values.
    !!
    !!  Return 0 if not found.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/25/2016
    !!
    !!  @param[in]  field_string   Character string identifying the desired variable
    !!
    !--------------------------------------------------------------------------------------------
    function get_adjoint_field_index(self,field_string) result(field_index)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: field_string

        integer(ik) :: field_index, ifield, nvar, ivar

        field_index = 0

        do ifield = 1,self%nadjoint_fields()
            if (field_string == self%adjoint_fields(ifield)%name) then
                field_index = self%adjoint_fields(ifield)%get_cache_index()
                exit
            end if
        end do

    end function get_adjoint_field_index
    !*******************************************************************************************





    !>  Search for a field string in the self%io_fields list. If found, return index of the field.
    !!  A set of fields could be stored in any order. So, when a field is initialized, it
    !!  is initialized with an index indicating its location in the set. That index is used to 
    !!  access the correct solution data values.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/6/2017
    !!
    !!  @param[in]  field_string   Character string identifying the desired variable
    !!
    !------------------------------------------------------------------------------------------
    function get_io_field_index(self,field_string) result(field_index)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: field_string

        integer(ik) :: field_index, ifield
        logical     :: found = .false.

        field_index = 0

        ! Search for character string in self%auxiliary_fields array. If found set index
        do ifield = 1,self%nio_fields()
            if (field_string == self%io_fields(ifield)%name) then
                field_index = ifield
                found = .true.
                exit
            end if
        end do

    end function get_io_field_index
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





    !>  Return number of adjoint fields that have been added.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/4/2017
    !!
    !------------------------------------------------------------------------------------------
    function nadjoint_fields(self) result(nfields)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%adjoint_fields)) then
            nfields = size(self%adjoint_fields)
        else
            nfields = 0
        end if

    end function nadjoint_fields
    !******************************************************************************************





    !>  Return number of io fields that have been added.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !------------------------------------------------------------------------------------------
    function nio_fields(self) result(nfields)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%io_fields)) then
            nfields = size(self%io_fields)
        else
            nfields = 0
        end if

    end function nio_fields
    !******************************************************************************************





    !>  Clear any IO fields that have been added.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/16/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear_io_fields(self)
        class(properties_t),    intent(inout)   :: self

        if (allocated(self%io_fields)) then
            deallocate(self%io_fields)
        end if

    end subroutine clear_io_fields
    !******************************************************************************************



    !>  Clear primary fields that have been added, only used in post-processing.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/31/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear_primary_fields(self)
        class(properties_t),    intent(inout)   :: self

        if (allocated(self%primary_fields)) then
            deallocate(self%primary_fields)
        end if

    end subroutine clear_primary_fields
    !******************************************************************************************



    !>  Clear model fields that have been added, only used in post-processing.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/31/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear_model_fields(self)
        class(properties_t),    intent(inout)   :: self

        if (allocated(self%model_fields)) then
            deallocate(self%model_fields)
        end if

    end subroutine clear_model_fields
    !******************************************************************************************



    !>  Clear adjoint fields that have been added.
    !!  Used in post-process and mesh-sensitivies computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/9/2018
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear_adjoint_fields(self)
        class(properties_t),    intent(inout)   :: self

        if (allocated(self%adjoint_fields)) then
            deallocate(self%adjoint_fields)
        end if

    end subroutine clear_adjoint_fields
    !******************************************************************************************







end module type_properties
