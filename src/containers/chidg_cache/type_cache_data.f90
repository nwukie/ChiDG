module type_cache_data
#include <messenger.h>
    use mod_kinds,                  only: ik, rk
    use mod_constants,              only: INTERIOR, BOUNDARY, CHIMERA, ZERO
    use type_cache_data_field,      only: cache_data_field_t
    use type_mesh,                  only: mesh_t
    use type_properties,            only: properties_t
    use type_seed,                  only: seed_t
    use DNAD_D
    implicit none





    !>  Cache data. Includes data differentiated with respect to all dependent elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    type, public :: cache_data_t

        type(cache_data_field_t),   allocatable :: fields(:)

    contains

        procedure   :: resize
        procedure   :: set_data
        procedure   :: get_data

        procedure   :: get_field_index

    end type cache_data_t
    !*******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    subroutine resize(self,cache_component,mesh,prop,idomain_l,ielement_l,iface)
        class(cache_data_t),    intent(inout)           :: self
        character(*),           intent(in)              :: cache_component
        type(mesh_t),           intent(in)              :: mesh(:)
        type(properties_t),     intent(in)              :: prop(:)
        integer(ik),            intent(in)              :: idomain_l
        integer(ik),            intent(in)              :: ielement_l
        integer(ik),            intent(in), optional    :: iface

        integer(ik)                 :: nprimary_fields,nmodel_fields, ntotal_fields, ierr, &
                                       ChiID, donor_idomain, ifield, iprimary_field, imodel_field, &
                                       nauxiliary_fields, iauxiliary_field
        logical                     :: interior_face, chimera_face, boundary_face, reallocate
        character(:),   allocatable :: field, user_msg


        !
        ! Check if iface was provided for face-type caches
        !
        if ((trim(cache_component) == 'face interior') .or. &
            (trim(cache_component) == 'face exterior')) then
            if (.not. present(iface)) then
                user_msg = "chidg_data%resize: Tried to resize face cache, but &
                            no face index was specified. Try providing iface to the call"
                call chidg_signal(FATAL,user_msg)
            end if
        end if
        


        select case (trim(cache_component))
            case ('element')
                nprimary_fields   = prop(idomain_l)%nprimary_fields()
                nauxiliary_fields = prop(idomain_l)%nauxiliary_fields()
                nmodel_fields     = prop(idomain_l)%nmodel_fields()
            case ('face interior')
                nprimary_fields   = prop(idomain_l)%nprimary_fields()
                nauxiliary_fields = prop(idomain_l)%nauxiliary_fields()
                nmodel_fields     = prop(idomain_l)%nmodel_fields()


            case ('face exterior')

                interior_face = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
                chimera_face  = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
                boundary_face = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)

                if (interior_face .or. boundary_face) then
                    nprimary_fields   = prop(idomain_l)%nprimary_fields()
                    nauxiliary_fields = prop(idomain_l)%nauxiliary_fields()
                    nmodel_fields     = prop(idomain_l)%nmodel_fields()
                else if (chimera_face) then
                    ChiID = mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                    ! To handle different fields on either side of the chimera boundary,
                    ! we will have to store the properties_t of the chimera donors so
                    ! we can query them here. For now, just use the source domain
                    ! and assume they have the same fields.
                    !neqns = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_neqns%at(1)
                    nprimary_fields   = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_neqns%at(1)
                    nauxiliary_fields = prop(idomain_l)%nauxiliary_fields()
                    nmodel_fields     = prop(idomain_l)%nmodel_fields()
                else
                    call chidg_signal(FATAL,"cache_data: Face type wasn't recognized")
                end if



            case default
                user_msg = "chidg_data%resize: An invalid value for the parameter 'type' &
                            from cache_info was returned in the accept call. cache_type &
                            is supposed to be either 'element', 'face interior', or &
                            'face exterior',  to indicate the cache type where the data &
                            is to be stored."
                call chidg_signal_one(FATAL,user_msg,cache_component)

        end select

        !
        ! Compute total number of fields to store
        !
        ntotal_fields = nprimary_fields + nauxiliary_fields + nmodel_fields


        !
        ! Re/Allocate equations
        !
        if (allocated(self%fields)) then
            reallocate = size(self%fields) /= ntotal_fields
            if (reallocate) deallocate(self%fields)
        end if

        if (.not. allocated(self%fields)) then
            allocate(self%fields(ntotal_fields), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if


        ! Resize primary fields
        ifield = 1
        do iprimary_field = 1,prop(idomain_l)%nprimary_fields()
            field = prop(idomain_l)%get_primary_field_name(iprimary_field)
            call self%fields(ifield)%resize(field,cache_component,mesh,prop,idomain_l,ielement_l,iface)
            ifield = ifield + 1
        end do


        ! Resize model fields
        do imodel_field = 1,prop(idomain_l)%nmodel_fields()
            field = prop(idomain_l)%get_model_field_name(imodel_field)
            call self%fields(ifield)%resize(field,cache_component,mesh,prop,idomain_l,ielement_l,iface)
            ifield = ifield + 1
        end do


        ! Resize auxiliary fields
        do iauxiliary_field = 1,prop(idomain_l)%nauxiliary_fields()
            field = prop(idomain_l)%get_auxiliary_field_name(iauxiliary_field)
            call self%fields(ifield)%resize(field,cache_component,mesh,prop,idomain_l,ielement_l,iface)
            ifield = ifield + 1
        end do


    end subroutine resize
    !*******************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    subroutine set_data(self,field,cache_data,data_type,idirection,seed)
        class(cache_data_t),    intent(inout)   :: self
        character(*),           intent(in)      :: field
        type(AD_D),             intent(in)      :: cache_data(:)
        character(*),           intent(in)      :: data_type
        integer(ik),            intent(in)      :: idirection
        type(seed_t),           intent(in)      :: seed

        integer(ik) :: field_index

        ! Get field index
        field_index = self%get_field_index(field)

        call self%fields(field_index)%set_data(cache_data,data_type,idirection,seed)


    end subroutine set_data
    !*******************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_data(self,field,cache_type,idirection,seed) result(cache_data)
        class(cache_data_t),    intent(inout)           :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: cache_type
        integer(ik),            intent(in)              :: idirection
        type(seed_t),           intent(in)              :: seed

        
        type(AD_D), allocatable, dimension(:) :: cache_data
        
        integer(ik)                 :: igq, iseed, field_index
        logical                     :: seed_found, has_seed
        character(:),   allocatable :: user_msg


        ! Get field index
        field_index = self%get_field_index(field)

        !
        ! Try to find data for field_index that was linearized wrt seed
        !
        select case (trim(cache_type))
            case('value')

                ! Try to find value differentiated wrt seed.
                seed_found = .false.
                do iseed = 1,size(self%fields(field_index)%value_seeds)
                    
                    has_seed = (self%fields(field_index)%value_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(field_index)%value_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(field_index)%value(:,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! values and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(field_index)%value(:,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if



            case('derivative')

                ! Try to find derivative differentiated wrt seed.
                seed_found = .false.
                do iseed = 1,size(self%fields(field_index)%derivative_seeds)
                    
                    has_seed = (self%fields(field_index)%derivative_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(field_index)%derivative_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(field_index)%derivative(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! derivative and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(field_index)%derivative(:,idirection,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if



            case('lift face')

                ! Try to find lift differentiated wrt seed.
                seed_found = .false.
                do iseed = 1,size(self%fields(field_index)%lift_seeds)
                    
                    has_seed = (self%fields(field_index)%lift_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(field_index)%lift_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(field_index)%lift_face(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! lift and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(field_index)%lift_face(:,idirection,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if


            case('lift element')

                ! Try to find lift differentiated wrt seed.
                seed_found = .false.
                do iseed = 1,size(self%fields(field_index)%lift_seeds)
                    
                    has_seed = (self%fields(field_index)%lift_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(field_index)%lift_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(field_index)%lift_element(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! lift and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(field_index)%lift_element(:,idirection,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if

            case default

                user_msg = "cache_data%get_data: Invalid data type for getting data. Options are &
                            'value', 'derivative', 'lift face', or 'lift element'"
                call chidg_signal_one(FATAL,user_msg,trim(cache_type))


        end select



        if (.not. allocated(cache_data(1)%xp_ad_) ) call chidg_signal(FATAL,"get_data: derivatives not allocated")


    end function get_data
    !*********************************************************************************






    !>  Given a field name, find the index of the field in the cache.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/2/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    function get_field_index(self,field) result(field_index)
        class(cache_data_t),    intent(in)  :: self
        character(*),           intent(in)  :: field

        integer(ik) :: ifield, field_index

        field_index = 0
        do ifield = 1,size(self%fields)
            if (trim(self%fields(ifield)%name) == trim(field)) then
                field_index = ifield
                exit
            end if
        end do

    end function get_field_index
    !*********************************************************************************











end module type_cache_data
