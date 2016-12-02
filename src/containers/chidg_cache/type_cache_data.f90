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

        integer(ik)                 :: neqns, ierr, ChiID, donor_idomain, ifield
        logical                     :: interior_face, chimera_face, boundary_face
        character(:),   allocatable :: user_msg


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
                neqns = mesh(idomain_l)%elems(ielement_l)%neqns


            case ('face interior')
                neqns = mesh(idomain_l)%elems(ielement_l)%neqns


            case ('face exterior')

                interior_face = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
                chimera_face  = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
                boundary_face = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)

                if (interior_face) then
                    neqns = mesh(idomain_l)%neqns
                else if (boundary_face) then
                    neqns = mesh(idomain_l)%neqns
                else if (chimera_face) then
                    ChiID = mesh(idomain_l)%faces(ielement_l,iface)%ChiID
                    neqns = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_neqns%at(1)
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
        ! Re/Allocate equations
        !
        if (allocated(self%fields)) then

            if (size(self%fields) /= neqns) then
                deallocate(self%fields)
                allocate(self%fields(neqns), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%fields(neqns), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if



        !
        ! Resize each equation storage
        !
        do ifield = 1,size(self%fields)
            call self%fields(ifield)%resize(cache_component,mesh,prop,idomain_l,ielement_l,iface)
        end do !ifield



    end subroutine resize
    !*******************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    subroutine set_data(self,cache_data,data_type,idirection,seed,ifield)
        class(cache_data_t),    intent(inout)   :: self
        type(AD_D),             intent(in)      :: cache_data(:)
        character(*),           intent(in)      :: data_type
        integer(ik),            intent(in)      :: idirection
        type(seed_t),           intent(in)      :: seed
        integer(ik),            intent(in)      :: ifield


        call self%fields(ifield)%set_data(cache_data,data_type,idirection,seed)


    end subroutine set_data
    !*******************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_data(self,cache_type,idirection,seed,ifield) result(cache_data)
        class(cache_data_t),    intent(inout)           :: self
        character(*),           intent(in)              :: cache_type
        integer(ik),            intent(in)              :: idirection
        type(seed_t),           intent(in)              :: seed
        integer(ik),            intent(in)              :: ifield

        
        type(AD_D), allocatable, dimension(:) :: cache_data
        
        integer(ik)                 :: igq, iseed
        logical                     :: seed_found, has_seed
        character(:),   allocatable :: user_msg



        !
        ! Try to find data for ifield that was linearized wrt seed
        !
        select case (trim(cache_type))
            case('value')

                ! Try to find value differentiated wrt seed.
                !
                seed_found = .false.
                do iseed = 1,size(self%fields(ifield)%value_seeds)
                    
                    has_seed = (self%fields(ifield)%value_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(ifield)%value_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(ifield)%value(:,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! values and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(ifield)%value(:,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if



            case('derivative')

                ! Try to find derivative differentiated wrt seed.
                !
                seed_found = .false.
                do iseed = 1,size(self%fields(ifield)%derivative_seeds)
                    
                    has_seed = (self%fields(ifield)%derivative_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(ifield)%derivative_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(ifield)%derivative(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! derivative and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(ifield)%derivative(:,idirection,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if



            case('lift face')

                ! Try to find lift differentiated wrt seed.
                !
                seed_found = .false.
                do iseed = 1,size(self%fields(ifield)%lift_seeds)
                    
                    has_seed = (self%fields(ifield)%lift_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(ifield)%lift_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(ifield)%lift_face(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! lift and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(ifield)%lift_face(:,idirection,1)

                    do igq = 1,size(cache_data)
                        cache_data(igq)%xp_ad_ = ZERO
                    end do
                end if


            case('lift element')

                ! Try to find lift differentiated wrt seed.
                !
                seed_found = .false.
                do iseed = 1,size(self%fields(ifield)%lift_seeds)
                    
                    has_seed = (self%fields(ifield)%lift_seeds(iseed)%idomain_g  == seed%idomain_g) .and. &
                               (self%fields(ifield)%lift_seeds(iseed)%ielement_g == seed%ielement_g)

                    if (has_seed) then
                        cache_data = self%fields(ifield)%lift_element(:,idirection,iseed)
                        seed_found = .true.
                        exit
                    end if

                end do

                ! If the current component doesn't have a linearization wrt seed, just take any
                ! lift and zero out autodiff
                if (.not. seed_found) then
                    cache_data = self%fields(ifield)%lift_element(:,idirection,1)

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



















end module type_cache_data
