module type_chidg_cache
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: NFACES, CACHE_FACE_INTERIOR, CACHE_FACE_EXTERIOR, NO_ID
    use type_mesh,      only: mesh_t
    use type_properties,    only: properties_t
    use type_seed,          only: seed_t
    use DNAD_D

    use type_cache_data,    only: cache_data_t
    implicit none


    !>  A cache container that precomputes data at the quadrature node sets that
    !!  can then be used by all the element/face operators for a given element.
    !!
    !!  Before any of the element/face operators are called for a particular element, 
    !!  a chidg_cache instance is initialized for the element. This involved computing
    !!  the interpolation of solution data from basis modes to the quadrature nodes
    !!  and also the linearization of those processes. It also includes calculation
    !!  of the diffusion lifting operators and interpolation of the lifting modes
    !!  to the quadrature node sets.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    type, public :: chidg_cache_t

        !type(cache_data_t)              :: element
        type(cache_data_t), allocatable :: element
        type(cache_data_t), allocatable :: faces(:,:)      ! (nfaces, side). side=1(INTERIOR), size=2(EXTERIOR)
        !type(cache_data_t)             :: faces(NFACES,2) ! Causes segfault due to incomplete compiler finalization support.

        logical :: lift

    contains

        procedure   :: resize

        procedure   :: set_data
        procedure   :: get_data

        procedure   :: check_field_exists

    end type chidg_cache_t
    !**************************************************************************************


contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine resize(self,mesh,prop,idomain_l,ielement_l,differentiate,lift)
        class(chidg_cache_t),   intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh
        type(properties_t),     intent(in)      :: prop(:)
        integer(ik),            intent(in)      :: idomain_l
        integer(ik),            intent(in)      :: ielement_l
        logical,                intent(in)      :: differentiate
        logical,                intent(in)      :: lift

        integer(ik) :: iface, ierr

        if (allocated(self%element)) deallocate(self%element)
        if (allocated(self%faces)) deallocate(self%faces)


        !
        ! Store logical indicating if lift is stored
        !
        self%lift = lift



        !
        ! Allocate face cache's. Fixes SegFault that occurs when these are declared as static 
        ! arrays at compile time (faces(NFACES,2)) because gfortran doesn't have complete
        ! finalization procedures implemented yet.
        !
        if (.not. allocated(self%faces)) then
            allocate(self%element, self%faces(NFACES,2), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if


        !
        ! Allocate storage for element cache
        !
        call self%element%resize('element',mesh,prop,idomain_l,ielement_l,differentiate=differentiate)



        !
        ! Allocate storage for faces cache
        !
        do iface = 1,size(self%faces,1)

            call self%faces(iface,CACHE_FACE_INTERIOR)%resize('face interior',mesh,prop,idomain_l,ielement_l,iface,differentiate=differentiate)
            call self%faces(iface,CACHE_FACE_EXTERIOR)%resize('face exterior',mesh,prop,idomain_l,ielement_l,iface,differentiate=differentiate)

        end do


    end subroutine resize
    !**************************************************************************************





    !>  Accept externally computed data to be stored in the cache.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_data(self,field,cache_component,cache_data,data_type,idirection,seed,iface)
        class(chidg_cache_t),   intent(inout)           :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: cache_component
        type(AD_D),             intent(in)              :: cache_data(:)
        character(*),           intent(in)              :: data_type
        integer(ik),            intent(in)              :: idirection
        type(seed_t),           intent(in)              :: seed
        integer(ik),            intent(in), optional    :: iface


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




        !
        ! Call accept routine on correct data component
        !
        select case(cache_component)
            case('element')
                call self%element%set_data(field,cache_data,data_type,idirection,seed)

            case('face interior')
                call self%faces(iface,1)%set_data(field,cache_data,data_type,idirection,seed)

            case('face exterior')
                call self%faces(iface,2)%set_data(field,cache_data,data_type,idirection,seed)

            case default 
                user_msg = "chidg_cache%set_data: An invalid value for the cache_component incoming parameter  &
                            Valid values are either 'element', 'face interior', or 'face exterior' to indicate &
                            the cache type where the data is to be stored."
                call chidg_signal_one(FATAL,user_msg,cache_component)

        end select




    end subroutine set_data
    !****************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_data(self,field,cache_component,cache_type,idirection,seed,iface) result(cache_data)
        class(chidg_cache_t),   intent(inout)           :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: cache_component
        character(*),           intent(in)              :: cache_type
        integer(ik),            intent(in)              :: idirection
        type(seed_t),           intent(in)              :: seed
        integer(ik),            intent(in), optional    :: iface


        type(AD_D), allocatable, dimension(:) :: cache_data



        select case (trim(cache_component))
            case ('element')
                cache_data = self%element%get_data(field,cache_type,idirection,seed)

            case ('face interior')
                cache_data = self%faces(iface,1)%get_data(field,cache_type,idirection,seed)

            case ('face exterior')
                cache_data = self%faces(iface,2)%get_data(field,cache_type,idirection,seed)

            case default
                call chidg_signal(FATAL,'chidg_cache%get_data: Error in cache_component string')
        end select



    end function get_data
    !*****************************************************************************************







    !>  Check if field exists in cache.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!  TODO: test
    !!
    !----------------------------------------------------------------------------------------
    function check_field_exists(self,field,cache_component,iface) result(exists)
        class(chidg_cache_t),   intent(inout)           :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: cache_component
        integer(ik),            intent(in), optional    :: iface

        logical :: exists

        select case (trim(cache_component))
            case ('element')
                exists = (self%element%get_field_index(field) /= NO_ID)
            case ('face interior')
                exists = (self%faces(iface,1)%get_field_index(field) /= NO_ID)
            case ('face exterior')
                exists = (self%faces(iface,2)%get_field_index(field) /= NO_ID)
            case default
                call chidg_signal(FATAL,'chidg_cache%check_field_exists: Error in cache_component string')
        end select



    end function check_field_exists
    !*****************************************************************************************












end module type_chidg_cache
