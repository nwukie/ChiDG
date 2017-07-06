module type_cache_data_field
#include <messenger.h>
    use mod_kinds,              only: ik
    use mod_constants,          only: INTERIOR, BOUNDARY, CHIMERA, NFACES, ZERO
    use type_mesh,              only: mesh_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t
    use DNAD_D
    implicit none





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, public :: cache_data_field_t

        character(:),   allocatable :: name

        type(AD_D),     allocatable :: value(:,:)           ! (nnodes, ndepend_value)
        type(seed_t),   allocatable :: value_seeds(:)


        type(AD_D),     allocatable :: gradient(:,:,:)      ! (nnodes, ndimension, ndepend_deriv)
        type(seed_t),   allocatable :: gradient_seeds(:)


        type(AD_D),     allocatable :: lift_face(:,:,:)     ! (nnodes_face, ndimension, ndepend_deriv)
        type(AD_D),     allocatable :: lift_element(:,:,:)  ! (nnodes_vol,  ndimension, ndepend_deriv)
        type(seed_t),   allocatable :: lift_seeds(:)

    contains

        procedure   :: resize
        procedure   :: set_data
        procedure   :: clear
        
        procedure   :: get_ndepend_face_exterior

    end type cache_data_field_t
    !*************************************************************************************




contains








    !>  Resize the cache storage for a given field.
    !!
    !!  Note: For cache_types = 'face interior' and 'face exterior', iface must be provided
    !!        in the subroutine call.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine resize(self,field,cache_component,mesh,prop,idomain_l,ielement_l,iface,differentiate)
        class(cache_data_field_t),  intent(inout)           :: self
        character(*),               intent(in)              :: field
        character(*),               intent(in)              :: cache_component
        type(mesh_t),               intent(in)              :: mesh
        type(properties_t),         intent(in)              :: prop(:)
        integer(ik),                intent(in)              :: idomain_l
        integer(ik),                intent(in)              :: ielement_l
        integer(ik),                intent(in), optional    :: iface
        logical,                    intent(in), optional    :: differentiate

        character(:),   allocatable :: user_msg
        integer(ik)                 :: nnodes, nnodes_vol, nnodes_face, ndepend_value, &
                                       ndepend_deriv, ierr, iface_loop, iseed
        logical                     :: reallocate


        !
        ! Set name
        !
        self%name = field



        !
        ! Get number of nodes(nnodes - face/element), number of elements that the function 
        ! value depends on(ndepend_value), number of elements that the function derivative 
        ! depends on(ndepend_deriv).
        !
        select case (trim(cache_component))


            case('element')
                nnodes      = mesh%domain(idomain_l)%elems(ielement_l)%basis_s%nnodes_ie()
                nnodes_vol  = mesh%domain(idomain_l)%elems(ielement_l)%basis_s%nnodes_ie()
                nnodes_face = mesh%domain(idomain_l)%elems(ielement_l)%basis_s%nnodes_if()

                ! Interior element
                !ndepend_value = 1
                ! The potential here is that model field values depend on exterior elements.
                ndepend_value = 1
                do iface_loop = 1,NFACES
                    ndepend_value = ndepend_value + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface_loop)
                end do

                ! Interior element + All Exterior Elements
                ndepend_deriv = 1
                do iface_loop = 1,NFACES
                    ndepend_deriv = ndepend_deriv + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface_loop)
                end do



            case('face interior')

                nnodes      = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_if()
                nnodes_vol  = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_ie()
                nnodes_face = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_if()

                ! Interior element + Face Exterior Elements
                !ndepend_value = 1
                ! The potential here is that model field values depend on exterior elements.
                ndepend_value = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)

                ! Interior element + Face Exterior Elements
                ndepend_deriv = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)


            case('face exterior')

                nnodes      = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_if()
                nnodes_vol  = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_ie()
                nnodes_face = mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%nnodes_if()

                ! Exterior Elements
                ndepend_value = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)

                ! Interior element + Face Exterior Elements
                ndepend_deriv = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)


            case default
                user_msg = "cache_data_field%resize: cache_type string wasn't an acceptable value. &
                            Make sure it is 'element', 'face interior', or 'face exterior'."
                call chidg_signal_one(FATAL,user_msg, cache_component)
        end select


        !
        ! Override ndepend_value, ndepend_deriv if we are not differentiating.
        ! Default if 'differentiate' is not present, continue with differentiation.
        !
        if (present(differentiate)) then
            if (.not. differentiate) then
                ndepend_value = 1
                ndepend_deriv = 1
            end if
        end if



        !
        ! Re/Allocate 'value' component
        !
        if (allocated(self%value)) then
            reallocate = (size(self%value,1) /= nnodes) .or. &
                         (size(self%value,2) /= ndepend_value)

            if (reallocate) deallocate(self%value, self%value_seeds)
        end if

        if (.not. allocated(self%value)) then
            allocate(self%value(nnodes,ndepend_value), &
                     self%value_seeds(ndepend_value), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if



        !
        ! Re/Allocate 'gradient', 'lift' components
        !
        if (allocated(self%gradient)) then
            reallocate = ( (size(self%gradient,1) /= nnodes)       .or. &
                           (size(self%gradient,2) /= 3)            .or. &
                           (size(self%gradient,3) /= ndepend_deriv) )

            if (reallocate) deallocate(self%gradient, self%gradient_seeds, &
                                       self%lift_face, self%lift_element, self%lift_seeds)
        end if


        if ( .not. allocated(self%gradient) ) then
            allocate(self%gradient(nnodes,3,ndepend_deriv),       &
                     self%gradient_seeds(ndepend_deriv),          &
                     self%lift_face(nnodes_face,3,ndepend_deriv),   &
                     self%lift_element(nnodes_vol,3,ndepend_deriv), &
                     self%lift_seeds(ndepend_deriv),    stat=ierr)
            if (ierr /= 0) call AllocationError
        end if



        !
        ! Clear entries
        !
        call self%clear()



    end subroutine resize
    !*************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine set_data(self,cache_data,data_type,idirection,seed)
        class(cache_data_field_t),  intent(inout)   :: self
        type(AD_D),                 intent(in)      :: cache_data(:)
        character(*),               intent(in)      :: data_type
        integer(ik),                intent(in)      :: idirection
        type(seed_t),               intent(in)      :: seed

        character(:), allocatable   :: user_msg
        logical                     :: seed_found, empty_seed
        integer(ik)                 :: iseed, seed_location


        

        select case (trim(data_type))

            !
            ! Set variable 'value' data
            !
            case('value')
                ! Search to see if a value differentiated wrt seed already exists
                seed_location = 0
                seed_found = .false.
                do iseed = 1,size(self%value_seeds)
                    seed_found = ( (seed%idomain_g  == self%value_seeds(iseed)%idomain_g) .and. &
                                   (seed%ielement_g == self%value_seeds(iseed)%ielement_g) )

                    if (seed_found) then
                        seed_location = iseed
                        exit
                    end if
                end do


                ! If matching seed was not found, find first empty seed location and place there
                if (.not. seed_found) then
                    do iseed = 1,size(self%value_seeds)
                        empty_seed = (self%value_seeds(iseed)%idomain_g == 0)

                        if (empty_seed) then
                            seed_location = iseed
                            exit
                        end if

                    end do
                end if

                
                user_msg = "cache_data_field%set_data: Did not find a location to put the data."
                if (seed_location == 0) call chidg_signal(FATAL,user_msg)

                ! Store data
                self%value(:,seed_location) = cache_data
                self%value_seeds(seed_location) = seed



            !
            ! Set variable 'gradient' data
            !
            case('gradient')
                ! Search to see if a value differentiated wrt seed already exists
                seed_location = 0
                seed_found = .false.
                do iseed = 1,size(self%gradient_seeds)
                    seed_found = ( (seed%idomain_g  == self%gradient_seeds(iseed)%idomain_g) .and. &
                                   (seed%ielement_g == self%gradient_seeds(iseed)%ielement_g) )

                    if (seed_found) then
                        seed_location = iseed
                        exit
                    end if
                end do


                ! If matching seed was not found, find first empty seed location and place there
                if (.not. seed_found) then
                    do iseed = 1,size(self%gradient_seeds)
                        empty_seed = (self%gradient_seeds(iseed)%idomain_g == 0)

                        if (empty_seed) then
                            seed_location = iseed
                            exit
                        end if

                    end do
                end if

                
                user_msg = "cache_data_field%set_data: Did not find a location to put the data."
                if (seed_location == 0) call chidg_signal(FATAL,user_msg)


                ! Store data
                self%gradient(:,idirection,seed_location) = cache_data
                self%gradient_seeds(seed_location) = seed


            !
            ! Set variable 'lift' data
            !
            case('lift face')
                ! Search to see if a value differentiated wrt seed already exists
                seed_location = 0
                seed_found = .false.
                do iseed = 1,size(self%lift_seeds)
                    seed_found = ( (seed%idomain_g  == self%lift_seeds(iseed)%idomain_g) .and. &
                                   (seed%ielement_g == self%lift_seeds(iseed)%ielement_g) )

                    if (seed_found) then
                        seed_location = iseed
                        exit
                    end if
                end do



                ! If matching seed was not found, find first empty seed location and place there
                if (.not. seed_found) then
                    do iseed = 1,size(self%lift_seeds)
                        empty_seed = (self%lift_seeds(iseed)%idomain_g == 0)

                        if (empty_seed) then
                            seed_location = iseed
                            exit
                        end if

                    end do
                end if


                user_msg = "cache_data_field%set_data: Did not find a location to put the data."
                if (seed_location == 0) call chidg_signal(FATAL,user_msg)
                
                ! Store data
                self%lift_face(:,idirection,seed_location) = cache_data
                self%lift_seeds(seed_location) = seed



            case('lift element')
                ! Search to see if a value differentiated wrt seed already exists
                seed_location = 0
                seed_found = .false.
                do iseed = 1,size(self%lift_seeds)
                    seed_found = ( (seed%idomain_g  == self%lift_seeds(iseed)%idomain_g) .and. &
                                   (seed%ielement_g == self%lift_seeds(iseed)%ielement_g) )

                    if (seed_found) then
                        seed_location = iseed
                        exit
                    end if
                end do



                ! If matching seed was not found, find first empty seed location and place there
                if (.not. seed_found) then
                    do iseed = 1,size(self%lift_seeds)
                        empty_seed = (self%lift_seeds(iseed)%idomain_g == 0)

                        if (empty_seed) then
                            seed_location = iseed
                            exit
                        end if

                    end do
                end if


                user_msg = "cache_data_field%set_data: Did not find a location to put the data."
                if (seed_location == 0) call chidg_signal(FATAL,user_msg)
                
                ! Store data
                self%lift_element(:,idirection,seed_location) = cache_data
                self%lift_seeds(seed_location) = seed





            case default
                user_msg = "cache_data_field%store: The incoming variable data_type did &
                            not have an valid value. Acceptable entries are 'value', &
                            'gradient', 'lift face', or 'lift element'"
                call chidg_signal_one(FATAL,user_msg,data_type)

        end select


    end subroutine set_data
    !************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    function get_ndepend_face_exterior(self,mesh,idomain_l,ielement_l,iface) result(ndepend)
        class(cache_data_field_t),  intent(in)  :: self
        type(mesh_t),               intent(in)  :: mesh
        integer(ik),                intent(in)  :: idomain_l
        integer(ik),                intent(in)  :: ielement_l
        integer(ik),                intent(in)  :: iface

        character(:),   allocatable :: user_msg
        integer(ik)                 :: ndepend, ChiID, group_ID, patch_ID, face_ID
        logical                     :: conforming_face, chimera_face, boundary_face


        conforming_face = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
        chimera_face    = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
        boundary_face   = (mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)

        if (conforming_face) then
            ndepend = 1   ! Exterior element


        else if (chimera_face) then
            ChiID   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID
            ndepend = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%ndonors()


        else if (boundary_face) then
            group_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
            patch_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
            face_ID  = mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID
            ndepend  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

        else
            user_msg = "cache_data_field%get_ndepend_face_exterior: Invalid face type detected."
            call chidg_signal(FATAL,user_msg)

        end if



    end function get_ndepend_face_exterior
    !************************************************************************************



    
    !>  Clear cache data field entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/3/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine clear(self)
        class(cache_data_field_t),  intent(inout)   :: self

        integer(ik) :: iseed


!        !
!        ! Zero values, gradients, lift
!        !
!        if (allocated(self%value)) then
!            self%value = ZERO
!        end if
!
!        if (allocated(self%gradient)) then
!            self%gradient = ZERO
!        end if
!
!        if (allocated(self%lift_face)) then
!            self%lift_face = ZERO
!        end if
!
!        if (allocated(self%lift_element)) then
!            self%lift_element = ZERO
!        end if



        !
        ! Clear seed contents, zero values
        !
        do iseed = 1,size(self%value_seeds)
            call self%value_seeds(iseed)%clear()
        end do

        do iseed = 1,size(self%gradient_seeds)
            call self%gradient_seeds(iseed)%clear()
        end do

        do iseed = 1,size(self%lift_seeds)
            call self%lift_seeds(iseed)%clear()
        end do



    end subroutine clear
    !***********************************************************************************










end module type_cache_data_field
