module type_cache_data_equation
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_constants,      only: INTERIOR, BOUNDARY, CHIMERA, NFACES
    use type_mesh,          only: mesh_t
    use type_seed,          only: seed_t
    use DNAD_D
    implicit none





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/6/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, public :: cache_data_equation_t

        type(AD_D),     allocatable :: value(:,:)        ! (nnodes, ndepend_value)
        type(seed_t),   allocatable :: value_seeds(:)


        type(AD_D),     allocatable :: derivative(:,:,:) ! (nnodes, nderiv_dimension, ndepend_deriv)
        type(seed_t),   allocatable :: derivative_seeds(:)


        type(AD_D),     allocatable :: lift(:,:,:)        ! (nnodes, nderiv_dimension, ndepend_deriv)
        type(seed_t),   allocatable :: lift_seeds(:)

    contains

        procedure   :: resize
        procedure   :: set_data
        
        procedure   :: get_ndepend_face_exterior

    end type cache_data_equation_t
    !*************************************************************************************




contains








    !>  Resize the cache storage for a given equation.
    !!
    !!  Note: For cache_types = 'face interior' and 'face exterior', iface must be provided
    !!        in the subroutine call.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/7/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine resize(self,cache_component,mesh,idomain_l,ielement_l,iface)
        class(cache_data_equation_t),   intent(inout)           :: self
        character(len=*),               intent(in)              :: cache_component
        type(mesh_t),                   intent(in)              :: mesh(:)
        integer(ik),                    intent(in)              :: idomain_l
        integer(ik),                    intent(in)              :: ielement_l
        integer(ik),                    intent(in), optional    :: iface

        integer(ik) :: nnodes, ndepend_value, ndepend_deriv, ierr, iface_loop, iseed
        logical     :: reallocate



        !
        ! Get number of nodes(nnodes - face/element), number of elements that the function value depends on(ndepend_value), number of 
        ! elements that the function derivative depends on(ndepend_deriv).
        !
        select case (trim(cache_component))


            case('element')
                nnodes = mesh(idomain_l)%elems(ielement_l)%gq%vol%nnodes

                ! Interior element
                ndepend_value = 1

                ! Interior element + All Exterior Elements
                ndepend_deriv = 1
                do iface_loop = 1,NFACES
                    ndepend_deriv = ndepend_deriv + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface_loop)
                end do



            case('face interior')

                nnodes = mesh(idomain_l)%faces(ielement_l,iface)%gq%face%nnodes

                ! Interior element
                ndepend_value = 1

                ! Interior element + Exterior Elements
                ndepend_deriv = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)

                ! Interior element + Exterior Elements


            case('face exterior')
                nnodes = mesh(idomain_l)%faces(ielement_l,iface)%gq%face%nnodes

                ! Exterior Elements
                ndepend_value = self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)

                ! Interior element + Exterior Elements
                ndepend_deriv = 1 + self%get_ndepend_face_exterior(mesh,idomain_l,ielement_l,iface)


            case default
                call chidg_signal_one(FATAL,"cache_data_equation%resize: cache_type string wasn't an acceptable value. &
                                                make sure it is either 'element', 'face interior', or 'face exterior'.", cache_component)
        end select



        !
        ! Re/Allocate 'value' component
        !
        if (allocated(self%value)) then

            reallocate = (size(self%value,1) /= nnodes) .or. (size(self%value,2) /= ndepend_value)
            if (reallocate) then
                deallocate(self%value, self%value_seeds)
                allocate(self%value(nnodes,ndepend_value), self%value_seeds(ndepend_value), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%value(nnodes,ndepend_value), self%value_seeds(ndepend_value), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if





        !
        ! Re/Allocate 'derivative', 'lift' components
        !
        if (allocated(self%derivative)) then

            reallocate = ( (size(self%derivative,1) /= nnodes)       .or. &
                           (size(self%derivative,2) /= 3)            .or. &
                           (size(self%derivative,3) /= ndepend_deriv) )

            if (reallocate) then
                deallocate(self%derivative, self%lift, self%derivative_seeds, self%lift_seeds)
                allocate(self%derivative(nnodes,3,ndepend_deriv), &
                         self%lift(      nnodes,3,ndepend_deriv), &
                         self%derivative_seeds(ndepend_deriv),    &
                         self%lift_seeds(      ndepend_deriv), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%derivative(nnodes,3,ndepend_deriv), &
                     self%lift(      nnodes,3,ndepend_deriv), &
                     self%derivative_seeds(ndepend_deriv),    &
                     self%lift_seeds(ndepend_deriv), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if





        !
        ! Clear seed contents
        !
        do iseed = 1,size(self%value_seeds)
            call self%value_seeds(iseed)%clear()
        end do

        do iseed = 1,size(self%derivative_seeds)
            call self%derivative_seeds(iseed)%clear()
        end do

        do iseed = 1,size(self%lift_seeds)
            call self%lift_seeds(iseed)%clear()
        end do


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
        class(cache_data_equation_t),   intent(inout)   :: self
        type(AD_D),                     intent(in)      :: cache_data(:)
        character(len=*),               intent(in)      :: data_type
        integer(ik),                    intent(in)      :: idirection
        type(seed_t),                   intent(in)      :: seed

        character(len=:), allocatable   :: msg
        logical     :: seed_found, empty_seed
        integer(ik) :: iseed, seed_location


        

        select case (trim(data_type))

            !
            ! Set variable 'value' data
            !
            case('value')
                ! Search to see if a value differentiated wrt seed already exist
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

                if (seed_location == 0) call chidg_signal(FATAL,"cache_data_equation%set_data: Did not find a location to put the data")

                ! Store data
                self%value(:,seed_location) = cache_data
                self%value_seeds(seed_location) = seed



            !
            ! Set variable 'derivative' data
            !
            case('derivative')
                ! Search to see if a value differentiated wrt seed already exists
                seed_location = 0
                seed_found = .false.
                do iseed = 1,size(self%derivative_seeds)
                    seed_found = ( (seed%idomain_g  == self%derivative_seeds(iseed)%idomain_g) .and. &
                                   (seed%ielement_g == self%derivative_seeds(iseed)%ielement_g) )

                    if (seed_found) then
                        seed_location = iseed
                        exit
                    end if
                end do


                ! If matching seed was not found, find first empty seed location and place there
                if (.not. seed_found) then
                    do iseed = 1,size(self%derivative_seeds)
                        empty_seed = (self%derivative_seeds(iseed)%idomain_g == 0)

                        if (empty_seed) then
                            seed_location = iseed
                            exit
                        end if

                    end do
                end if

                
                if (seed_location == 0) call chidg_signal(FATAL,"cache_data_equation%set_data: Did not find a location to put the data")


                ! Store data
                self%derivative(:,idirection,seed_location) = cache_data
                self%derivative_seeds(seed_location) = seed


            !
            ! Set variable 'lift' data
            !
            case('lift')
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


                if (seed_location == 0) call chidg_signal(FATAL,"cache_data_equation%set_data: Did not find a location to put the data")
                
                ! Store data
                self%lift(:,idirection,seed_location) = cache_data
                self%lift_seeds(seed_location) = seed





            case default
                msg = "cache_data_equation%store: The incoming variable data_type did not have an &
                             valid value. Acceptable entries are 'value', 'derivative', or 'lift'"
                call chidg_signal_one(FATAL,msg,data_type)

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
        class(cache_data_equation_t),   intent(in)  :: self
        type(mesh_t),                   intent(in)  :: mesh(:)
        integer(ik),                    intent(in)  :: idomain_l
        integer(ik),                    intent(in)  :: ielement_l
        integer(ik),                    intent(in)  :: iface

        integer(ik) :: ndepend, ChiID
        logical     :: conforming_face, chimera_face, boundary_face


        conforming_face = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR)
        chimera_face    = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA )
        boundary_face   = (mesh(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY)

        if (conforming_face) then
            ndepend = 1   ! Exterior element


        else if (chimera_face) then
            ChiID   = mesh(idomain_l)%faces(ielement_l,iface)%ChiID
            ndepend = mesh(idomain_l)%chimera%recv%data(ChiID)%ndonors()


        else if (boundary_face) then
            ndepend = mesh(idomain_l)%faces(ielement_l,iface)%BC_ndepend

        end if



    end function get_ndepend_face_exterior
    !************************************************************************************














end module type_cache_data_equation
