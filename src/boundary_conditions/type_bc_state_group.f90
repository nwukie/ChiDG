module type_bc_state_group
#include <messenger.h>
    use mod_constants,          only: NO_ID
    use type_bcvector,          only: bcvector_t
    use type_bc_state,          only: bc_state_t
    use type_bc_state_wrapper,  only: bc_state_wrapper_t
    use type_mesh,              only: mesh_t

    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use mpi_f08,                only: mpi_comm, MPI_LOGICAL, MPI_Comm_split, MPI_Allgather, MPI_UNDEFINED
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: bc_state_group_t
        
        character(:),               allocatable :: name         ! boundary state group name
        character(:),               allocatable :: family       ! boundary state group family

        class(bc_state_wrapper_t),  allocatable :: bc_state(:)  ! boundary state functions

        type(mpi_comm)                          :: bc_COMM      ! MPI communicator for bc procs
        integer(ik)                             :: bc_ID        ! bc state group identifier

    contains

        procedure   :: set_name
        procedure   :: get_name
        procedure   :: set_family
        procedure   :: get_family


        procedure   :: add_bc_state
        procedure   :: new_bc_state
        procedure   :: remove_states
        procedure   :: nbc_states


        procedure   :: init_comm
        procedure   :: init_coupling
        procedure   :: init_specialized
        procedure   :: propagate_coupling



    end type bc_state_group_t
    !*****************************************************************************************



contains



    !>  Set the bc_state_group name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_name(self,name)
        class(bc_state_group_t),    intent(inout)   :: self
        character(*),               intent(in)      :: name

        self%name = name

    end subroutine set_name
    !******************************************************************************************





    !>  Return bc_state_group name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    function get_name(self) result(name)
        class(bc_state_group_t),  intent(inout)   :: self

        character(:),   allocatable :: name, user_msg

        if (allocated(self%name)) then
            name = self%name
        else
            user_msg = "bc_state_group%get_name: It looks like the boundary condition group &
                        name was never set. Make sure bc_group%set_name gets called in the &
                        boundary condition group initialization routine"
            call chidg_signal(FATAL,user_msg)
        end if

    end function get_name
    !******************************************************************************************





    !>  Set the bc_group family.
    !!
    !!  bc_family may be:
    !!      - Wall
    !!      - Inlet
    !!      - Outlet
    !!      - Symmetry
    !!      - Periodic
    !!      - Farfield
    !!      - Scalar
    !!      - Extrapolation
    !!      - Empty
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_family(self,family)
        class(bc_state_group_t),  intent(inout)   :: self
        character(*),       intent(in)      :: family

        character(:),   allocatable :: user_msg

        if ( (trim(family) == 'Wall'    )       .or. &
             (trim(family) == 'Inlet'   )       .or. &
             (trim(family) == 'Outlet'  )       .or. &
             (trim(family) == 'Symmetry')       .or. &
             (trim(family) == 'Periodic')       .or. &
             (trim(family) == 'Farfield')       .or. &
             (trim(family) == 'Scalar'  )       .or. &
             (trim(family) == 'Extrapolation')  .or. &
             (trim(family) == 'Empty'   ) ) then

            self%family = family

        else
            user_msg = "bc_state_group%set_family: The string passed in to set the boundary &
                        condition family did not match any of valid boundary condition &
                        families. These include: 'Wall', 'Inlet', 'Outlet', 'Symmetry', &
                        'Periodic', 'Farfield', 'Scalar', 'Extrapolation'"
            call chidg_signal_one(FATAL,user_msg,family)
        end if

    end subroutine set_family
    !******************************************************************************************






    !>  Return the bc_state_group family.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_family(self) result(family)
        class(bc_state_group_t),  intent(inout)   :: self

        character(:),   allocatable :: family, user_msg

        if (allocated(self%family)) then
            family = self%family
        else
            family = ' '
        end if

    end function get_family
    !******************************************************************************************







    !>  Add a bc_state object to the bc_group.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc_state(self,bc_state)
        class(bc_state_group_t),  intent(inout)   :: self
        class(bc_state_t),        intent(in)      :: bc_state 

        character(:),   allocatable :: group_family, state_family, user_msg
        integer(ik)                 :: ierr, state_ID
        logical                     :: add_state
        

        !
        ! Check bc_state family conforms to any already added.
        !
        group_family = self%get_family()
        state_family = bc_state%get_family()
        add_state = (group_family == ' ') .or. (trim(group_family) == trim(state_family))


        !
        ! Add to vector of bc_states on the group.
        !
        if (add_state) then
            call self%set_family(trim(state_family))
            state_ID = self%new_bc_state()
            allocate(self%bc_state(state_ID)%state, source=bc_state, stat=ierr)
            if (ierr /= 0) call AllocationError
        else
            user_msg = "bc_state_group%add_bc_state: An attempt was made to add a bc_state &
                        object to a bc_state_group with dissimilar family. As a rule, &
                        bc_state_group objects may only contain bc_state objects of a &
                        single family."
            call chidg_signal_one(FATAL,user_msg, bc_state%get_name())
        end if

    end subroutine add_bc_state
    !******************************************************************************************









    !>  Extend the storage for bc_state functions on the group. Return an index location
    !!  for the new bc_state object in self%bc_states(:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function new_bc_state(self) result(state_ID)
        class(bc_state_group_t),    intent(inout)   :: self

        integer(ik)                             :: state_ID, ierr
        type(bc_state_wrapper_t),   allocatable :: temp_states(:)


        !
        ! Resize array storage
        !
        allocate(temp_states(self%nbc_states() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%nbc_states() > 0) then
            temp_states(1:size(self%bc_state)) = self%bc_state(1:size(self%bc_state))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_states,self%bc_state)
        


        !
        ! Set patch identifier of newly allocated patch that will be returned
        !
        state_ID = self%nbc_states()



    end function new_bc_state
    !**************************************************************************************








    !>  Remove bc_states that have been allocated.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/7/2017
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine remove_states(self)  
        class(bc_state_group_t),    intent(inout)   :: self


        if (allocated(self%bc_state)) deallocate(self%bc_state)


    end subroutine remove_states
    !**************************************************************************************








    !>  Return the number of bc_state objects in the group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !--------------------------------------------------------------------------------------
    function nbc_states(self) result(n)
        class(bc_state_group_t),    intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%bc_state)) then
            n = size(self%bc_state)
        else
            n = 0
        end if

    end function nbc_states
    !***************************************************************************************







    !>  Initialize parallel communicators on the bc_state group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_comm(self,mesh)
        class(bc_state_group_t),    intent(in)  :: self
        type(mesh_t),               intent(in)  :: mesh

        logical                     :: irank_has_geometry, ranks_have_geometry(NRANK)
        integer(ik)                 :: ierr, color
        character(:),   allocatable :: user_msg



!        !
!        ! Check if current processor contains geometry associated with the bc_t
!        !
!        irank_has_geometry = allocated(self%bc_patch)
!
!
!
!        !
!        ! Send this single information to all and receive back ordered information from all
!        !
!        call MPI_Allgather(irank_has_geometry,1,MPI_LOGICAL,ranks_have_geometry,1,MPI_LOGICAL,ChiDG_COMM,ierr)
!        user_msg = "bc%init_bc_comm: Error in collective MPI_Allgather for determining which &
!                    MPI ranks contain portions of a boundary condition bc_patch."
!        if (ierr /= 0) call chidg_signal(FATAL,user_msg)
!
!
!        !
!        ! Create a new MPI communicator for the current boundary condition 
!        ! that includes only those processors with bc_patch data that
!        ! has been allocated; indicating they contain a portion of the 
!        ! bc geometry.
!        !
!        if (irank_has_geometry) then
!            color = 1
!        else
!            color = MPI_UNDEFINED
!        end if
!
!
!        call MPI_Comm_split(ChiDG_COMM, color, IRANK, self%bc_COMM, ierr)
!        user_msg = "bc%init_bc_comm: Error in collective MPI_Comm_split when trying &
!                    to create a communicator for exchanging boundary condition data &
!                    between processors."
!        if (ierr /= 0) call chidg_signal(FATAL,user_msg)


    end subroutine init_comm
    !***************************************************************************************









    !>  Implementation specific routine for bc_state objects.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_specialized(self,mesh)
        class(bc_state_group_t),    intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh

        integer(ik)                 :: iop, group_ID


        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        if (allocated(self%bc_state)) then

            group_ID = mesh%get_bc_patch_group_id(self%name)
            if (group_ID /= NO_ID) then
                if (mesh%bc_patch_group(group_ID)%npatches() > 0) then

                    do iop = 1,size(self%bc_state)
                        call self%bc_state(iop)%state%init_bc_specialized(mesh, group_ID, self%bc_COMM)
                    end do !iop

                end if !bc_patch
            end if

        end if !bc_state




    end subroutine init_specialized
    !****************************************************************************************








    !>  Initialize boundary condition coupling.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_coupling(self,mesh)
        class(bc_state_group_t),    intent(inout)  :: self
        type(mesh_t),               intent(inout)  :: mesh


        integer(ik) :: iop, group_ID

        !
        ! Have bc_operators initialize the boundary condition coupling
        !
        if (allocated(self%bc_state)) then

            group_ID = mesh%get_bc_patch_group_id(self%name)
            if (group_ID /= NO_ID) then
                if (mesh%bc_patch_group(group_ID)%npatches() > 0) then

                    do iop = 1,size(self%bc_state)
                        call self%bc_state(iop)%state%init_bc_coupling(mesh,group_ID)
                    end do !iop

                end if !bc_patch
            end if !NO_ID

        end if !bc_state




    end subroutine init_coupling
    !****************************************************************************************







    !>  Propagate boundary condition coupling to mesh data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine propagate_coupling(self,mesh)
        class(bc_state_group_t),    intent(in)      :: self
        type(mesh_t),               intent(inout)   :: mesh

        integer(ik) :: group_ID, patch_ID, face_ID, idom, ielem, iface


        !
        ! set ncoupled elements back to mesh face
        !
        if (allocated(self%bc_state)) then

            group_ID = mesh%get_bc_patch_group_id(self%name)
            if (group_ID /= NO_ID) then
                if (mesh%bc_patch_group(group_ID)%npatches() > 0) then

                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                        do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                            idom  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l(face_ID)
                            ielem = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                            iface = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                            mesh%domain(idom)%faces(ielem,iface)%bc_ndepend = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

                        end do !face_ID
                    end do !patch_ID

                end if !bc_patch
            end if !NO_ID

        end if !bc_state


    end subroutine propagate_coupling
    !****************************************************************************************















end module type_bc_state_group
