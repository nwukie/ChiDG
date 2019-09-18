module type_rbf_mesh_motion
#include <messenger.h>
    use mod_constants
    use mod_kinds,                              only: rk, ik
    use type_mesh_motion,                       only: mesh_motion_t
    use type_rbf_info,                          only: rbf_info_t
    use type_rbf_source_vector,                 only: rbf_source_vector_t
    use type_rbf_interpolation,                 only: rbf_interpolation_t
    use type_rbf_mm_driver_wrapper,             only: rbf_mm_driver_wrapper_t
    use type_rbf_mm_driver,                     only: rbf_mm_driver_t
    use pmm_rbf_mm_driver,                      only: rbf_mm_driver_pmm
    use type_chidg_worker,                      only: chidg_worker_t
    use type_mesh,                            only: mesh_t
    use type_svector,                           only: svector_t
    use mod_string
    implicit none


    !>  Mesh motion class. Abstract interface to enable treatment of rbf mesh motion.
    !!
    !!  Data flow:
    !!  gridfile :-- readrbfmeshmotion_hdf --> pmm_group, pmm_domain_data ...
    !!  :-- init_pmm_group, init_pmm_domain --> pmm(:), mesh
    !!
    !!  @author Eric Wolf
    !!  @date   3/16/2017
    !!
    !---------------------------------------------------------------------------------------------------------------


    type, public, extends(mesh_motion_t)     :: rbf_mesh_motion_t
        
        !
        ! Initialization Data Structures: These are populated from the HDF5 gridfile,
        ! then used to initialize the other data structures
        !

        ! rbf_info contains the types of base and explicit RBFs, the base radius,
        ! and the base fraction (fraction of source nodes used as base nodes)
        class(rbf_info_t), allocatable                          :: rbf_info

        class(rbf_source_vector_t), allocatable                 :: rbf_sources


        ! rbfi is the actual RBF interpolation object
        class(rbf_interpolation_t),     allocatable             :: rbfi

        ! dnodes and vnodes are the values of grid displacement and velocity at source nodes
        ! Their ordering should conform to the ordering of nodes in 
        ! rbfi%node_patch%nodes(:,:), of size (nnodes,3)
        real(rk),                       allocatable             :: dnodes(:,:)
        real(rk),                       allocatable             :: vnodes(:,:)

        real(rk),                       allocatable             :: dcoeff(:,:)
        real(rk),                       allocatable             :: vcoeff(:,:)

        ! rbf_mm_driver_ID(nnodes) gives the driver index (from the drivers object below)
        ! used to update the source values for each node
        integer(ik),                    allocatable             :: rbf_mm_driver_ID(:)

        ! rbf_mm_driver(ndrivers) holds the RBF driver objects used to provide
        ! source node values for displacement and velocity
        ! For now, these are computed according to PMMFs, but in the future
        ! these drivers could be part of an interface to a structural solver

        ! Question: should this be rbf_mm_driver_vector_t instead?
        class(rbf_mm_driver_wrapper_t), allocatable             :: rbf_mm_driver(:)

    contains
        ! mesh_motion_t deferred procedures
        procedure   :: init_name
        procedure   :: init
        procedure   :: update
        procedure   :: apply
   

        ! Specialized procedures
        procedure   :: add_rbf_info
        procedure   :: add_rbf_mm_driver
        procedure   :: new_rbf_mm_driver
        procedure   :: remove_states
        procedure   :: nrbf_mm_drivers

              

    end type rbf_mesh_motion_t
    
contains

    !
    ! Deferred Procedures
    !
    !>
    !! Description: Initializes name for factory method.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine init_name(self)
        class(rbf_mesh_motion_t),   intent(inout)   :: self

        call self%set_family_name("RBF")
    end subroutine init_name
    !********************************************************************************
    !>
    !! Description: Initializes the RBFMM object.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self, mesh)
        class(rbf_mesh_motion_t),   intent(inout)   :: self
        type(mesh_t),            intent(inout)      :: mesh 

        type(svector_t)     :: patch_names
        type(string_t)      :: patch_name
        integer(ik)         :: ipatch, npatches, inode, nnodes, nnodes_base, ierr
        real(rk) :: radius_base(3)

        call self%init_name()
        ! Add the RBFs to the RBFI

        call self%rbfi%add_rbf(self%rbf_info%rbf_name)
        call self%rbfi%add_rbf_explicit(self%rbf_info%rbf_name_explicit)

        ! Form the source node patch and register
        npatches = self%rbf_sources%size()

        do ipatch = 1, npatches
            ! Get the patch name
            call patch_name%set(self%rbf_sources%data(ipatch)%patch_name)

            call patch_names%push_back(patch_name)

            ! Add the MM Driver
            call self%add_rbf_mm_driver(self%rbf_sources%data(ipatch)%driver)

        end do

        call self%rbfi%assemble_rbf_patch(mesh,patch_names)

        nnodes = self%rbfi%rbf_node_patch%nnodes_patch
        allocate(self%dnodes(nnodes,3), &
                 self%vnodes(nnodes,3), &
                 self%dcoeff(nnodes,3), &
                 self%vcoeff(nnodes,3), &
                self%rbf_mm_driver_ID(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Assign each node to a driver
        do inode = 1, nnodes
            self%rbf_mm_driver_ID(inode) = self%rbfi%rbf_node_patch%patch_ID(inode)
        end do


        ! Construct the RBFI
        nnodes_base = nint(nnodes*self%rbf_info%base_fraction, ik)
        radius_base = self%rbf_info%radius
        call self%rbfi%construct_rbf_interpolation(nnodes_base, radius_base)

    end subroutine init

    !>
    !! Description: Updates the RBF source nodes according to their respective drivers.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine update(self,mesh,time)
        class(rbf_mesh_motion_t),   intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        real(rk),                   intent(in)      :: time

        integer(ik) :: inode, idir, driver_ID
        real(rk) :: node(3), val(3)

        ! Loop over RBF source nodes and apply the driver
        do inode = 1, self%rbfi%rbf_node_patch%nnodes_patch

            driver_ID = self%rbf_mm_driver_ID(inode)
            node(1:3) = self%rbfi%rbf_node_patch%nodes(inode,1:3)

! internal compiler error when using fbounds=check
!            !self%dnodes(inode,1:3) = self%rbf_mm_driver(driver_ID)%driver%compute_disp(time, node)
!            !self%vnodes(inode,1:3) = self%rbf_mm_driver(driver_ID)%driver%compute_vel(time, node)
!            val = self%rbf_mm_driver(driver_ID)%driver%compute_disp(time, node)
!            self%dnodes(inode,:) = val!self%rbf_mm_driver(driver_ID)%driver%compute_disp(time, node)
!            val = self%rbf_mm_driver(driver_ID)%driver%compute_vel(time, node)
!            self%vnodes(inode,:) = val!self%rbf_mm_driver(driver_ID)%driver%compute_vel(time, node)


print*, 'WARNING: uncomment in type_rbf_mesh_motion.f90'
!            val = self%rbf_mm_driver(driver_ID)%driver%compute_disp(time,node)
!            self%dnodes(inode,:) = val
!
!            val = self%rbf_mm_driver(driver_ID)%driver%compute_vel(time,node)
!            self%vnodes(inode,:) = val
            

        end do


        ! Solve for coefficients
        do idir = 1,3
            self%dcoeff(:,idir) = self%rbfi%solve(self%dnodes(:,idir))
            self%vcoeff(:,idir) = self%rbfi%solve(self%vnodes(:,idir))
        end do

    end subroutine update
    !*******************************************************************************




    !>
    !! Description: Applies the RBFMM to the mesh nodes
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine apply(self,mesh,time)
        class(rbf_mesh_motion_t),   intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh 
        real(rk),                       intent(in)      :: time

        integer(ik) :: idir, inode, nnodes, idom, mm_ID

        do idom = 1,mesh%ndomains()
            mm_ID = mesh%domain(idom)%mm_ID
            if (mm_ID == self%mm_ID) then
                do idir = 1,3
                    do inode = 1, size(mesh%domain(idom)%nodes,1)
                        mesh%domain(idom)%dnodes(inode,idir) = self%rbfi%evaluate(self%dcoeff(:,idir),mesh%domain(idom)%nodes(inode,:))
                        mesh%domain(idom)%vnodes(inode,idir) = self%rbfi%evaluate(self%vcoeff(:,idir),mesh%domain(idom)%nodes(inode,:))
                    end do
                end do
            end if
        end do

    end subroutine apply

    !
    ! Specialized Procedures
    !

    !>
    !! Description: Adds an input rbf_info object to the RBFMM object. Used in initialization.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine add_rbf_info(self,rbf_info)
        class(rbf_mesh_motion_t),  intent(inout)   :: self
        class(rbf_info_t),        intent(in)      :: rbf_info

        integer(ik)                                     :: ierr
        if (allocated(self%rbf_info)) deallocate(self%rbf_info)

        allocate(self%rbf_info, source = rbf_info, stat = ierr)
        if (ierr /= 0) call AllocationError

    end subroutine add_rbf_info

    !>
    !! Description: Adds an RBFMM driver object to the RBFMM object.
    !!
    !! @author: Eric M. Wolf
    !! @date:   12/29/2017 
    !!
    !--------------------------------------------------------------------------------
    subroutine add_rbf_mm_driver(self,rbf_mm_driver)
        class(rbf_mesh_motion_t),  intent(inout)   :: self
        class(rbf_mm_driver_t),        intent(in)      :: rbf_mm_driver 

        character(:),   allocatable :: group_family, state_family, user_msg
        integer(ik)                 :: ierr, state_ID
        logical                     :: add_state
        

        !
        ! Check rbf_mm_driver family conforms to any already added.
        !
        !group_family = self%get_family()
        !state_family = rbf_mm_driver%get_family()
        !add_state = (trim(group_family) == ''                ) .or. &
        !            (trim(group_family) == 'Empty'           ) .or. &
        !            (trim(group_family) == trim(state_family))


        !
        ! Add to vector of rbf_mm_drivers on the group.
        !
       ! if (add_state) then
            !call self%set_family(trim(state_family))
        state_ID = self%new_rbf_mm_driver()
        allocate(self%rbf_mm_driver(state_ID)%driver, source=rbf_mm_driver, stat=ierr)
        if (ierr /= 0) call AllocationError
        !else
        !    user_msg = "rbf_mm_driver_group%add_rbf_mm_driver: An attempt was made to add a rbf_mm_driver &
        !                object to a rbf_mm_driver_group with dissimilar family. As a rule, &
        !                rbf_mm_driver_group objects may only contain rbf_mm_driver objects of a &
        !                single family."
        !    call chidg_signal(FATAL,user_msg)
        !end if

    end subroutine add_rbf_mm_driver
    !******************************************************************************************









    !>  Extend the storage for rbf_mm_driver functions on the group. Return an index location
    !!  for the new rbf_mm_driver object in self%rbf_mm_drivers(:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function new_rbf_mm_driver(self) result(state_ID)
        class(rbf_mesh_motion_t),    intent(inout)   :: self

        integer(ik)                             :: state_ID, ierr
        type(rbf_mm_driver_wrapper_t),   allocatable :: temp_states(:)


        !
        ! Resize array storage
        !
        allocate(temp_states(self%nrbf_mm_drivers() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%nrbf_mm_drivers() > 0) then
            temp_states(1:size(self%rbf_mm_driver)) = self%rbf_mm_driver(1:size(self%rbf_mm_driver))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_states,self%rbf_mm_driver)
        


        !
        ! Set patch identifier of newly allocated patch that will be returned
        !
        state_ID = self%nrbf_mm_drivers()



    end function new_rbf_mm_driver
    !**************************************************************************************








    !>  Remove rbf_mm_drivers that have been allocated.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/7/2017
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine remove_states(self)  
        class(rbf_mesh_motion_t),    intent(inout)   :: self

        if (allocated(self%rbf_mm_driver)) deallocate(self%rbf_mm_driver)

    end subroutine remove_states
    !**************************************************************************************








    !>  Return the number of rbf_mm_driver objects in the group.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !--------------------------------------------------------------------------------------
    function nrbf_mm_drivers(self) result(n)
        class(rbf_mesh_motion_t),    intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%rbf_mm_driver)) then
            n = size(self%rbf_mm_driver)
        else
            n = 0
        end if

    end function nrbf_mm_drivers
    !***************************************************************************************








end module type_rbf_mesh_motion
