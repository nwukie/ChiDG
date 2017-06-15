module type_mesh
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NO_ID
    use type_point,                 only: point_t
    use type_domain,                only: domain_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use type_bc_patch,              only: bc_patch_t
    use type_bc_patch_group,        only: bc_patch_group_t
    use type_ivector,               only: ivector_t
    implicit none
    private




    !> Data type for mesh information
    !!      - contains array of elements, array of faces for each element
    !!      - calls initialization procedure for elements and faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !----------------------------------------------------------------------------------
    type, public :: mesh_t

        integer(ik)                         :: ntime_

        type(domain_t),         allocatable :: domain(:)
        type(bc_patch_group_t), allocatable :: bc_patch_group(:)

    contains

        ! Domain procedures
        procedure           :: add_domain
        procedure           :: get_domain_id
        procedure, private  :: new_domain
        procedure           :: ndomains


        ! Boundary patch procedures
        procedure           :: add_bc_patch
        procedure           :: get_bc_patch_group_id
        procedure, private  :: new_bc_patch_group
        procedure           :: nbc_patch_groups


        ! Resouce management
        procedure           :: release


        ! Parallel communication pattern
        procedure           :: get_recv_procs
        procedure           :: get_send_procs


        ! Extra routines for testing private procedures
        procedure           :: stub_new_domain

    end type mesh_t
    !**********************************************************************************





contains




    !>  Add a domain to the mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine add_domain(self, name, nodes, connectivity, nelements_g, coord_system, eqn_ID)
        class(mesh_t),                  intent(inout)   :: self
        character(*),                   intent(in)      :: name
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: nelements_g
        character(*),                   intent(in)      :: coord_system
        integer(ik),                    intent(in)      :: eqn_ID

        integer(ik) :: idomain_l


        !
        ! Create new domain
        !
        idomain_l = self%new_domain()


        !
        ! Set name
        !
        self%domain(idomain_l)%name = trim(name)


        !
        ! Initialize new domain
        !
        call self%domain(idomain_l)%init_geom(idomain_l,    &
                                              nelements_g,  &
                                              nodes,        &
                                              connectivity, &
                                              coord_system )

        
        !
        ! Set equation set identifier
        !
        call self%domain(idomain_l)%init_eqn(eqn_ID)

    end subroutine add_domain
    !***********************************************************************************










    !>  Extend domain storage by one and return the domain identifier of the newly 
    !!  allocated domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!
    !-----------------------------------------------------------------------------------
    function new_domain(self) result(idomain_l)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)                 :: idomain_l, ierr
        type(domain_t), allocatable :: temp_domains(:)



        !
        ! Resize array storage
        !
        allocate(temp_domains(self%ndomains() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%ndomains() > 0) then
            temp_domains(1:size(self%domain)) = self%domain(1:size(self%domain))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_domains,self%domain)
        


        !
        ! Set domain identifier of newly allocated domain that will be returned
        !
        idomain_l = self%ndomains()


    end function new_domain
    !*********************************************************************************






    !>  Return the number of domains in the chidg_data_t instance.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function ndomains(self) result(ndomains_)
        class(mesh_t),  intent(in)      :: self

        integer :: ndomains_

        if (allocated(self%domain)) then
            ndomains_ = size(self%domain)
        else
            ndomains_ = 0
        end if

    end function ndomains
    !********************************************************************************








    !>  Given a domain name, return the domain identifier so that it can be 
    !!  indexed in self%domains(:).
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_domain_id(self,domain_name) result(domain_index)
        class(mesh_t),  intent(in)  :: self
        character(*),   intent(in)  :: domain_name

        character(:),   allocatable :: user_msg
        integer(ik)  :: idom
        integer(ik)  :: domain_index
        
        domain_index = 0

        do idom = 1,self%ndomains()
            if ( trim(domain_name) == trim(self%domain(idom)%name) ) then
                domain_index = idom
                exit
            end if
        end do


        user_msg = "mesh%get_domain_index: No domain was found with a name &
                    matching the incoming string"
        if (domain_index == 0) call chidg_signal_one(FATAL,user_msg,domain_name)


    end function get_domain_id
    !********************************************************************************






    !>  Add a bc_patch to the specified group.
    !!
    !!  Notes:
    !!      - the bc_patch_group will always be created here
    !!      - the bc_patch may or may not get created. It depends on if the current
    !!        processor contains geometry associated with the boundary patch. If
    !!        geometry is found on the local processor, the patch will be created.
    !!        If geometry is not found, the patch will not be created.
    !!
    !!
    !!  @param[in]  domain_name     Name of the domain associated with the patch
    !!  @param[in]  group_name      Name of the boundary condition group
    !!  @param[in]  bc_connectivity Boundary condition face connectivity to initialize
    !!  @param[in]  bc_ID           Identifier for the associated bc_state_group
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine add_bc_patch(self,domain_name, group_name, bc_connectivity, bc_ID)
        class(mesh_t),                  intent(inout)   :: self
        character(*),                   intent(in)      :: domain_name
        character(*),                   intent(in)      :: group_name
        type(boundary_connectivity_t),  intent(in)      :: bc_connectivity
        integer(ik),                    intent(in)      :: bc_ID

        integer(ik) :: group_ID, idomain

        !
        ! Find patch group identifier. If none, call new_bc_patch_group.
        !   - don't create patches for 'Empty' groups. These should be left ORPHAN
        !     so that they get picked up by Chimera.
        !
        if ( (trim(group_name) /= 'Empty') .and. &
             (trim(group_name) /= 'empty') ) then


            group_ID = self%get_bc_patch_group_id(group_name)

            if (group_ID == NO_ID) then
                group_ID = self%new_bc_patch_group()
                self%bc_patch_group(group_ID)%name     = trim(group_name)
                self%bc_patch_group(group_ID)%group_ID = group_ID
            end if


            !
            ! Find domain identifier
            !
            idomain = self%get_domain_id(domain_name)


            !
            ! Add bc_patch
            !
            call self%bc_patch_group(group_ID)%add_bc_patch(self%domain(idomain), bc_connectivity, bc_ID)


        end if



    end subroutine add_bc_patch
    !*********************************************************************************







    !>  Extend the allocation of bc_patch_group's by one and return index of new
    !!  bc_patch_group object in self%bc_patch_group(:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------
    function new_bc_patch_group(self) result(group_ID)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)                         :: group_ID, ierr
        type(bc_patch_group_t), allocatable :: temp_groups(:)



        !
        ! Resize array storage
        !
        allocate(temp_groups(self%nbc_patch_groups() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%nbc_patch_groups() > 0) then
            temp_groups(1:size(self%bc_patch_group)) = self%bc_patch_group(1:size(self%bc_patch_group))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp_groups,self%bc_patch_group)
        


        !
        ! Set domain identifier of newly allocated domain that will be returned
        !
        group_ID = self%nbc_patch_groups()


    end function new_bc_patch_group
    !*********************************************************************************







    !>  Return the number of bc_patch_group's on the mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    function nbc_patch_groups(self) result(n)
        class(mesh_t),  intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%bc_patch_group)) then
            n = size(self%bc_patch_group)
        else
            n = 0
        end if

    end function nbc_patch_groups
    !*********************************************************************************







    !>  Given a patch group name, return its index on self%bc_patch_group(:).
    !!
    !!  If no group is found, return NO_ID.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    function get_bc_patch_group_id(self,group_name_in) result(group_ID)
        class(mesh_t),  intent(in)  :: self
        character(*),   intent(in)  :: group_name_in

        integer(ik)                 :: igroup, group_ID
        character(:),   allocatable :: group_name
        logical                     :: found_group

        group_ID = NO_ID
        do igroup = 1,self%nbc_patch_groups()

            found_group = trim(group_name_in) == trim(self%bc_patch_group(igroup)%name)
            
            if (found_group) group_ID = igroup
            if (found_group) exit

        end do


    end function get_bc_patch_group_id
    !********************************************************************************







    !>  Release any allocated resource on the object.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/5/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine release(self)
        class(mesh_t),  intent(inout)   :: self


        if (allocated(self%domain)) deallocate(self%domain)


    end subroutine release
    !*********************************************************************************





    !>  Return the MPI ranks that the current rank is receiving information from.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    function get_recv_procs(self) result(recv_procs_array)
        class(mesh_t),  intent(in)  :: self

        integer(ik)                 :: idom, igroup, iproc
        integer(ik),    allocatable :: recv_procs_dom(:), recv_procs_bc(:), recv_procs_array(:)
        type(ivector_t)             :: recv_procs


        !
        ! Accumulate processor ranks that we are receiving from: domains
        !
        do idom = 1,self%ndomains()
            recv_procs_dom = self%domain(idom)%get_recv_procs()

            do iproc = 1,size(recv_procs_dom)
                call recv_procs%push_back_unique(recv_procs_dom(iproc))
            end do ! iproc

        end do ! idom



        !
        ! Accumulate processor ranks that we are receiving from: bc_patch_groups
        !
        do igroup = 1,self%nbc_patch_groups()
            recv_procs_bc = self%bc_patch_group(igroup)%get_recv_procs()

            do iproc = 1,size(recv_procs_bc)
                call recv_procs%push_back_unique(recv_procs_bc(iproc))
            end do ! iproc

        end do ! igroup



        !
        ! Return as integer array
        !
        recv_procs_array = recv_procs%data()


    end function get_recv_procs
    !*********************************************************************************







    !>  Return the MPI ranks that the current rank is sending information to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    function get_send_procs(self) result(send_procs_array)
        class(mesh_t),  intent(in)  :: self

        integer(ik)                 :: idom, iproc, igroup
        integer(ik),    allocatable :: send_procs_dom(:), send_procs_bc(:), send_procs_array(:)
        type(ivector_t)             :: send_procs


        !
        ! Accumulate processor ranks that we are receiving from: domains
        !
        do idom = 1,self%ndomains()
            send_procs_dom = self%domain(idom)%get_send_procs()

            do iproc = 1,size(send_procs_dom)
                call send_procs%push_back_unique(send_procs_dom(iproc))
            end do ! iproc

        end do ! idom



        !
        ! Accumulate processor ranks that we are receiving from: bc_patch_groups
        !
        do igroup = 1,self%nbc_patch_groups()
            send_procs_bc = self%bc_patch_group(igroup)%get_send_procs()

            do iproc = 1,size(send_procs_bc)
                call send_procs%push_back_unique(send_procs_bc(iproc))
            end do ! iproc

        end do ! igroup



        !
        ! Return as integer array
        !
        send_procs_array = send_procs%data()


    end function get_send_procs
    !*********************************************************************************










    !>  A public interface for calling 'new_domain'. 
    !!
    !!  Used for tests.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !--------------------------------------------------------------------------------
    function stub_new_domain(self) result(idomain)
        class(mesh_t),  intent(inout)   :: self

        integer(ik) :: idomain

        idomain = self%new_domain()

    end function stub_new_domain
    !********************************************************************************




end module type_mesh
