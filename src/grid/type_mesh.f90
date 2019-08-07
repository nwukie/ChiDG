module type_mesh
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NO_ID, INTERIOR, NFACES, CARTESIAN, CYLINDRICAL, NO_PROC
    use mod_grid,                   only: NFACE_CORNERS
    use mod_chidg_mpi,              only: NRANK, IRANK
    use type_element,               only: element_t
    use type_domain,                only: domain_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use type_bc_patch,              only: bc_patch_t
    use type_bc_patch_group,        only: bc_patch_group_t
    use type_ivector,               only: ivector_t
    use type_element_info,          only: element_info_t, element_info
    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08,                    only: mpi_isend, mpi_recv, mpi_integer4, mpi_real8, &
                                          mpi_waitall, mpi_request, mpi_status_ignore,  &
                                          mpi_statuses_ignore, mpi_character, mpi_sum,  &
                                          mpi_max, mpi_logical, mpi_lor, mpi_comm
    use type_octree,                only: octree_t
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

        integer(ik)                         :: ntime_         = NO_ID
        integer(ik)                         :: mesh_dof_start = NO_ID ! Based on Fortran 1-indexing

        ! Local data
        type(domain_t),         allocatable :: domain(:)
        type(bc_patch_group_t), allocatable :: bc_patch_group(:)

        ! Parallel data
        integer(ik),            allocatable :: nelements_per_proc(:)
        type(element_t),        allocatable :: parallel_element(:)
        type(mpi_request_vector_t)          :: comm_requests

        ! Tree data
        real(rk), allocatable               :: global_nodes(:,:)
        type(octree_t)                      :: octree

        ! Interpolation node set ('Uniform', 'Quadrature')
        character(:),           allocatable :: interpolation

    contains

        ! Mesh procedures
        procedure           :: nelements
        procedure           :: get_dof_start
        procedure           :: get_element_info

        ! Domain procedures
        procedure           :: add_domain
        procedure           :: get_domain_id
        procedure, private  :: new_domain
        procedure           :: ndomains

        ! Boundary patch group procedures
        procedure, private  :: new_bc_patch_group
        procedure           :: get_bc_patch_group_id
        procedure           :: nbc_patch_groups
        procedure           :: add_bc_patch

        ! Chimera
        procedure           :: assemble_chimera_data

        ! Octree
        procedure           :: set_global_nodes

        ! Resouce management
        procedure           :: release

        ! Parallel communication pattern
        procedure           :: get_recv_procs
        procedure           :: get_send_procs

        procedure           :: get_proc_ninterior_neighbors
        procedure           :: get_proc_nchimera_donors
        procedure           :: get_nelements_recv

        procedure           :: find_parallel_element
        procedure           :: new_parallel_element
        procedure           :: nparallel_elements
        procedure           :: get_parallel_dofs


        ! Extra routines for testing private procedures
        procedure           :: stub_new_domain

        ! Communication 
        procedure           :: init_comm_local
        procedure           :: init_comm_global

        procedure           :: comm_nelements
        procedure           :: comm_domain_procs
        procedure           :: comm_send
        procedure           :: comm_recv
        procedure           :: comm_wait

    end type mesh_t
    !**********************************************************************************





contains




    !>  Add a domain to the mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2017
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine add_domain(self, name, nodes, dnodes, vnodes, connectivity, nelements_g, coord_system, eqn_ID)
        class(mesh_t),                  intent(inout)   :: self
        character(*),                   intent(in)      :: name
        real(rk),                       intent(in)      :: nodes(:,:)
        real(rk),                       intent(in)      :: dnodes(:,:)
        real(rk),                       intent(inout)   :: vnodes(:,:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: nelements_g
        character(*),                   intent(in)      :: coord_system
        integer(ik),                    intent(in)      :: eqn_ID

        integer(ik) :: idomain_l

        integer(ik) :: inode, msg1, msg2, funit
        real(rk)    :: R1_velocity_1, R1_velocity_2, R1_velocity_3, &
                       R2_velocity_1, R2_velocity_2, R2_velocity_3
        integer(ik) :: R1_domain_min, R1_domain_max, &
                       R2_domain_min, R2_domain_max
        logical     :: file_exists

        namelist /region_one/ R1_velocity_1, R1_velocity_2, R1_velocity_3, R1_domain_min, R1_domain_max
        namelist /region_two/ R2_velocity_1, R2_velocity_2, R2_velocity_3, R2_domain_min, R2_domain_max


        !
        ! Create new domain, initialize
        !
        idomain_l = self%new_domain()
        self%domain(idomain_l)%name = trim(name)


        !
        ! Initialize reference domain
        !
        call self%domain(idomain_l)%init_geom(idomain_l,    &
                                              nelements_g,  &
                                              nodes,        &
                                              connectivity, &
                                              coord_system )
        



        !
        ! Check for grid velocity specification in grid_velocity.nml
        !
        inquire(file='grid_velocity.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=funit,form='formatted',file='grid_velocity.nml')
            read(funit,nml=region_one,iostat=msg1)
            read(funit,nml=region_two,iostat=msg2)
            close(funit)

            ! Region 1 
            if (msg1 == 0 .and. self%domain(idomain_l)%idomain_g >= R1_domain_min .and. self%domain(idomain_l)%idomain_g <= R1_domain_max) then
                print*, 'setting velocity: ', R1_velocity_1, R1_velocity_2, R1_velocity_3, ' on domain: ', self%domain(idomain_l)%idomain_g
                do inode = 1,size(self%domain(idomain_l)%vnodes,1)
                    select case(trim(coord_system))
                        case('Cartesian')
                            vnodes(inode,1:3) = [R1_velocity_1,R1_velocity_2,R1_velocity_3]
                        case('Cylindrical')
                            vnodes(inode,1:3) = [R1_velocity_1,nodes(inode,1)*R1_velocity_2,R1_velocity_3]
                        case default
                            call chidg_signal_one(FATAL,"mesh%add_domain: Invalid coordinate system.",trim(coord_system))
                    end select
                end do
            end if

            ! Region 2
            if (msg2 == 0 .and. self%domain(idomain_l)%idomain_g >= R2_domain_min .and. self%domain(idomain_l)%idomain_g <= R2_domain_max) then
                print*, 'setting velocity: ', R2_velocity_1, R2_velocity_2, R2_velocity_3, ' on domain: ', self%domain(idomain_l)%idomain_g
                do inode = 1,size(self%domain(idomain_l)%vnodes,1)
                    select case(trim(coord_system))
                        case('Cartesian')
                            vnodes(inode,1:3) = [R2_velocity_1,R2_velocity_2,R2_velocity_3]
                        case('Cylindrical')
                            vnodes(inode,1:3) = [R2_velocity_1,nodes(inode,1)*R2_velocity_2,R2_velocity_3]
                        case default
                            call chidg_signal_one(FATAL,"domain%init_geom: Invalid coordinate system.",trim(coord_system))
                    end select
                end do
            end if

        end if






        !
        ! Initialize ALE 
        !
        call self%domain(idomain_l)%set_displacements_velocities(dnodes,vnodes)

        
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
    subroutine add_bc_patch(self,domain_name, group_name, patch_name, bc_connectivity, bc_ID)
        class(mesh_t),                  intent(inout)   :: self
        character(*),                   intent(in)      :: domain_name
        character(*),                   intent(in)      :: group_name
        character(*),                   intent(in)      :: patch_name
        type(boundary_connectivity_t),  intent(in)      :: bc_connectivity
        integer(ik),                    intent(in)      :: bc_ID

        integer(ik) :: group_ID, idomain

        !
        ! Add new group, if not already in existence
        !
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
        call self%bc_patch_group(group_ID)%add_bc_patch(self%domain(idomain), patch_name, bc_connectivity, bc_ID)


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






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/26/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine assemble_chimera_data(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)             :: idom, ChiID, idonor, donor_proc, npts, ipt,      &
                                   idomain_g, ielement_g, idomain_l, ielement_l,    &
                                   pelem_ID, igq
        real(rk),   allocatable :: ref_coords(:,:)
        real(rk)                :: ale_g_pt, ale_g_grad_pt(3), &
                                   ale_Dinv_pt(3,3), interp_coords_vel_pt(3), xi, eta, zeta
        logical                 :: parallel_donor

 
        !
        ! Construct receiver interpolations
        !
        !
        do idom = 1,self%ndomains()
           do ChiID = 1,self%domain(idom)%chimera%nreceivers()
              do idonor = 1,self%domain(idom)%chimera%recv(ChiID)%ndonors()
               
        
                  ! Get donor location
                  donor_proc = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%iproc
                  ref_coords = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%coords
                  npts = size(ref_coords,1)
        
                  parallel_donor = (donor_proc /= IRANK)
        
                  if (parallel_donor) then
                      idomain_g   = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_g
                      ielement_g  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_g
                      pelem_ID    = self%find_parallel_element(idomain_g,ielement_g)
                      if (pelem_ID == NO_ID) call chidg_signal_three(FATAL,"mesh%assemble_chimera_data: did not find parallel chimera element.", IRANK, idomain_g, ielement_g)

                      do ipt = 1,npts
                          igq  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%node_index(ipt)
                          xi   = ref_coords(ipt,1)
                          eta  = ref_coords(ipt,2)
                          zeta = ref_coords(ipt,3)
                          call self%parallel_element(pelem_ID)%ale_point(xi, eta, zeta,ale_g_pt,ale_g_grad_pt,ale_Dinv_pt,interp_coords_vel_pt)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g(igq)               = ale_g_pt
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad1(igq)         = ale_g_grad_pt(1)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad2(igq)         = ale_g_grad_pt(2)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad3(igq)         = ale_g_grad_pt(3)
                          self%domain(idom)%chimera%recv(ChiID)%ale_Dinv(:,:,igq)        = ale_Dinv_pt
                          self%domain(idom)%chimera%recv(ChiID)%interp_coords_vel(igq,:) = interp_coords_vel_pt
                      end do
        
                   else
        
                      idomain_l   = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_l
                      ielement_l  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_l
        
                      do ipt = 1, npts
                          igq  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%node_index(ipt)
                          xi   = ref_coords(ipt,1)
                          eta  = ref_coords(ipt,2)
                          zeta = ref_coords(ipt,3)
                          call self%domain(idomain_l)%elems(ielement_l)%ale_point(xi, eta, zeta,ale_g_pt,ale_g_grad_pt,ale_Dinv_pt,interp_coords_vel_pt)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g(igq)               = ale_g_pt
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad1(igq)         = ale_g_grad_pt(1)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad2(igq)         = ale_g_grad_pt(2)
                          self%domain(idom)%chimera%recv(ChiID)%ale_g_grad3(igq)         = ale_g_grad_pt(3)
                          self%domain(idom)%chimera%recv(ChiID)%ale_Dinv(:,:,igq)        = ale_Dinv_pt
                          self%domain(idom)%chimera%recv(ChiID)%interp_coords_vel(igq,:) = interp_coords_vel_pt
                      end do
                   end if

         
              end do !idonor
           end do !ChiID
        end do !idom




    end subroutine assemble_chimera_data
    !********************************************************************************



    !>
    !!
    !! @author  Eric M. Wolf
    !! @date    08/30/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine set_global_nodes(self,global_nodes)
        class(mesh_t),    intent(inout)           :: self
        real(rk),           intent(in)              :: global_nodes(:,:)

        if (allocated(self%global_nodes)) deallocate(self%global_nodes)
        self%global_nodes = global_nodes

    end subroutine set_global_nodes
    !*****************************************************************************************








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

        if (allocated(self%domain))             deallocate(self%domain)
        if (allocated(self%bc_patch_group))     deallocate(self%bc_patch_group)
        if (allocated(self%nelements_per_proc)) deallocate(self%nelements_per_proc)
        if (allocated(self%parallel_element))   deallocate(self%parallel_element)
        if (allocated(self%global_nodes))       deallocate(self%global_nodes)

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









    !>  Return the number of interior neighbor elements that live on processor, iproc.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/19/2017
    !!
    !---------------------------------------------------------------------------------
    function get_proc_ninterior_neighbors(self,iproc) result(n)
        class(mesh_t),  intent(in)  :: self
        integer(ik),    intent(in)  :: iproc

        integer(ik) :: idom, ielem, iface, n


        !
        ! Accumulate number of INTERIOR neighbors that live on iproc
        !
        n = 0
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%get_nelements_local()
                do iface = 1,NFACES
                    if (self%domain(idom)%faces(ielem,iface)%ineighbor_proc == iproc) n = n + 1
                end do !iface
            end do !ielem
        end do !idom
        


    end function get_proc_ninterior_neighbors
    !**********************************************************************************





    !>  Return the number of chimera donors that live on processor, iproc.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/27/2017
    !!
    !---------------------------------------------------------------------------------
    function get_proc_nchimera_donors(self,iproc) result(n)
        class(mesh_t),  intent(in)  :: self
        integer(ik),    intent(in)  :: iproc

        integer(ik)     :: idom, ChiID, idonor, n, donor_domain_g, donor_element_g, i
        type(ivector_t) :: donor_domains_g, donor_elements_g
        logical         :: already_counted, proc_has_donor


        !
        ! Accumulate number of CHIMERA donors that live on iproc
        !
        do idom = 1,self%ndomains()
            do ChiID = 1,self%domain(idom)%chimera%nreceivers()
                do idonor = 1,self%domain(idom)%chimera%recv(ChiID)%ndonors()

                    proc_has_donor = (self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%iproc == iproc)
                    if (proc_has_donor) then

                        donor_domain_g  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_g
                        donor_element_g = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_g

                        ! Check if donor was already counted
                        already_counted = .false.
                        do i = 1,donor_elements_g%size()
                            already_counted = (donor_domains_g%at(i) == donor_domain_g .and. donor_elements_g%at(i) == donor_element_g)
                            if (already_counted) exit
                        end do

                        
                        if (.not. already_counted) then
                            call donor_domains_g%push_back(donor_domain_g)
                            call donor_elements_g%push_back(donor_element_g)
                        end if

                    end if
                    

                end do !idonor
            end do !ChiID
        end do !idom
        

        !
        ! Get number of counted donors coming from iproc
        !
        n = donor_elements_g%size()

    end function get_proc_nchimera_donors
    !**********************************************************************************









    !>  Return the number of elements being received into self%parallel_elements
    !!  by the mesh.
    !!
    !!  NOTE: currently, we are only receiving off-processor INTERIOR and CHIMERA elements
    !!        into self%parallel_elements. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/26/2017
    !!
    !-----------------------------------------------------------------------------------
    function get_nelements_recv(self) result(n)
        class(mesh_t),  intent(in)  :: self

        integer(ik)     :: idom, ielem, iface, irecv, idonor, donor_iproc, idomain_g, ielement_g, &
                           idomain_g_list, ielement_g_list, ientry, n
        logical         :: parallel_donor, already_counted
        type(ivector_t) :: recv_domain_g, recv_element_g


        ! Accumulate parallel INTERIOR neighbors
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%nelements()
                do iface = 1,NFACES

                    if (self%domain(idom)%faces(ielem,iface)%ftype == INTERIOR .and. self%domain(idom)%faces(ielem,iface)%ineighbor_proc /= IRANK) then
                        idomain_g  = self%domain(idom)%faces(ielem,iface)%ineighbor_domain_g
                        ielement_g = self%domain(idom)%faces(ielem,iface)%ineighbor_element_g

                        already_counted = .false.
                        do ientry = 1,recv_element_g%size()
                            idomain_g_list  = recv_domain_g%at(ientry)
                            ielement_g_list = recv_element_g%at(ientry)
                            if ( (idomain_g_list  == idomain_g ) .and. &
                                 (ielement_g_list == ielement_g) ) already_counted = .true.
                            if (already_counted) exit
                        end do

                        if (.not. already_counted) then
                            call recv_domain_g%push_back(idomain_g)
                            call recv_element_g%push_back(ielement_g)
                        end if

                    end if

                end do
            end do
        end do


        ! Accumulate parallel CHIMERA neighbors
        do idom = 1,self%ndomains()
            do irecv = 1,self%domain(idom)%chimera%nreceivers()
                do idonor = 1,self%domain(idom)%chimera%recv(irecv)%ndonors()

                    donor_iproc = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%elem_info%iproc

                    parallel_donor = (donor_iproc /= IRANK)
                    if (parallel_donor) then
                        idomain_g  = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%elem_info%idomain_g
                        ielement_g = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%elem_info%ielement_g

                        already_counted = .false.
                        do ientry = 1,recv_element_g%size()
                            idomain_g_list  = recv_domain_g%at(ientry)
                            ielement_g_list = recv_element_g%at(ientry)
                            if ( (idomain_g_list  == idomain_g ) .and. &
                                 (ielement_g_list == ielement_g) ) already_counted = .true.
                            if (already_counted) exit
                        end do

                        if (.not. already_counted) then
                            call recv_domain_g%push_back(idomain_g)
                            call recv_element_g%push_back(ielement_g)
                        end if
                    end if


                end do !idonor
            end do !irecv
        end do !idom




        !
        ! Get size of accumulated elements being received
        !
        n = recv_element_g%size()

    end function get_nelements_recv
    !***********************************************************************************










    !>  Try to locate an element in self%parallel_element(:). If found, return valid
    !!  pelem_ID for self%parallel_element(pelem_ID). If not found, return NO_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !-----------------------------------------------------------------------------------
    function find_parallel_element(self,idomain_g, ielement_g) result(pelem_ID)
        class(mesh_t),  intent(in)  :: self
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: ielement_g

        integer(ik) :: ielem, pelem_ID

        pelem_ID = NO_ID
        do ielem = 1,self%nparallel_elements()

            if ( (self%parallel_element(ielem)%idomain_g == idomain_g) .and.    &
                 (self%parallel_element(ielem)%ielement_g == ielement_g) ) then
                pelem_ID = ielem
                exit
            end if

        end do !ielem

    end function find_parallel_element
    !***********************************************************************************





    !>  Extend allocation for parallel_element and return new identifier, pelem_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !----------------------------------------------------------------------------------
    function new_parallel_element(self) result(pelem_ID)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)                     :: pelem_ID, ierr
        type(element_t),    allocatable :: temp(:)

        
        ! Resize array storage
        allocate(temp(self%nparallel_elements() + 1), stat=ierr)


        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%nparallel_elements() > 0) then
            temp(1:size(self%parallel_element)) = self%parallel_element(1:size(self%parallel_element))
        end if


        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        call move_alloc(temp,self%parallel_element)
        

        ! Set domain identifier of newly allocated domain that will be returned
        pelem_ID = self%nparallel_elements()


    end function new_parallel_element
    !**********************************************************************************





    !>  Return number of parallel elements on the mesh.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !---------------------------------------------------------------------------------
    function nparallel_elements(self) result(n)
        class(mesh_t),  intent(in)  :: self

        integer(ik) :: n

        if (allocated(self%parallel_element)) then
            n = size(self%parallel_element)
        else
            n = 0
        end if

    end function nparallel_elements
    !**********************************************************************************







    !>  Return the global indices of all parallel degrees-of-freedom that need
    !!  accessed on the local process.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/11/2019
    !!
    !!  Returned indices based on Fortran 1-indexing
    !!
    !!
    !----------------------------------------------------------------------------------
    function get_parallel_dofs(self) result(parallel_dofs)
        class(mesh_t),  intent(in)  :: self

        integer(ik) :: dof_start, nterms_s, ntime, nfields, pelem_ID, idof

        type(ivector_t)             :: parallel_dofs_vector
        integer(ik),    allocatable :: parallel_dofs(:)


        ! Loop through parallel_elements and accumulate dof indices
        do pelem_ID = 1,self%nparallel_elements()

            dof_start = self%parallel_element(pelem_ID)%dof_start
            nfields   = self%parallel_element(pelem_ID)%nfields
            nterms_s  = self%parallel_element(pelem_ID)%nterms_s
            ntime     = self%parallel_element(pelem_ID)%ntime

            ! Store each dof index
            do idof = dof_start, dof_start + (nfields*nterms_s*ntime - 1)
                call parallel_dofs_vector%push_back(idof)
            end do !idof

        end do !pelem_ID

        ! Return integer array
        parallel_dofs = parallel_dofs_vector%data()

    end function get_parallel_dofs
    !**********************************************************************************











    !>  Communicate number of elements on each process.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/24/2019
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_nelements(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)                 :: ierr
        integer(ik), allocatable    :: buffer(:)

        ! Allocate storage for number of elements on each process
        if (allocated(self%nelements_per_proc)) deallocate(self%nelements_per_proc)
        allocate(self%nelements_per_proc(NRANK), stat=ierr)
        if (ierr /= 0) call AllocationError
        self%nelements_per_proc = 0
        buffer = self%nelements_per_proc


        ! Fill current RANK entry. Index from 1, so add 1 since IRANK is referenced from 0
        self%nelements_per_proc(IRANK+1) = self%nelements()


        ! Reduce results
        call MPI_AllReduce(self%nelements_per_proc,buffer,NRANK,MPI_INTEGER4,MPI_SUM,ChiDG_COMM,ierr)
        self%nelements_per_proc = buffer


    end subroutine comm_nelements
    !*********************************************************************************




    !>  Communicate domain partition neighbors.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/22/2019
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_domain_procs(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)     :: ierr, idomain_g_local_max, idomain_g_global_max, iproc, idom, idomain_g, idomain_l_index
        type(ivector_t) :: procs
        logical         :: ranks_with_domain(NRANK), ranks_with_domain_reduced(NRANK)

        ! Find largest domain index, proc-local
        idomain_g_local_max = 0
        do idom = 1,self%ndomains()
            if (idomain_g_local_max < self%domain(idom)%idomain_g) idomain_g_local_max = self%domain(idom)%idomain_g
        end do


        ! Get global idomain_g max 
        call MPI_AllReduce(idomain_g_local_max, idomain_g_global_max, 1, MPI_INTEGER4, MPI_MAX, ChiDG_COMM, ierr)


        ! For each domain, check if local proc has and register if so. Then reduce registration
        ! for domain and collect participating processes into list to store for domain.
        do idomain_g = 1,idomain_g_global_max

            ranks_with_domain = .false.
            do idom = 1,self%ndomains()
                if (self%domain(idom)%idomain_g == idomain_g) then
                    ! Indexing from 1, so add 1 to IRANK since it is 0-based.
                    ranks_with_domain(IRANK+1) = .true.
                    idomain_l_index = idom
                    exit
                end if
            end do


            ! Reduce results
            ranks_with_domain_reduced = .false.
            call MPI_AllReduce(ranks_with_domain,ranks_with_domain_reduced,NRANK,MPI_LOGICAL,MPI_LOR,ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'mesh%comm_neighbors: error calling MPI_AllReduce.')


            ! If current proc has domain, register all participating procs 
            if (ranks_with_domain(IRANK+1)) then

                ! Collect participating procs
                call procs%clear()
                do iproc = 1,NRANK
                    if (ranks_with_domain_reduced(iproc)) then
                        call procs%push_back(iproc-1)
                    end if
                end do

                ! Register
                self%domain(idomain_l_index)%procs = procs%data()
            end if

        end do


    end subroutine comm_domain_procs
    !*********************************************************************************














    !>  Initialiate parallel communication of mesh quantities.
    !!
    !!  Currently communicates just ALE quantities.
    !!
    !!  Order of operations:
    !!      1: mesh%comm_send()
    !!      2: mesh%comm_recv()
    !!      3: mesh%comm_wait()
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_send(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)         :: idom, ielem, iface, isend, isend_proc, iproc, ierr,  &
                               send_size_a, send_size_b, send_size_c, send_size_d,  &
                               idomain_l, ielement_l
        type(mpi_request)   :: request(7)
        logical             :: interior_face, parallel_neighbor


        !
        ! Send interior face neighbor DATA: important because we need face data, not just element data.
        !
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%nelements()
                do iface = 1,NFACES

                    interior_face     = (self%domain(idom)%faces(ielem,iface)%ftype == INTERIOR)
                    parallel_neighbor = (self%domain(idom)%faces(ielem,iface)%ineighbor_proc /= IRANK)

                    !
                    ! For INTERIOR faces that have off-processor neighbors we need to communicate ALE data.
                    !
                    if ( interior_face .and. parallel_neighbor ) then

                        associate ( face => self%domain(idom)%faces(ielem,iface) ) 

                        send_size_a = size(face%neighbor_location)
                        send_size_b = size(face%interp_coords_vel)
                        send_size_c = size(face%ale_g)
                        send_size_d = size(face%ale_Dinv)

                        ! First, send neighbor location. This way, the receiving processor knows where to put the data.
                        ! Next, send all ALE information
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%neighbor_location, send_size_a, mpi_integer4, face%ineighbor_proc, 0, ChiDG_COMM, request(1), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%interp_coords_vel, send_size_b, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(2), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g,             send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(3), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad1,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(4), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad2,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(5), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad3,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(6), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_Dinv,          send_size_d, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(7), ierr)

                        call self%comm_requests%push_back(request(1))
                        call self%comm_requests%push_back(request(2))
                        call self%comm_requests%push_back(request(3))
                        call self%comm_requests%push_back(request(4))
                        call self%comm_requests%push_back(request(5))
                        call self%comm_requests%push_back(request(6))
                        call self%comm_requests%push_back(request(7))

                        end associate

                    end if


                end do !iface
            end do !ielem
        end do !idom





        !
        ! Send interior face neighbor elements: these will be initialized in target proc mesh%parallel_elements.
        ! That way, we can query parallel_elements to decide order of access in parallel vector storage, for example.
        !
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%nelements()
                do iface = 1,NFACES

                    interior_face     = (self%domain(idom)%faces(ielem,iface)%ftype == INTERIOR)
                    parallel_neighbor = (self%domain(idom)%faces(ielem,iface)%ineighbor_proc /= IRANK)

                    !
                    ! For INTERIOR faces that have off-processor neighbors we 
                    ! need to send that element to the neighbor processor.
                    !
                    if ( interior_face .and. parallel_neighbor ) then

                        associate ( face => self%domain(idom)%faces(ielem,iface) ) 

                        idomain_l  = self%domain(idom)%elems(ielem)%idomain_l
                        ielement_l = self%domain(idom)%elems(ielem)%ielement_l

                        send_size_a = size(self%domain(idomain_l)%elems(ielement_l)%connectivity)
                        send_size_b = size(self%domain(idomain_l)%elems(ielement_l)%node_coords)
                        send_size_c = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_def)
                        send_size_d = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel)
                
                        ! Send element to off-processor interior neighbor
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_location,             5, mpi_integer4, face%ineighbor_proc, 0, ChiDG_COMM, request(1), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 9, mpi_integer4, face%ineighbor_proc, 0, ChiDG_COMM, request(2), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords,        send_size_b, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(3), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_def,    send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(4), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel,    send_size_d, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(5), ierr)

                        call self%comm_requests%push_back(request(1))
                        call self%comm_requests%push_back(request(2))
                        call self%comm_requests%push_back(request(3))
                        call self%comm_requests%push_back(request(4))
                        call self%comm_requests%push_back(request(5))

                        end associate

                    end if


                end do !iface
            end do !ielem
        end do !idom



        !
        ! Send chimera donor elements: these will be initialized in target proc mesh%parallel_elements
        !
        do idom = 1,self%ndomains()
            ! For each element that needs sent
            do isend = 1,self%domain(idom)%chimera%nsend()
                ! For each processor the current element gets sent to
                do isend_proc = 1,self%domain(idom)%chimera%send(isend)%nsend_procs()

                    ! Get the processor rank we are sending to
                    iproc = self%domain(idom)%chimera%send(isend)%send_procs%at(isend_proc)

                    ! If receiver is off-processor, send reference and physical nodes/velocities
                    if (iproc /= IRANK) then

                        idomain_l  = self%domain(idom)%chimera%send(isend)%idomain_l
                        ielement_l = self%domain(idom)%chimera%send(isend)%ielement_l

                        send_size_a = size(self%domain(idomain_l)%elems(ielement_l)%connectivity)
                        send_size_b = size(self%domain(idomain_l)%elems(ielement_l)%node_coords)
                        send_size_c = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_def)
                        send_size_d = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel)
                
                        ! First send location of donor
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_location,             5, mpi_integer4, iproc, 0, ChiDG_COMM, request(1), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 9, mpi_integer4, iproc, 0, ChiDG_COMM, request(2), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords,        send_size_b, mpi_real8,    iproc, 0, ChiDG_COMM, request(3), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_def,    send_size_c, mpi_real8,    iproc, 0, ChiDG_COMM, request(4), ierr)
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel,    send_size_d, mpi_real8,    iproc, 0, ChiDG_COMM, request(5), ierr)


                        call self%comm_requests%push_back(request(1))
                        call self%comm_requests%push_back(request(2))
                        call self%comm_requests%push_back(request(3))
                        call self%comm_requests%push_back(request(4))
                        call self%comm_requests%push_back(request(5))

                    end if 

                end do !isend_procs
            end do !isend
        end do !idom















! PREVIOUS


!        !
!        ! Send interior face data
!        !
!        do idom = 1,self%ndomains()
!            do ielem = 1,self%domain(idom)%nelem
!                do iface = 1,NFACES
!
!                    interior_face     = (self%domain(idom)%faces(ielem,iface)%ftype == INTERIOR)
!                    parallel_neighbor = (self%domain(idom)%faces(ielem,iface)%ineighbor_proc /= IRANK)
!
!                    !
!                    ! For INTERIOR faces that have off-processor neighbors we need to communicate ALE data.
!                    !
!                    if ( interior_face .and. parallel_neighbor ) then
!
!                        associate ( face => self%domain(idom)%faces(ielem,iface) ) 
!
!                        send_size_a = size(face%neighbor_location)
!                        send_size_b = size(face%interp_coords_vel)
!                        send_size_c = size(face%ale_g)
!                        send_size_d = size(face%ale_Dinv)
!
!                        ! First, send neighbor location. This way, the receiving processor knows where to put the data.
!                        ! Next, send all ALE information
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%neighbor_location, send_size_a, mpi_integer4, face%ineighbor_proc, 0, ChiDG_COMM, request(1), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%interp_coords_vel, send_size_b, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(2), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g,             send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(3), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad1,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(4), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad2,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(5), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_g_grad3,       send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(6), ierr)
!                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%ale_Dinv,          send_size_d, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(7), ierr)
!
!                        call self%comm_requests%push_back(request(1))
!                        call self%comm_requests%push_back(request(2))
!                        call self%comm_requests%push_back(request(3))
!                        call self%comm_requests%push_back(request(4))
!                        call self%comm_requests%push_back(request(5))
!                        call self%comm_requests%push_back(request(6))
!                        call self%comm_requests%push_back(request(7))
!
!                        end associate
!
!                    end if
!
!
!                end do !iface
!            end do !ielem
!        end do !idom
!
!
!
!        !
!        ! Send chimera donors
!        !
!        do idom = 1,self%ndomains()
!            ! For each element that needs sent
!            do isend = 1,self%domain(idom)%chimera%nsend()
!                ! For each processor the current element gets sent to
!                do isend_proc = 1,self%domain(idom)%chimera%send(isend)%nsend_procs()
!
!                    ! Get the processor rank we are sending to
!                    iproc = self%domain(idom)%chimera%send(isend)%send_procs%at(isend_proc)
!
!                    ! If receiver is off-processor, send reference and physical nodes/velocities
!                    if (iproc /= IRANK) then
!
!                        idomain_l  = self%domain(idom)%chimera%send(isend)%idomain_l
!                        ielement_l = self%domain(idom)%chimera%send(isend)%ielement_l
!
!                        send_size_a = size(self%domain(idomain_l)%elems(ielement_l)%connectivity)
!                        send_size_b = size(self%domain(idomain_l)%elems(ielement_l)%node_coords)
!                        send_size_c = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_def)
!                        send_size_d = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel)
!                
!                        ! First send location of donor
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_location,             5, mpi_integer4, iproc, 0, ChiDG_COMM, request(1), ierr)
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 9, mpi_integer4, iproc, 0, ChiDG_COMM, request(2), ierr)
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords,        send_size_b, mpi_real8,    iproc, 0, ChiDG_COMM, request(3), ierr)
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_def,    send_size_c, mpi_real8,    iproc, 0, ChiDG_COMM, request(4), ierr)
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel,    send_size_d, mpi_real8,    iproc, 0, ChiDG_COMM, request(5), ierr)
!
!
!                        call self%comm_requests%push_back(request(1))
!                        call self%comm_requests%push_back(request(2))
!                        call self%comm_requests%push_back(request(3))
!                        call self%comm_requests%push_back(request(4))
!                        call self%comm_requests%push_back(request(5))
!
!                    end if 
!
!                end do !isend_procs
!            end do !isend
!        end do !idom
!
!
!!        do group_ID = 1,self%nbc_patch_groups()
!!            do patch_ID = 1,self%bc_patch_group(group_ID)%npatches()
!!                if (self%bc_patch_group(group_ID)%patch(patch_ID)%spatial_coupling = 'Global') then
!!
!!                    ! Send each element in the patch to each other processor in bc_comm
!!                    call MPI_Comm_Size(bc_comm, bc_NRANK, ierr)
!!                    call MPI_Comm_Rank(bc_comm, bc_IRANK, ierr)
!!
!!                    do iproc = 0,
!!
!!                end if !global spatial coupling
!!            end do !patch_ID
!!        end do !group_ID
!
!
!!        !
!!        ! Send coupled bc elements
!!        !
!!        do idom = 1,self%ndomains()
!!            ! For each element that needs sent
!!            do isend = 1,self%domain(idom)%chimera%nsend()
!!                ! For each processor the current element gets sent to
!!                do isend_proc = 1,self%domain(idom)%chimera%send(isend)%nsend_procs()
!!
!!                    ! Get the processor rank we are sending to
!!                    iproc = self%domain(idom)%chimera%send(isend)%send_procs%at(isend_proc)
!!
!!                    ! If receiver is off-processor, send reference and physical nodes/velocities
!!                    if (iproc /= IRANK) then
!!
!!                        idomain_l  = self%domain(idom)%chimera%send(isend)%idomain_l
!!                        ielement_l = self%domain(idom)%chimera%send(isend)%ielement_l
!!
!!                        send_size_a = size(self%domain(idomain_l)%elems(ielement_l)%connectivity)
!!                        send_size_b = size(self%domain(idomain_l)%elems(ielement_l)%node_coords)
!!                        send_size_c = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_def)
!!                        send_size_d = size(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel)
!!                
!!                        ! First send location of donor
!!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_location,             5, mpi_integer4, iproc, 0, ChiDG_COMM, request(1), ierr)
!!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 9, mpi_integer4, iproc, 0, ChiDG_COMM, request(2), ierr)
!!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords,        send_size_b, mpi_real8,    iproc, 0, ChiDG_COMM, request(3), ierr)
!!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_def,    send_size_c, mpi_real8,    iproc, 0, ChiDG_COMM, request(4), ierr)
!!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%node_coords_vel,    send_size_d, mpi_real8,    iproc, 0, ChiDG_COMM, request(5), ierr)
!!
!!
!!                        call self%comm_requests%push_back(request(1))
!!                        call self%comm_requests%push_back(request(2))
!!                        call self%comm_requests%push_back(request(3))
!!                        call self%comm_requests%push_back(request(4))
!!                        call self%comm_requests%push_back(request(5))
!!
!!                    end if 
!!
!!                end do !isend_procs
!!            end do !isend
!!        end do !idom





    end subroutine comm_send
    !*********************************************************************************





    !>  Initiates parallel recieve of mesh quantities.
    !!
    !!  Order of operations:
    !!      1: mesh%comm_send()
    !!      2: mesh%comm_recv()
    !!      3: mesh%comm_wait()
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(mesh_t),  intent(inout)   :: self

        character(:),   allocatable :: coord_system
        real(rk),       allocatable :: nodes(:,:), nodes_def(:,:), nodes_vel(:,:), nodes_disp(:,:)
        integer(ik),    allocatable :: recv_procs(:), connectivity(:)
        integer(ik)                 :: idom, ielem, iface, iproc, irecv, ierr,      &
                                       etype, nnodes, nterms_s, nterms_c, nfields,  &
                                       ntime, pelem_ID, interpolation_level,        &
                                       idomain_g, ielement_g, coordinate_system,    & 
                                       recv_size_a, recv_size_b, recv_size_c,       &
                                       neighbor_location(7), element_location(5),       &
                                       element_data(9), spacedim, inode, dof_start, &
                                       ineighbor_domain_g, ineighbor_element_g,     &
                                       recv_dof, ChiID, idonor, idonor_domain_g,    &
                                       idonor_element_g
        logical :: interior_face, parallel_neighbor, parallel_donor


        !
        ! Receive interior face neighbor DATA
        !
        recv_procs = self%get_recv_procs()
        do iproc = 1,size(recv_procs)
            do irecv = 1,self%get_proc_ninterior_neighbors(recv_procs(iproc))

                ! The sending proc sent its neighbor location, which is a face on our local processor 
                ! here where we will store the incoming ALE data
                ! face_location = [idomain_g, idomain_l, ielement_g, ielement_l, iface, dof_start]
                call mpi_recv(neighbor_location, size(neighbor_location), mpi_integer4, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                idom  = neighbor_location(2)
                ielem = neighbor_location(4)
                iface = neighbor_location(5)

                recv_size_a = size(self%domain(idom)%faces(ielem,iface)%neighbor_interp_coords_vel)
                recv_size_b = size(self%domain(idom)%faces(ielem,iface)%neighbor_ale_g)
                recv_size_c = size(self%domain(idom)%faces(ielem,iface)%neighbor_ale_Dinv)

                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_interp_coords_vel, recv_size_a, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_ale_g,             recv_size_b, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_ale_g_grad1,       recv_size_b, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_ale_g_grad2,       recv_size_b, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_ale_g_grad3,       recv_size_b, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_ale_Dinv,          recv_size_c, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)

            end do !irecv
        end do !iproc




        !
        ! Receive/construct parallel neighbors/chimera donors
        !
        do iproc = 1,size(recv_procs)
            do irecv = 1,(self%get_proc_ninterior_neighbors(recv_procs(iproc)) + self%get_proc_nchimera_donors(recv_procs(iproc)))


                ! element_location = [idomain_g, idomain_l, ielement_g, ielement_l, iproc]
                call mpi_recv(element_location, 5, mpi_integer4,  recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                idomain_g  = element_location(1)
                ielement_g = element_location(3)



                ! element_data = [element_type, spacedim, coordinate_system, nfields, nterms_s, nterms_c, ntime, interpolation_level]
                call mpi_recv(element_data, 9, mpi_integer4,  recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                etype               = element_data(1)
                spacedim            = element_data(2)
                coordinate_system   = element_data(3)
                nfields             = element_data(4)
                nterms_s            = element_data(5)
                nterms_c            = element_data(6)
                ntime               = element_data(7)
                interpolation_level = element_data(8)
                dof_start           = element_data(9)
                nnodes = (etype+1)*(etype+1)*(etype+1)
                

                ! Allocate buffers and receive: nodes, displacements, and velocities. 
                ! These quantities are located at the element support nodes, not interpolation
                ! nodes.
                if (allocated(nodes)) deallocate(nodes, nodes_def, nodes_vel, connectivity)
                allocate(nodes(       nnodes,3), &
                         nodes_def(   nnodes,3), &
                         nodes_vel(   nnodes,3), &
                         connectivity(nnodes  ), stat=ierr)
                if (ierr /= 0) call AllocationError


                call mpi_recv(nodes,     nnodes*3, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(nodes_def, nnodes*3, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(nodes_vel, nnodes*3, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)

                
                !
                ! Compute node displacements
                !
                nodes_disp = nodes_def - nodes

                !
                ! Build local connectivity
                !   : we construct the parallel element using just a local ordering
                !   : so connectivity starts at 1 and goes to the number of nodes in 
                !   : the element, nnodes.
                !   :
                !   :   connectivity = [1, 2, 3, 4, 5, 6, 7, 8 ...]
                !   :
                !   : We assume here that the displacements and velocities are ordered
                !   : appropriately.
                !
                do inode = 1,nnodes
                    connectivity(inode) = inode
                end do



                !
                ! Check for existing parallel element. If one does not
                ! exist, get an identifier for a new parallel element.
                !
                pelem_ID = self%find_parallel_element(idomain_g,ielement_g)

                if (pelem_ID == NO_ID) pelem_ID = self%new_parallel_element()


                ! Initialize element geometry
                select case(coordinate_system)
                    case(CARTESIAN)
                        coord_system = 'Cartesian'
                    case(CYLINDRICAL)
                        coord_system = 'Cylindrical'
                    case default
                        call chidg_signal(FATAL,"element%comm_recv: invalid coordinate system.")
                end select

                
                ! Construct/initialize/reinitialize parallel element
                if (.not. self%parallel_element(pelem_ID)%geom_initialized) then
                    call self%parallel_element(pelem_ID)%init_geom(nodes,connectivity,etype,element_location,trim(coord_system))
                end if


                ! NOTE: we want to be able to reinitialize parallel_element solution data-space (e.g. for wall_distance calculation)
                call self%parallel_element(pelem_ID)%init_sol(self%interpolation,interpolation_level,nterms_s,nfields,ntime,dof_start,dof_local_start=NO_ID)


                call self%parallel_element(pelem_ID)%set_displacements_velocities(nodes_disp,nodes_vel)
                call self%parallel_element(pelem_ID)%update_interpolations_ale()



            end do !irecv
        end do !iproc



        !
        ! Update recv_dof indices for all parallel elements. 
        ! NOTE: this should be done after the previous loop structure to allow parallel elements
        ! to be reinitialized (for example if we are reinitializing nterms_s).
        !
        recv_dof = 1    ! DOF indices are based on Fortran 1-indexing
        do pelem_ID = 1,self%nparallel_elements()

            ! Set current element
            self%parallel_element(pelem_ID)%recv_dof = recv_dof

            ! Update start of potential next element
            nterms_s = self%parallel_element(pelem_ID)%nterms_s 
            nfields  = self%parallel_element(pelem_ID)%nfields
            ntime    = self%parallel_element(pelem_ID)%ntime

            recv_dof = recv_dof + nterms_s*nfields*ntime

        end do



        !
        ! Initialize interior parallel neighbor element ID and recv_dof
        !
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%nelements()
                do iface = 1,NFACES

                    interior_face     = (self%domain(idom)%faces(ielem,iface)%ftype == INTERIOR)
                    parallel_neighbor = (self%domain(idom)%faces(ielem,iface)%ineighbor_proc /= IRANK)

                    ! For INTERIOR faces that have off-processor neighbors we need to find the parallel 
                    ! element ID in the local parallel element list
                    if ( interior_face .and. parallel_neighbor ) then

                        ineighbor_domain_g   = self%domain(idom)%faces(ielem,iface)%ineighbor_domain_g
                        ineighbor_element_g  = self%domain(idom)%faces(ielem,iface)%ineighbor_element_g
                        pelem_ID = self%find_parallel_element(ineighbor_domain_g,ineighbor_element_g)
                        if (pelem_ID == NO_ID) call chidg_signal(FATAL,'mesh%comm_recv: could not find interior neighbor parallel element.')

                        self%domain(idom)%faces(ielem,iface)%ineighbor_pelem_ID = pelem_ID
                        self%domain(idom)%faces(ielem,iface)%recv_dof           = self%parallel_element(pelem_ID)%recv_dof

                    end if

                end do !iface
            end do !ielem
        end do !idom



        !
        ! Initialize chimera parallel neighbor element ID and recv_dof
        !
        do idom = 1,self%ndomains()
            do ChiID = 1,self%domain(idom)%chimera%nreceivers()
                do idonor = 1,self%domain(idom)%chimera%recv(ChiID)%ndonors()

                    parallel_donor = (self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%iproc /= IRANK)
                    if (parallel_donor) then
                        idonor_domain_g  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_g
                        idonor_element_g = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_g
                        pelem_ID = self%find_parallel_element(idonor_domain_g,idonor_element_g)
                        if (pelem_ID == NO_ID) call chidg_signal(FATAL,'mesh%comm_recv: could not find overset donor parallel element.')

                        self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%pelem_ID = pelem_ID
                        self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%recv_dof = self%parallel_element(pelem_ID)%recv_dof
                    end if

                end do !idonor
            end do !ChiID
        end do !idom



        !
        ! Assemble overset interpolations on receiver exterior state
        !
        call self%assemble_chimera_data()


    end subroutine comm_recv
    !*********************************************************************************





    !>  Wait until all outstanding requests initiated by mesh%comm_send() have
    !!  been completed.
    !!
    !!  Order of operations:
    !!      1: mesh%comm_send()
    !!      2: mesh%comm_recv()
    !!      3: mesh%comm_wait()
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_wait(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik) :: nwait, ierr

        nwait = self%comm_requests%size()
        call mpi_waitall(nwait, self%comm_requests%data(1:nwait), mpi_statuses_ignore, ierr)
        call self%comm_requests%clear()

    end subroutine comm_wait
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



    !>  Return number of elements on the processor-local mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/21/2018
    !!
    !--------------------------------------------------------------------------------
    function nelements(self) result(nelem)
        class(mesh_t),  intent(in)  :: self

        integer(ik) :: nelem, idom

        nelem = 0
        do idom = 1,self%ndomains()
            nelem = nelem + self%domain(idom)%nelements()
        end do

    end function nelements
    !********************************************************************************





    !>  Return the global dof starting index for a given element.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/24/2019
    !!
    !--------------------------------------------------------------------------------
    function get_dof_start(self,idomain_l,ielement_l) result(dof_start)
        class(mesh_t),  intent(in)  :: self
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_l


        integer(ik) :: dof_start, ndof_prev, idom, end_elem, ielem

        ! Sum DOF on ranks 0 to IRANK-1. Indexing goes to IRANK since Fortran starts at 1.
        ! WARNING: ASSUMING ALL ELEMENTS ON OTHER PROCS HAVE SAME nterms_s
        ndof_prev = self%mesh_dof_start - 1


        ! Accumulate dofs on current RANK until idomain_l,ielement_l
        do idom = 1,idomain_l

            if (idom < idomain_l) then
                end_elem = self%domain(idom)%nelements()
            else if (idom == idomain_l) then
                end_elem = ielement_l-1
            end if

            do ielem = 1,end_elem
                ndof_prev = ndof_prev + self%domain(idom)%elems(ielem)%nterms_s*self%domain(idom)%elems(ielem)%nfields*self%domain(idom)%elems(ielem)%ntime
            end do

        end do


        ! dof index for requested element starts on dof after all previous
        dof_start = ndof_prev + 1


    end function get_dof_start
    !********************************************************************************






    !>  Return element info object.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   2/27/2019
    !!
    !--------------------------------------------------------------------------------
    function get_element_info(self,idomain_l,ielement_l) result(elem_info)
        class(mesh_t),  intent(in)  :: self
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_l

        type(element_info_t)    :: elem_info

        elem_info = element_info(idomain_g       = self%domain(idomain_l)%elems(ielement_l)%idomain_g,       &
                                 idomain_l       = self%domain(idomain_l)%elems(ielement_l)%idomain_l,       &
                                 ielement_g      = self%domain(idomain_l)%elems(ielement_l)%ielement_g,      &
                                 ielement_l      = self%domain(idomain_l)%elems(ielement_l)%ielement_l,      &
                                 iproc           = self%domain(idomain_l)%elems(ielement_l)%iproc,           &
                                 pelem_ID        = NO_ID,                                                    &
                                 eqn_ID          = self%domain(idomain_l)%elems(ielement_l)%eqn_ID,          &
                                 nfields         = self%domain(idomain_l)%elems(ielement_l)%nfields,         &
                                 ntime           = self%domain(idomain_l)%elems(ielement_l)%ntime,           &
                                 nterms_s        = self%domain(idomain_l)%elems(ielement_l)%nterms_s,        &
                                 nterms_c        = self%domain(idomain_l)%elems(ielement_l)%nterms_c,        &
                                 dof_start       = self%domain(idomain_l)%elems(ielement_l)%dof_start,       &
                                 dof_local_start = self%domain(idomain_l)%elems(ielement_l)%dof_local_start, &
                                 recv_comm       = self%domain(idomain_l)%elems(ielement_l)%recv_comm,       &
                                 recv_domain     = self%domain(idomain_l)%elems(ielement_l)%recv_domain,     &
                                 recv_element    = self%domain(idomain_l)%elems(ielement_l)%recv_element,    &
                                 recv_dof        = self%domain(idomain_l)%elems(ielement_l)%recv_dof)


    end function get_element_info
    !********************************************************************************




    !>
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine init_comm_local(self)
        class(mesh_t),  intent(inout)  :: self

        integer(ik) :: idom

        ! All ranks initialize local communication
        call write_line("Initialize: proc-local neighbor communication...", io_proc=GLOBAL_MASTER)
        do idom = 1,self%ndomains()
            call self%domain(idom)%init_comm_local()
        end do

    end subroutine init_comm_local
    !********************************************************************************





    !>
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine init_comm_global(self,ChiDG_COMM)
        class(mesh_t),  intent(inout)   :: self
        type(mpi_comm), intent(in)      :: ChiDG_COMM

        integer(ik) :: iproc, idom, idomain_g, nfaces_search, ierr

        integer(ik), allocatable :: face_search_corners(:,:), face_owner_rank(:), face_owner_rank_reduced(:)
        logical     :: has_domain, searching_mesh, searching

        call write_line("Initialize: proc-global neighbor communication...", io_proc=GLOBAL_MASTER)
        do iproc = 0,NRANK-1

            ! For current rank, send out request for neighbors
            if (iproc == IRANK) then

                ! Search
                do idom = 1,self%ndomains()

                    ! Broadcast that a mesh from iproc is searching for face neighbors
                    searching_mesh=.true.
                    call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                    ! send global domain index of mesh being searched
                    call MPI_BCast(self%domain(idom)%idomain_g,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)

                    ! Trigger comm for current domain
                    call self%domain(idom)%init_comm_global(ChiDG_COMM)

                end do

                ! Broadcast that a mesh from iproc is searching for face neighbors
                searching_mesh=.false.
                call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

            end if !iproc


            if (iproc /= IRANK) then

                    ! Does iproc have any domains in mesh
                    searching_mesh=.true.
                    call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                    do while (searching_mesh)

                        ! What domain is being searched
                        call MPI_BCast(idomain_g,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)

                        ! Get information about faces being searched for
                        call MPI_BCast(nfaces_search,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)


                        ! Allocate buffer to receive face search data
                        if (allocated(face_search_corners)) deallocate(face_search_corners)
                        if (allocated(face_owner_rank)) deallocate(face_owner_rank)
                        if (allocated(face_owner_rank_reduced)) deallocate(face_owner_rank_reduced)
                        allocate(face_search_corners(nfaces_search,NFACE_CORNERS), &
                                 face_owner_rank(nfaces_search), &
                                 face_owner_rank_reduced(nfaces_search), stat=ierr)
                        if (ierr /= 0) call AllocationError

                        ! Get faces being searched for
                        call MPI_BCast(face_search_corners,nfaces_search*NFACE_CORNERS,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)


                        ! Search through local domains and check if we have part of domain being searched
                        has_domain = .false.
                        do idom = 1,self%ndomains()
                            has_domain = ( idomain_g == self%domain(idom)%idomain_g )
                            if ( has_domain ) exit
                        end do
                        
                        ! Initialize face_owner_rank, and try and find owner on local process
                        face_owner_rank = NO_PROC
                        if (has_domain) call self%domain(idom)%find_face_owner(face_search_corners,face_owner_rank)


                        ! Synchronize result across all ranks
                        call MPI_Reduce(face_owner_rank,face_owner_rank_reduced,nfaces_search,MPI_INTEGER4,MPI_MAX,iproc,ChiDG_COMM,ierr)
                        if (ierr /= 0) call chidg_signal(FATAL,'mesh%init_comm_global: error reducing face owners.')

                        ! For each face we found, send information
                        if (has_domain) call self%domain(idom)%transmit_face_info(face_search_corners,face_owner_rank,iproc,ChiDG_COMM)

                        ! Is there another domain to be searched
                        call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                    end do !searching_mesh

            end if !iproc


        end do !iproc
        

    end subroutine init_comm_global
    !********************************************************************************






end module type_mesh
