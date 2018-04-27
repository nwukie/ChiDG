module type_mesh
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NO_ID, INTERIOR, NFACES, CARTESIAN, CYLINDRICAL
    use type_element,               only: element_t
    use type_domain,                only: domain_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use type_bc_patch,              only: bc_patch_t
    use type_bc_patch_group,        only: bc_patch_group_t
    use type_ivector,               only: ivector_t
    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08,                    only: mpi_isend, mpi_recv, mpi_integer4, mpi_real8, &
                                          mpi_waitall, mpi_request, mpi_status_ignore,  &
                                          mpi_statuses_ignore, mpi_character
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

        ! Local data
        type(domain_t),         allocatable :: domain(:)
        type(bc_patch_group_t), allocatable :: bc_patch_group(:)

        ! Parallel data
        type(element_t),        allocatable :: parallel_element(:)
        type(mpi_request_vector_t)          :: comm_requests

    contains

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


        ! Resouce management
        procedure           :: release


        ! Parallel communication pattern
        procedure           :: get_recv_procs
        procedure           :: get_send_procs

        procedure           :: get_proc_ninterior_neighbors
        procedure           :: get_proc_nchimera_donors
        !procedure           :: get_proc_nchimera_receivers
        procedure           :: get_nelements_recv

        procedure           :: find_parallel_element
        procedure           :: new_parallel_element
        procedure           :: nparallel_elements


        ! Extra routines for testing private procedures
        procedure           :: stub_new_domain

        ! Communication 
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
        real(rk),                       intent(in)      :: vnodes(:,:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        integer(ik),                    intent(in)      :: nelements_g
        character(*),                   intent(in)      :: coord_system
        integer(ik),                    intent(in)      :: eqn_ID

        integer(ik) :: idomain_l


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
                  donor_proc = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%iproc
                  ref_coords = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%coords
                  npts = size(ref_coords,1)
        
                  parallel_donor = (donor_proc /= IRANK)
        
                  if (parallel_donor) then
                      idomain_g   = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_g
                      ielement_g  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_g
                      pelem_ID    = self%find_parallel_element(idomain_g,ielement_g)
                      if (pelem_ID == NO_ID) call chidg_signal_three(FATAL,"mesh%assemble_chimera_data: did not find parallel chimera element.", IRANK, idomain_g, ielement_g)
        
                      do ipt = 1, npts
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
        
                      idomain_l   = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_l
                      ielement_l  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_l
        
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

                    proc_has_donor = (self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%iproc == iproc)
                    if (proc_has_donor) then

                        donor_domain_g  = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_g
                        donor_element_g = self%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_g

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
    !!  NOTE: currently, we are only receiving off-processor chimera donors into 
    !!        self%parallel_elements. So, even though we may have interior neighbors
    !!        that are off-processor, those elements are not included in this count
    !!        atthe moment, because they are not stored in self%parallel_elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/26/2017
    !!
    !!
    !-----------------------------------------------------------------------------------
    function get_nelements_recv(self) result(n)
        class(mesh_t),  intent(in)  :: self

        integer(ik)     :: idom, irecv, idonor, donor_iproc, idomain_g, ielement_g, &
                           idomain_g_list, ielement_g_list, ientry, n
        logical         :: parallel_donor, already_counted
        type(ivector_t) :: donor_domain_g, donor_element_g


        do idom = 1,self%ndomains()
            do irecv = 1,self%domain(idom)%chimera%nreceivers()
                do idonor = 1,self%domain(idom)%chimera%recv(irecv)%ndonors()

                    donor_iproc = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%iproc

                    parallel_donor = (donor_iproc /= IRANK)
                    if (parallel_donor) then
                        idomain_g  = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%idomain_g
                        ielement_g = self%domain(idom)%chimera%recv(irecv)%donor(idonor)%ielement_g

                        already_counted = .false.
                        do ientry = 1,donor_element_g%size()
                            idomain_g_list  = donor_domain_g%at(ientry)
                            ielement_g_list = donor_element_g%at(ientry)
                            if ( (idomain_g_list  == idomain_g ) .and. &
                                 (ielement_g_list == ielement_g) ) already_counted = .true.
                            if (already_counted) exit
                        end do

                        if (.not. already_counted) then
                            call donor_domain_g%push_back(idomain_g)
                            call donor_element_g%push_back(ielement_g)
                        end if
                    end if


                end do !idonor
            end do !irecv
        end do !idom



        !
        ! Get size of accumulated elements being received
        !
        n = donor_element_g%size()

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

        
        !
        ! Resize array storage
        !
        allocate(temp(self%nparallel_elements() + 1), stat=ierr)



        ! Copy previously initialized instances to new array. Be careful about pointers 
        ! components here. For example, a pointer from a face to an element would no 
        ! longer be valid in the new array.
        if (self%nparallel_elements() > 0) then
            temp(1:size(self%parallel_element)) = self%parallel_element(1:size(self%parallel_element))
        end if



        !
        ! Move resized temp allocation back to mesh container. 
        ! Be careful about pointer components here! Their location in memory has changed.
        !
        call move_alloc(temp,self%parallel_element)
        


        !
        ! Set domain identifier of newly allocated domain that will be returned
        !
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
        ! Send interior face data
        !
        do idom = 1,self%ndomains()
            do ielem = 1,self%domain(idom)%nelem
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
        ! Send chimera donors
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
                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 8, mpi_integer4, iproc, 0, ChiDG_COMM, request(2), ierr)
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


!        do group_ID = 1,self%nbc_patch_groups()
!            do patch_ID = 1,self%bc_patch_group(group_ID)%npatches()
!                if (self%bc_patch_group(group_ID)%patch(patch_ID)%spatial_coupling = 'Global') then
!
!                    ! Send each element in the patch to each other processor in bc_comm
!                    call MPI_Comm_Size(bc_comm, bc_NRANK, ierr)
!                    call MPI_Comm_Rank(bc_comm, bc_IRANK, ierr)
!
!                    do iproc = 0,
!
!                end if !global spatial coupling
!            end do !patch_ID
!        end do !group_ID


!        !
!        ! Send coupled bc elements
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
!                        call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_data,                 8, mpi_integer4, iproc, 0, ChiDG_COMM, request(2), ierr)
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
                                       face_location(5), element_location(5),       &
                                       element_data(8), spacedim, inode

        logical :: parallel_donor

        real(rk), allocatable :: ref_coords(:,:)
        real(rk) :: xi, eta, zeta, det_jacobian_grid_pt, det_jacobian_grid_grad_pt(3),&
            inv_jacobian_grid_pt(3,3), grid_vel_pt(3)
        integer(ik) :: npts, ipt, donor_proc, idonor, idomain_l, ielement_l,  igq, ChiID


        !
        ! Receive interior face data
        !
        recv_procs = self%get_recv_procs()
        do iproc = 1,size(recv_procs)
            do irecv = 1,self%get_proc_ninterior_neighbors(recv_procs(iproc))

                ! The sending proc sent its neighbor location, which is a face on our local processor 
                ! here where we will store the incoming ALE data
                ! face_location = [idomain_g, idomain_l, ielement_g, ielement_l, iface]
                call mpi_recv(face_location, 5, mpi_integer4, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                idom  = face_location(2)
                ielem = face_location(4)
                iface = face_location(5)

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
        ! Receive/construct parallel chimera donors
        !
        do iproc = 1,size(recv_procs)
            do irecv = 1,self%get_proc_nchimera_donors(recv_procs(iproc))



                ! element_location = [idomain_g, idomain_l, ielement_g, ielement_l, iproc]
                call mpi_recv(element_location, 5, mpi_integer4,  recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                idomain_g  = element_location(1)
                ielement_g = element_location(3)


                ! element_data = [element_type, spacedim, coordinate_system, nfields, nterms_s, nterms_c, ntime, interpolation_level]
                call mpi_recv(element_data,     8, mpi_integer4,  recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                etype               = element_data(1)
                spacedim            = element_data(2)
                coordinate_system   = element_data(3)
                nfields             = element_data(4)
                nterms_s            = element_data(5)
                nterms_c            = element_data(6)
                ntime               = element_data(7)
                interpolation_level = element_data(8)
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


                !
                ! Initialize element geometry
                !
                select case(coordinate_system)
                    case(CARTESIAN)
                        coord_system = 'Cartesian'
                    case(CYLINDRICAL)
                        coord_system = 'Cylindrical'
                    case default
                        call chidg_signal(FATAL,"element%comm_recv: invalid coordinate system.")
                end select

                
                !
                ! Construct/initialize/reinitialize parallel element
                !
                if (.not. self%parallel_element(pelem_ID)%geom_initialized) then
                    call self%parallel_element(pelem_ID)%init_geom(nodes,connectivity,etype,element_location,trim(coord_system))
                end if


                call self%parallel_element(pelem_ID)%init_sol('Quadrature',interpolation_level,nterms_s,nfields,ntime)
                call self%parallel_element(pelem_ID)%set_displacements_velocities(nodes_disp,nodes_vel)
                call self%parallel_element(pelem_ID)%update_interpolations_ale()


            end do !irecv
        end do !iproc



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




end module type_mesh
