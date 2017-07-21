module type_mesh
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NO_ID, INTERIOR, NFACES
    use type_domain,                only: domain_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    use type_bc_patch,              only: bc_patch_t
    use type_bc_patch_group,        only: bc_patch_group_t
    use type_ivector,               only: ivector_t
    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08,                    only: mpi_isend, mpi_recv, mpi_integer4, mpi_real8, &
                                          mpi_waitall, mpi_request, mpi_status_ignore,  &
                                          mpi_statuses_ignore
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

        ! Boundary patch procedures
        procedure           :: add_bc_patch

        ! Resouce management
        procedure           :: release


        ! Parallel communication pattern
        procedure           :: get_recv_procs
        procedure           :: get_send_procs

        procedure           :: get_proc_ninterior_neighbors
        !procedure           :: get_proc_nchimera_donors
        !procedure           :: get_proc_nchimera_receivers


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
        call self%domain(idomain_l)%init_ale(dnodes,vnodes)


        
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
        ! Find patch group identifier. If none, call new_bc_patch_group.
        !   - don't create patches for 'Empty' groups. These should be left ORPHAN
        !     so that they get picked up by Chimera.
        !
        !if ( (trim(group_name) /= 'Empty') .and. &
        !     (trim(group_name) /= 'empty') ) then


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


        !end if



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











    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_send(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik)         :: idom, ielem, iface, idonor, iproc, ierr,    &
                               send_size_a, send_size_b, send_size_c, send_size_d
        type(mpi_request)   :: request(4)
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
                        send_size_b = size(face%grid_vel)
                        send_size_c = size(face%det_jacobian_grid)
                        send_size_d = size(face%inv_jacobian_grid)

                        ! First, send neighbor location. This way, the receiving processor knows where to put the data.
                        ! Next, send all ALE information
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%neighbor_location, send_size_a, mpi_integer4, face%ineighbor_proc, 0, ChiDG_COMM, request(1), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%grid_vel,          send_size_b, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(2), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%det_jacobian_grid, send_size_c, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(3), ierr)
                        call mpi_isend(self%domain(idom)%faces(ielem,iface)%inv_jacobian_grid, send_size_d, mpi_real8,    face%ineighbor_proc, 0, ChiDG_COMM, request(4), ierr)

                        call self%comm_requests%push_back(request(1))
                        call self%comm_requests%push_back(request(2))
                        call self%comm_requests%push_back(request(3))
                        call self%comm_requests%push_back(request(4))

                        end associate

                    end if


                end do !iface
            end do !ielem
        end do !idom



!        !
!        ! Send chimera donors
!        !
!        do idom = 1,self%ndomains()
!            do idonor = 1,self%domain(idom)%chimera%send%ndonors()
!
!                iproc = self%domain(idom)%chimera%send%receiver_proc%at(idonor)
!
!                ! If receiver is off-processor, send reference and physical nodes/velocities
!                if (iproc /= IRANK) then
!
!                    idomain_l  = self%domain(idom)%chimera%send%donors(idonor)%idomain_l
!                    ielement_l = self%domain(idom)%chimera%send%donors(idonor)%ielement_l
!
!                    send_size_a = size(self%domain(idomain_l)%elems(ielement_l)%elem_pts)
!                    send_size_b = size(self%domain(idomain_l)%elems(ielement_l)%dnodes_l)
!                    send_size_c = size(self%domain(idomain_l)%elems(ielement_l)%vnodes_l)
!            
!                    ! First send location of donor
!                    call mpi_isend(receiver_location                                                    4, mpi_integer4, iproc, 0, ChiDG_COMM, request(1), ierr)
!                    call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%element_location,           4, mpi_integer4, iproc, 0, ChiDG_COMM, request(1), ierr)
!                    call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%elem_pts,         send_size_a, mpi_real8,    iproc, 0, ChiDG_COMM, request(2), ierr)
!                    call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%dnodes_l,         send_size_b, mpi_real8,    iproc, 0, ChiDG_COMM, request(3), ierr)
!                    call mpi_isend(self%domain(idomain_l)%elems(ielement_l)%vnodes_l,         send_size_c, mpi_real8,    iproc, 0, ChiDG_COMM, request(4), ierr)
!
!                    call self%comm_requests%push_back(request(1))
!                    call self%comm_requests%push_back(request(2))
!                    call self%comm_requests%push_back(request(3))
!                    call self%comm_requests%push_back(request(4))
!
!                end if 
!
!            end do !idonor
!        end do !idom



    end subroutine comm_send
    !*********************************************************************************





    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(mesh_t),  intent(inout)   :: self

        integer(ik),    allocatable :: recv_procs(:)
        integer(ik)                 :: idom, ielem, iface, iproc, irecv, ierr, &
                                       recv_size_a, recv_size_b, recv_size_c, face_location(5)


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

                recv_size_a = size(self%domain(idom)%faces(ielem,iface)%neighbor_grid_vel)
                recv_size_b = size(self%domain(idom)%faces(ielem,iface)%neighbor_det_jacobian_grid)
                recv_size_c = size(self%domain(idom)%faces(ielem,iface)%neighbor_inv_jacobian_grid)

                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_grid_vel,          recv_size_a, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_det_jacobian_grid, recv_size_b, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)
                call mpi_recv(self%domain(idom)%faces(ielem,iface)%neighbor_inv_jacobian_grid, recv_size_c, mpi_real8, recv_procs(iproc), 0, ChiDG_COMM, mpi_status_ignore, ierr)

            end do !irecv
        end do !iproc





        !
        ! Receive chimera donors
        !





    end subroutine comm_recv
    !*********************************************************************************





    !>  Wait until all outstanding requests initiated by mesh%comm_send() have
    !!  been completed.
    !!
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
