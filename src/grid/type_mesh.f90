module type_mesh
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NFACES,XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                          ORPHAN, INTERIOR, BOUNDARY, CHIMERA, TWO_DIM, THREE_DIM, NO_NEIGHBOR_FOUND, NEIGHBOR_FOUND
    use mod_grid,                   only: FACE_CORNERS
    use mod_chidg_mpi,              only: IRANK,NRANK
    use mpi_f08

    use type_point,                 only: point_t
    use type_element,               only: element_t
    use type_face,                  only: face_t
    use type_chimera,               only: chimera_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_element_connectivity,  only: element_connectivity_t
    implicit none
    private


    !> Data type for mesh information
    !!      - contains array of elements, array of faces for each element
    !!      - calls initialization procedure for elements and faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/27/2016
    !!
    !------------------------------------------------------------------------------------------------------------
    type, public :: mesh_t

        !
        ! Integer parameters
        !
        integer(ik)                     :: spacedim   = 0               !< Number of spatial dimensions
        integer(ik)                     :: neqns      = 0               !< Number of equations being solved
        integer(ik)                     :: nterms_s   = 0               !< Number of terms in the solution expansion
        integer(ik)                     :: nterms_c   = 0               !< Number of terms in the grid coordinate expansion
        integer(ik)                     :: nelem      = 0               !< Number of total elements
        integer(ik)                     :: ntime      = 0               !< Number of time instances

        !
        ! Grid data
        !
        integer(ik)                     :: idomain_g
        integer(ik)                     :: idomain_l
        
        type(point_t),    allocatable   :: nodes(:)                     !< Original node points for the domain
        type(element_t),  allocatable   :: elems(:)                     !< Element storage (1:nelem)
        type(face_t),     allocatable   :: faces(:,:)                   !< Face storage    (1:nelem,1:nfaces)
        type(chimera_t)                 :: chimera                      !< Chimera interface data

        !
        ! Initialization flags
        !
        logical                         :: geomInitialized = .false.    !< Status of geometry initialization
        logical                         :: solInitialized  = .false.    !< Status of numerics initialization

    contains

        procedure           :: init_geom
        procedure           :: init_sol

        procedure, private  :: init_elems_geom
        procedure, private  :: init_elems_sol
        procedure, private  :: init_faces_geom
        procedure, private  :: init_faces_sol

        procedure           :: init_comm_local
        procedure           :: init_comm_global

        ! Utilities
        procedure, private  :: find_neighbor_local
        procedure, private  :: find_neighbor_global
        procedure           :: handle_neighbor_request

        final               :: destructor

    end type mesh_t
    !************************************************************************************************************





contains






    !>  Mesh geometry initialization procedure
    !!
    !!  Sets number of terms in coordinate expansion for the entire domain
    !!  and calls sub-initialization routines for individual element and face geometry
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  nterms_c    Number of terms in the coordinate expansion
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine init_geom(self,idomain_l,spacedim,nterms_c, nodes, connectivity)
        class(mesh_t),                  intent(inout), target   :: self
        integer(ik),                    intent(in)              :: idomain_l
        integer(ik),                    intent(in)              :: spacedim
        integer(ik),                    intent(in)              :: nterms_c
        type(point_t),                  intent(in)              :: nodes(:)
        type(domain_connectivity_t),    intent(in)              :: connectivity


        !
        ! Store number of terms in coordinate expansion and domain index
        !
        self%spacedim  = spacedim
        self%nterms_c  = nterms_c
        self%idomain_g = connectivity%get_domain_index()
        self%idomain_l = idomain_l
        self%nodes     = nodes


        !
        ! Call geometry initialization for elements and faces
        !
        call self%init_elems_geom(spacedim,nodes,connectivity)
        call self%init_faces_geom(spacedim,nodes,connectivity)


        !
        ! Confirm initialization
        !
        self%geomInitialized = .true.


    end subroutine init_geom
    !************************************************************************************************************











    !>  Mesh numerics initialization procedure
    !!
    !!  Sets number of equations being solved, number of terms in the solution expansion and
    !!  calls sub-initialization routines for individual element and face numerics
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  neqns       Number of equations being solved in the current domain
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine init_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s


        !
        ! Store number of equations and number of terms in solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s


        !
        ! Call numerics initialization for elements and faces
        !
        call self%init_elems_sol(neqns,nterms_s) 
        call self%init_faces_sol()               

        
        !
        ! Confirm initialization
        !
        self%solInitialized = .true.


    end subroutine init_sol
    !************************************************************************************************************












    !>  Mesh - element initialization procedure
    !!
    !!  Computes the number of elements based on the element mapping selected and
    !!  calls the element initialization procedure on individual elements.
    !!
    !!  TODO: Generalize for non-block structured ness. Eliminate dependence on, xi, eta, zeta directions.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine init_elems_geom(self,spacedim,nodes,connectivity)
        class(mesh_t),                  intent(inout)   :: self
        integer(ik),                    intent(in)      :: spacedim
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity


        type(point_t),  allocatable     :: points_l(:)
        type(element_connectivity_t)    :: element_connectivity

        integer(ik)                ::   ierr,     ipt,       ielem_l,           &
                                        ipt_xi,   ipt_eta,   ipt_zeta,          &
                                        ixi,      ieta,      izeta,             &
                                        xi_start, eta_start, zeta_start,        &
                                        nelem_xi, nelem_eta, nelem_zeta, nelem, &
                                        neqns,    nterms_s,  nnodes, nterms_c,  &
                                        npts_1d, mapping, idomain_l, inode


        !
        ! Compute number of 1d points for a single element
        !
        npts_1d = 0
        
        if ( spacedim == THREE_DIM ) then
            do while (npts_1d*npts_1d*npts_1d < self%nterms_c)
                npts_1d = npts_1d + 1       ! really just computing the cubed root of nterms_c, the number of terms in the coordinate expansion
            end do

        else if ( spacedim == TWO_DIM ) then
            do while (npts_1d*npts_1d < self%nterms_c)
                npts_1d = npts_1d + 1       ! really just computing the cubed root of nterms_c, the number of terms in the coordinate expansion
            end do

        end if


        !
        ! Store number of elements in each direction along with total number of elements
        !
        nelem           = connectivity%get_nelements()
        self%nelem      = nelem
        mapping         = (npts_1d - 1)     !> 1 - linear, 2 - quadratic, 3 - cubic, etc.

        !
        ! Allocate element storage
        !
        allocate(self%elems(nelem), points_l(self%nterms_c), stat=ierr)
        if(ierr /= 0) stop "Memory allocation error: init_elements"

        !
        ! Call geometry initialization for each element
        !
        idomain_l = self%idomain_l
        do ielem_l = 1,nelem

            ! Element geometry initialization
            element_connectivity = connectivity%get_element_connectivity(ielem_l)
            call self%elems(ielem_l)%init_geom(spacedim,nodes,element_connectivity,idomain_l,ielem_l)

        end do ! ielem


    end subroutine init_elems_geom
    !**************************************************************************************************************















    !>  Mesh - element solution data initialization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  neqns       Number of equations in the domain equation set
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine init_elems_sol(self,neqns,nterms_s)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: neqns
        integer(ik),    intent(in)      :: nterms_s
        integer(ik) :: ielem


        !
        ! Store number of equations and number of terms in the solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s


        !
        ! Call the numerics initialization procedure for each element
        !
        do ielem = 1,self%nelem
            call self%elems(ielem)%init_sol(self%neqns,self%nterms_s)
        end do


    end subroutine init_elems_sol
    !***************************************************************************************************************















    !>  Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_faces_geom(self,spacedim,nodes,connectivity)
        class(mesh_t),                  intent(inout)   :: self
        integer(ik),                    intent(in)      :: spacedim
        type(point_t),                  intent(in)      :: nodes(:)
        type(domain_connectivity_t),    intent(in)      :: connectivity

        integer(ik)     :: ielem, iface, ierr

        !
        ! Allocate face storage array
        !
        allocate(self%faces(self%nelem,NFACES),stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"mesh%init_faces_geom -- face allocation error")


        !
        ! Loop through each element and call initialization for each face
        !
        do ielem = 1,self%nelem
            do iface = 1,NFACES

                ! Call face geometry initialization
                call self%faces(ielem,iface)%init_geom(iface,self%elems(ielem))


            end do !iface
        end do !ielem




    end subroutine init_faces_geom
    !**************************************************************************************************************
















    !>  Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine init_faces_sol(self)
        class(mesh_t), intent(inout)  :: self

        integer(ik) :: ielem, iface

        !
        ! Loop through elements
        !
        do ielem = 1,self%nelem

            !
            ! Loop through faces and call numerics initialization routine
            !
            do iface = 1,NFACES

                call self%faces(ielem,iface)%init_sol(self%elems(ielem))

            end do ! iface

        end do ! ielem


    end subroutine init_faces_sol
    !***************************************************************************************************************















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine init_comm_local(self)
        class(mesh_t),                  intent(inout)   :: self


        integer(ik)             :: ixi,ieta,izeta,iface,ftype,ielem,ierr, ielem_neighbor
        integer(ik)             :: corner_one, corner_two, corner_three, corner_four
        integer(ik)             :: node_indices(4)
        logical                 :: boundary_face = .false.
        logical                 :: includes_node_one, includes_node_two, includes_node_three, includes_node_four
        logical                 :: neighbor_element

        integer(ik)             :: ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l, ineighbor_proc, neighbor_status

        !
        ! Loop through each local element and call initialization for each face
        !
        do ielem = 1,self%nelem
            do iface = 1,NFACES
                neighbor_status = NO_NEIGHBOR_FOUND

                !
                ! Check if face has neighbor on local partition
                !
                if ( self%faces(ielem,iface)%ftype == ORPHAN ) then
                    call self%find_neighbor_local(ielem,iface,  &
                                                  ineighbor_domain_g,        &
                                                  ineighbor_domain_l,        &
                                                  ineighbor_element_g,       &
                                                  ineighbor_element_l,       &
                                                  ineighbor_proc,            &
                                                  neighbor_status)


                    
                    !
                    ! If no neighbor found, either boundary condition face or chimera face
                    !
                    if ( neighbor_status == NEIGHBOR_FOUND ) then
                        ! Neighbor data should already be set, from previous routines. Set face type.
                        ftype = INTERIOR

                    else
                        ! Default ftype to ORPHAN face and clear neighbor index data.
                        ftype = ORPHAN      ! This should be processed later; either by a boundary condition(ftype=1), or a chimera boundary(ftype=2)
                        ineighbor_domain_g  = 0
                        ineighbor_domain_l  = 0
                        ineighbor_element_g = 0
                        ineighbor_element_l = 0
                        ineighbor_proc      = 0

                    end if


                    !
                    ! Call face neighbor initialization routine
                    !
                    call self%faces(ielem,iface)%init_neighbor(ftype,ineighbor_domain_g,ineighbor_domain_l,ineighbor_element_g,ineighbor_element_l,ineighbor_proc)

                end if


            end do !iface
        end do !ielem



    end subroutine init_comm_local
    !**************************************************************************************************************




















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine init_comm_global(self,ChiDG_COMM)
        class(mesh_t),                  intent(inout)   :: self
        type(mpi_comm),                 intent(in)      :: ChiDG_COMM


        integer(ik)             :: ixi,ieta,izeta,iface,ftype,ielem,ierr, ielem_neighbor
        integer(ik)             :: corner_one, corner_two, corner_three, corner_four
        integer(ik)             :: node_indices(4)
        logical                 :: boundary_face = .false.
        logical                 :: includes_node_one, includes_node_two, includes_node_three, includes_node_four
        logical                 :: neighbor_element, searching

        integer(ik)             :: ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l, ineighbor_proc, neighbor_status




        do ielem = 1,self%nelem
            do iface = 1,NFACES
                neighbor_status = NO_NEIGHBOR_FOUND

                !
                ! Check if face has neighbor on another MPI rank
                !
                if ( self%faces(ielem,iface)%ftype == ORPHAN ) then
                    ! send search request for neighbor face among global MPI ranks.
                    searching = .true.
                    call MPI_Bcast(searching,1,MPI_LOGICAL,IRANK,ChiDG_COMM,ierr)

                    call self%find_neighbor_global(ielem,iface,             &
                                                   ineighbor_domain_g,      &
                                                   ineighbor_domain_l,      &
                                                   ineighbor_element_g,     &
                                                   ineighbor_element_l,     &
                                                   ineighbor_proc,          &
                                                   neighbor_status,         &
                                                   ChiDG_COMM)
                            
                
                    !
                    ! If no neighbor found, either boundary condition face or chimera face
                    !
                    if ( neighbor_status == NEIGHBOR_FOUND ) then
                        ! Neighbor data should already be set, from previous routines. Set face type.
                        ftype = INTERIOR

                    else
                        ! Default ftype to ORPHAN face and clear neighbor index data.
                        ftype = ORPHAN      ! This should be processed later; either by a boundary condition(ftype=1), or a chimera boundary(ftype=2)
                        ineighbor_domain_g  = 0
                        ineighbor_domain_l  = 0
                        ineighbor_element_g = 0
                        ineighbor_element_l = 0
                        ineighbor_proc      = 0

                    end if


                    !
                    ! Call face neighbor initialization routine
                    !
                    call self%faces(ielem,iface)%init_neighbor(ftype,ineighbor_domain_g,ineighbor_domain_l,ineighbor_element_g,ineighbor_element_l,ineighbor_proc)

                end if


            end do !iface
        end do !ielem


        ! End search for global faces
        searching = .false.
        call MPI_Bcast(searching,1,MPI_LOGICAL,IRANK,ChiDG_COMM,ierr)





    end subroutine init_comm_global
    !**************************************************************************************************************















    !>  Outside of this subroutine, it should have already been determined that a neighbor request was initiated
    !!  from another processor and the current processor contains part of the domain of interest. This routine
    !!  receives corner indices from the requesting processor and tries to find a match in the current mesh.
    !!  The status of the element match is sent back. If a match was 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine handle_neighbor_request(self,iproc,ChiDG_COMM)
        class(mesh_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: iproc
        type(mpi_comm), intent(in)      :: ChiDG_COMM

        integer(ik) :: idomain_g, ielem_l, imesh
        integer(ik) :: ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l
        integer(ik) :: data(4), corner_indices(4)
        integer     :: ierr
        logical     :: includes_corner_one, includes_corner_two, includes_corner_three, includes_corner_four
        logical     :: searching, has_domain, neighbor_element



        ! Receive corner indices of face to be matched
        call MPI_Recv(corner_indices,4,MPI_INTEGER4,iproc,2,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)

        ! Loop through local domain and try to find a match
        ! Test the incoming face nodes against local elements, if all face nodes are also contained in an element, then they are neighbors.
        neighbor_element = .false.
        do ielem_l = 1,self%nelem
            includes_corner_one   = any( self%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(1) )
            includes_corner_two   = any( self%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(2) )
            includes_corner_three = any( self%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(3) )
            includes_corner_four  = any( self%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(4) )
            neighbor_element = ( includes_corner_one .and. includes_corner_two .and. includes_corner_three .and. includes_corner_four )

            if ( neighbor_element ) then
                ineighbor_domain_g  = self%elems(ielem_l)%connectivity%get_domain_index()
                ineighbor_domain_l  = self%elems(ielem_l)%idomain_l
                ineighbor_element_g = self%elems(ielem_l)%connectivity%get_element_index()
                ineighbor_element_l = self%elems(ielem_l)%ielement_l

                data = [ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l]
                exit
            end if

        end do ! ielem_l



        ! Send element-found status. If found, send element index information.
        call MPI_Send(neighbor_element,1,MPI_LOGICAL,iproc,3,ChiDG_COMM,ierr)

        if ( neighbor_element ) then
            call MPI_Send(data,4,MPI_INTEGER4,iproc,4,ChiDG_COMM,ierr)
        end if




    end subroutine handle_neighbor_request
    !**************************************************************************************************************




















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine find_neighbor_local(self,ielem_l,iface,ineighbor_domain_g,ineighbor_domain_l,ineighbor_element_g,ineighbor_element_l,ineighbor_proc,neighbor_status)
        class(mesh_t),                  intent(inout)   :: self
        integer(ik),                    intent(in)      :: ielem_l
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(inout)   :: ineighbor_domain_g
        integer(ik),                    intent(inout)   :: ineighbor_domain_l
        integer(ik),                    intent(inout)   :: ineighbor_element_g
        integer(ik),                    intent(inout)   :: ineighbor_element_l
        integer(ik),                    intent(inout)   :: ineighbor_proc
        integer(ik),                    intent(inout)   :: neighbor_status

        integer(ik) :: corner_one, corner_two, corner_three, corner_four
        integer(ik) :: corner_indices(4), ielem_neighbor, mapping
        logical     :: includes_corner_one, includes_corner_two, includes_corner_three, includes_corner_four
        logical     :: neighbor_element

        neighbor_status = NO_NEIGHBOR_FOUND

        ! Get the indices of the corner nodes that correspond to the current face in an element connectivity list
        mapping = self%elems(ielem_l)%connectivity%get_element_mapping()
        corner_one   = FACE_CORNERS(iface,1,mapping)
        corner_two   = FACE_CORNERS(iface,2,mapping)
        corner_three = FACE_CORNERS(iface,3,mapping)
        corner_four  = FACE_CORNERS(iface,4,mapping)


        ! For the current face, get the indices of the coordinate nodes for the corners
        corner_indices(1) = self%elems(ielem_l)%connectivity%get_element_node(corner_one)
        corner_indices(2) = self%elems(ielem_l)%connectivity%get_element_node(corner_two)
        corner_indices(3) = self%elems(ielem_l)%connectivity%get_element_node(corner_three)
        corner_indices(4) = self%elems(ielem_l)%connectivity%get_element_node(corner_four)

        
        ! Test the face nodes against other elements, if all face nodes are also contained in another element, then they are neighbors.
        neighbor_element = .false.
        do ielem_neighbor = 1,self%nelem
            if (ielem_neighbor /= ielem_l ) then
                includes_corner_one   = any( self%elems(ielem_neighbor)%connectivity%get_element_nodes() == corner_indices(1) )
                includes_corner_two   = any( self%elems(ielem_neighbor)%connectivity%get_element_nodes() == corner_indices(2) )
                includes_corner_three = any( self%elems(ielem_neighbor)%connectivity%get_element_nodes() == corner_indices(3) )
                includes_corner_four  = any( self%elems(ielem_neighbor)%connectivity%get_element_nodes() == corner_indices(4) )

                neighbor_element = ( includes_corner_one .and. includes_corner_two .and. includes_corner_three .and. includes_corner_four )

                if ( neighbor_element ) then
                    ineighbor_domain_g  = self%elems(ielem_neighbor)%connectivity%get_domain_index()
                    ineighbor_domain_l  = self%idomain_l
                    ineighbor_element_g = self%elems(ielem_neighbor)%connectivity%get_element_index()
                    ineighbor_element_l = ielem_neighbor
                    ineighbor_proc      = self%elems(ielem_neighbor)%connectivity%get_element_partition()
                    neighbor_status     = NEIGHBOR_FOUND
                    exit
                end if

            end if
        end do


    end subroutine find_neighbor_local
    !***************************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------
    subroutine find_neighbor_global(self,ielem_l,iface,ineighbor_domain_g,ineighbor_domain_l,ineighbor_element_g,ineighbor_element_l,ineighbor_proc,neighbor_status,ChiDG_COMM)
        class(mesh_t),                  intent(inout)   :: self
        integer(ik),                    intent(in)      :: ielem_l
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(inout)   :: ineighbor_domain_g
        integer(ik),                    intent(inout)   :: ineighbor_domain_l
        integer(ik),                    intent(inout)   :: ineighbor_element_g
        integer(ik),                    intent(inout)   :: ineighbor_element_l
        integer(ik),                    intent(inout)   :: ineighbor_proc
        integer(ik),                    intent(inout)   :: neighbor_status
        type(mpi_comm),                 intent(in)      :: ChiDG_COMM

        integer(ik) :: corner_one, corner_two, corner_three, corner_four
        integer(ik) :: corner_indices(4), data(4), ielem_neighbor, mapping, iproc, idomain_g, ierr
        logical     :: includes_corner_one, includes_corner_two, includes_corner_three, includes_corner_four
        logical     :: neighbor_element, has_domain

        neighbor_status = NO_NEIGHBOR_FOUND

        ! Get the indices of the corner nodes that correspond to the current face in an element connectivity list
        mapping      = self%elems(ielem_l)%connectivity%get_element_mapping()
        corner_one   = FACE_CORNERS(iface,1,mapping)
        corner_two   = FACE_CORNERS(iface,2,mapping)
        corner_three = FACE_CORNERS(iface,3,mapping)
        corner_four  = FACE_CORNERS(iface,4,mapping)


        ! For the current face, get the indices of the coordinate nodes for the corners
        corner_indices(1) = self%elems(ielem_l)%connectivity%get_element_node(corner_one)
        corner_indices(2) = self%elems(ielem_l)%connectivity%get_element_node(corner_two)
        corner_indices(3) = self%elems(ielem_l)%connectivity%get_element_node(corner_three)
        corner_indices(4) = self%elems(ielem_l)%connectivity%get_element_node(corner_four)

        
        ! Test the face nodes against other elements, if all face nodes are also contained in another element, then they are neighbors.
        neighbor_element = .false.




        do iproc = 0,NRANK-1
            if ( iproc /= IRANK ) then

                ! send global domain index of mesh being searched
                idomain_g = self%elems(ielem_l)%connectivity%get_domain_index()
                call MPI_Send(idomain_g,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,ierr)


                ! Check if other MPI rank has domain 
                call MPI_Recv(has_domain,1,MPI_LOGICAL,iproc,1,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)


                ! If so, send corner information
                if ( has_domain ) then
                    ! Send corner indices
                    call MPI_Send(corner_indices,4,MPI_INTEGER4,iproc,2,ChiDG_COMM,ierr)

                    ! Get status from proc on if it has neighbor element
                    call MPI_Recv(neighbor_element,1,MPI_LOGICAL,iproc,3,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)

                    if (neighbor_element) then
                        call MPI_Recv(data,4,MPI_INTEGER4,iproc,4,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        ineighbor_domain_g  = data(1)
                        ineighbor_domain_l  = data(2)
                        ineighbor_element_g = data(3)
                        ineighbor_element_l = data(4)
                        ineighbor_proc      = iproc
                        neighbor_status     = NEIGHBOR_FOUND
                    end if
                end if


            end if
        end do




    end subroutine find_neighbor_global
    !***************************************************************************************************************













    subroutine destructor(self)
        type(mesh_t), intent(inout) :: self

    
    end subroutine


end module type_mesh
