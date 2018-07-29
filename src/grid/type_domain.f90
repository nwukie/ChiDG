module type_domain
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                          ORPHAN, INTERIOR, BOUNDARY, CHIMERA, TWO_DIM, &
                                          THREE_DIM, NO_NEIGHBOR_FOUND, NEIGHBOR_FOUND, &
                                          NO_PROC, NFACES, ZERO, NO_PMM_ASSIGNED
    use mod_grid,                   only: FACE_CORNERS
    use mod_chidg_mpi,              only: IRANK, NRANK, GLOBAL_MASTER
    use mpi_f08

    use type_element,               only: element_t
    use type_face,                  only: face_t
    use type_ivector,               only: ivector_t
    use type_chimera,               only: chimera_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_element_connectivity,  only: element_connectivity_t
    implicit none
    private


    !>  Domain data type.
    !!
    !!  A domain_t contains arrays of elements and faces that define the geometry.
    !!  It also contains information about Chimera interfaces. 
    !!
    !!
    !!  For each element in a domain, there is an entry in domain%elems(:). In this
    !!  way, elements can be accessed by element index as:
    !!      domain%elems(ielem)
    !!
    !!  For each element, there are six faces. So, a face(iface) for a given element(ielem) 
    !!  can be accessed as:
    !!      domain%faces(ielem,iface)
    !!
    !!  Information on any Chimera interfaces for the domain are contains in domain%chimera
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/27/2016
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/4/2016
    !!  @note   restructure, create domain_t from previously mesh_t
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: domain_t

        character(:),   allocatable     :: name

        !
        ! Integer parameters
        !
        integer(ik)                     :: idomain_g
        integer(ik)                     :: idomain_l

        integer(ik)                     :: neqns       = 0     ! N-equations being solved
        integer(ik)                     :: nterms_s    = 0     ! N-terms in solution expansion
        integer(ik)                     :: nelements_g = 0     ! Number of elements in the global domain
        integer(ik)                     :: nelem       = 0     ! Number of total elements
        integer(ik)                     :: ntime       = 0     ! Number of time instances
        !integer(ik)                     :: eqn_ID      = NO_EQUATION_SET
        integer(ik)                     :: pmm_ID      = NO_PMM_ASSIGNED
        character(:),   allocatable     :: coordinate_system   ! 'Cartesian' or 'Cylindrical'

        
        !
        ! domain data
        !
        real(rk),           allocatable :: nodes(:,:)      ! Nodes of the reference domain.                 Proc-global. (nnodes, 3-coords)
        real(rk),           allocatable :: dnodes(:,:)     ! Node displacements: node_ale = nodes + dnodes. Proc-global. (nnodes, 3-coords)
        real(rk),           allocatable :: vnodes(:,:)     ! Node velocities:                               Proc-global. (nnodes, 3-coords)
        type(element_t),    allocatable :: elems(:)        ! Element storage (1:nelem)
        type(face_t),       allocatable :: faces(:,:)      ! Face storage (1:nelem,1:nfaces)
        

        ! chimera interfaces container
        type(chimera_t)                 :: chimera  


        !
        ! Initialization flags
        !
        logical   :: geomInitialized          = .false. ! Status of geometry initialization
        logical   :: solInitialized           = .false. ! Status of numerics initialization
        logical   :: local_comm_initialized   = .false. ! Status of processor-local comm init
        logical   :: global_comm_initialized  = .false. ! Status of processor-global comm init

    contains

        procedure           :: init_geom                ! geometry init for elements and faces 
        procedure           :: init_sol                 ! init data depending on solution order for elements and faces
        procedure           :: init_eqn                 ! initialize the equation set identifier on the mesh

        procedure, private  :: init_elems_geom          ! Loop through elements init geometry
        procedure, private  :: init_elems_sol           ! Loop through elements init data depending on the solution order
        procedure, private  :: init_faces_geom          ! Loop through faces init geometry
        procedure, private  :: init_faces_sol           ! Loop through faces init data depending on the solution order

        procedure           :: init_comm_local          ! For faces, find proc-local neighbors, initialize face neighbor indices 
        procedure           :: init_comm_global         ! For faces, find neighbors across procs, initialize face neighbor indices

        ! ALE
        procedure, public   :: set_displacements_velocities
        procedure           :: update_interpolations_ale

        ! Utilities
        procedure, private  :: find_neighbor_local      ! Try to find a neighbor for a particular face on the local processor
        procedure, private  :: find_neighbor_global     ! Try to find a neighbor for a particular face across processors
        procedure           :: handle_neighbor_request  ! When a neighbor request from another processor comes in, 
                                                        ! check if current processor contains neighbor


        procedure,  public  :: get_recv_procs           ! Return proc ranks receiving from (neighbor+chimera)
        procedure,  public  :: get_recv_procs_local     ! Return proc ranks receiving neighbor data from
        procedure,  public  :: get_recv_procs_chimera   ! Return proc ranks receiving chimera data from 

        procedure,  public  :: get_send_procs           ! Return proc ranks sending to (neighbor+chimera)
        procedure,  public  :: get_send_procs_local     ! Return proc ranks sending neighbor data to
        procedure,  public  :: get_send_procs_chimera   ! Return proc ranks sending chimera data to

        procedure,  public  :: get_nelements_global             ! Return number of elements in the global domain
        procedure,  public  :: get_nelements_local              ! Return number of elements in the processor-local domain
        procedure,  public  :: nelements => get_nelements_local ! Included for framework consistency

        final               :: destructor

    end type domain_t
    !*****************************************************************************************





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
    !!  @param[in]  idomain_l       Proc-local domain index.
    !!  @param[in]  nelements_g     Proc-global number of elements in the domain.
    !!  @param[in]  nodes           Proc-global node list.                  (nnodes, 3-coords)
    !!  @param[in]  dnodes          Proc-global node coordinate delta list. (nnodes, 3-coords)
    !!  @param[in]  connectivity    Proc-local connectivities.
    !!  @param[in]  coord_system    Coordinate system of the nodal coordinates.
    !!  
    !!  TODO: test dnodes initialization
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_geom(self,idomain_l,nelements_g,nodes,connectivity,coord_system)
        class(domain_t),                intent(inout)   :: self
        integer(ik),                    intent(in)      :: idomain_l
        integer(ik),                    intent(in)      :: nelements_g
        real(rk),                       intent(in)      :: nodes(:,:)
        type(domain_connectivity_t),    intent(in)      :: connectivity
        character(*),                   intent(in)      :: coord_system

        integer(ik) :: inode, msg1, msg2, funit
        real(rk)    :: R1_velocity_1, R1_velocity_2, R1_velocity_3, &
                       R2_velocity_1, R2_velocity_2, R2_velocity_3
        integer(ik) :: R1_domain_min, R1_domain_max, &
                       R2_domain_min, R2_domain_max
        logical     :: file_exists

        namelist /region_one/ R1_velocity_1, R1_velocity_2, R1_velocity_3, R1_domain_min, R1_domain_max
        namelist /region_two/ R2_velocity_1, R2_velocity_2, R2_velocity_3, R2_domain_min, R2_domain_max

        !
        ! Store number of terms in coordinate expansion and domain index
        !
        self%idomain_g    = connectivity%get_domain_index()
        self%idomain_l    = idomain_l
        self%nelements_g  = nelements_g


        !
        ! Initialize nodes:
        !   Reference nodes = nodes
        !   Default node coordinate deltas = zero
        !
        self%nodes  = nodes
        self%dnodes = nodes
        self%vnodes = nodes
        self%dnodes = ZERO
        self%vnodes = ZERO


        !
        ! Call geometry initialization for elements and faces
        !
        call self%init_elems_geom(nodes,connectivity,coord_system)
        call self%init_faces_geom()


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
            if (msg1 == 0 .and. self%idomain_g >= R1_domain_min .and. self%idomain_g <= R1_domain_max) then
                print*, 'setting velocity: ', R1_velocity_1, R1_velocity_2, R1_velocity_3, ' on domain: ', self%idomain_g
                do inode = 1,size(self%vnodes,1)
                    select case(trim(coord_system))
                        case('Cartesian')
                            self%vnodes(inode,1:3) = [R1_velocity_1,R1_velocity_2,R1_velocity_3]
                        case('Cylindrical')
                            self%vnodes(inode,1:3) = [R1_velocity_1,self%nodes(inode,1)*R1_velocity_2,R1_velocity_3]
                        case default
                            call chidg_signal_one(FATAL,"domain%init_geom: Invalid coordinate system.",trim(coord_system))
                    end select
                end do
                call self%set_displacements_velocities(self%dnodes,self%vnodes)
            end if

            ! Region 2
            if (msg2 == 0 .and. self%idomain_g >= R2_domain_min .and. self%idomain_g <= R2_domain_max) then
                print*, 'setting velocity: ', R2_velocity_1, R2_velocity_2, R2_velocity_3, ' on domain: ', self%idomain_g
                do inode = 1,size(self%vnodes,1)
                    select case(trim(coord_system))
                        case('Cartesian')
                            self%vnodes(inode,1:3) = [R2_velocity_1,R2_velocity_2,R2_velocity_3]
                        case('Cylindrical')
                            self%vnodes(inode,1:3) = [R2_velocity_1,self%nodes(inode,1)*R2_velocity_2,R2_velocity_3]
                        case default
                            call chidg_signal_one(FATAL,"domain%init_geom: Invalid coordinate system.",trim(coord_system))
                    end select
                end do
                call self%set_displacements_velocities(self%dnodes,self%vnodes)
            end if

        end if




        !
        ! Set coordinate system and confirm initialization 
        !
        self%coordinate_system = coord_system
        self%geomInitialized = .true.


    end subroutine init_geom
    !*****************************************************************************************








    !>  Initialize ALE data from node displacement data. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !!
    !!  TODO: Test
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_displacements_velocities(self,dnodes,vnodes)
        class(domain_t),        intent(inout)   :: self
        real(rk),               intent(in)      :: dnodes(:,:)
        real(rk),               intent(in)      :: vnodes(:,:)

        integer(ik) :: ielem, iface

        do ielem = 1,self%nelem
            call self%elems(ielem)%set_displacements_velocities(dnodes,vnodes)
            do iface = 1,NFACES
                call self%faces(ielem,iface)%set_displacements_velocities(self%elems(ielem))
            end do !iface
        end do !ielem


    end subroutine set_displacements_velocities
    !*****************************************************************************************


    !>  Initialize ALE data from node displacement data. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2017
    !!
    !!
    !!  TODO: Test
    !!
    !----------------------------------------------------------------------------------------
    subroutine update_interpolations_ale(self)
        class(domain_t),        intent(inout)   :: self

        integer(ik) :: ielem, iface

        do ielem = 1,self%nelem
            call self%elems(ielem)%update_interpolations_ale()
            do iface = 1,NFACES
                call self%faces(ielem,iface)%update_interpolations_ale(self%elems(ielem))
            end do !iface
        end do !ielem


    end subroutine update_interpolations_ale
    !*****************************************************************************************







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
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/9/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_sol(self,interpolation,level,nterms_s,neqns,ntime)
        class(domain_t),    intent(inout)   :: self
        character(*),       intent(in)      :: interpolation
        integer(ik),        intent(in)      :: level
        integer(ik),        intent(in)      :: nterms_s
        integer(ik),        intent(in)      :: neqns
        integer(ik),        intent(in)      :: ntime

        !
        ! Store number of equations and number of terms in solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s
        self%ntime    = ntime

        !
        ! Call numerics initialization for elements and faces
        !
        call self%init_elems_sol(interpolation,level,nterms_s,neqns,ntime)
        call self%init_faces_sol()               

        call self%update_interpolations_ale()
        !
        ! Confirm initialization
        !
        self%solInitialized = .true.

    end subroutine init_sol
    !*****************************************************************************************








    !>  Initialize the equation set identifier on the mesh.
    !!
    !!  Sets the equation set identifier self%eqn_ID that can be used to acces
    !!  the equation_set_t object on chidg_data as chidg_data%eqnset(eqn_ID)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_eqn(self,eqn_ID)
        class(domain_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: eqn_ID

        integer(ik) :: ielem

        !
        ! Store number of equations and number of terms in solution expansion
        !
        !self%eqn_ID = eqn_ID

        
        !
        ! Assign all elements in the domain to the equation set identifier.
        !
        do ielem = 1,self%nelements()
            call self%elems(ielem)%init_eqn(eqn_ID)
        end do

    end subroutine init_eqn
    !*****************************************************************************************










    !>  Mesh - element initialization procedure
    !!
    !!  Computes the number of elements based on the element mapping selected and
    !!  calls the element initialization procedure on individual elements.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  points_g    Rank-3 matrix of coordinate points defining a block mesh
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_elems_geom(self,nodes,domain_connectivity,coord_system)
        class(domain_t),                intent(inout)   :: self
        real(rk),                       intent(in)      :: nodes(:,:)
        type(domain_connectivity_t),    intent(in)      :: domain_connectivity
        character(*),                   intent(in)      :: coord_system


        type(element_connectivity_t)    :: element_connectivity
        integer(ik)                     :: ierr, nelem, location(5), etype,             &
                                           idomain_g, ielement_g, idomain_l, ielement_l
        integer(ik),    allocatable     :: connectivity(:)


        !
        ! Store total number of elements
        !
        nelem      = domain_connectivity%get_nelements()
        self%nelem = nelem


        !
        ! Allocate element storage
        !
        allocate(self%elems(nelem), stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"mesh%init_elems_geom: Memory allocation error: init_elements")


        !
        ! Call geometry initialization for each element
        !
        idomain_l = self%idomain_l
        do ielement_l = 1,nelem

            element_connectivity = domain_connectivity%get_element_connectivity(ielement_l)
            connectivity = element_connectivity%get_element_nodes()
            idomain_g    = element_connectivity%get_domain_index()
            ielement_g   = element_connectivity%get_element_index()
            !location     = [idomain_g, idomain_l, ielement_g, ielement_l]
            location     = [idomain_g, idomain_l, ielement_g, ielement_l, IRANK]
            etype        = element_connectivity%get_element_mapping()

            call self%elems(ielement_l)%init_geom(nodes,connectivity,etype,location,coord_system)

        end do ! ielem


    end subroutine init_elems_geom
    !*****************************************************************************************








    !>  Mesh - element solution data initialization
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  neqns       Number of equations in the domain equation set
    !!  @param[in]  nterms_s    Number of terms in the solution expansion
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   11/5/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_elems_sol(self,interpolation,level,nterms_s,neqns,ntime)
        class(domain_t),    intent(inout)   :: self
        character(*),       intent(in)      :: interpolation
        integer(ik),        intent(in)      :: level
        integer(ik),        intent(in)      :: nterms_s
        integer(ik),        intent(in)      :: neqns
        integer(ik),        intent(in)      :: ntime
        integer(ik)                         :: ielem


        !
        ! Store number of equations and number of terms in the solution expansion
        !
        self%neqns    = neqns
        self%nterms_s = nterms_s
        self%ntime    = ntime

        !
        ! Call the numerics initialization procedure for each element
        !
        do ielem = 1,self%nelem
            call self%elems(ielem)%init_sol(interpolation,level,self%nterms_s,self%neqns,ntime) 
        end do


    end subroutine init_elems_sol
    !*****************************************************************************************










    !>  Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_faces_geom(self)
        class(domain_t),                intent(inout)   :: self

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

                call self%faces(ielem,iface)%init_geom(iface,self%elems(ielem))

            end do !iface
        end do !ielem


    end subroutine init_faces_geom
    !*****************************************************************************************









    !>  Mesh - face initialization procedure
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_faces_sol(self)
        class(domain_t), intent(inout)  :: self

        integer(ik) :: ielem, iface

        !
        ! Loop through elements, faces and call initialization that depends on 
        ! the solution basis.
        !
        do ielem = 1,self%nelem
            do iface = 1,NFACES

                call self%faces(ielem,iface)%init_sol(self%elems(ielem))

            end do ! iface
        end do ! ielem


    end subroutine init_faces_sol
    !******************************************************************************************















    !>  Initialize processor-local, interior neighbor communication.
    !!
    !!  For each face without an interior neighbor, search the current mesh for a 
    !!  potential neighbor element/face by trying to match the corner indices of the elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_comm_local(self)
        class(domain_t),                  intent(inout)   :: self


        integer(ik)             :: iface,ftype,ielem,ierr, ielem_neighbor,              &
                                   corner_one, corner_two, corner_three, corner_four,   &
                                   node_indices(4),                                     &
                                   ineighbor_domain_g,  ineighbor_domain_l,             &
                                   ineighbor_element_g, ineighbor_element_l,            &
                                   ineighbor_face,      ineighbor_proc,                 &
                                   ineighbor_neqns,     ineighbor_nterms_s,             &
                                   neighbor_status, idomain_g, idomain_l, ielement_g,   &
                                   ielement_l, nterms_s, neqns

        logical                 :: boundary_face = .false.
        logical                 :: includes_node_one, includes_node_two,    &
                                   includes_node_three, includes_node_four, &
                                   neighbor_element

        !
        ! Loop through each local element and call initialization for each face
        !
        do ielem = 1,self%nelem

            do iface = 1,NFACES

                !
                ! Check if face has neighbor on local partition.
                !   - ORPHAN means the exterior state is empty and we want to try and find a connection
                !
                if ( self%faces(ielem,iface)%ftype == ORPHAN ) then
                    call self%find_neighbor_local(ielem,iface,              &
                                                  ineighbor_domain_g,       &
                                                  ineighbor_domain_l,       &
                                                  ineighbor_element_g,      &
                                                  ineighbor_element_l,      &
                                                  ineighbor_face,           &
                                                  ineighbor_proc,           &
                                                  neighbor_status)

                !
                !   - INTERIOR means the neighbor is already connected, but maybe we want to reinitialize
                !     the info, for example nterms_s if the order has been increased. So here,
                !     we just access the location that is already initialized.
                !
                else if ( self%faces(ielem,iface)%ftype == INTERIOR )  then
                    
                    ineighbor_domain_g  = self%faces(ielem,iface)%ineighbor_domain_g
                    ineighbor_domain_l  = self%faces(ielem,iface)%ineighbor_domain_l
                    ineighbor_element_g = self%faces(ielem,iface)%ineighbor_element_g
                    ineighbor_element_l = self%faces(ielem,iface)%ineighbor_element_l
                    ineighbor_face      = self%faces(ielem,iface)%ineighbor_face
                    ineighbor_proc      = self%faces(ielem,iface)%ineighbor_proc
                    neighbor_status     = NEIGHBOR_FOUND

                end if


                    
                !
                ! If no neighbor found, either boundary condition face or chimera face
                !
                if ( (self%faces(ielem,iface)%ftype == ORPHAN) .or. &
                     (self%faces(ielem,iface)%ftype == INTERIOR) ) then

                    if ( neighbor_status == NEIGHBOR_FOUND ) then

                        ftype               = INTERIOR
                        ineighbor_neqns     = self%elems(ineighbor_element_l)%neqns
                        ineighbor_nterms_s  = self%elems(ineighbor_element_l)%nterms_s


                    else
                        ! Default ftype to ORPHAN face and clear neighbor index data.
                        ! ftype should be processed later; either by a boundary conditions (ftype=1), 
                        ! or a chimera boundary (ftype = 2)
                        ftype = ORPHAN
                        ineighbor_domain_g  = 0
                        ineighbor_domain_l  = 0
                        ineighbor_element_g = 0
                        ineighbor_element_l = 0
                        ineighbor_face      = 0
                        ineighbor_neqns     = 0
                        ineighbor_nterms_s  = 0
                        ineighbor_proc      = NO_PROC

                    end if


                    !
                    ! Call face neighbor initialization routine
                    !
                    call self%faces(ielem,iface)%set_neighbor(ftype,ineighbor_domain_g, ineighbor_domain_l,    &
                                                                    ineighbor_element_g,ineighbor_element_l,   &
                                                                    ineighbor_face,     ineighbor_neqns,       &
                                                                    ineighbor_nterms_s, ineighbor_proc)

                    !
                    ! Also, initialize neighbor face at the same time so we don't
                    ! have to do the search again. 
                    !
                    ! Only can initialize opposite neighbor if opposite element is on-proc.
                    !
                    if ( (neighbor_status == NEIGHBOR_FOUND) .and. (ineighbor_proc == IRANK) ) then
                        idomain_g  = self%elems(ielem)%idomain_g
                        idomain_l  = self%elems(ielem)%idomain_l
                        ielement_g = self%elems(ielem)%ielement_g
                        ielement_l = self%elems(ielem)%ielement_l
                        neqns      = self%elems(ielem)%neqns
                        nterms_s   = self%elems(ielem)%nterms_s
                        call self%faces(ineighbor_element_l,ineighbor_face)%set_neighbor(ftype,idomain_g,  idomain_l,  &
                                                                                               ielement_g, ielement_l, &
                                                                                               iface,      neqns,      &
                                                                                               nterms_s,   IRANK)
                    end if

                end if

            end do !iface
        end do !ielem


        ! Set initialized
        self%local_comm_initialized = .true.

    end subroutine init_comm_local
    !*****************************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/17/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_comm_global(self,ChiDG_COMM)
        class(domain_t),                  intent(inout)   :: self
        type(mpi_comm),                 intent(in)      :: ChiDG_COMM


        integer(ik)             :: iface,ftype,ielem,ierr, ielem_neighbor,              &
                                   corner_one, corner_two, corner_three, corner_four,   &
                                   node_indices(4),                                     &
                                   ineighbor_domain_g,  ineighbor_domain_l,             &
                                   ineighbor_element_g, ineighbor_element_l,            &
                                   ineighbor_face,      ineighbor_proc,                 &
                                   ineighbor_neqns,     ineighbor_nterms_s, neighbor_status

        logical                 :: includes_node_one, includes_node_two,    &
                                   includes_node_three, includes_node_four, &
                                   neighbor_element, searching


        real(rk)                                :: neighbor_h(3)
        real(rk), allocatable, dimension(:,:)   :: neighbor_grad1,   neighbor_grad2,    &
                                                   neighbor_grad3,   neighbor_br2_face, &
                                                   neighbor_br2_vol, neighbor_invmass



        do ielem = 1,self%nelem
            do iface = 1,NFACES

                !
                ! Check if face has neighbor on another MPI rank.
                !
                !   Do this for ORPHAN faces, that are looking for a potential neighbor
                !   Do this also for INTERIOR faces with off-processor neighbors, in case 
                !   this is being called as a reinitialization routine, so that 
                !   element-specific information gets updated, such as neighbor_grad1, 
                !   etc. because these could have changed if the order of the solution changed
                !
                if ( (self%faces(ielem,iface)%ftype == ORPHAN) .or.         &
                     ( (self%faces(ielem,iface)%ftype == INTERIOR) .and.    &
                       (self%faces(ielem,iface)%ineighbor_proc /= IRANK) )  &
                   ) then

                    ! send search request for neighbor face among global MPI ranks.
                    searching = .true.
                    call MPI_Bcast(searching,1,MPI_LOGICAL,IRANK,ChiDG_COMM,ierr)

                    call self%find_neighbor_global(ielem,iface,             &
                                                   ineighbor_domain_g,      &
                                                   ineighbor_domain_l,      &
                                                   ineighbor_element_g,     &
                                                   ineighbor_element_l,     &
                                                   ineighbor_face,          &
                                                   ineighbor_neqns,         &
                                                   ineighbor_nterms_s,      &
                                                   ineighbor_proc,          &
                                                   neighbor_grad1,          &
                                                   neighbor_grad2,          &
                                                   neighbor_grad3,          &
                                                   neighbor_br2_face,       &
                                                   neighbor_br2_vol,        &
                                                   neighbor_invmass,        &
                                                   neighbor_h,              &
                                                   neighbor_status,         &
                                                   ChiDG_COMM)
                            
                
                    !
                    ! If no neighbor found, either boundary condition face or chimera face
                    !
                    if ( neighbor_status == NEIGHBOR_FOUND ) then
                        ! Neighbor data should already be set, from previous routines. 
                        ! Set face type.
                        ftype = INTERIOR

                        !
                        ! Set neighbor data
                        !
                        self%faces(ielem,iface)%neighbor_h        = neighbor_h
                        self%faces(ielem,iface)%neighbor_grad1    = neighbor_grad1
                        self%faces(ielem,iface)%neighbor_grad2    = neighbor_grad2
                        self%faces(ielem,iface)%neighbor_grad3    = neighbor_grad3
                        self%faces(ielem,iface)%neighbor_br2_face = neighbor_br2_face
                        self%faces(ielem,iface)%neighbor_br2_vol  = neighbor_br2_vol
                        self%faces(ielem,iface)%neighbor_invmass  = neighbor_invmass

                    else
                        ! Default ftype to ORPHAN face and clear neighbor index data.
                        ! ftype should be processed later; either by a boundary 
                        ! condition(ftype=1), or a chimera boundary(ftype=2)
                        ! 
                        ftype = ORPHAN
                        ineighbor_domain_g  = 0
                        ineighbor_domain_l  = 0
                        ineighbor_element_g = 0
                        ineighbor_element_l = 0
                        ineighbor_face      = 0
                        ineighbor_neqns     = 0
                        ineighbor_nterms_s  = 0
                        ineighbor_proc      = NO_PROC

                    end if


                    !
                    ! Call face neighbor initialization routine
                    !
                    call self%faces(ielem,iface)%set_neighbor(ftype,ineighbor_domain_g,  ineighbor_domain_l,   &
                                                                    ineighbor_element_g, ineighbor_element_l,  &
                                                                    ineighbor_face,      ineighbor_neqns,      &
                                                                    ineighbor_nterms_s,  ineighbor_proc)


                end if


            end do !iface
        end do !ielem


        ! End search for global faces
        searching = .false.
        call MPI_Bcast(searching,1,MPI_LOGICAL,IRANK,ChiDG_COMM,ierr)


        ! Set initialized
        self%global_comm_initialized = .true.


    end subroutine init_comm_global
    !*****************************************************************************************








    !>  Outside of this subroutine, it should have already been determined that a 
    !!  neighbor request was initiated from another processor and the current processor 
    !!  contains part of the domain of interest. This routine receives corner indices from 
    !!  the requesting processor and tries to find a match in the current mesh. The status 
    !!  of the element match is sent back. If a match was 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine handle_neighbor_request(self,iproc,ChiDG_COMM)
        class(domain_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: iproc
        type(mpi_comm), intent(in)      :: ChiDG_COMM

        integer(ik) :: ielem_l, iface, ierr,                                &
                       ineighbor_domain_g, ineighbor_domain_l,              &
                       ineighbor_element_g, ineighbor_element_l,            &
                       ineighbor_face, ineighbor_neqns, ineighbor_nterms_s, &
                       data(7), corner_indices(4), grad_size(2),            &
                       invmass_size(2), br2_face_size(2), br2_vol_size(2), size_data(8)
        logical     :: includes_corner_one, includes_corner_two, &
                       includes_corner_three, includes_corner_four, neighbor_element


        ! Receive corner indices of face to be matched
        call MPI_Recv(corner_indices,4,MPI_INTEGER4,iproc,2,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)


        ! Loop through local domain and try to find a match
        ! Test the incoming face nodes against local elements, if all face nodes are 
        ! also contained in an element, then they are neighbors.
        neighbor_element = .false.
        do ielem_l = 1,self%nelem
            includes_corner_one   = any( self%elems(ielem_l)%connectivity == corner_indices(1) )
            includes_corner_two   = any( self%elems(ielem_l)%connectivity == corner_indices(2) )
            includes_corner_three = any( self%elems(ielem_l)%connectivity == corner_indices(3) )
            includes_corner_four  = any( self%elems(ielem_l)%connectivity == corner_indices(4) )
            neighbor_element = ( includes_corner_one   .and. &
                                 includes_corner_two   .and. &
                                 includes_corner_three .and. &
                                 includes_corner_four )

            if ( neighbor_element ) then
                !
                ! Get indices for neighbor element
                !
                ineighbor_domain_g  = self%elems(ielem_l)%idomain_g
                ineighbor_domain_l  = self%elems(ielem_l)%idomain_l
                ineighbor_element_g = self%elems(ielem_l)%ielement_g
                ineighbor_element_l = self%elems(ielem_l)%ielement_l
                ineighbor_neqns     = self%elems(ielem_l)%neqns
                ineighbor_nterms_s  = self%elems(ielem_l)%nterms_s

                
                !
                ! Get face index connected to the requesting element
                !
                iface = self%elems(ielem_l)%get_face_from_corners(corner_indices)
                ineighbor_face = iface

                data = [ineighbor_domain_g,  ineighbor_domain_l,    &
                        ineighbor_element_g, ineighbor_element_l,   &
                        ineighbor_face,      ineighbor_neqns,       &
                        ineighbor_nterms_s]

                exit
            end if

        end do ! ielem_l



        ! Send element-found status. If found, send element index information.
        call MPI_Send(neighbor_element,1,MPI_LOGICAL,iproc,3,ChiDG_COMM,ierr)

        if ( neighbor_element ) then
            !
            ! Send Indices
            !
            call MPI_Send(data,7,MPI_INTEGER4,iproc,4,ChiDG_COMM,ierr)

            !
            ! Send Element Data
            !
            !grad_size(1)     = size(self%faces(ielem_l,iface)%grad1,1)
            !grad_size(2)     = size(self%faces(ielem_l,iface)%grad1,2)
            !br2_face_size(1) = size(self%faces(ielem_l,iface)%br2_face,1)
            !br2_face_size(2) = size(self%faces(ielem_l,iface)%br2_face,2)
            !br2_vol_size(1)  = size(self%faces(ielem_l,iface)%br2_vol,1)
            !br2_vol_size(2)  = size(self%faces(ielem_l,iface)%br2_vol,2)
            !invmass_size(1)  = size(self%elems(ielem_l)%invmass,1)
            !invmass_size(2)  = size(self%elems(ielem_l)%invmass,2)
            size_data(1) = size(self%faces(ielem_l,iface)%grad1,1)
            size_data(2) = size(self%faces(ielem_l,iface)%grad1,2)
            size_data(3) = size(self%faces(ielem_l,iface)%br2_face,1)
            size_data(4) = size(self%faces(ielem_l,iface)%br2_face,2)
            size_data(5) = size(self%faces(ielem_l,iface)%br2_vol,1)
            size_data(6) = size(self%faces(ielem_l,iface)%br2_vol,2)
            size_data(7) = size(self%elems(ielem_l)%invmass,1)
            size_data(8) = size(self%elems(ielem_l)%invmass,2)

            call MPI_Send(size_data,8,MPI_INTEGER4,iproc,5,ChiDG_COMM,ierr)
            !call MPI_Send(grad_size,    2,MPI_INTEGER4,iproc,5,ChiDG_COMM,ierr)
            !call MPI_Send(br2_face_size,2,MPI_INTEGER4,iproc,6,ChiDG_COMM,ierr)
            !call MPI_Send(br2_vol_size, 2,MPI_INTEGER4,iproc,7,ChiDG_COMM,ierr)
            !call MPI_Send(invmass_size, 2,MPI_INTEGER4,iproc,8,ChiDG_COMM,ierr)

            !call MPI_Send(self%faces(ielem_l,iface)%grad1,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc, 9,ChiDG_COMM,ierr)
            !call MPI_Send(self%faces(ielem_l,iface)%grad2,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc,10,ChiDG_COMM,ierr)
            !call MPI_Send(self%faces(ielem_l,iface)%grad3,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc,11,ChiDG_COMM,ierr)
            !call MPI_Send(self%faces(ielem_l,iface)%br2_face,   br2_face_size(1)*br2_face_size(2),  MPI_REAL8,iproc,12,ChiDG_COMM,ierr)
            !call MPI_Send(self%faces(ielem_l,iface)%br2_vol,    br2_vol_size(1)*br2_vol_size(2),    MPI_REAL8,iproc,13,ChiDG_COMM,ierr)
            !call MPI_Send(self%elems(ielem_l)%invmass,          invmass_size(1)*invmass_size(2),    MPI_REAL8,iproc,14,ChiDG_COMM,ierr)
            !call MPI_Send(self%elems(ielem_l)%h,                3,                                  MPI_REAL8,iproc,15,ChiDG_COMM,ierr)

            call MPI_Send(self%faces(ielem_l,iface)%grad1,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc, 6,ChiDG_COMM,ierr)
            call MPI_Send(self%faces(ielem_l,iface)%grad2,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc, 7,ChiDG_COMM,ierr)
            call MPI_Send(self%faces(ielem_l,iface)%grad3,      grad_size(1)*grad_size(2),          MPI_REAL8,iproc, 8,ChiDG_COMM,ierr)
            call MPI_Send(self%faces(ielem_l,iface)%br2_face,   br2_face_size(1)*br2_face_size(2),  MPI_REAL8,iproc, 9,ChiDG_COMM,ierr)
            call MPI_Send(self%faces(ielem_l,iface)%br2_vol,    br2_vol_size(1)*br2_vol_size(2),    MPI_REAL8,iproc,10,ChiDG_COMM,ierr)
            call MPI_Send(self%elems(ielem_l)%invmass,          invmass_size(1)*invmass_size(2),    MPI_REAL8,iproc,11,ChiDG_COMM,ierr)
            call MPI_Send(self%elems(ielem_l)%h,                3,                                  MPI_REAL8,iproc,12,ChiDG_COMM,ierr)
        end if


    end subroutine handle_neighbor_request
    !*****************************************************************************************










    !>  For given element/face indices, try to find a potential interior neighbor. That is, 
    !!  a matching element within the current domain and on the current processor(local).
    !!
    !!  If found, return neighbor info.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine find_neighbor_local(self,ielem_l,iface,ineighbor_domain_g, ineighbor_domain_l,   &
                                                      ineighbor_element_g,ineighbor_element_l,  &
                                                      ineighbor_face,     ineighbor_proc,       &
                                                      neighbor_status)
        class(domain_t),                  intent(inout)   :: self
        integer(ik),                    intent(in)      :: ielem_l
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(inout)   :: ineighbor_domain_g
        integer(ik),                    intent(inout)   :: ineighbor_domain_l
        integer(ik),                    intent(inout)   :: ineighbor_element_g
        integer(ik),                    intent(inout)   :: ineighbor_element_l
        integer(ik),                    intent(inout)   :: ineighbor_face
        integer(ik),                    intent(inout)   :: ineighbor_proc
        integer(ik),                    intent(inout)   :: neighbor_status

        integer(ik),    allocatable :: element_nodes(:)
        integer(ik) :: corner_one, corner_two, corner_three, corner_four,   &
                       corner_indices(4), ielem_neighbor, mapping
        logical     :: includes_corner_one, includes_corner_two, &
                       includes_corner_three, includes_corner_four, neighbor_element

        neighbor_status = NO_NEIGHBOR_FOUND

        !
        ! Get the element-local node indices of the corner nodes that correspond 
        ! to the current face in an element connectivity list
        !
        !mapping = self%elems(ielem_l)%connectivity%get_element_mapping()
        mapping = self%elems(ielem_l)%element_type
        corner_one   = FACE_CORNERS(iface,1,mapping)
        corner_two   = FACE_CORNERS(iface,2,mapping)
        corner_three = FACE_CORNERS(iface,3,mapping)
        corner_four  = FACE_CORNERS(iface,4,mapping)

        
        !
        ! For the current face, get the global-indices of the coordinate nodes 
        ! for the corners
        !
        !corner_indices(1) = self%elems(ielem_l)%connectivity%get_element_node(corner_one)
        !corner_indices(2) = self%elems(ielem_l)%connectivity%get_element_node(corner_two)
        !corner_indices(3) = self%elems(ielem_l)%connectivity%get_element_node(corner_three)
        !corner_indices(4) = self%elems(ielem_l)%connectivity%get_element_node(corner_four)
        corner_indices(1) = self%elems(ielem_l)%connectivity(corner_one)
        corner_indices(2) = self%elems(ielem_l)%connectivity(corner_two)
        corner_indices(3) = self%elems(ielem_l)%connectivity(corner_three)
        corner_indices(4) = self%elems(ielem_l)%connectivity(corner_four)

        
        !
        ! Test the global face node indices against other elements. If all face nodes 
        ! are also contained in another element, then they are neighbors.
        !
        neighbor_element = .false.
        do ielem_neighbor = 1,self%nelem
            if (ielem_neighbor /= ielem_l) then

                !element_nodes = self%elems(ielem_neighbor)%connectivity%get_element_nodes()
                element_nodes = self%elems(ielem_neighbor)%connectivity
                includes_corner_one   = any( element_nodes == corner_indices(1) )
                includes_corner_two   = any( element_nodes == corner_indices(2) )
                includes_corner_three = any( element_nodes == corner_indices(3) )
                includes_corner_four  = any( element_nodes == corner_indices(4) )

                neighbor_element = ( includes_corner_one   .and. &
                                     includes_corner_two   .and. &
                                     includes_corner_three .and. &
                                     includes_corner_four )

                if ( neighbor_element ) then
                    ineighbor_domain_g  = self%elems(ielem_neighbor)%idomain_g
                    ineighbor_domain_l  = self%elems(ielem_neighbor)%idomain_l
                    ineighbor_element_g = self%elems(ielem_neighbor)%ielement_g
                    ineighbor_element_l = self%elems(ielem_neighbor)%ielement_l
                    ineighbor_face      = self%elems(ielem_neighbor)%get_face_from_corners(corner_indices)
                    !ineighbor_proc      = self%elems(ielem_neighbor)%connectivity%get_element_partition()
                    ineighbor_proc      = IRANK
                    neighbor_status     = NEIGHBOR_FOUND
                    exit
                end if

            end if
        end do


    end subroutine find_neighbor_local
    !*****************************************************************************************








    !>  Search for an interior neighbor element across all processors.
    !!
    !!  Pass the corner indices for matching on a potential neighbor element so they can 
    !!  be checked by elements on another processor. If element is found, get the grad1, 
    !!  grad2, grad3, invmass etc. information that is element specific to that neighbor, 
    !!  which is located on another processor. This allows us to compute cartesian 
    !!  derivatives as they would be computed in the neighbor element, without sending 
    !!  the entire element representation across.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine find_neighbor_global(self,ielem_l,iface,ineighbor_domain_g,  ineighbor_domain_l,  &
                                                       ineighbor_element_g, ineighbor_element_l, &
                                                       ineighbor_face,      ineighbor_neqns,     &
                                                       ineighbor_nterms_s,  ineighbor_proc,      &
                                                       neighbor_grad1, neighbor_grad2, neighbor_grad3, &
                                                       neighbor_br2_face, neighbor_br2_vol,      &
                                                       neighbor_invmass,    &
                                                       neighbor_h,          &
                                                       neighbor_status,     &
                                                       ChiDG_COMM)
        class(domain_t),                intent(inout)   :: self
        integer(ik),                    intent(in)      :: ielem_l
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(inout)   :: ineighbor_domain_g
        integer(ik),                    intent(inout)   :: ineighbor_domain_l
        integer(ik),                    intent(inout)   :: ineighbor_element_g
        integer(ik),                    intent(inout)   :: ineighbor_element_l
        integer(ik),                    intent(inout)   :: ineighbor_face
        integer(ik),                    intent(inout)   :: ineighbor_neqns
        integer(ik),                    intent(inout)   :: ineighbor_nterms_s
        integer(ik),                    intent(inout)   :: ineighbor_proc
        real(rk),   allocatable,        intent(inout)   :: neighbor_grad1(:,:)
        real(rk),   allocatable,        intent(inout)   :: neighbor_grad2(:,:)
        real(rk),   allocatable,        intent(inout)   :: neighbor_grad3(:,:)
        real(rk),   allocatable,        intent(inout)   :: neighbor_br2_face(:,:)
        real(rk),   allocatable,        intent(inout)   :: neighbor_br2_vol(:,:)
        real(rk),   allocatable,        intent(inout)   :: neighbor_invmass(:,:)
        real(rk),                       intent(inout)   :: neighbor_h(3)
        integer(ik),                    intent(inout)   :: neighbor_status
        type(mpi_comm),                 intent(in)      :: ChiDG_COMM

        integer(ik) :: corner_one, corner_two, corner_three, corner_four,           &
                       corner_indices(4), data(7), mapping, iproc, idomain_g, ierr, &
                       grad_size(2), invmass_size(2), br2_face_size(2), br2_vol_size(2), size_data(8)
        logical     :: neighbor_element, has_domain


        neighbor_status = NO_NEIGHBOR_FOUND

        ! Get the indices of the corner nodes that correspond to the current face 
        ! in an element connectivity list.
        mapping      = self%elems(ielem_l)%element_type
        corner_one   = FACE_CORNERS(iface,1,mapping)
        corner_two   = FACE_CORNERS(iface,2,mapping)
        corner_three = FACE_CORNERS(iface,3,mapping)
        corner_four  = FACE_CORNERS(iface,4,mapping)


        ! For the current face, get the indices of the coordinate nodes for 
        ! the corners defining a face
        corner_indices(1) = self%elems(ielem_l)%connectivity(corner_one)
        corner_indices(2) = self%elems(ielem_l)%connectivity(corner_two)
        corner_indices(3) = self%elems(ielem_l)%connectivity(corner_three)
        corner_indices(4) = self%elems(ielem_l)%connectivity(corner_four)

        
        ! Test the face nodes against other elements, if all face nodes are also 
        ! contained in another element, then they are neighbors.
        neighbor_element = .false.



        do iproc = 0,NRANK-1
            if ( iproc /= IRANK ) then

                ! send global domain index of mesh being searched
                idomain_g = self%elems(ielem_l)%idomain_g
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
                        call MPI_Recv(data,7,MPI_INTEGER4,iproc,4,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        ineighbor_domain_g  = data(1)
                        ineighbor_domain_l  = data(2)
                        ineighbor_element_g = data(3)
                        ineighbor_element_l = data(4)
                        ineighbor_face      = data(5)
                        ineighbor_neqns     = data(6)
                        ineighbor_nterms_s  = data(7)
                        ineighbor_proc      = iproc

                        call MPI_Recv(size_data,8,MPI_INTEGER4,iproc,5,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(grad_size,    2,MPI_INTEGER4,iproc,5,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(br2_face_size,2,MPI_INTEGER4,iproc,6,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(br2_vol_size, 2,MPI_INTEGER4,iproc,7,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(invmass_size, 2,MPI_INTEGER4,iproc,8,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                        grad_size     = size_data(1:2)
                        br2_face_size = size_data(3:4)
                        br2_vol_size  = size_data(5:6)
                        invmass_size  = size_data(7:8)

                        if (allocated(neighbor_grad1)) deallocate(neighbor_grad1,neighbor_grad2,neighbor_grad3, &
                                                                neighbor_br2_face, neighbor_br2_vol, neighbor_invmass)
                        allocate(neighbor_grad1(grad_size(1),grad_size(2)), &
                                 neighbor_grad2(grad_size(1),grad_size(2)), &
                                 neighbor_grad3(grad_size(1),grad_size(2)), &
                                 neighbor_br2_face(br2_face_size(1),br2_face_size(2)), &
                                 neighbor_br2_vol(br2_vol_size(1),br2_vol_size(2)),    &
                                 neighbor_invmass(invmass_size(1),invmass_size(2)),  stat=ierr)
                        if (ierr /= 0) call AllocationError

                        !call MPI_Recv(neighbor_grad1,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc,  9, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_grad2,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc, 10, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_grad3,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc, 11, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_br2_face,  br2_face_size(1)*br2_face_size(2),  MPI_REAL8, iproc, 12, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_br2_vol,   br2_vol_size(1)*br2_vol_size(2),    MPI_REAL8, iproc, 13, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_invmass,   invmass_size(1)*invmass_size(2),    MPI_REAL8, iproc, 14, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        !call MPI_Recv(neighbor_h,         3,                                  MPI_REAL8, iproc, 15, ChiDG_COMM, MPI_STATUS_IGNORE,ierr) 
                        call MPI_Recv(neighbor_grad1,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc,  6, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_grad2,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc,  7, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_grad3,     grad_size(1)*grad_size(2),          MPI_REAL8, iproc,  8, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_br2_face,  br2_face_size(1)*br2_face_size(2),  MPI_REAL8, iproc,  9, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_br2_vol,   br2_vol_size(1)*br2_vol_size(2),    MPI_REAL8, iproc, 10, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_invmass,   invmass_size(1)*invmass_size(2),    MPI_REAL8, iproc, 11, ChiDG_COMM, MPI_STATUS_IGNORE,ierr)
                        call MPI_Recv(neighbor_h,         3,                                  MPI_REAL8, iproc, 12, ChiDG_COMM, MPI_STATUS_IGNORE,ierr) 

                        neighbor_status = NEIGHBOR_FOUND

                    end if
                end if


            end if
        end do




    end subroutine find_neighbor_global
    !*****************************************************************************************







    !>  Return the processor ranks that the current mesh is receiving from.
    !!
    !!  This includes interior neighbor elements located on another processor and also chimera
    !!  donor elements located on another processor.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_procs(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        type(ivector_t)             :: comm_procs_vector
        integer(ik),    allocatable :: comm_procs(:), comm_procs_local(:), comm_procs_chimera(:)
        integer(ik)                 :: iproc, proc

        !
        ! Test if global communication has been initialized
        !
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,"mesh%get_recv_procs: mesh global communication not initialized")

        

        !
        ! Get procs we are receiving neighbor data from
        !
        comm_procs_local = self%get_recv_procs_local()
        do iproc = 1,size(comm_procs_local)
            proc = comm_procs_local(iproc)
            call comm_procs_vector%push_back_unique(proc)
        end do !iproc



        !
        ! Get procs we are receiving chimera donors from
        !
        comm_procs_chimera = self%get_recv_procs_chimera()
        do iproc = 1,size(comm_procs_chimera)
            proc = comm_procs_chimera(iproc)
            call comm_procs_vector%push_back_unique(proc)
        end do !iproc



        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()



    end function get_recv_procs
    !*****************************************************************************************











    !>  Return the processor ranks that the current mesh is receiving interior 
    !!  neighbor elements from.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_procs_local(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        type(ivector_t)             :: comm_procs_vector
        integer(ik),    allocatable :: comm_procs(:)
        integer(ik)                 :: myrank, neighbor_rank, ielem, iface
        logical                     :: has_neighbor, comm_neighbor
        character(:),   allocatable :: user_msg

        !
        ! Test if global communication has been initialized
        !
        user_msg = "mesh%get_comm_procs: mesh global communication not initialized."
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,user_msg)


        !
        ! Get current processor rank
        !
        myrank = IRANK

        do ielem = 1,self%nelem
            do iface = 1,size(self%faces,2)

                ! Get face properties
                has_neighbor = ( self%faces(ielem,iface)%ftype == INTERIOR )

                ! For interior neighbor
                if ( has_neighbor ) then

                    ! Get neighbor processor rank. If off-processor, add to list uniquely
                    neighbor_rank = self%faces(ielem,iface)%ineighbor_proc
                    comm_neighbor = ( myrank /= neighbor_rank )
                    if ( comm_neighbor ) call comm_procs_vector%push_back_unique(neighbor_rank)

                end if

            end do !iface
        end do !ielem


        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()

    end function get_recv_procs_local
    !*****************************************************************************************











    !>  Return the processor ranks that the current mesh is receiving chimera donor 
    !!  elements from.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_procs_chimera(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        character(:),   allocatable :: user_msg
        integer(ik),    allocatable :: comm_procs(:)
        integer(ik)                 :: myrank, ielem, iface, ChiID, idonor, donor_rank
        logical                     :: is_chimera, comm_donor
        type(ivector_t)             :: comm_procs_vector

        !
        ! Test if global communication has been initialized
        !
        user_msg = "mesh%get_recv_procs_chimera: mesh global communication not initialized."
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,user_msg)


        !
        ! Get current processor rank
        !
        myrank = IRANK

        do ielem = 1,self%nelem
            do iface = 1,size(self%faces,2)

                ! Get face properties
                is_chimera   = ( self%faces(ielem,iface)%ftype == CHIMERA  )

                if ( is_chimera ) then

                    ! Loop through donor elements. If off-processor, add to list uniquely
                    ChiID = self%faces(ielem,iface)%ChiID
                    do idonor = 1,self%chimera%recv(ChiID)%ndonors()
                        !donor_rank = self%chimera%recv(ChiID)%donor_proc%at(idonor)
                        donor_rank = self%chimera%recv(ChiID)%donor(idonor)%iproc
                        comm_donor = ( myrank /= donor_rank )
                        if ( comm_donor ) call comm_procs_vector%push_back_unique(donor_rank)
                    end do !idonor

                end if !is_chimera


            end do !iface
        end do !ielem


        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()

    end function get_recv_procs_chimera
    !*****************************************************************************************










    !>  Return the processor ranks that the current mesh is sending to.
    !!
    !!  This includes processors that are being sent interior neighbor elements and also 
    !!  chimera donor elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_send_procs(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        type(ivector_t)             :: comm_procs_vector
        integer(ik),    allocatable :: comm_procs(:), comm_procs_local(:), comm_procs_chimera(:)
        integer(ik)                 :: iproc, proc
        character(:),   allocatable :: user_msg



        !
        ! Test if global communication has been initialized
        !
        user_msg = "mesh%get_send_procs: mesh global communication not initialized."
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,user_msg)



        !
        ! Get procs we are receiving neighbor data from
        !
        comm_procs_local = self%get_send_procs_local()
        do iproc = 1,size(comm_procs_local)
            proc = comm_procs_local(iproc)
            call comm_procs_vector%push_back_unique(proc)
        end do !iproc



        !
        ! Get procs we are receiving chimera donors from
        !
        comm_procs_chimera = self%get_send_procs_chimera()
        do iproc = 1,size(comm_procs_chimera)
            proc = comm_procs_chimera(iproc)
            call comm_procs_vector%push_back_unique(proc)
        end do !iproc



        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()




    end function get_send_procs
    !*****************************************************************************************








    !>  Return the processor ranks that the current mesh is sending to.
    !!
    !!  This includes processors that are being sent interior neighbor elements.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_send_procs_local(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        type(ivector_t)             :: comm_procs_vector
        integer(ik),    allocatable :: comm_procs(:)
        integer(ik)                 :: myrank, neighbor_rank, ielem, iface
        logical                     :: has_neighbor, comm_neighbor
        character(:),   allocatable :: user_msg

        !
        ! Test if global communication has been initialized
        !
        user_msg = "mesh%get_send_procs_local: mesh global communication not initialized."
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,user_msg)


        !
        ! Get current processor rank
        !
        myrank = IRANK

        do ielem = 1,self%nelem
            do iface = 1,size(self%faces,2)

                ! Get face properties
                has_neighbor = ( self%faces(ielem,iface)%ftype == INTERIOR )

                ! For interior neighbor
                if ( has_neighbor ) then
                    ! Get neighbor processor rank. If off-processor add to list uniquely
                    neighbor_rank = self%faces(ielem,iface)%ineighbor_proc
                    comm_neighbor = ( myrank /= neighbor_rank )
                    if ( comm_neighbor ) call comm_procs_vector%push_back_unique(neighbor_rank)
                end if
                

            end do !iface
        end do !ielem


        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()

    end function get_send_procs_local
    !*****************************************************************************************











    !>  Return the processor ranks that the current mesh is sending to.
    !!
    !!  This includes processors that are being sent chimera donor elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_send_procs_chimera(self) result(comm_procs)
        class(domain_t),   intent(in)  :: self

        type(ivector_t)             :: comm_procs_vector
        integer(ik),    allocatable :: comm_procs(:)
        integer(ik)                 :: myrank, isend_elem, isend_proc, send_rank
        logical                     :: comm_donor
        character(:),   allocatable :: user_msg

        !
        ! Test if global communication has been initialized
        !
        user_msg = "mesh%get_send_procs_chimera: mesh global communication not initialized."
        if ( .not. self%global_comm_initialized) call chidg_signal(WARN,user_msg)


        !
        ! Get current processor rank
        !
        myrank = IRANK


        !
        ! Collect processors that we are sending chimera donor elements to
        !
        do isend_elem = 1,self%chimera%nsend()
            do isend_proc = 1,self%chimera%send(isend_elem)%nsend_procs()

                ! Get donor rank. If off-processor, add to list uniquely.
                send_rank = self%chimera%send(isend_elem)%send_procs%at(isend_proc)
                comm_donor = (myrank /= send_rank)
                if ( comm_donor ) call comm_procs_vector%push_back_unique(send_rank)

            end do !isend_proc
        end do !isend_elem





        !
        ! Set vector data to array to be returned.
        !
        comm_procs = comm_procs_vector%data()


    end function get_send_procs_chimera
    !*****************************************************************************************







    !>  Return the number of elements in the original unpartitioned domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nelements_global(self) result(nelements)
        class(domain_t),   intent(in)  :: self

        integer(ik) :: nelements

        nelements = self%nelements_g

    end function get_nelements_global
    !****************************************************************************************






    !>  Return the number of elements in the local, possibly partitioned, domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_nelements_local(self) result(nelements)
        class(domain_t),   intent(in)  :: self

        integer(ik) :: nelements

        nelements = self%nelem

    end function get_nelements_local
    !****************************************************************************************





    !>
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(domain_t), intent(inout) :: self

    
    end subroutine destructor
    !*************************************************************************************


end module type_domain
