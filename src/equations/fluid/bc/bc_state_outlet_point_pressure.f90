module bc_state_outlet_point_pressure
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO, NO_PROC

    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mod_chidg_mpi,          only: IRANK
    use mpi_f08,                only: mpi_comm, MPI_REAL8, MPI_Allgather, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_NULL, MPI_Comm_compare
    use DNAD_D
    implicit none


    !>  Name: Outlet - Point Pressure
    !!      : Set static pressure at a single point
    !!      : Extrapolate everything at all other points. Essentially, d/dxi = 0
    !!
    !!  Options:
    !!      : Static Pressure
    !!      : Coordinate-1
    !!      : Coordinate-2
    !!      : Coordinate-3
    !!
    !!  Behavior:
    !!      : Boundary conditions only compute things at quadrature nodes. 
    !!      : This boundary condition will search for the quadrature node closest
    !!      : to the user-specified coordinate. The input pressure will be set at this 
    !!      : node only, which is probably slightly different than the user-specified
    !!      : coordinate. As order of accuracy is increased, more nodes exist in the node
    !!      : set and the node set point should get closer to the user-specified value.
    !!      
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_point_pressure_t

        !
        ! Data initialized in init_bc_specialized that defines at which particular
        ! quadrature node the user-specified static pressure will be set.
        !
        integer(ik) :: node_idomain_g  = 0
        integer(ik) :: node_ielement_g = 0
        integer(ik) :: node_iface      = 0
        integer(ik) :: node_index      = 0

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_specialized  ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type outlet_point_pressure_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_point_pressure_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name('Outlet - Point Pressure')
        call self%set_family('Outlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('Static Pressure','Required')
        call self%bcproperties%add('Coordinate-1',   'Required')
        call self%bcproperties%add('Coordinate-2',   'Required')
        call self%bcproperties%add('Coordinate-3',   'Required')


    end subroutine init
    !********************************************************************************







    !>  Specialized initialization routine.
    !!
    !!  Find the quadrature node across the boundary condition patch that is 
    !!  closest to the user-specified value that we would like to set a pressure
    !!  at.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/21/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_specialized(self,mesh,bc_patch,bc_COMM)
        class(outlet_point_pressure_t), intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(bc_patch_t),               intent(in)      :: bc_patch(:)
        type(mpi_comm),                 intent(in)      :: bc_comm

        type(point_t)                               :: point
        type(point_t),  allocatable, dimension(:)   :: coords, quad_pts
        real(rk),       allocatable, dimension(:)   :: node1, node2, node3, dist_face
        integer(ik)                                 :: iface_bc, iface, ielement, ipatch, &
                                                       node_idomain_g, node_ielement_g, node_iface, &
                                                       node_index, minloc_dist_face,     &
                                                       idomain_l, ielement_l, rank_min, bc_IRANK, bc_NRANK, ierr
        real(rk)                                    :: time, dist_bc
        real(rk),       allocatable                 :: dist_bc_global(:)
        
        
        print*, IRANK, 'point - 1'

        !
        ! Get user-specified node. Need coord,time here because thats the
        ! function interface. They aren't really doing anything.
        !
        point%c1_ = 0._rk
        point%c2_ = 0._rk
        point%c3_ = 0._rk
        coords = [point]
        time   = 0._rk
        node1  = self%bcproperties%compute('Coordinate-1',time,coords)
        node2  = self%bcproperties%compute('Coordinate-2',time,coords)
        node3  = self%bcproperties%compute('Coordinate-3',time,coords)


        print*, IRANK, 'point - 2'


        !
        ! Loop through proc-local patches, find quadrature node closest to user-specified node.
        !
        dist_bc = huge(ONE)  !initialize closest distance on the proc-local bc


        !
        ! Search 
        !
        do ipatch = 1,size(bc_patch)
            do iface_bc = 1,bc_patch(ipatch)%nfaces()

        print*, IRANK, 'point - 3'
                ! get face location
                idomain_l  = bc_patch(ipatch)%idomain_l_%at(iface_bc)
                ielement_l = bc_patch(ipatch)%ielement_l_%at(iface_bc)
                iface      = bc_patch(ipatch)%iface_%at(iface_bc)
                
                ! get points at quadrature nodes for current face
                quad_pts = mesh(idomain_l)%faces(ielement_l,iface)%quad_pts

                ! compute distance from quadrature nodes to user node
                dist_face = sqrt( (quad_pts(:)%c1_ - node1(1))**TWO + &
                                  (quad_pts(:)%c2_ - node2(1))**TWO + &
                                  (quad_pts(:)%c3_ - node3(1))**TWO ) 

                ! Get location of gq minimum distance node on current face
                minloc_dist_face = minloc(dist_face,1)
            
        print*, IRANK, 'point - 4'

                ! If new minimum is found, record location, update local dist_bc
                if (dist_face(minloc_dist_face) < dist_bc) then
                    node_idomain_g  = mesh(idomain_l)%elems(ielement_l)%idomain_g
                    node_ielement_g = mesh(idomain_l)%elems(ielement_l)%ielement_g
                    node_iface      = iface
                    node_index      = minloc_dist_face

                    dist_bc = dist_face(minloc_dist_face)
                end if
                
                ! compute distance from quadrature nodes to user node
            
            end do !iface_bc
        end do !ipatch




        print*, IRANK, 'point - 5'
        !
        ! Get rank of current proc in bc communicator and total number of ranks
        !
        call MPI_Comm_rank(bc_COMM,bc_IRANK)
        call MPI_Comm_size(bc_COMM,bc_NRANK)

        print*, IRANK, 'point - 6'

        !
        ! Have communicator root gather minimum values from all members
        !
        allocate(dist_bc_global(bc_NRANK), stat=ierr)
        if (ierr /= 0) call AllocationError
        print*, 'bc_NRANK: ', bc_NRANK
        print*, (bc_COMM == MPI_COMM_NULL)
        print*, MPI_Comm_compare(bc_COMM, MPI_COMM_NULL, result)
        call MPI_Allgather(dist_bc, 1, MPI_REAL8, dist_bc_global, bc_NRANK, MPI_REAL8, bc_COMM)
!
!        print*, IRANK, 'point - 7'
!
!        !
!        ! Determine which processor has the minimum distance.
!        ! Minus 1 because ranks are 0-based.
!        !
!        rank_min = minloc(dist_bc_global,1) - 1
!
!        print*, IRANK, 'point - 8'
!
!        !
!        ! If minimum distance is on current proc, store node values
!        !
!        if (rank_min == bc_IRANK) then
!            self%node_idomain_g  = node_idomain_g
!            self%node_ielement_g = node_ielement_g
!            self%node_iface      = node_iface
!            self%node_index      = node_index
!        end if
!
!        print*, IRANK, 'point - 9'


    end subroutine init_bc_specialized
    !********************************************************************************











    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(outlet_point_pressure_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            u_bc,   v_bc,    w_bc,  H_bc, p_m, p_bc


        logical                                     :: face_has_node
        real(rk)                                    :: time, gam_m
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            p_user, norm_1, norm_2, norm_3


        !
        ! Get back pressure from function.
        !
        coords = worker%coords()
        time   = worker%time()
        p_user = self%bcproperties%compute('Static Pressure',time,coords)



        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')



        drho_dx_m  = worker%get_primary_field_face('Density'   , 'grad1', 'face interior')
        drho_dy_m  = worker%get_primary_field_face('Density'   , 'grad2', 'face interior')
        drho_dz_m  = worker%get_primary_field_face('Density'   , 'grad3', 'face interior')

        drhou_dx_m = worker%get_primary_field_face('Momentum-1', 'grad1', 'face interior')
        drhou_dy_m = worker%get_primary_field_face('Momentum-1', 'grad2', 'face interior')
        drhou_dz_m = worker%get_primary_field_face('Momentum-1', 'grad3', 'face interior')

        drhov_dx_m = worker%get_primary_field_face('Momentum-2', 'grad1', 'face interior')
        drhov_dy_m = worker%get_primary_field_face('Momentum-2', 'grad2', 'face interior')
        drhov_dz_m = worker%get_primary_field_face('Momentum-2', 'grad3', 'face interior')

        drhow_dx_m = worker%get_primary_field_face('Momentum-3', 'grad1', 'face interior')
        drhow_dy_m = worker%get_primary_field_face('Momentum-3', 'grad2', 'face interior')
        drhow_dz_m = worker%get_primary_field_face('Momentum-3', 'grad3', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face('Energy'    , 'grad1', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face('Energy'    , 'grad2', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face('Energy'    , 'grad3', 'face interior')




        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
        end if


        !
        ! Compute gamma
        !
        gam_m = 1.4_rk


        !
        ! Extrapolate density and momentum
        !
        density_bc = density_m
        mom1_bc    = mom1_m
        mom2_bc    = mom2_m
        mom3_bc    = mom3_m


        !
        ! Compute velocities
        !
        u_bc = mom1_bc/density_bc
        v_bc = mom2_bc/density_bc
        w_bc = mom3_bc/density_bc


        !
        ! Compute boundary condition energy
        !
        p_m = (gam_m-ONE)*(energy_m - HALF*( (mom1_m*mom1_m) + (mom2_m*mom2_m) + (mom3_m*mom3_m) )/density_m )

        face_has_node = (worker%element_info%idomain_g  == self%node_idomain_g)  .and. &
                        (worker%element_info%ielement_g == self%node_ielement_g) .and. &
                        (worker%iface                   == self%node_iface)

        p_bc = p_m
        if (face_has_node) then

            !
            ! Face contains user-specified node. 
            ! Set pressure at particular node, keep extrapolated values for all others
            !
            p_bc(self%node_index) = p_user(1)

        end if !face_has_node


        energy_bc = p_bc/(gam_m - ONE)  + (density_bc*HALF)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * worker%coordinate('1','boundary')
        end if



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')





        call worker%store_bc_state('Density'   , drho_dx_m,  'grad1')
        call worker%store_bc_state('Density'   , drho_dy_m,  'grad2')
        call worker%store_bc_state('Density'   , drho_dz_m,  'grad3')
                                                
        call worker%store_bc_state('Momentum-1', drhou_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-1', drhou_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-1', drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-2', drhov_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-2', drhov_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-2', drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-3', drhow_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-3', drhow_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-3', drhow_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Energy'    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state('Energy'    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state('Energy'    , drhoE_dz_m, 'grad3')







    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_outlet_point_pressure
