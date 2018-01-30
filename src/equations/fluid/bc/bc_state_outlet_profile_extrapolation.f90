module bc_state_outlet_profile_extrapolation
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO
    use mod_fluid,              only: gam, Rgas

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use ieee_arithmetic
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_profile_extrapolation_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
!        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! Boundary condition function implementation

    end type outlet_profile_extrapolation_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_profile_extrapolation_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("Outlet - Profile Extrapolation")
        call self%set_family("Outlet")


        !
        ! Add functions
        !
        call self%bcproperties%add('Static Pressure','Required')


    end subroutine init
    !********************************************************************************





!    !>  Initialize boundary group coupling.
!    !!
!    !!  For this profile extrapolation outlet, each patch face is coupled with the 
!    !!  same one face in the bc_group. This coupling occurs because each face uses 
!    !!  a pressure offset pressure that is computed from the same face. 
!    !!
!    !!  Coupling initialization:
!    !!      1: each process loops through its local faces, initializes coupling
!    !!         of all local faces with all other local faces.
!    !!
!    !!      2: loop through ranks in bc_COMM
!    !!          a: iproc broadcasts information about its coupling to bc_COMM
!    !!          b: all other procs receive from iproc and initialize parallel coupling
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   12/28/2017
!    !!
!    !--------------------------------------------------------------------------------
!    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
!        class(outlet_profile_extrapolation_t),  intent(inout)   :: self
!        type(mesh_t),                           intent(inout)   :: mesh
!        integer(ik),                            intent(in)      :: group_ID
!        type(mpi_comm),                         intent(in)      :: bc_COMM
!
!        integer(ik) :: patch_ID, face_ID, elem_ID, patch_ID_coupled, face_ID_coupled,   &
!                       idomain_g, idomain_l, ielement_g, ielement_l, iface,             &
!                       bc_IRANK, bc_NRANK, ierr, iproc, nbc_elements,     &
!                       ielem, neqns, nterms_s, ngq, ibc
!
!        integer(ik) :: idomain_g_coupled, idomain_l_coupled, ielement_g_coupled, ielement_l_coupled, &
!                       iface_coupled, proc_coupled
!
!        real(rk),       allocatable :: interp_coords_def(:,:)
!        real(rk),       allocatable :: areas(:)
!        real(rk)                    :: total_area
!
!
!
!        !
!        ! For each face, initialize coupling with all faces on the current processor.
!        !
!        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
!            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!                
!                !
!                ! Loop through, initialize coupling
!                !
!                do patch_ID_coupled = 1,mesh%bc_patch_group(group_ID)%npatches()
!                    do face_ID_coupled = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!
!                        !
!                        ! Get block-element index of current face_ID_coupled
!                        !
!                        idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_g()
!                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_l()
!                        ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_g(face_ID_coupled)
!                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_l(face_ID_coupled)
!                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%iface(     face_ID_coupled)
!
!
!                        neqns             = mesh%domain(idomain_l)%faces(ielement_l,iface)%neqns
!                        nterms_s          = mesh%domain(idomain_l)%faces(ielement_l,iface)%nterms_s
!                        total_area        = mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area
!                        areas             = mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas
!                        interp_coords_def = mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def
!
!
!
!                        !
!                        ! For the face (patch_ID,face_ID) add the element on (patch_ID_coupled,face_ID_coupled)
!                        !
!                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g,  &
!                                                                                                        idomain_l,  &
!                                                                                                        ielement_g, &
!                                                                                                        ielement_l, &
!                                                                                                        iface,      &
!                                                                                                        IRANK)
!
!                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g,     &
!                                                                                                             ielement_g,    &
!                                                                                                             neqns,         &
!                                                                                                             nterms_s,      &
!                                                                                                             total_area,    &
!                                                                                                             areas,         &
!                                                                                                             interp_coords_def)
!
!
!                    end do ! face_ID_couple
!                end do ! patch_ID_couple
!
!            end do ! face_ID
!        end do ! patch_ID
!
!
!
!
!
!
!
!        !
!        ! Get bc_NRANK, bc_IRANK from bc_COMM
!        !
!        call MPI_Comm_Size(bc_COMM, bc_NRANK, ierr)
!        call MPI_Comm_Rank(bc_COMM, bc_IRANK, ierr)
!
!
!
!
!
!        !
!        ! Initialize coupling with faces on other processors
!        !
!        do iproc = 0,bc_NRANK-1
!
!
!
!            !
!            ! Send local elements out
!            !
!            if (iproc == bc_IRANK) then
!
!
!                nbc_elements = mesh%bc_patch_group(group_ID)%nfaces()
!                call MPI_Bcast(IRANK,        1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!
!                do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
!                    do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
!                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
!                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
!                        
!                        ! Broadcast element for coupling
!                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g(),         1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l(),         1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID), 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID), 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID),      1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!
!                        ! Broadcast auxiliary data
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%neqns,      1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%nterms_s,   1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!                        ngq = size(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def,1)
!                        call MPI_Bcast(ngq,                                                                          1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas,          ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,1),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,2),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,3),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c1_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c2_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c3_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!                    end do ! face_ID
!                end do ! patch_ID
!            
!
!
!
!
!
!
!
!
!            !
!            ! All other processors recieve
!            !
!            else
!
!
!                call MPI_Bcast(proc_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!
!
!                !
!                ! For the face (patch_ID,face_ID) add each element from the sending proc
!                !
!                do ielem = 1,nbc_elements
!
!                    ! Receive coupled element
!                    call MPI_BCast(idomain_g_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(idomain_l_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(ielement_g_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(ielement_l_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(iface_coupled,      1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!
!                    ! Receive auxiliary data
!                    call MPI_BCast(neqns,     1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(nterms_s,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    call MPI_BCast(total_area,1, MPI_INTEGER, iproc, bc_COMM, ierr)
!
!
!                    call MPI_BCast(ngq, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
!                    if (allocated(areas) ) deallocate(areas, interp_coords_def)
!                    allocate(areas(ngq), interp_coords_def(ngq,3), stat=ierr)
!                    if (ierr /= 0) call AllocationError
!
!
!                    call MPI_BCast(areas,           ngq, MPI_REAL8, iproc, bc_COMM, ierr)
!                    call MPI_BCast(interp_coords_def(:,1), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
!                    call MPI_BCast(interp_coords_def(:,2), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
!                    call MPI_BCast(interp_coords_def(:,3), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
!
!
!                    !
!                    ! Each face on the current proc adds the off-processor element to their list 
!                    ! of coupled elems
!                    !
!                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
!                        do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()
!
!                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g_coupled,     &
!                                                                                                            idomain_l_coupled,     &
!                                                                                                            ielement_g_coupled,    &
!                                                                                                            ielement_l_coupled,    &
!                                                                                                            iface_coupled,         &
!                                                                                                            proc_coupled)
!
!                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g_coupled,     &
!                                                                                                                 ielement_g_coupled,    &
!                                                                                                                 neqns,                 &
!                                                                                                                 nterms_s,              &
!                                                                                                                 total_area,            &
!                                                                                                                 areas,                 &
!                                                                                                                 interp_coords_def)
!
!
!
!
!
!                        end do ! face_ID
!                    end do ! patch_ID
!
!                end do !ielem
!
!
!
!
!            end if
!
!
!
!
!            call MPI_Barrier(bc_COMM,ierr)
!        end do
!
!
!
!
!
!
!
!    end subroutine init_bc_coupling
!    !******************************************************************************************




















    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_profile_extrapolation_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_bc,   v_bc,    w_bc,  T_m, T_bc, p_m, p_bc, p_element
            

        real(rk),   allocatable, dimension(:) :: r
        real(rk),   allocatable, dimension(:) :: p_input
        !type(AD_D), allocatable, dimension(:) :: p_bc

        type(AD_D)  :: p_diff



        !
        ! Get back pressure from function.
        !
        p_input = self%bcproperties%compute('Static Pressure',worker%time(),worker%coords())
        !p_bc = worker%get_field('Pressure_TEMP', 'value', 'face interior')


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density_m = worker%get_field('Density'    , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1' , 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2' , 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3' , 'value', 'face interior')
        energy_m  = worker%get_field('Energy'     , 'value', 'face interior')
        p_m       = worker%get_field('Pressure'   , 'value', 'face interior')
        T_m       = worker%get_field('Temperature', 'value', 'face interior')
        p_element = worker%get_field('Pressure'   , 'value', 'element')
        print*, 'nnodes 2d: ', size(p_m)
        print*, 'nnodes 3d: ', size(p_element)
        ! P2
        !p_m = p_element(101:125)
        ! P3
        !p_m = p_element(181:216)
        ! P4
        p_m = p_element(449:512)
        ! P5
        !p_m = p_element(649:729)
        ! P6
        !p_m = p_element(1211:1331)



        grad1_density_m = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m    = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field('Momentum-3', 'grad3', 'face interior')
        
        grad1_energy_m  = worker%get_field('Energy'    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_field('Energy'    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_field('Energy'    , 'grad3', 'face interior')




        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if


        !
        ! Extrapolate temperature and velocity
        !
        T_bc = T_m
        u_bc = mom1_m/density_m
        v_bc = mom2_m/density_m
        w_bc = mom3_m/density_m


        !
        ! Extrapolate pressure, adjust by dp for a point
        !
        ! Confirmed, signs are correct
        p_diff = (p_m(1) - p_input(1))
        p_bc = p_m - p_diff

        

        !
        ! Compute density, momentum, energy
        !
        density_bc = p_bc/(Rgas*T_bc)
        mom1_bc    = u_bc*density_bc
        mom2_bc    = v_bc*density_bc
        mom3_bc    = w_bc*density_bc
        energy_bc  = p_bc/(gam - ONE) + (density_bc*HALF)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')





        call worker%store_bc_state('Density'   , grad1_density_m, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density_m, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', grad1_mom1_m,    'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1_m,    'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-2', grad1_mom2_m,    'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2_m,    'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-3', grad1_mom3_m,    'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3_m,    'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3_m,    'grad3')
                                                
        call worker%store_bc_state('Energy'    , grad1_energy_m,  'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy_m,  'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy_m,  'grad3')







    end subroutine compute_bc_state
    !**************************************************************************************






end module bc_state_outlet_profile_extrapolation
