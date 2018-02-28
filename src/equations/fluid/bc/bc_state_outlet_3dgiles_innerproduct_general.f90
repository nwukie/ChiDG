module bc_state_outlet_3dgiles_innerproduct_general
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_fluid,              only: gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast, MPI_MIN, MPI_MAX
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none





    !>  Name: Outlet - 3D Giles
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!  Behavior:
    !!      
    !!  References:
    !!              
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_3dgiles_innerproduct_general_t

        real(rk),   allocatable :: r(:)
        real(rk),   allocatable :: theta(:)
        real(rk)                :: theta_ref

        type(element_info_t),   allocatable :: donor(:,:)
        real(rk),               allocatable :: donor_coord(:,:,:)
        

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Implement coupling pattern
        procedure   :: init_bc_postcomm     ! Implement specialized initialization
        procedure   :: compute_bc_state     ! boundary condition function implementation

        procedure   :: compute_averages
        procedure   :: compute_fourier_decomposition
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

    end type outlet_3dgiles_innerproduct_general_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name('Outlet - 3D Giles Innerproduct General')
        call self%set_family('Outlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('Average Pressure','Required')
        call self%bcproperties%add('Pitch',           'Required')


    end subroutine init
    !********************************************************************************







    !>  Initialize boundary group coupling.
    !!
    !!  For this LODI-based outlet, each patch face is coupled with every other
    !!  face in the bc_group. This coupling occurs because each face uses an
    !!  average pressure that is computed over the group. The average pressure
    !!  calculation couples every element on the group. This coupling is initialized
    !!  here.
    !!
    !!  Coupling initialization:
    !!      1: each process loops through its local faces, initializes coupling
    !!         of all local faces with all other local faces.
    !!
    !!      2: loop through ranks in bc_COMM
    !!          a: iproc broadcasts information about its coupling to bc_COMM
    !!          b: all other procs receive from iproc and initialize parallel coupling
    !!
    !!  @author Nathan A. average_pressure
    !!  @date   4/18/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
        class(outlet_3dgiles_innerproduct_general_t),    intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(in)      :: group_ID
        type(mpi_comm),             intent(in)      :: bc_COMM

        integer(ik) :: patch_ID, face_ID, elem_ID, patch_ID_coupled, face_ID_coupled,   &
                       idomain_g, idomain_l, ielement_g, ielement_l, iface,             &
                       bc_IRANK, bc_NRANK, ierr, iproc, nbc_elements,     &
                       ielem, neqns, nterms_s, ngq, ibc

        integer(ik) :: idomain_g_coupled, idomain_l_coupled, ielement_g_coupled, ielement_l_coupled, &
                       iface_coupled, proc_coupled

        real(rk),       allocatable :: interp_coords_def(:,:)
        real(rk),       allocatable :: areas(:)
        real(rk)                    :: total_area



        !
        ! For each face, initialize coupling with all faces on the current processor.
        !
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                
                !
                ! Loop through, initialize coupling
                !
                do patch_ID_coupled = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID_coupled = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()


                        !
                        ! Get block-element index of current face_ID_coupled
                        !
                        idomain_g  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_g()
                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%idomain_l()
                        ielement_g = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_g(face_ID_coupled)
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%ielement_l(face_ID_coupled)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID_coupled)%iface(     face_ID_coupled)


                        neqns      = mesh%domain(idomain_l)%faces(ielement_l,iface)%neqns
                        nterms_s   = mesh%domain(idomain_l)%faces(ielement_l,iface)%nterms_s
                        total_area = mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area
                        areas      = mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas
                        interp_coords_def = mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def



                        !
                        ! For the face (patch_ID,face_ID) add the element on (patch_ID_coupled,face_ID_coupled)
                        !
                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g,  &
                                                                                                        idomain_l,  &
                                                                                                        ielement_g, &
                                                                                                        ielement_l, &
                                                                                                        iface,      &
                                                                                                        IRANK)

                        call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g,     &
                                                                                                             ielement_g,    &
                                                                                                             neqns,         &
                                                                                                             nterms_s,      &
                                                                                                             total_area,    &
                                                                                                             areas,         &
                                                                                                             interp_coords_def)


                    end do ! face_ID_couple
                end do ! patch_ID_couple

            end do ! face_ID
        end do ! patch_ID







        !
        ! Get bc_NRANK, bc_IRANK from bc_COMM
        !
        call MPI_Comm_Size(bc_COMM, bc_NRANK, ierr)
        call MPI_Comm_Rank(bc_COMM, bc_IRANK, ierr)





        !
        ! Initialize coupling with faces on other processors
        !
        do iproc = 0,bc_NRANK-1



            !
            ! Send local elements out
            !
            if (iproc == bc_IRANK) then


                nbc_elements = mesh%bc_patch_group(group_ID)%nfaces()
                call MPI_Bcast(IRANK,        1, MPI_INTEGER, iproc, bc_COMM, ierr)
                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)


                do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                    do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                        idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                        ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                        iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
                        
                        ! Broadcast element for coupling
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g(),         1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l(),         1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID), 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID), 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID),      1, MPI_INTEGER, iproc, bc_COMM, ierr)


                        ! Broadcast auxiliary data
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%neqns,      1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%nterms_s,   1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%total_area, 1, MPI_INTEGER, iproc, bc_COMM, ierr)

                        ngq = size(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def,1)
                        call MPI_Bcast(ngq,                                                                          1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas,          ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,1),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,2),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:,3),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)

                    end do ! face_ID
                end do ! patch_ID
            








            !
            ! All other processors recieve
            !
            else


                call MPI_Bcast(proc_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                call MPI_Bcast(nbc_elements, 1, MPI_INTEGER, iproc, bc_COMM, ierr)



                !
                ! For the face (patch_ID,face_ID) add each element from the sending proc
                !
                do ielem = 1,nbc_elements

                    ! Receive coupled element
                    call MPI_BCast(idomain_g_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(idomain_l_coupled,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(ielement_g_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(ielement_l_coupled, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(iface_coupled,      1, MPI_INTEGER, iproc, bc_COMM, ierr)


                    ! Receive auxiliary data
                    call MPI_BCast(neqns,     1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(nterms_s,  1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    call MPI_BCast(total_area,1, MPI_INTEGER, iproc, bc_COMM, ierr)


                    call MPI_BCast(ngq, 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                    if (allocated(areas) ) deallocate(areas, interp_coords_def)
                    allocate(areas(ngq), interp_coords_def(ngq,3), stat=ierr)
                    if (ierr /= 0) call AllocationError


                    call MPI_BCast(areas,                  ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,1), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,2), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(interp_coords_def(:,3), ngq, MPI_REAL8, iproc, bc_COMM, ierr)


                    !
                    ! Each face on the current proc adds the off-processor element to their list 
                    ! of coupled elems
                    !
                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                        do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%add_coupled_element(face_ID, idomain_g_coupled,     &
                                                                                                            idomain_l_coupled,     &
                                                                                                            ielement_g_coupled,    &
                                                                                                            ielement_l_coupled,    &
                                                                                                            iface_coupled,         &
                                                                                                            proc_coupled)

                            call mesh%bc_patch_group(group_ID)%patch(patch_ID)%set_coupled_element_data(face_ID, idomain_g_coupled,     &
                                                                                                                 ielement_g_coupled,    &
                                                                                                                 neqns,                 &
                                                                                                                 nterms_s,              &
                                                                                                                 total_area,            &
                                                                                                                 areas,                 &
                                                                                                                 interp_coords_def)


                        end do ! face_ID
                    end do ! patch_ID

                end do !ielem

            end if

            call MPI_Barrier(bc_COMM,ierr)
        end do


    end subroutine init_bc_coupling
    !*************************************************************************************








    !>  Default specialized initialization procedure. This is called from the base bc%init procedure
    !!  and can be overwritten by derived types to implement specialized initiailization details.
    !!
    !!  By default, this routine does nothing. However, a particular bc_state_t could reimplement
    !!  this routine to perform some specialized initialization calculations during initialization.
    !!
    !!  For example, a point pressure outlet boundary condition may want to find a particular 
    !!  quadrature node to set pressure at. init_bc_specialized could be defined for that
    !!  bc_state_t implementation to search the quadrature nodes over all the bc_patch faces
    !!  to find the correct node to set the pressure at.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2018
    !!
    !----------------------------------------------------------------------------------------------
    subroutine init_bc_postcomm(self,mesh,group_ID,bc_comm)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(mesh_t),                                   intent(inout)   :: mesh
        integer(ik),                                    intent(in)      :: group_ID
        type(mpi_comm),                                 intent(in)      :: bc_comm


        call self%analyze_bc_geometry(mesh,group_ID,bc_comm)


        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm)


    end subroutine init_bc_postcomm
    !**********************************************************************************************









!
!    !>  Update the area-averaged pressure for the boundary condition.
!    !!
!    !!  @author Nathan A. average_pressure
!    !!  @date   3/31/2017
!    !!
!    !!
!    !-------------------------------------------------------------------------------------
!    subroutine compute_averages(self,worker,bc_COMM, vel1_avg, vel2_avg, vel3_avg, density_avg, p_avg)
!        class(outlet_3dgiles_innerproduct_general_t),    intent(inout)   :: self
!        type(chidg_worker_t),       intent(inout)   :: worker
!        type(mpi_comm),             intent(in)      :: bc_COMM
!        type(AD_D),                 intent(inout)   :: vel1_avg
!        type(AD_D),                 intent(inout)   :: vel2_avg
!        type(AD_D),                 intent(inout)   :: vel3_avg
!        type(AD_D),                 intent(inout)   :: density_avg
!        type(AD_D),                 intent(inout)   :: p_avg
!
!        type(face_info_t)   :: face_info
!
!        type(AD_D), allocatable,    dimension(:)    ::  &
!            density, mom1, mom2, mom3, energy, p, u, v, w
!
!        type(AD_D)  :: p_integral, u_integral, v_integral, w_integral, density_integral,    &
!                       face_density, face_u, face_v, face_w, face_p
!
!
!        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
!                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
!                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
!                       ielement_l_coupled, iface_coupled
!
!        real(rk),   allocatable,    dimension(:)    :: weights, areas, r
!        real(rk)    :: face_area, total_area
!
!
!
!        !
!        ! Zero integrated quantities
!        !
!        total_area = ZERO
!
!
!        ! Get location on domain
!        idomain_l  = worker%element_info%idomain_l
!        ielement_l = worker%element_info%ielement_l
!        iface      = worker%iface
!
!        ! Get location on bc_patch_group
!        group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
!        patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
!        face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID
!
!
!
!
!
!        !
!        ! Loop through coupled faces and compute their contribution to the average pressure
!        !
!        do icoupled = 1,worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)
!
!            !
!            ! Get solution
!            !
!            idensity = 1
!            imom1    = 2
!            imom2    = 3
!            imom3    = 4
!            ienergy  = 5
!            itime    = 1
!
!
!            !
!            ! Get face info from coupled element we want to interpolate from
!            !
!            idomain_g_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g( icoupled)
!            idomain_l_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l( icoupled)
!            ielement_g_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(icoupled)
!            ielement_l_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(icoupled)
!            iface_coupled      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%iface(     icoupled)
!
!            face_info%idomain_g  = idomain_g_coupled
!            face_info%idomain_l  = idomain_l_coupled
!            face_info%ielement_g = ielement_g_coupled
!            face_info%ielement_l = ielement_l_coupled
!            face_info%iface      = iface_coupled
!
!            
!            !
!            ! Interpolate coupled element solution on face of coupled element
!            !
!            density = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, idensity, itime, 'value', ME)
!            mom1    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom1,    itime, 'value', ME)
!            mom2    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom2,    itime, 'value', ME)
!            mom3    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom3,    itime, 'value', ME)
!            energy  = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, ienergy,  itime, 'value', ME)
!
!            !r = worker%coordinate('1','boundary')
!            !if (worker%coordinate_system() == 'Cylindrical') then
!            !    mom2 = mom2 / r
!            !end if
!            if (worker%coordinate_system() == 'Cylindrical') then
!                mom2 = mom2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%interp_coords_def(:,1)
!            end if
!
!
!            
!            !
!            ! Compute quantities for averaging
!            !
!            u = mom1 / density
!            v = mom2 / density
!            w = mom3 / density
!            p = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)
!
!
!            !
!            ! Get weights + areas
!            !
!            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
!            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
!            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area
!
!
!
!            !
!            ! Integrate and contribute to average
!            !
!            face_density = sum(density * areas * weights)
!            face_u       = sum(u       * areas * weights)
!            face_v       = sum(v       * areas * weights)
!            face_w       = sum(w       * areas * weights)
!            face_p       = sum(p       * areas * weights)
!
!
!            !
!            ! Allocate derivatives and clear integral for first face.
!            !
!            if (icoupled == 1) then
!                density_integral = face_u
!                u_integral       = face_u
!                v_integral       = face_u
!                w_integral       = face_u
!                p_integral       = face_u
!                density_integral = ZERO
!                u_integral       = ZERO
!                v_integral       = ZERO
!                w_integral       = ZERO
!                p_integral       = ZERO
!            end if
!
!
!            !
!            ! Accumulate face contribution.
!            !
!            density_integral = density_integral + face_density
!            u_integral       = u_integral       + face_u
!            v_integral       = v_integral       + face_v
!            w_integral       = w_integral       + face_w
!            p_integral       = p_integral       + face_p
!
!            
!            ! Accumulate surface area
!            total_area = total_area + face_area
!
!
!        end do !icoupled
!
!
!
!                                                      
!        !                                             
!        ! Compute average pressure:
!        !   area-weighted pressure integral over the total area
!        !   
!        !
!        vel1_avg    = u_integral       / total_area
!        vel2_avg    = v_integral       / total_area
!        vel3_avg    = w_integral       / total_area
!        density_avg = density_integral / total_area
!        p_avg       = p_integral       / total_area
!
!
!
!    end subroutine compute_averages
!    !************************************************************************************



    !>  Compute averaged quantities over the face. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, u_avg, v_avg, w_avg, density_avg, p_avg)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        type(mpi_comm),                                 intent(in)      :: bc_COMM
        type(AD_D),                                     intent(inout)   :: u_avg
        type(AD_D),                                     intent(inout)   :: v_avg
        type(AD_D),                                     intent(inout)   :: w_avg
        type(AD_D),                                     intent(inout)   :: density_avg
        type(AD_D),                                     intent(inout)   :: p_avg

        type(AD_D)          :: p_integral, u_integral, v_integral, w_integral, density_integral, &
                               face_density, face_u, face_v, face_w, face_p
        type(face_info_t)   :: face_info

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3, energy, p,       &
            u, v, w

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r

        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
                       ielement_l_coupled, iface_coupled
        real(rk)    :: face_area, total_area



        !
        ! Zero integrated quantities
        !
        total_area = ZERO


        ! Get location on domain
        idomain_l  = worker%element_info%idomain_l
        ielement_l = worker%element_info%ielement_l
        iface      = worker%iface

        ! Get location on bc_patch_group
        group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
        patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
        face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID


        !
        ! Loop through coupled faces and compute their contribution to the average pressure
        !
        do icoupled = 1,worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

            !
            ! Get solution
            !
            idensity = 1
            imom1    = 2
            imom2    = 3
            imom3    = 4
            ienergy  = 5
            itime    = 1


            !
            ! Get face info from coupled element we want to interpolate from
            !
            idomain_g_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g( icoupled)
            idomain_l_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l( icoupled)
            ielement_g_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(icoupled)
            ielement_l_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(icoupled)
            iface_coupled      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%iface(     icoupled)

            face_info%idomain_g  = idomain_g_coupled
            face_info%idomain_l  = idomain_l_coupled
            face_info%ielement_g = ielement_g_coupled
            face_info%ielement_l = ielement_l_coupled
            face_info%iface      = iface_coupled

            
            !
            ! Interpolate coupled element solution on face of coupled element
            !
            density = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, idensity, itime, 'value', ME)
            mom1    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom1,    itime, 'value', ME)
            mom2    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom2,    itime, 'value', ME)
            mom3    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom3,    itime, 'value', ME)
            energy  = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, ienergy,  itime, 'value', ME)

            if (worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%coordinate_system == CYLINDRICAL) then
                mom2 = mom2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%interp_coords_def(:,1)
            end if
            !r = worker%coordinate('1','boundary')
            !if (worker%coordinate_system() == 'Cylindrical') then
            !    mom2 = mom2 / r
            !end if


            
            !
            ! Compute velocities and pressure
            !
            u = mom1 / density
            v = mom2 / density
            w = mom3 / density
            p = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )



            !
            ! Get weights + areas
            !
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area



            !
            ! Integrate and contribute to average
            !
            face_u       = sum(u       * areas * weights)
            face_v       = sum(v       * areas * weights)
            face_w       = sum(w       * areas * weights)
            face_density = sum(density * areas * weights)
            face_p       = sum(p       * areas * weights)



            if (allocated(u_integral%xp_ad_)) then
                u_integral = u_integral + face_u
            else
                u_integral = face_u
            end if

            if (allocated(v_integral%xp_ad_)) then
                v_integral = v_integral + face_v
            else
                v_integral = face_v
            end if

            if (allocated(w_integral%xp_ad_)) then
                w_integral = w_integral + face_w
            else
                w_integral = face_w
            end if

            if (allocated(p_integral%xp_ad_)) then
                p_integral = p_integral + face_p
            else
                p_integral = face_p
            end if

            if (allocated(density_integral%xp_ad_)) then
                density_integral = density_integral + face_density
            else
                density_integral = face_density
            end if


            total_area = total_area + face_area


        end do !icoupled



        !
        ! Compute average pressure:
        !   area-weighted pressure integral over the total area
        !   
        !
        u_avg       = u_integral       / total_area
        v_avg       = v_integral       / total_area
        w_avg       = w_integral       / total_area
        density_avg = density_integral / total_area
        p_avg       = p_integral       / total_area



    end subroutine compute_averages
    !*******************************************************************************************







    !>  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop
        type(mpi_comm),                                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,                            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,                           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            vel1_bc, vel2_bc, vel3_bc, pressure_bc,                                     &
            vel1_m,  vel2_m,  vel3_m,  pressure_m,                                      &
            c1,    c2,    c3,    c4,    c5,                                             &
            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar,             &
            ddensity, dvel1, dvel2, dvel3, dpressure

        type(AD_D), allocatable, dimension(:,:) ::                                              &
            density_hat_real, vel1_hat_real, vel2_hat_real, vel3_hat_real, pressure_hat_real,   &
            density_hat_imag, vel1_hat_imag, vel2_hat_imag, vel3_hat_imag, pressure_hat_imag,   &
            c1_hat_real,      c2_hat_real,   c3_hat_real,   c4_hat_real,   c5_hat_real,         &
            c1_hat_imag,      c2_hat_imag,   c3_hat_imag,   c4_hat_imag,   c5_hat_imag,         &
            c5_hat_real_gq,   c5_hat_imag_gq


        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg, c_avg,              &
                       ddensity_mean, dvel1_mean, dvel2_mean, dvel3_mean, dpressure_mean,           &
                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r,  &
                       A3_real, A3_imag, A4_real, A4_imag, beta

        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset
        type(point_t),  allocatable                 :: coords(:)
        integer                                     :: i, ngq, ivec, imode, iradius, nmodes, ierr, igq



        !
        ! Get back pressure from function.
        !
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',           worker%time(),worker%coords())



        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')



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
        ! Store boundary gradient state. Grad(Q_bc). Do this here, before we
        ! compute any transformations for cylindrical.
        !
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




        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)
        end if

        !
        ! Compute velocity and pressure
        !
        vel1_m = mom1_m/density_m
        vel2_m = mom2_m/density_m
        vel3_m = mom3_m/density_m
        pressure_m = worker%get_field('Pressure', 'value', 'face interior')


        !
        ! Update average pressure
        !
        call self%compute_averages(worker,bc_COMM,vel1_avg, vel2_avg, vel3_avg, density_avg, pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)




        !
        ! Compute Fourier decomposition at set of radial stations: 
        !   : U_hat(nmodes,nradius)
        !
        call self%compute_fourier_decomposition(worker,bc_COMM,                        &
                                                density_hat_real,  density_hat_imag,   &
                                                vel1_hat_real,     vel1_hat_imag,      &
                                                vel2_hat_real,     vel2_hat_imag,      &
                                                vel3_hat_real,     vel3_hat_imag,      &
                                                pressure_hat_real, pressure_hat_imag,  &
                                                c1_hat_real,       c1_hat_imag,        &
                                                c2_hat_real,       c2_hat_imag,        &
                                                c3_hat_real,       c3_hat_imag,        &
                                                c4_hat_real,       c4_hat_imag,        &
                                                c5_hat_real,       c5_hat_imag)




        !
        ! Solve for c5 using nonreflecting condition
        !
        nmodes = size(density_hat_real,1)
        do iradius = 1,size(self%r)

            !
            ! Get average parts
            !
            density_bar_r  = density_hat_real( 1,iradius)
            vel1_bar_r     = vel1_hat_real(    1,iradius)
            vel2_bar_r     = vel2_hat_real(    1,iradius)
            vel3_bar_r     = vel3_hat_real(    1,iradius)
            pressure_bar_r = pressure_hat_real(1,iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)


            !
            ! The imaginary part of beta has already been accounted for in
            ! the expressions for A2 and A3
            !
            beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
            A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
            A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)

            A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
            A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)



            !
            ! Compute c5 according to nonreflecting condition
            !
            !   hat{c5} = A3*hat{c3}  -  A4*hat{c4}
            !
            !do imode = 2,(nmodes-1)/2 ! -1 here because the first mode is treated with 1D characteristics
            do imode = 2,nmodes ! -1 here because the first mode is treated with 1D characteristics

                c5_hat_real(imode,iradius) = (A3_real*c3_hat_real(imode,iradius) - A3_imag*c3_hat_imag(imode,iradius))  &   ! A3*c3 (real)
                                           - (A4_real*c4_hat_real(imode,iradius) - A4_imag*c4_hat_imag(imode,iradius))      ! A4*c4 (real)
                c5_hat_imag(imode,iradius) = (A3_imag*c3_hat_real(imode,iradius) + A3_real*c3_hat_imag(imode,iradius))  &   ! A3*c3 (imag)
                                           - (A4_imag*c4_hat_real(imode,iradius) + A4_real*c4_hat_imag(imode,iradius))      ! A4*c4 (imag)

            end do !imode


        end do !iradius


        !
        ! Interpolate c5 to the correct radial stations for quadrature nodes
        !
        coords = worker%coords()
        allocate(c5_hat_real_gq(size(c5_hat_real,1),size(coords)), &
                 c5_hat_imag_gq(size(c5_hat_imag,1),size(coords)), stat=ierr)
        if (ierr /= 0) call AllocationError
        c5_hat_real_gq = c5_hat_real(1,1)
        c5_hat_imag_gq = c5_hat_real(1,1)
        c5_hat_real_gq = ZERO   
        c5_hat_imag_gq = ZERO
        do igq = 1,size(coords)
            do imode = 2,nmodes ! not interpolating mode1, so it remains zero and isn't present in idft
                c5_hat_real_gq(imode,igq) = interpolate_linear_ad(self%r,c5_hat_real(imode,:),coords(igq)%c1_)
                c5_hat_imag_gq(imode,igq) = interpolate_linear_ad(self%r,c5_hat_imag(imode,:),coords(igq)%c1_)
            end do
        end do



        !
        ! Evaluate c5 at radius to correct theta
        !
        c5_3d = c5_hat_real_gq(1,:)
        c5_3d = ZERO
        do igq = 1,size(coords)
            theta_offset = coords(igq)%c2_ - self%theta_ref
            ! We include all modes here for generality, but we already set mode1 to zero
            ! so we are only getting the perturbation part.
            c5_3d(igq:igq) = idft_eval(c5_hat_real_gq(:,igq),c5_hat_imag_gq(:,igq),[theta_offset],pitch(1))
        end do






        !
        ! Handle perturbation from local radial mean (m /= 0)
        !
        c1_3d = c5_3d
        c2_3d = c5_3d
        c3_3d = c5_3d
        c4_3d = c5_3d
        c1_3d = ZERO
        c2_3d = ZERO
        c3_3d = ZERO
        c4_3d = ZERO
        density_bar  = c5_3d
        vel1_bar     = c5_3d
        vel2_bar     = c5_3d
        vel3_bar     = c5_3d
        pressure_bar = c5_3d
        density_bar  = ZERO
        vel1_bar     = ZERO
        vel2_bar     = ZERO
        vel3_bar     = ZERO
        pressure_bar = ZERO
        do igq = 1,size(coords)
            density_bar(igq)  = interpolate_linear_ad(self%r, density_hat_real( 1,:), coords(igq)%c1_)
            vel1_bar(igq)     = interpolate_linear_ad(self%r, vel1_hat_real(    1,:), coords(igq)%c1_)
            vel2_bar(igq)     = interpolate_linear_ad(self%r, vel2_hat_real(    1,:), coords(igq)%c1_)
            vel3_bar(igq)     = interpolate_linear_ad(self%r, vel3_hat_real(    1,:), coords(igq)%c1_)
            pressure_bar(igq) = interpolate_linear_ad(self%r, pressure_hat_real(1,:), coords(igq)%c1_)
        end do
        c_bar = sqrt(gam*pressure_bar/density_bar)





        !
        ! Compute perturbation from radius-local mean
        !
        ddensity  = density_m  - density_bar
        dvel1     = vel1_m     - vel1_bar
        dvel2     = vel2_m     - vel2_bar
        dvel3     = vel3_m     - vel3_bar
        dpressure = pressure_m - pressure_bar



        !
        ! Compute characteristics 1-4 associated with local perturbation. 
        ! 5 was already handled from nonreflecting condition on Fourier modes.
        !
        c1_3d = (-c_bar*c_bar)*ddensity    +  (ONE)*dpressure
        c2_3d = (density_bar*c_bar)*dvel1
        c3_3d = (density_bar*c_bar)*dvel2
        c4_3d = (density_bar*c_bar)*dvel3  +  (ONE)*dpressure


        !
        ! Handle m=0 perturbation
        !
        c1_1d = density_m
        c2_1d = density_m
        c3_1d = density_m
        c4_1d = density_m
        c5_1d = density_m
        c1_1d = ZERO
        c2_1d = ZERO
        c3_1d = ZERO
        c4_1d = ZERO
        c5_1d = ZERO



        !
        ! Compute 1-4 characteristics from extrapolation
        !
        ddensity  = density_bar  - density_avg 
        dvel1     = vel1_bar     - vel1_avg
        dvel2     = vel2_bar     - vel2_avg
        dvel3     = vel3_bar     - vel3_avg
        dpressure = pressure_bar - pressure_avg
        do igq = 1,size(ddensity)
            c1_1d(igq) = -c_avg*c_avg*ddensity(igq)    +  dpressure(igq)
            c2_1d(igq) = density_avg*c_avg*dvel1(igq)
            c3_1d(igq) = density_avg*c_avg*dvel2(igq)
            c4_1d(igq) = density_avg*c_avg*dvel3(igq)  +  dpressure(igq)
        end do


        !
        ! Compute characteristic 5 to achieve average pressure
        !
        c5_1d          = -TWO*(pressure_avg - p_user(1))
        ddensity_mean  =  c5_1d(1)/(TWO*c_avg*c_avg)
        dvel3_mean     = -c5_1d(1)/(TWO*density_avg*c_avg)
        dpressure_mean = HALF*c5_1d(1)

        dvel1_mean = dvel3_mean
        dvel2_mean = dvel3_mean
        dvel1_mean = ZERO
        dvel2_mean = ZERO



        !
        ! Contribute average part to boundary state
        !
        density_bc  = density_m
        vel1_bc     = density_m
        vel2_bc     = density_m
        vel3_bc     = density_m
        pressure_bc = density_m

        
        !
        ! Compose boundary state beginning with average
        !
        density_bc  = density_avg
        vel1_bc     = vel1_avg
        vel2_bc     = vel2_avg
        vel3_bc     = vel3_avg
        pressure_bc = pressure_avg

        
        !
        ! Contribute perturbation from boundary-global 1D characteristic update
        !
        do igq = 1,size(c1_1d)
            density_bc(igq)  = density_bc(igq)   +  (-ONE/(c_avg*c_avg))*c1_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c4_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c5_1d(igq)
            vel1_bc(igq)     = vel1_bc(igq)      +  (ONE/(density_avg*c_avg))*c2_1d(igq)
            vel2_bc(igq)     = vel2_bc(igq)      +  (ONE/(density_avg*c_avg))*c3_1d(igq)
            vel3_bc(igq)     = vel3_bc(igq)      +  (ONE/(TWO*density_avg*c_avg))*c4_1d(igq)  -  (ONE/(TWO*density_avg*c_avg))*c5_1d(igq)
            pressure_bc(igq) = pressure_bc(igq)  +  HALF*c4_1d(igq)  +  HALF*c5_1d(igq)
        end do


        !
        ! Contribute perturbation from perturbation about radius-local mean from 
        ! quasi-3d nrbc.
        !
        density_bc  = density_bc   +  (-ONE/(c_bar*c_bar))*c1_3d  +  (ONE/(TWO*c_bar*c_bar))*c4_3d  +  (ONE/(TWO*c_bar*c_bar))*c5_3d
        vel1_bc     = vel1_bc      +  (ONE/(density_bar*c_bar))*c2_3d
        vel2_bc     = vel2_bc      +  (ONE/(density_bar*c_bar))*c3_3d
        vel3_bc     = vel3_bc      +  (ONE/(TWO*density_bar*c_bar))*c4_3d  -  (ONE/(TWO*density_bar*c_bar))*c5_3d
        pressure_bc = pressure_bc  +  HALF*c4_3d  +  HALF*c5_3d


        !
        ! Form conserved variables
        !
        density_bc = density_bc
        mom1_bc    = density_bc*vel1_bc
        mom2_bc    = density_bc*vel2_bc
        mom3_bc    = density_bc*vel3_bc
        energy_bc  = pressure_bc/(gam - ONE)  + (density_bc*HALF)*(vel1_bc*vel1_bc + vel2_bc*vel2_bc + vel3_bc*vel3_bc)


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        !
        ! Store boundary condition state. Q_bc
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')




    end subroutine compute_bc_state
    !*********************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_fourier_decomposition(self,worker,bc_comm,                   &
                                             density_hat_real,  density_hat_imag,   &
                                             vel1_hat_real,     vel1_hat_imag,      &
                                             vel2_hat_real,     vel2_hat_imag,      &
                                             vel3_hat_real,     vel3_hat_imag,      &
                                             pressure_hat_real, pressure_hat_imag,  &
                                             c1_hat_real,       c1_hat_imag,        &
                                             c2_hat_real,       c2_hat_imag,        &
                                             c3_hat_real,       c3_hat_imag,        &
                                             c4_hat_real,       c4_hat_imag,        &
                                             c5_hat_real,       c5_hat_imag)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        type(mpi_comm),                                 intent(in)      :: bc_comm
        type(AD_D),     allocatable,                    intent(inout)   :: density_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: density_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel1_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel1_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel2_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel2_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel3_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: vel3_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: pressure_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: pressure_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c1_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c1_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c2_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c2_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c3_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c3_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c4_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c4_hat_imag(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c5_hat_real(:,:)
        type(AD_D),     allocatable,                    intent(inout)   :: c5_hat_imag(:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp,   &
            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
            c1,         c2,         c3,         c4,         c5,                                 &
            ddensity,   dvel1,      dvel2,      dvel3,      dpressure

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik)             :: nmodes, nradius, ntheta, iradius, itheta, ncoeff, ierr
        real(rk)                :: dtheta, dtheta_n, theta, z, midpoint(3)
        real(rk),   allocatable :: pitch(:), physical_nodes(:,:)



        !
        ! Determine z-location of current face and assume entire boundary is constant-z
        !
        if (worker%iface == XI_MIN) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([-ONE,ZERO,ZERO],'Deformed')
        else if (worker%iface == XI_MAX) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([ ONE,ZERO,ZERO],'Deformed')
        else if (worker%iface == ETA_MIN) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([ZERO,-ONE,ZERO],'Deformed')
        else if (worker%iface == ETA_MAX) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([ZERO, ONE,ZERO],'Deformed')
        else if (worker%iface == ZETA_MIN) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([ZERO,ZERO,-ONE],'Deformed')
        else if (worker%iface == ZETA_MAX) then
            midpoint = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%physical_point([ZERO,ZERO, ONE],'Deformed')
        end if
        z = midpoint(3)


        !
        ! Define Fourier discretization
        !
        nmodes  = 5
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff
        

        !
        ! Allocate interpolation nodes
        !
        allocate(physical_nodes(ntheta,3), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Allocate storage in result
        !
        allocate(density_hat_real( ncoeff,nradius), density_hat_imag( ncoeff,nradius),  &
                 vel1_hat_real(    ncoeff,nradius), vel1_hat_imag(    ncoeff,nradius),  &
                 vel2_hat_real(    ncoeff,nradius), vel2_hat_imag(    ncoeff,nradius),  &
                 vel3_hat_real(    ncoeff,nradius), vel3_hat_imag(    ncoeff,nradius),  &
                 pressure_hat_real(ncoeff,nradius), pressure_hat_imag(ncoeff,nradius),  &
                 c1_hat_real(      ncoeff,nradius), c1_hat_imag(      ncoeff,nradius),  &
                 c2_hat_real(      ncoeff,nradius), c2_hat_imag(      ncoeff,nradius),  &
                 c3_hat_real(      ncoeff,nradius), c3_hat_imag(      ncoeff,nradius),  &
                 c4_hat_real(      ncoeff,nradius), c4_hat_imag(      ncoeff,nradius),  &
                 c5_hat_real(      ncoeff,nradius), c5_hat_imag(      ncoeff,nradius), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Initialize theta discretization parameters
        !
        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta



        !
        ! Perform Fourier decomposition at each radial station.
        !
        do iradius = 1,nradius

            !
            ! Construct theta discretization
            !
            do itheta = 1,ntheta
                theta = self%theta_ref + (itheta-1)*dtheta_n
                physical_nodes(itheta,:) = [self%r(iradius),theta,z]
            end do


            !
            ! Interpolate solution to physical_nodes at current radial station
            !
            !density = worker%interpolate_field_general('Density',    physical_nodes, try_offset=[ZERO,-pitch(1),ZERO])
            !mom1    = worker%interpolate_field_general('Momentum-1', physical_nodes, try_offset=[ZERO,-pitch(1),ZERO])
            !mom2    = worker%interpolate_field_general('Momentum-2', physical_nodes, try_offset=[ZERO,-pitch(1),ZERO])
            !mom3    = worker%interpolate_field_general('Momentum-3', physical_nodes, try_offset=[ZERO,-pitch(1),ZERO])
            !energy  = worker%interpolate_field_general('Energy',     physical_nodes, try_offset=[ZERO,-pitch(1),ZERO])
            density = worker%interpolate_field_general('Density',    physical_nodes, donors=self%donor(iradius,:), donor_coords=self%donor_coord(iradius,:,:))
            mom1    = worker%interpolate_field_general('Momentum-1', physical_nodes, donors=self%donor(iradius,:), donor_coords=self%donor_coord(iradius,:,:))
            mom2    = worker%interpolate_field_general('Momentum-2', physical_nodes, donors=self%donor(iradius,:), donor_coords=self%donor_coord(iradius,:,:))
            mom3    = worker%interpolate_field_general('Momentum-3', physical_nodes, donors=self%donor(iradius,:), donor_coords=self%donor_coord(iradius,:,:))
            energy  = worker%interpolate_field_general('Energy',     physical_nodes, donors=self%donor(iradius,:), donor_coords=self%donor_coord(iradius,:,:))

            if (worker%coordinate_system() == 'Cylindrical') then
                mom2 = mom2/self%r(iradius)  ! convert to tangential momentum
            end if



            !
            ! Compute velocities and pressure
            !
            vel1 = mom1/density
            vel2 = mom2/density
            vel3 = mom3/density
            pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )



            !
            ! Compute Fourier transform
            !
            call dft(density,  density_real_tmp,  density_imag_tmp )
            call dft(vel1,     vel1_real_tmp,     vel1_imag_tmp    )
            call dft(vel2,     vel2_real_tmp,     vel2_imag_tmp    )
            call dft(vel3,     vel3_real_tmp,     vel3_imag_tmp    )
            call dft(pressure, pressure_real_tmp, pressure_imag_tmp)




            density_hat_real( :,iradius) = density_real_tmp
            density_hat_imag( :,iradius) = density_imag_tmp
            vel1_hat_real(    :,iradius) = vel1_real_tmp
            vel1_hat_imag(    :,iradius) = vel1_imag_tmp
            vel2_hat_real(    :,iradius) = vel2_real_tmp
            vel2_hat_imag(    :,iradius) = vel2_imag_tmp
            vel3_hat_real(    :,iradius) = vel3_real_tmp
            vel3_hat_imag(    :,iradius) = vel3_imag_tmp
            pressure_hat_real(:,iradius) = pressure_real_tmp
            pressure_hat_imag(:,iradius) = pressure_imag_tmp

            
            !
            ! Get average term
            !
            density_bar  = density_hat_real( 1,iradius)
            vel1_bar     = vel1_hat_real(    1,iradius)
            vel2_bar     = vel2_hat_real(    1,iradius)
            vel3_bar     = vel3_hat_real(    1,iradius)
            pressure_bar = pressure_hat_real(1,iradius)
            c_bar = sqrt(gam*pressure_bar/density_bar)


            !
            ! Compute perturbation
            !
            ddensity  = density  - density_bar
            dvel1     = vel1     - vel1_bar
            dvel2     = vel2     - vel2_bar
            dvel3     = vel3     - vel3_bar
            dpressure = pressure - pressure_bar


            !
            ! Convert perturbation to 1D characteristics
            !
            c1 = ddensity
            c2 = ddensity
            c3 = ddensity
            c4 = ddensity
            c5 = ddensity
            do itheta = 1,size(ddensity)
                c1(itheta) = -(c_bar*c_bar)*ddensity(itheta)    +  (ONE)*dpressure(itheta)
                c2(itheta) = (density_bar*c_bar)*dvel1(itheta)
                c3(itheta) = (density_bar*c_bar)*dvel2(itheta)
                c4(itheta) = (density_bar*c_bar)*dvel3(itheta)  +  (ONE)*dpressure(itheta)
                c5(itheta) = -(density_bar*c_bar)*dvel3(itheta) +  (ONE)*dpressure(itheta)
            end do


            !
            ! Compute Fourier transform of characteristic variables
            !
            call dft(c1, c1_real_tmp, c1_imag_tmp)
            call dft(c2, c2_real_tmp, c2_imag_tmp)
            call dft(c3, c3_real_tmp, c3_imag_tmp)
            call dft(c4, c4_real_tmp, c4_imag_tmp)
            call dft(c5, c5_real_tmp, c5_imag_tmp)


            c1_hat_real(:,iradius) = c1_real_tmp
            c2_hat_real(:,iradius) = c2_real_tmp
            c3_hat_real(:,iradius) = c3_real_tmp
            c4_hat_real(:,iradius) = c4_real_tmp
            c5_hat_real(:,iradius) = c5_real_tmp

            c1_hat_imag(:,iradius) = c1_imag_tmp
            c2_hat_imag(:,iradius) = c2_imag_tmp
            c3_hat_imag(:,iradius) = c3_imag_tmp
            c4_hat_imag(:,iradius) = c4_imag_tmp
            c5_hat_imag(:,iradius) = c5_imag_tmp

        end do !iradius



    end subroutine compute_fourier_decomposition
    !*********************************************************************************








    !>  Determine rmin, rmax, and average theta.
    !!
    !!  Initialize radial stations and reference theta for Fourier transform.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine analyze_bc_geometry(self,mesh,group_ID,bc_comm)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(mesh_t),                                   intent(in)      :: mesh
        integer(ik),                                    intent(in)      :: group_ID
        type(mpi_comm),                                 intent(in)      :: bc_comm

        integer(ik) :: bc_idomain_l, bc_ielement_l, patch_ID, face_ID, iface, ierr, inode
        real(rk)    :: face_rmin, face_rmax, local_rmin, local_rmax, global_rmin, global_rmax,  &
                       face_thetamin, face_thetamax, local_thetamin, local_thetamax, global_thetamin, global_thetamax
        real(rk)    :: ref_nodes(8,3), physical_nodes(8,3)


        !
        ! Search for min/max radius on local processor
        !
        local_rmin     =  HUGE(1._rk) ! any radius will be smaller than this, so it is guarunteed to be reset.
        local_rmax     = -HUGE(1._rk) ! any radius will be larger than this, so it is guarunteed to be reset.
        local_thetamin =  HUGE(1._rk) ! any theta will be smaller than this, so it is guarunteed to be reset.
        local_thetamax = -HUGE(1._rk) ! any theta will be larger than this, so it is guarunteed to be reset.
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                iface = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                ! Pick points to evaluate coordinates
                if (iface == XI_MIN) then
                    ref_nodes(1,:) = [-ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [-ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [-ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [-ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [-ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [-ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [-ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [-ONE, ZERO, -ONE]
                else if (iface == XI_MAX) then
                    ref_nodes(1,:) = [ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [ONE, ZERO, -ONE]

                else if (iface == ETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, -ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, -ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, -ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, ONE, -ONE]

                else if (iface == ZETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, -ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, -ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, -ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, -ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, -ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ZETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, -ONE, ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, ONE]

                else
                    call chidg_signal(FATAL,"outlet_3dgiles_innerproduct_general: analyze_bc_geometry, invalid face indec.")
                end if



                ! Evaluate physical coordinates on face edges
                bc_idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                bc_ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                do inode = 1,size(ref_nodes,1)
                    physical_nodes(inode,:) = mesh%domain(bc_idomain_l)%elems(bc_ielement_l)%physical_point(ref_nodes(inode,:),'Deformed')
                end do !inode

                ! Get face min/max radius
                face_rmin     = minval(physical_nodes(:,1))
                face_rmax     = maxval(physical_nodes(:,1))
                face_thetamin = minval(physical_nodes(:,2))
                face_thetamax = maxval(physical_nodes(:,2))

                ! Update processor-local value if new min/max values were found on the face
                if (face_rmin < local_rmin)         local_rmin     = face_rmin
                if (face_rmax > local_rmax)         local_rmax     = face_rmax
                if (face_thetamin < local_thetamin) local_thetamin = face_thetamin
                if (face_thetamax > local_thetamax) local_thetamax = face_thetamax

            end do !face_ID
        end do !patch_ID


        ! Reduce processor local values to determine boundary-global min/max values
        call MPI_AllReduce(local_rmin,    global_rmin,    1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_rmax,    global_rmax,    1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        call MPI_AllReduce(local_thetamin,global_thetamin,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_thetamax,global_thetamax,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        

        ! Create radial stations
        self%r = linspace(global_rmin,global_rmax,10)

        ! Compute theta_ref
        self%theta_ref = (global_thetamin + global_thetamax)/TWO

    end subroutine analyze_bc_geometry
    !********************************************************************************







    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine initialize_fourier_discretization(self,mesh,group_ID,bc_comm)
        class(outlet_3dgiles_innerproduct_general_t),   intent(inout)   :: self
        type(mesh_t),                                   intent(in)      :: mesh
        integer(ik),                                    intent(in)      :: group_ID
        type(mpi_comm),                                 intent(in)      :: bc_comm

        integer(ik)                 :: nmodes, ncoeff, nradius, ntheta, idomain_l, ielement_l, iface, iradius, itheta, ierr
        real(rk)                    :: dtheta, dtheta_n, midpoint(3), try_offset(3), node(3), z
        real(rk),       allocatable :: pitch(:)
        character(:),   allocatable :: user_msg
        logical                     :: donor_found



        !
        ! Determine z-location of some face on the boundary and assume 
        ! entire boundary is constant-z
        !
        idomain_l  = mesh%bc_patch_group(group_ID)%patch(1)%idomain_l()
        ielement_l = mesh%bc_patch_group(group_ID)%patch(1)%ielement_l(1)
        iface      = mesh%bc_patch_group(group_ID)%patch(1)%iface(1)
        if (iface == XI_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([-ONE,ZERO,ZERO],'Deformed')
        else if (iface == XI_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ ONE,ZERO,ZERO],'Deformed')
        else if (iface == ETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,-ONE,ZERO],'Deformed')
        else if (iface == ETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO, ONE,ZERO],'Deformed')
        else if (iface == ZETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO,-ONE],'Deformed')
        else if (iface == ZETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO, ONE],'Deformed')
        end if
        z = midpoint(3)





        !
        ! Define Fourier discretization
        !
        nmodes  = 5
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff
        


        !
        ! Initialize theta discretization parameters
        !
        !pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta


        !
        ! Construct theta discretization at each radius
        !
        allocate(self%theta(ntheta), stat=ierr)
        if (ierr /= 0) call AllocationError
        do itheta = 1,ntheta
            self%theta(itheta) = self%theta_ref + (itheta-1)*dtheta_n
        end do


        !
        ! Donor search offset, if needed
        !
        try_offset = [ZERO, -pitch(1), ZERO]


        !
        ! For each radial station, initialized donor for each node in theta grid
        !
        allocate(self%donor(nradius,ntheta), self%donor_coord(nradius,ntheta,3), stat=ierr)
        if (ierr /= 0) call AllocationError
        do iradius = 1,size(self%r)
            do itheta = 1,size(self%theta)

                node = [self%r(iradius), self%theta(itheta), z]

                !
                ! Try processor-LOCAL elements
                !
                call find_gq_donor(mesh,                                    &
                                   node,                                    &
                                   [ZERO,ZERO,ZERO],                        &
                                   face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                   self%donor(iradius,itheta),              &
                                   self%donor_coord(iradius,itheta,1:3),    &
                                   donor_found)

                !
                ! Try LOCAL elements with try_offset if still not found 
                !
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                    &
                                       node,                                    &
                                       try_offset,                              &
                                       face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                       self%donor(iradius,itheta),              &
                                       self%donor_coord(iradius,itheta,1:3),    &
                                       donor_found)

                end if



                !
                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                !
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                [ZERO,ZERO,ZERO],                       &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_coord(iradius,itheta,1:3),   &
                                                donor_found)
                end if

                
                !
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                !
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                try_offset,                             &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_coord(iradius,itheta,1:3),   &
                                                donor_found)
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_state_outlet_3dgiles_innerproduct_general%initialize_fourier_discretization: &
                            no donor element found for Fourier discretization node."
                if (.not. donor_found) call chidg_signal(FATAL,user_msg)


            end do !itheta
        end do !iradius



    end subroutine initialize_fourier_discretization
    !**************************************************************************************












end module bc_state_outlet_3dgiles_innerproduct_general
