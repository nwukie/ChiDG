module bc_state_outlet_LODI_pressure
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO, NO_PROC, ME

    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use type_ivector,           only: ivector_t
    use type_face_info,         only: face_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_SUM, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast
    use DNAD_D
    implicit none





    !>  Name: Outlet - LODI Pressure
    !!      : Update average pressure using LODI with transverse terms
    !!      : Extrapolate other characteristics
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!  Behavior:
    !!      
    !!  References:
    !!      [1] Koupper et al."Compatibility of Characteristic Boundary Conditions wth 
    !!                         Radial Equilibrium in Turbomachinery Simulation."
    !!                         AIAA Journal, Vol. 52, No. 12, December 2014.
    !!
    !!      [2] Granet et al. "Comparison of Nonreflecting Outlet Boundary Conditions for 
    !!                         Compressible Solvers on Unstructured Grids."
    !!                         AIAA Journal, Vol. 48, No. 10, October 2010.
    !!
    !!      [3] Yoo et al. "Characteristic boundary conditions for direct simulations of
    !!                      turbulent counterflow flames." 
    !!                      Combustion Theory and Modelling, Vol. 9, No. 4, November 2005, 
    !!                      pp. 617-646.
    !!              
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   4/20/2017
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_LODI_pressure_t


    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! boundary condition function implementation

        procedure   :: compute_averages

    end type outlet_LODI_pressure_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_LODI_pressure_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name('Outlet - LODI Pressure')
        call self%set_family('Outlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('Average Pressure','Required')


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
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
        class(outlet_LODI_pressure_t),  intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        integer(ik),                    intent(in)      :: group_ID
        type(mpi_comm),                 intent(in)      :: bc_COMM

        integer(ik) :: patch_ID, face_ID, elem_ID, patch_ID_coupled, face_ID_coupled,   &
                       idomain_g, idomain_l, ielement_g, ielement_l, iface,             &
                       bc_IRANK, bc_NRANK, ierr, iproc, nbc_elements,     &
                       ielem, neqns, nterms_s, ngq, ibc

        integer(ik) :: idomain_g_coupled, idomain_l_coupled, ielement_g_coupled, ielement_l_coupled, &
                       iface_coupled, proc_coupled

        !type(point_t),  allocatable :: quad_pts(:)
        real(rk),       allocatable :: quad_pts(:,:)
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
                        quad_pts   = mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts



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
                                                                                                             quad_pts)


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

                        ngq = size(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts,1)
                        call MPI_Bcast(ngq,                                                                 1, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%differential_areas, ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:,1),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:,2),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:,3),      ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:)%c1_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:)%c2_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%quad_pts(:)%c3_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)

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
                    if (allocated(areas) ) deallocate(areas, quad_pts)
                    allocate(areas(ngq), quad_pts(ngq,3), stat=ierr)
                    if (ierr /= 0) call AllocationError


                    call MPI_BCast(areas,         ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(quad_pts(:,1), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(quad_pts(:,2), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    call MPI_BCast(quad_pts(:,3), ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    !call MPI_BCast(quad_pts(:)%c1_, ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    !call MPI_BCast(quad_pts(:)%c2_, ngq, MPI_REAL8, iproc, bc_COMM, ierr)
                    !call MPI_BCast(quad_pts(:)%c3_, ngq, MPI_REAL8, iproc, bc_COMM, ierr)


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
                                                                                                                 quad_pts)





                        end do ! face_ID
                    end do ! patch_ID

                end do !ielem




            end if




            call MPI_Barrier(bc_COMM,ierr)
        end do







    end subroutine init_bc_coupling
    !******************************************************************************************






    !>  Update the area-averaged pressure for the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, u_avg, v_avg, w_avg, density_avg, p_avg)
        class(outlet_LODI_pressure_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_COMM
        type(AD_D),                     intent(inout)   :: u_avg
        type(AD_D),                     intent(inout)   :: v_avg
        type(AD_D),                     intent(inout)   :: w_avg
        type(AD_D),                     intent(inout)   :: density_avg
        type(AD_D),                     intent(inout)   :: p_avg

        type(AD_D)          :: face_p, face_u, face_v, face_w, face_density, p_integral, u_integral, v_integral, w_integral, density_integral
        type(face_info_t)   :: face_info

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom_1, mom_2, mom_3, energy, p,    &
            u, v, w, c, M, vmag

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r

        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
                       ielement_l_coupled, iface_coupled
        real(rk)    :: face_area, total_area, gam


        gam = 1.4_rk

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
            mom_1   = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom1,    itime, 'value', ME)
            mom_2   = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom2,    itime, 'value', ME)
            mom_3   = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, imom3,    itime, 'value', ME)
            energy  = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, ienergy,  itime, 'value', ME)

            if (worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%coordinate_system == 'Cylindrical') then
                !mom_2 = mom_2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%quad_pts(:)%c1_
                mom_2 = mom_2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%quad_pts(:,1)
            end if


            
            !
            ! Compute velocity
            !
            u = mom_1 / density
            v = mom_2 / density
            w = mom_3 / density


            !
            ! Compute pressure over the face
            !
            p = (gam - ONE)*(energy - HALF*density*(u*u + v*v + w*w))


            !
            ! Compute speed of sound and Mach
            !
            c = sqrt(gam * p / density)

            
            !
            ! Compute Mach number
            !
            vmag = sqrt(u*u + v*v + w*w)
            M = vmag/c


            !
            ! Get weights + areas
            !
            !weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%gq%face%weights(:,iface_coupled)
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%ref_s%weights(iface_coupled)
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

            if (allocated(p_integral%xp_ad_)) then
                p_integral = p_integral + face_p
            else
                p_integral = face_p
            end if


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







    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_LODI_pressure_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,                            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,                           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_bc, v_bc, w_bc, p_bc,                                                     &
            u_m,  v_m,  w_m,  p_m, invdensity,                                          &
            du,   dv,   dw,   dp,  ddensity,  c,                                        &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,                      &
            du_ddensity, du_dmom1,                                                      &
            dv_ddensity, dv_dmom2,                                                      &
            dw_ddensity, dw_dmom3,                                                      &
            grad1_p, grad2_p, grad3_p,                                                  &
            grad1_u, grad2_u, grad3_u,                                                  &
            grad1_v, grad2_v, grad3_v,                                                  &
            grad1_w, grad2_w, grad3_w,                                                  &
            lamda_1, lamda_234, lamda_5,                                                &
            L1, L2, L3, L4, L5,                                                         &
            T1, T2, T3, T4, T5, M_avg_array, c1_p, c2_p, c3_p, c4_p, c5_p, u_p, v_p, w_p, density_p, p_p, &
            c_bar_gq, p_bar_gq, u_bar_gq, v_bar_gq, w_bar_gq, density_bar_gq

        type(AD_D)  :: c_bar, p_bar, u_bar, v_bar, w_bar, density_bar, c1_bar, c2_bar, c3_bar, c4_bar, c5_bar


        logical                                     :: face_has_node
        real(rk)                                    :: gam, beta, k
        real(rk),       allocatable, dimension(:)   ::  &
            p_user, norm_1, norm_2, norm_3, r


        !
        ! Get back pressure from function.
        !
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())

        !
        ! Compute gamma
        !
        gam = 1.4_rk


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')



        grad1_density_m = worker%get_primary_field_face('Density'   , 'grad1', 'face interior')
        grad2_density_m = worker%get_primary_field_face('Density'   , 'grad2', 'face interior')
        grad3_density_m = worker%get_primary_field_face('Density'   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_primary_field_face('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_primary_field_face('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_primary_field_face('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_primary_field_face('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_primary_field_face('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m    = worker%get_primary_field_face('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_primary_field_face('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_primary_field_face('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_primary_field_face('Momentum-3', 'grad3', 'face interior')
        
        grad1_energy_m  = worker%get_primary_field_face('Energy'    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_primary_field_face('Energy'    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_primary_field_face('Energy'    , 'grad3', 'face interior')



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if



        !
        ! Update average pressure
        !
        !p_bar = self%compute_average_pressure(worker,bc_COMM)
        call self%compute_averages(worker,bc_COMM,u_bar, v_bar, w_bar, density_bar, p_bar)


        print*, 'Average pressure: ', p_bar%x_ad_
!        print*, 'Average Mach: ', M_bar%x_ad_

        
        !
        ! Compute updates due to average 
        !
        !vmag_bar = sqrt(u_bar*u_bar + v_bar*v_bar + w_bar*w_bar)
        c_bar = sqrt(gam * p_bar / density_bar)



        !
        ! Compute characteristics due to average quantities
        !
        c1_bar = -(c_bar*c_bar)*density_bar    +  p_bar
        c2_bar =  (density_bar*c_bar) * u_bar  +  p_bar
        c3_bar =  (density_bar*c_bar) * v_bar
        c4_bar =  (density_bar*c_bar) * w_bar
        c5_bar = -(density_bar*c_bar) * u_bar  +  p_bar


        !
        ! Recompute incoming average characteristic
        !
        c5_bar = -TWO*(p_bar - p_user(1))




        !
        ! Compute boundary primitive variables from modified average characteristics
        !   1: initialize arrays at all quadrature nodes
        !   2: set to average values
        !
        density_bc = density_m
        u_bc       = density_m
        v_bc       = density_m
        w_bc       = density_m
        p_bc       = density_m

        density_bc = -(ONE/(c_bar*c_bar))*c1_bar  +  (ONE/(TWO*c_bar*c_bar))*c2_bar  +  (ONE/(TWO*c_bar*c_bar))*c5_bar
        u_bc       = (ONE/(TWO*density_bar*c_bar))*c2_bar  -  (ONE/(TWO*density_bar*c_bar))*c5_bar
        v_bc       = (ONE/(density_bar*c_bar))*c3_bar
        w_bc       = (ONE/(density_bar*c_bar))*c4_bar
        p_bc       = HALF*c2_bar  +  HALF*c5_bar






        !
        ! Compute velocities
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m


        !
        ! Compute pressure from extrapolated data
        !
        p_m = (gam-ONE)*(energy_m - HALF*( (mom1_m*mom1_m) + (mom2_m*mom2_m) + (mom3_m*mom3_m) )/density_m )
    




        !
        ! Compute perturbation from mean
        !
        u_p       = u_m       - u_bar
        v_p       = v_m       - v_bar
        w_p       = w_m       - w_bar
        p_p       = p_m       - p_bar
        density_p = density_m - density_bar


        c_bar_gq       = density_m
        u_bar_gq       = density_m
        v_bar_gq       = density_m
        w_bar_gq       = density_m
        p_bar_gq       = density_m
        density_bar_gq = density_m

        c_bar_gq = c_bar
        u_bar_gq = u_bar
        v_bar_gq = v_bar
        w_bar_gq = w_bar
        p_bar_gq = p_bar
        density_bar_gq = density_bar


        !
        ! Compute characteristics due to perturbation quantities
        !
        c1_p = -(c_bar_gq*c_bar_gq)*density_p     +   p_p
        c2_p =  (density_bar_gq*c_bar_gq) * u_p   +   p_p
        c3_p =  (density_bar_gq*c_bar_gq) * v_p
        c4_p =  (density_bar_gq*c_bar_gq) * w_p
        c5_p = -(density_bar_gq*c_bar_gq) * u_p   +   p_p


!        ! Cancel incoming perturbations
!        c5_p = ZERO
!
!
!        density_bc = density_bc + (-ONE/(c_bar_gq*c_bar_gq))*c1_p  +  (ONE/(TWO*c_bar_gq*c_bar_gq))*c2_p  +  (ONE/(TWO*c_bar_gq*c_bar_gq))*c5_p
!        u_bc       = u_bc       + (ONE/(TWO*density_bar_gq*c_bar_gq))*c2_p  -  (ONE/(TWO*density_bar_gq*c_bar_gq))*c5_p
!        v_bc       = v_bc       + (ONE/(density_bar_gq*c_bar_gq))*c3_p
!        w_bc       = w_bc       + (ONE/(density_bar_gq*c_bar_gq))*c4_p
!        p_bc       = p_bc       + HALF*c2_p  +  HALF*c5_p







        !
        ! LODI+transverse updates
        !   1: compute grad(p)
        !   2: compute grad(u)
        !   3: compute L, T
        !   4: compute update in primitive variables dU
        !   5: compute new values U_new = U_old + dU
        !   6: convert to conservative variables and store
        !


        !
        ! Compute pressure jacobians
        !
        dp_ddensity =  (gam - ONE)*HALF*(u_m*u_m + v_m*v_m + w_m*w_m)
        dp_dmom1    = -(gam - ONE)*u_m
        dp_dmom2    = -(gam - ONE)*v_m
        dp_dmom3    = -(gam - ONE)*w_m
        dp_denergy  = dp_dmom3 ! initialize derivatives
        dp_denergy  =  (gam - ONE)


        !
        ! Compute pressure gradient via chain rule:
        !
        grad1_p = dp_ddensity*grad1_density_m  +  dp_dmom1*grad1_mom1_m  +  dp_dmom2*grad1_mom2_m  +  dp_dmom3*grad1_mom3_m  +  dp_denergy*grad1_energy_m
        grad2_p = dp_ddensity*grad2_density_m  +  dp_dmom1*grad2_mom1_m  +  dp_dmom2*grad2_mom2_m  +  dp_dmom3*grad2_mom3_m  +  dp_denergy*grad2_energy_m
        grad3_p = dp_ddensity*grad3_density_m  +  dp_dmom1*grad3_mom1_m  +  dp_dmom2*grad3_mom2_m  +  dp_dmom3*grad3_mom3_m  +  dp_denergy*grad3_energy_m


        !
        ! compute velocity jacobians
        !
        invdensity  = ONE/density_m
        du_ddensity = -invdensity*invdensity*mom1_m
        dv_ddensity = -invdensity*invdensity*mom2_m
        dw_ddensity = -invdensity*invdensity*mom3_m

        du_dmom1 = invdensity
        dv_dmom2 = invdensity
        dw_dmom3 = invdensity


        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = du_ddensity*grad1_density_m  +  du_dmom1*grad1_mom1_m
        grad2_u = du_ddensity*grad2_density_m  +  du_dmom1*grad2_mom1_m
        grad3_u = du_ddensity*grad3_density_m  +  du_dmom1*grad3_mom1_m

        grad1_v = dv_ddensity*grad1_density_m  +  dv_dmom2*grad1_mom2_m
        grad2_v = dv_ddensity*grad2_density_m  +  dv_dmom2*grad2_mom2_m
        grad3_v = dv_ddensity*grad3_density_m  +  dv_dmom2*grad3_mom2_m

        grad1_w = dw_ddensity*grad1_density_m  +  dw_dmom3*grad1_mom3_m
        grad2_w = dw_ddensity*grad2_density_m  +  dw_dmom3*grad2_mom3_m
        grad3_w = dw_ddensity*grad3_density_m  +  dw_dmom3*grad3_mom3_m




!        !
!        ! Compute wave speeds. u_m because assumption is that boundary is normal to x-axis
!        !
!        c = sqrt(gam * p_m / density_m)
!        lamda_1   = u_bar - c_bar
!        lamda_234 = u_bar
!        lamda_5   = u_bar + c_bar
!
!
!
!        !
!        ! Compute 
!        !
!
!
!
!
!        !
!        ! Compute amplitudes of axial waves
!        !
!        L2 =      lamda_234*(grad1_density_m  -  grad1_p/(c_bar*c_bar))
!        L3 =      lamda_234*(grad1_v)
!        L4 =      lamda_234*(grad1_w)
!        L5 = HALF*lamda_5*(grad1_p + density_m*c*grad1_u)
!
!
!        !
!        ! Compute amplitudes of transverse waves
!        !
!        T1 = ( (v_m*grad2_p  + w_m*grad3_p )   +  gam*p_m*(grad2_v + grad3_w)  -  density_m*c*(v_m*grad2_u  +  w_m*grad3_u) )/TWO
!        T2 =   (grad2_mom2_m + grad3_mom3_m)   - (gam*p_m*(grad2_v + grad3_w)  +              (v_m*grad2_p  +  w_m*grad3_p) )/(c*c)
!        T3 = ( (v_m*grad2_v  + w_m*grad3_v )   +  invdensity*grad2_p )
!        T4 = ( (v_m*grad2_w  + w_m*grad3_w )   +  invdensity*grad3_p )
!        T5 = ( (v_m*grad2_p  + w_m*grad3_p )   +  gam*p_m*(grad2_v + grad3_w)  +  density_m*c*(v_m*grad2_u  +  w_m*grad3_u) )/TWO
!
!        !T1 = ZERO
!        !T2 = ZERO
!        !T3 = ZERO
!        !T4 = ZERO
!        !T5 = ZERO
!
!
!
!                                           
!                                                                                               
!
!        !
!        ! Compute incoming axial wave
!        !
!        k    = 10._rk
!        beta = 0.2_rk   ! best investigation yet says this should be the mean Mach number across the face
!
!        !M_avg_array = density_m
!        !M_avg_array = M_avg
!        !L1 = k*(p_avg - p_user) + (beta - ONE)*T1
!        !L1 = k*(p_avg - p_user) + (M_avg_array - ONE)*T1
!        L1 = k*(p_bar - p_user) - T1
!
!        !print*, 'L1: ', L1(1)%xp_ad_
!        !if (any(T5(1)%xp_ad_ > 0.1)) then
!        !    print*, 'L5: ', L5(1)%xp_ad_
!        !    print*, 'T5: ', T5(1)%xp_ad_
!        !end if
!
!        !
!        ! Compute perturbation in primitive variables due to the characteristics
!        !
!        du       = -(L5 - L1)/(density_m*c)  -  (T5 - T1)/(density_m*c)
!        dv       = -(L3)                     -  (T3)
!        dw       = -(L4)                     -  (T4)
!        ddensity = -(L2 + (L5+L1)/(c*c))     -  (T2 + (T5+T1)/(c*c))
!        dp       = -(L5 + L1)                -  (T5 + T1)
!
!        !du       = -(L5 - L1)/(density_m*c)  
!        !dv       = -(L3)                     
!        !dw       = -(L4)                     
!        !ddensity = -(L2 + (L5+L1)/(c*c))     
!        !dp       = -(L5 + L1)                
!
!
!
!        !
!        ! Update primitive variables and form conservative boundary values
!        !
!        u_bc       = u_m       + du
!        v_bc       = v_m       + dv
!        w_bc       = w_m       + dw
!        density_bc = density_m + ddensity
!        p_bc       = p_m       + dp



        !
        ! Form conserved variables
        !
        density_bc = density_bc
        mom1_bc    = density_bc*u_bc
        mom2_bc    = density_bc*v_bc
        mom3_bc    = density_bc*w_bc
        energy_bc  = p_bc/(gam - ONE)  + (density_bc*HALF)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)




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
    !**********************************************************************************************






end module bc_state_outlet_LODI_pressure
