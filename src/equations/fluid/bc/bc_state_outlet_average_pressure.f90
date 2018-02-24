module bc_state_outlet_average_pressure
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL
    use mod_fluid,              only: gam

    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_SUM, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast
    use DNAD_D
    implicit none





    !>  Name: Outlet - Average Pressure
    !!
    !!  Handle perturbations from average using one-dimensional characteristic analysis.
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
    !!  @date   2/23/2018
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_average_pressure_t


    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! boundary condition function implementation
        procedure   :: compute_averages

    end type outlet_average_pressure_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_average_pressure_t),   intent(inout) :: self
        
        ! Set name, family
        call self%set_name('Outlet - Average Pressure')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure','Required')

    end subroutine init
    !********************************************************************************




    !>  Initialize boundary group coupling.
    !!
    !!  Each element is coupled with every other element that belongs to the boundary
    !!  condition. This coupling occurs because each face uses an
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
        class(outlet_average_pressure_t),  intent(inout)   :: self
        type(mesh_t),                      intent(inout)   :: mesh
        integer(ik),                       intent(in)      :: group_ID
        type(mpi_comm),                    intent(in)      :: bc_COMM

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
                        interp_coords_def   = mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def



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


                    call MPI_BCast(areas,           ngq, MPI_REAL8, iproc, bc_COMM, ierr)
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
    !******************************************************************************************






    !>  Compute averaged quantities over the face. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, u_avg, v_avg, w_avg, density_avg, p_avg)
        class(outlet_average_pressure_t),   intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        type(mpi_comm),                     intent(in)      :: bc_COMM
        type(AD_D),                         intent(inout)   :: u_avg
        type(AD_D),                         intent(inout)   :: v_avg
        type(AD_D),                         intent(inout)   :: w_avg
        type(AD_D),                         intent(inout)   :: density_avg
        type(AD_D),                         intent(inout)   :: p_avg

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







    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. average_pressure
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_average_pressure_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop
        type(mpi_comm),                     intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,                            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,                           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            u_bc, v_bc, w_bc, p_bc,                                                     &
            u_m,  v_m,  w_m,  p_m,                                                      &
            ddensity_c,    du_c,    dv_c,    dw_c,    dp_c,                             &
            c1, c2, c3, c4, ddensity, dp, du, dv, dw


        type(AD_D)  :: p_avg, u_avg, v_avg, w_avg, density_avg, M_avg, c_avg, c4_1d,    &
                       ddensity_mean, du_mean, dv_mean, dw_mean, dp_mean

        real(rk),       allocatable, dimension(:)   ::  p_user, r
        integer :: i


        !
        ! Get back pressure from function.
        !
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())



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
        ! Update average pressure
        !
        call self%compute_averages(worker,bc_COMM,u_avg, v_avg, w_avg, density_avg, p_avg)
        c_avg = sqrt(gam*p_avg/density_avg)


        !
        ! Compute velocities
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m


        !
        ! Compute pressure from extrapolated data
        !
        p_m = worker%get_field('Pressure', 'value', 'face interior')
    

        
!        !
!        ! Define boundary primitive quantities
!        !
!        density_bc = density_m
!        u_bc       = u_m
!        v_bc       = v_m
!        w_bc       = w_m
!        p_bc       = p_user  +  (p_m - p_avg)


        !
        ! Compute update for average quantities
        !
        !c4_1d         = -TWO*(p_avg - p_user(1))
        !ddensity_mean =  c4_1d/(TWO*c_avg*c_avg)
        !du_mean       = -c4_1d/(TWO*density_avg*c_avg)
        !dp_mean       =  HALF*c4_1d
        c4_1d         = -TWO*(p_avg - p_user(1))
        ddensity_mean =  c4_1d/(TWO*c_avg*c_avg)
        dw_mean       = -c4_1d/(TWO*density_avg*c_avg)
        dp_mean       =  HALF*c4_1d


        !
        ! Compute perturbation from mean
        !
        du       = u_m       - u_avg
        dv       = v_m       - v_avg
        dw       = w_m       - w_avg
        ddensity = density_m - density_avg
        dp       = p_m       - p_avg


        !
        ! Compute 1D characteristics 
        !
        allocate(c1(size(dp)), c2(size(dp)), c3(size(dp)), c4(size(dp)))
        do i = 1,size(dp)
            !c1(i) = -c_avg*c_avg*ddensity(i)  +  dp(i)
            !c2(i) =  density_avg*c_avg*dv(i)
            !c3(i) =  density_avg*c_avg*du(i)  +  dp(i)
            !c4(i) = -density_avg*c_avg*du(i)  +  dp(i)
            c1(i) = -c_avg*c_avg*ddensity(i)  +  dp(i)
            c2(i) =  density_avg*c_avg*dv(i)
            c3(i) =  density_avg*c_avg*dw(i)  +  dp(i)
            c4(i) = -density_avg*c_avg*dw(i)  +  dp(i)
        end do




        !
        ! Compute update from characteristics for perturbation quantities: No contribution from c4
        !
        allocate(ddensity_c(size(dp)), du_c(size(dp)), dv_c(size(dp)), dw_c(size(dp)), dp_c(size(dp)))
        do i = 1,size(dp)
            !ddensity_c(i) = -c1(i)/(c_avg*c_avg)  +  c3(i)/(TWO*c_avg*c_avg) 
            !du_c(i)       =  c3(i)/(TWO*density_avg*c_avg)
            !dv_c(i)       =  c2(i)/(density_avg*c_avg)
            !dp_c(i)       =  c3(i)/TWO
            ddensity_c(i) = -c1(i)/(c_avg*c_avg)  +  c3(i)/(TWO*c_avg*c_avg) 
            dw_c(i)       =  c3(i)/(TWO*density_avg*c_avg)
            dv_c(i)       =  c2(i)/(density_avg*c_avg)
            dp_c(i)       =  c3(i)/TWO
        end do



        !
        ! Construct boundary state from average and perturbations
        !
        density_bc = density_m
        u_bc = density_m
        v_bc = density_m
        w_bc = density_m
        p_bc = density_m
        do i = 1,size(dp)
            !density_bc(i) = density_avg  +  ddensity_mean  +  ddensity_c(i)
            !u_bc(i)       = u_avg        +  du_mean        +  du_c(i)
            !v_bc(i)       = v_avg                          +  dv_c(i)
            !w_bc(i)       = w_m(i)
            !p_bc(i)       = p_avg        +  dp_mean        +  dp_c(i)
            density_bc(i) = density_avg  +  ddensity_mean  +  ddensity_c(i)
            u_bc(i)       = u_m(i)
            v_bc(i)       = v_avg                          +  dv_c(i)
            w_bc(i)       = w_avg        +  dw_mean        +  dw_c(i)
            p_bc(i)       = p_avg        +  dp_mean        +  dp_c(i)
        end do





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
        ! Store boundary condition state. Q_bc
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')





    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_outlet_average_pressure
