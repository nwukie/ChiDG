module bc_state_graddemo_gradp_extrapolate_outer
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL
    use mod_interpolate,        only: interpolate_face_autodiff
    use type_mesh,              only: mesh_t
    use type_face_info,         only: face_info_t
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: MPI_REAL8, MPI_SUM, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: graddemo_gradp_extrapolate_outer_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure

        procedure   :: compute_averages

    end type graddemo_gradp_extrapolate_outer_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_gradp_extrapolate_outer_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("graddemo gradp extrapolate outer")
        call self%set_family("Wall")


        !
        ! Add parameters
        !
        call self%bcproperties%add('Average Pressure',      'Required')
        call self%bcproperties%add('Pressure Gradient - 1', 'Required')
        call self%bcproperties%add('Pressure Gradient - 2', 'Required')
        call self%bcproperties%add('Pressure Gradient - 3', 'Required')

    end subroutine init
    !********************************************************************************




    !>  Initialize boundary group coupling.
    !!
    !!  This coupling occurs because each face uses an
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
    !!  @author Nathan A. Wukie average_pressure
    !!  @date   4/18/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
        class(graddemo_gradp_extrapolate_outer_t),  intent(inout)   :: self
        type(mesh_t),                               intent(inout)   :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_COMM

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
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c1_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c2_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)
                        !call MPI_Bcast(mesh%domain(idomain_l)%faces(ielement_l,iface)%interp_coords_def(:)%c3_,    ngq, MPI_INTEGER, iproc, bc_COMM, ierr)

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






    !>  Update the area-averaged pressure for the boundary condition.
    !!
    !!  @author Nathan A. average_pressure
    !!  @date   3/31/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, p_avg)
        class(graddemo_gradp_extrapolate_outer_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_COMM
        type(AD_D),                                 intent(inout)   :: p_avg

        type(AD_D)          :: face_p, p_integral
        type(face_info_t)   :: face_info

        type(AD_D), allocatable,    dimension(:)    :: pressure

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r

        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
                       ipressure, eqn_ID, group_ID, patch_ID, face_ID, &
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
            ! Get solution
            !
            eqn_ID = worker%mesh%domain(idomain_l)%elems(ielement_l)%eqn_ID
            ipressure = worker%prop(eqn_ID)%get_primary_field_index('Pressure_TEMP')
            itime     = 1

            
            !
            ! Interpolate coupled element solution on face of coupled element
            !
            pressure = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,face_info,worker%function_info, ipressure, itime, 'value', ME)


            !
            ! Get weights + areas
            !
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area



            !
            ! Integrate and contribute to average
            !
            face_p = sum(pressure*areas*weights)


            if (allocated(p_integral%xp_ad_)) then
                p_integral = p_integral + face_p
            else
                p_integral = face_p
            end if


            total_area = total_area + face_area


        end do !icoupled



        !
        ! Compute average pressure:
        !   area-weighted pressure integral over the total area
        !   
        !
        p_avg = p_integral / total_area



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
    !---------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(graddemo_gradp_extrapolate_outer_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            pressure, grad1_pressure, grad2_pressure, grad3_pressure

        type(AD_D)  :: p_avg, delta_p

        real(rk),   allocatable :: p_avg_user(:), grad1_p_user(:), grad2_p_user(:), grad3_p_user(:)


        !
        ! Get user parameter settings
        !
        p_avg_user   = self%bcproperties%compute('Average Pressure',     worker%time(),worker%coords())
        grad1_p_user = self%bcproperties%compute('Pressure Gradient - 1',worker%time(),worker%coords())
        grad2_p_user = self%bcproperties%compute('Pressure Gradient - 2',worker%time(),worker%coords())
        grad3_p_user = self%bcproperties%compute('Pressure Gradient - 3',worker%time(),worker%coords())



        !
        ! Interpolate interior solution to face quadrature nodes
        !
        pressure       = worker%get_field('Pressure_TEMP', 'value', 'face interior')
        grad1_pressure = worker%get_field('Pressure_TEMP', 'grad1', 'face interior')
        grad2_pressure = worker%get_field('Pressure_TEMP', 'grad2', 'face interior')
        grad3_pressure = worker%get_field('Pressure_TEMP', 'grad3', 'face interior')


        !
        ! Compute average pressure
        !
        call self%compute_averages(worker,bc_COMM,p_avg)



        !
        ! Compute average update
        !
        delta_p = p_avg - p_avg_user(1)


        !
        ! boundary pressure: extrapolating interior and adding average update.
        ! boundary pressure gradient: set to user-specified values
        !
        pressure       = pressure + delta_p
        grad1_pressure = grad1_p_user
        grad2_pressure = grad2_p_user
        grad3_pressure = grad3_p_user


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Pressure_TEMP', pressure,       'value')
        call worker%store_bc_state('Pressure_TEMP', grad1_pressure, 'grad1')
        call worker%store_bc_state('Pressure_TEMP', grad2_pressure, 'grad2')
        call worker%store_bc_state('Pressure_TEMP', grad3_pressure, 'grad3')


    end subroutine compute_bc_state
    !*******************************************************************************






end module bc_state_graddemo_gradp_extrapolate_outer
