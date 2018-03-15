module bc_state_graddemo_gradp_extrapolate_outer
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL
    use mod_interpolate,        only: interpolate_face_autodiff
    use mod_fluid,              only: gam
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
    !!  Call global coupling routine to initialize implicit coupling between each
    !!  element with every other element on the boundary, a result of averaging
    !!  and Fourier transform operations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/18/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_COMM)
        class(graddemo_gradp_extrapolate_outer_t),  intent(inout)   :: self
        type(mesh_t),                               intent(inout)   :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_COMM

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

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

! 1D Char
!        type(AD_D), allocatable, dimension(:)   ::  &
!            density, invdensity, dvel1_dmom1, dvel1_ddensity, grad1_density, grad1_mom1, mom1, grad1_pressure_m, grad2_pressure_m, grad3_pressure_m, c, pressure_m, grad1_vel1

! Axial Momentum
!        type(AD_D), allocatable, dimension(:)   ::  &
!            density, mom1, mom2, mom3, invdensity, vel1, vel2, vel3, &
!            grad1_density, grad2_density, grad3_density, grad1_vel1, grad2_vel2, grad3_vel3,    &
!            dvel1_ddensity, dvel2_ddensity, dvel3_ddensity, dvel1_dmom1, dvel2_dmom2, dvel3_dmom3,  &
!            grad3_mom3, grad2_mom2, grad1_mom1, grad2_mom1, grad3_mom1

! Axial Momentum + 1D Char
        type(AD_D), allocatable, dimension(:)   ::  &
            density, mom1, mom2, mom3, invdensity, vel1, vel2, vel3, &
            grad1_density, grad2_density, grad3_density, grad1_vel1, grad2_vel2, grad3_vel3,    &
            dvel1_ddensity, dvel2_ddensity, dvel3_ddensity, dvel1_dmom1, dvel2_dmom2, dvel3_dmom3,  &
            grad3_mom3, grad2_mom2, grad1_mom1, grad2_mom1, grad3_mom1, pressure_m, grad1_pressure_m,   &
            grad2_pressure_m, grad3_pressure_m, c

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



!        !-----------------------------------------------------
!        !
!        !   Try computing normal gradient from characteristics
!        !
!        !-----------------------------------------------------
!        density          = worker%get_field('Density',    'value', 'face interior')
!        mom1             = worker%get_field('Momentum-1', 'value', 'face interior')
!        invdensity       = ONE/density
!        dvel1_dmom1      = invdensity
!        dvel1_ddensity   = -invdensity*invdensity*mom1
!        grad1_density    = worker%get_field('Density'   , 'grad1', 'face interior')
!        grad1_mom1       = worker%get_field('Momentum-1', 'grad1', 'face interior')
!        grad1_vel1       = dvel1_ddensity*grad1_density  +  dvel1_dmom1*grad1_mom1
!        grad1_pressure_m = worker%get_field('Pressure_TEMP', 'grad1', 'face interior')
!        grad2_pressure_m = worker%get_field('Pressure_TEMP', 'grad2', 'face interior')
!        grad3_pressure_m = worker%get_field('Pressure_TEMP', 'grad3', 'face interior')
!        pressure_m       = worker%get_field('Pressure_TEMP', 'value', 'face interior')
!        c = sqrt(gam*pressure_m/density)
!        ! 1D Char seems to work here
!        grad1_pressure = -HALF*(grad1_pressure_m  +  density*c*grad1_vel1)


!        !--------------------------------------------------------------
!        !
!        !   Try computing normal gradient from normal momentum equation
!        !
!        !--------------------------------------------------------------
!        density          = worker%get_field('Density',    'value', 'face interior')
!        mom1             = worker%get_field('Momentum-1', 'value', 'face interior')
!        mom2             = worker%get_field('Momentum-2', 'value', 'face interior')
!        mom3             = worker%get_field('Momentum-3', 'value', 'face interior')
!        invdensity       = ONE/density
!
!        grad1_density = worker%get_field('Density', 'grad1', 'face interior')
!        grad2_density = worker%get_field('Density', 'grad2', 'face interior')
!        grad3_density = worker%get_field('Density', 'grad3', 'face interior')
!
!        grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'face interior')
!        grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'face interior')
!        grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'face interior')
!
!        grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'face interior')
!        grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'face interior')
!
!
!        !
!        ! Compute velocities
!        !
!        vel1 = mom1/density
!        vel2 = mom2/density
!        vel3 = mom3/density
!
!        !
!        ! compute velocity jacobians
!        !
!        dvel1_ddensity  = -invdensity*invdensity*mom1
!        dvel2_ddensity  = -invdensity*invdensity*mom2
!        dvel3_ddensity  = -invdensity*invdensity*mom3
!
!        dvel1_dmom1 = invdensity
!        dvel2_dmom2 = invdensity
!        dvel3_dmom3 = invdensity
!
!        !
!        ! compute velocity gradients via chain rule:
!        !
!        !   u = f(rho,rhou)
!        !
!        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
!        !
!        grad1_vel1 = dvel1_ddensity*grad1_density  +  dvel1_dmom1*grad1_mom1
!        grad2_vel2 = dvel2_ddensity*grad2_density  +  dvel2_dmom2*grad2_mom2
!        grad3_vel3 = dvel3_ddensity*grad3_density  +  dvel3_dmom3*grad3_mom3
!
!        
!        grad1_pressure = -( (vel1*grad1_mom1  +  mom1*grad1_vel1)  +  &
!                            (vel2*grad2_mom1  +  mom1*grad2_vel2)  +  &
!                            (vel3*grad3_mom1  +  mom1*grad3_vel3) )
!
!
!        !---------------------------------------------------------------
!        !
!        ! Combination of 1D char and axial momentum
!        !
!        !---------------------------------------------------------------
!        density          = worker%get_field('Density',    'value', 'face interior')
!        mom1             = worker%get_field('Momentum-1', 'value', 'face interior')
!        mom2             = worker%get_field('Momentum-2', 'value', 'face interior')
!        mom3             = worker%get_field('Momentum-3', 'value', 'face interior')
!        invdensity       = ONE/density
!
!        grad1_density = worker%get_field('Density', 'grad1', 'face interior')
!        grad2_density = worker%get_field('Density', 'grad2', 'face interior')
!        grad3_density = worker%get_field('Density', 'grad3', 'face interior')
!
!        grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'face interior')
!        grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'face interior')
!        grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'face interior')
!
!        grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'face interior')
!        grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'face interior')
!
!
!        grad1_pressure_m = worker%get_field('Pressure_TEMP', 'grad1', 'face interior')
!        grad2_pressure_m = worker%get_field('Pressure_TEMP', 'grad2', 'face interior')
!        grad3_pressure_m = worker%get_field('Pressure_TEMP', 'grad3', 'face interior')
!        pressure_m       = worker%get_field('Pressure_TEMP', 'value', 'face interior')
!        c = sqrt(gam*pressure_m/density)
!
!
!        !
!        ! Compute velocities
!        !
!        vel1 = mom1/density
!        vel2 = mom2/density
!        vel3 = mom3/density
!
!        !
!        ! compute velocity jacobians
!        !
!        dvel1_ddensity  = -invdensity*invdensity*mom1
!        dvel2_ddensity  = -invdensity*invdensity*mom2
!        dvel3_ddensity  = -invdensity*invdensity*mom3
!
!        dvel1_dmom1 = invdensity
!        dvel2_dmom2 = invdensity
!        dvel3_dmom3 = invdensity
!
!        !
!        ! compute velocity gradients via chain rule:
!        !
!        !   u = f(rho,rhou)
!        !
!        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
!        !
!        grad1_vel1 = dvel1_ddensity*grad1_density  +  dvel1_dmom1*grad1_mom1
!        grad2_vel2 = dvel2_ddensity*grad2_density  +  dvel2_dmom2*grad2_mom2
!        grad3_vel3 = dvel3_ddensity*grad3_density  +  dvel3_dmom3*grad3_mom3
!
!        
!        !grad1_pressure = -( (vel1*grad1_mom1  +  mom1*grad1_vel1)  +  &
!        !                    (vel2*grad2_mom1  +  mom1*grad2_vel2)  +  &
!        !                    (vel3*grad3_mom1  +  mom1*grad3_vel3) )
!
!        ! Contribution from 1D characteristics
!        grad1_pressure = -HALF*(grad1_pressure_m  +  density*c*grad1_vel1)
!
!        ! Contribution from axial momentum equation
!        grad1_pressure = grad1_pressure - ( (vel2*grad2_mom1  +  mom1*grad2_vel2)  +  &
!                                            (vel3*grad3_mom1  +  mom1*grad3_vel3) )

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
