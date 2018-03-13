module bc_state_auxiliary_boundary
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
    type, public, extends(bc_state_t) :: auxiliary_boundary_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure

        procedure   :: compute_averages

    end type auxiliary_boundary_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(auxiliary_boundary_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("Auxiliary Gradient Boundary")
        call self%set_family("Wall")


        !
        ! Add parameters
        !
        call self%bcproperties%add('Average Pressure',      'Required')

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
        class(auxiliary_boundary_t),  intent(inout)   :: self
        type(mesh_t),                               intent(inout)   :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_COMM

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !******************************************************************************************





    !>  Compute averaged quantities over the face. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, u_avg, v_avg, w_avg, density_avg, p_avg, c_avg)
        class(auxiliary_boundary_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_COMM
        type(AD_D),                     intent(inout)   :: u_avg
        type(AD_D),                     intent(inout)   :: v_avg
        type(AD_D),                     intent(inout)   :: w_avg
        type(AD_D),                     intent(inout)   :: density_avg
        type(AD_D),                     intent(inout)   :: p_avg
        type(AD_D),                     intent(inout)   :: c_avg

        type(AD_D)          :: p_integral, u_integral, v_integral, w_integral, density_integral, c_integral, &
                               face_density, face_u, face_v, face_w, face_p, face_c
        type(face_info_t)   :: face_info

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3, energy, p, c,    &
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

            

            !
            ! Compute velocities and pressure
            !
            u = mom1 / density
            v = mom2 / density
            w = mom3 / density
            p = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )
            c = sqrt(gam*p/density)


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
            face_c       = sum(c       * areas * weights)



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

            if (allocated(c_integral%xp_ad_)) then
                c_integral = c_integral + face_c
            else
                c_integral = face_c
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
        c_avg       = c_integral       / total_area



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
        class(auxiliary_boundary_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density,    grad1_density,      grad2_density,      grad3_density,      &
            pressure,   grad1_pressure,     grad2_pressure,     grad3_pressure,     &
            vel1,       grad1_vel1,         grad2_vel1,         grad3_vel1,         &
                        grad1_density_m,    grad2_density_m,    grad3_density_m,    &
                        grad1_vel1_m,       grad2_vel1_m,       grad3_vel1_m,       &
                        grad1_vel2_m,       grad2_vel2_m,       grad3_vel2_m,       &
                        grad1_vel3_m,       grad2_vel3_m,       grad3_vel3_m,       &
                        grad1_pressure_m,   grad2_pressure_m,   grad3_pressure_m,   &
                        grad1_phi1,         grad1_phi2,         grad1_phi3,         &
            density_m, mom1_m, vel1_m


        type(AD_D)  :: vel1_avg, vel2_avg, vel3_avg, p_avg, c_avg, density_avg, delta_p

        real(rk),   allocatable :: p_avg_user(:)
        integer(ik) :: igq


        !
        ! Get user parameter settings
        !
        p_avg_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())



        !
        ! Interpolate Auxiliary Fields and Gradients
        !
        density        = worker%get_field('Density_TEMP',    'value', 'face interior')
        vel1           = worker%get_field('Velocity-1_TEMP', 'value', 'face interior')
        pressure       = worker%get_field('Pressure_TEMP',   'value', 'face interior')

        grad1_density  = worker%get_field('Density_TEMP',    'grad1', 'face interior')
        grad2_density  = worker%get_field('Density_TEMP',    'grad2', 'face interior')
        grad3_density  = worker%get_field('Density_TEMP',    'grad3', 'face interior')

        grad1_vel1     = worker%get_field('Velocity-1_TEMP', 'grad1', 'face interior')
        grad2_vel1     = worker%get_field('Velocity-1_TEMP', 'grad2', 'face interior')
        grad3_vel1     = worker%get_field('Velocity-1_TEMP', 'grad3', 'face interior')

        grad1_pressure = worker%get_field('Pressure_TEMP',   'grad1', 'face interior')
        grad2_pressure = worker%get_field('Pressure_TEMP',   'grad2', 'face interior')
        grad3_pressure = worker%get_field('Pressure_TEMP',   'grad3', 'face interior')



        call compute_pressure_gradient( worker, grad1_pressure_m, grad2_pressure_m, grad3_pressure_m)
        call compute_density_gradient(  worker, grad1_density_m,  grad2_density_m,  grad3_density_m)
        call compute_velocity_gradients(worker, grad1_vel1_m,     grad2_vel1_m,     grad3_vel1_m,   &
                                                grad1_vel2_m,     grad2_vel2_m,     grad3_vel2_m,   &
                                                grad1_vel3_m,     grad2_vel3_m,     grad3_vel3_m)



        !
        ! Compute average pressure
        !
        call self%compute_averages(worker,bc_comm,vel1_avg,vel2_avg,vel3_avg,density_avg,p_avg,c_avg)


        !
        ! Compute axial derivatives of 1D characteristics
        !
        grad1_phi1 = grad1_density_m
        grad1_phi2 = grad1_density_m
        grad1_phi3 = grad1_density_m
        do igq = 1,size(grad1_density_m)
            grad1_phi1(igq) = -c_avg*c_avg*grad1_density_m(igq)     +  grad1_pressure_m(igq)
            grad1_phi2(igq) =  density_avg*c_avg*grad1_vel1_m(igq)  +  grad1_pressure_m(igq)
            grad1_phi3(igq) = -density_avg*c_avg*grad1_vel1_m(igq)  +  grad1_pressure_m(igq)
        end do

        
        !
        ! Zero incoming characteristic
        !
        grad1_phi3 = ZERO


        !
        ! Reconstruct axial derivatives of primitive variables after zeroing 
        ! incoming characteristic
        !
        do igq = 1,size(grad1_phi1)
            !grad1_density(igq)  = -(ONE/(c_avg*c_avg))*grad1_phi1(igq)  +  (ONE/(TWO*c_avg*c_avg))*grad1_phi2(igq)        +  (ONE/(TWO*c_avg*c_avg))*grad1_phi3(igq)
            !grad1_vel1(igq)     =                                          (ONE/(TWO*density_avg*c_avg))*grad1_phi2(igq)  -  (ONE/(TWO*density_avg*c_avg))*grad1_phi3(igq)
            grad1_pressure(igq) =                                          HALF*grad1_phi2(igq)                           +  HALF*grad1_phi3(igq)
        end do


        !
        ! Reverse axial gradient
        !
        grad1_density  = grad1_density
        grad1_vel1     = grad1_vel1
        !grad1_pressure = grad1_pressure
        !grad1_pressure = -0.25_rk*grad1_pressure
        !grad1_pressure = -HALF*(grad1_pressure_m  +  density*c*grad1_vel1)
        do igq = 1,size(grad1_density_m)
            grad1_pressure(igq) = -HALF*(grad1_pressure_m(igq)  +  density_avg*c_avg*grad1_vel1_m(igq))
        end do



        !
        ! Extrapolate transverse derivatives
        !
        grad2_density  = grad2_density
        grad3_density  = grad3_density

        grad2_vel1     = grad2_vel1
        grad3_vel1     = grad3_vel1

        grad2_pressure = grad2_pressure
        grad3_pressure = grad3_pressure


!        !
!        ! Override characteristic gradient
!        !
!        grad1_density  = ZERO
!        grad1_vel1     = ZERO
!        grad1_pressure = ZERO

        !
        ! Compute average update
        !
        delta_p = p_avg - p_avg_user(1)

        
        !
        ! Compute velocity/density from original problem
        !
        density_m = worker%get_field('Density',    'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        vel1_m = mom1_m/density_m


        !
        ! Compute boundary values
        !
        density  = density_m
        vel1     = vel1_m
        pressure = pressure + delta_p


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density_TEMP',    density,        'value')
        call worker%store_bc_state('Density_TEMP',    grad1_density,  'grad1')
        call worker%store_bc_state('Density_TEMP',    grad2_density,  'grad2')
        call worker%store_bc_state('Density_TEMP',    grad3_density,  'grad3')

        call worker%store_bc_state('Velocity-1_TEMP', vel1,           'value')
        call worker%store_bc_state('Velocity-1_TEMP', grad1_vel1,     'grad1')
        call worker%store_bc_state('Velocity-1_TEMP', grad2_vel1,     'grad2')
        call worker%store_bc_state('Velocity-1_TEMP', grad3_vel1,     'grad3')

        call worker%store_bc_state('Pressure_TEMP',   pressure,       'value')
        call worker%store_bc_state('Pressure_TEMP',   grad1_pressure, 'grad1')
        call worker%store_bc_state('Pressure_TEMP',   grad2_pressure, 'grad2')
        call worker%store_bc_state('Pressure_TEMP',   grad3_pressure, 'grad3')



    end subroutine compute_bc_state
    !*******************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_pressure_gradient(worker,grad1_p, grad2_p, grad3_p) 
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_p 
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_p 
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_p 

        type(AD_D), allocatable, dimension(:)   ::                              &
            density,       mom1,       mom2,       mom3,       energy,          &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,    &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,    &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy,    &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy

        real(rk),   allocatable, dimension(:)   :: r

        !
        ! Interpolate solution to quadrature nodes
        !
        density       = worker%get_field('Density',    'value', 'face interior')
        mom1          = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2          = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3          = worker%get_field('Momentum-3', 'value', 'face interior')
        energy        = worker%get_field('Energy',     'value', 'face interior')

        grad1_density = worker%get_field('Density',    'grad1', 'face interior', override_lift=.true.)
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', 'face interior', override_lift=.true.)
        grad1_energy  = worker%get_field('Energy',     'grad1', 'face interior', override_lift=.true.)


        grad2_density = worker%get_field('Density',    'grad2', 'face interior', override_lift=.true.)
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'face interior', override_lift=.true.)
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'face interior', override_lift=.true.)
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', 'face interior', override_lift=.true.)
        grad2_energy  = worker%get_field('Energy',     'grad2', 'face interior', override_lift=.true.)


        grad3_density = worker%get_field('Density',    'grad3', 'face interior', override_lift=.true.)
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', 'face interior', override_lift=.true.)
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', 'face interior', override_lift=.true.)
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', 'face interior', override_lift=.true.)
        grad3_energy  = worker%get_field('Energy',     'grad3', 'face interior', override_lift=.true.)


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if



        !
        ! Compute pressure jacobians
        !
        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
        dp_dmom1    = -(gam-ONE)*mom1/density
        dp_dmom2    = -(gam-ONE)*mom2/density
        dp_dmom3    = -(gam-ONE)*mom3/density
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)


        !
        ! Compute pressure gradient
        !
        grad1_p = dp_ddensity * grad1_density  + &
                  dp_dmom1    * grad1_mom1     + &
                  dp_dmom2    * grad1_mom2     + &
                  dp_dmom3    * grad1_mom3     + &
                  dp_denergy  * grad1_energy

        grad2_p = dp_ddensity * grad2_density  + &
                  dp_dmom1    * grad2_mom1     + &
                  dp_dmom2    * grad2_mom2     + &
                  dp_dmom3    * grad2_mom3     + &
                  dp_denergy  * grad2_energy

        grad3_p = dp_ddensity * grad3_density  + &
                  dp_dmom1    * grad3_mom1     + &
                  dp_dmom2    * grad3_mom2     + &
                  dp_dmom3    * grad3_mom3     + &
                  dp_denergy  * grad3_energy


    end subroutine compute_pressure_gradient
    !******************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_velocity_gradients(worker,  &
                                          grad1_vel1, grad2_vel1, grad3_vel1,   &
                                          grad1_vel2, grad2_vel2, grad3_vel2,   &
                                          grad1_vel3, grad2_vel3, grad3_vel3)
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel1
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel2
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_vel3
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_vel3
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_vel3

        type(AD_D), allocatable, dimension(:)   ::              &
            density,       mom1,       mom2,       mom3,        &
                           vel1,       vel2,       vel3,        &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3,  &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3,  &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3

        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density       = worker%get_field('Density',    'value', 'face interior')
        mom1          = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2          = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3          = worker%get_field('Momentum-3', 'value', 'face interior')

        grad1_density = worker%get_field('Density',    'grad1', 'face interior', override_lift=.true.)
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'face interior', override_lift=.true.)
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', 'face interior', override_lift=.true.)


        grad2_density = worker%get_field('Density',    'grad2', 'face interior',override_lift=.true.)
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'face interior',override_lift=.true.)
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', 'face interior',override_lift=.true.)


        grad3_density = worker%get_field('Density',    'grad3', 'face interior',override_lift=.true.)
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', 'face interior',override_lift=.true.)
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', 'face interior',override_lift=.true.)

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if

        
        !
        ! Compute velocities
        !
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density


        !
        ! Compute velocity gradient
        !
        grad1_vel1 = (ONE/density)*grad1_mom1  -  (vel1/density)*grad1_density
        grad2_vel1 = (ONE/density)*grad2_mom1  -  (vel1/density)*grad2_density
        grad3_vel1 = (ONE/density)*grad3_mom1  -  (vel1/density)*grad3_density

        grad1_vel2 = (ONE/density)*grad1_mom2  -  (vel2/density)*grad1_density
        grad2_vel2 = (ONE/density)*grad2_mom2  -  (vel2/density)*grad2_density
        grad3_vel2 = (ONE/density)*grad3_mom2  -  (vel2/density)*grad3_density

        grad1_vel3 = (ONE/density)*grad1_mom3  -  (vel3/density)*grad1_density
        grad2_vel3 = (ONE/density)*grad2_mom3  -  (vel3/density)*grad2_density
        grad3_vel3 = (ONE/density)*grad3_mom3  -  (vel3/density)*grad3_density


    end subroutine compute_velocity_gradients
    !*******************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2018
    !!
    !------------------------------------------------------------------------------
    subroutine compute_density_gradient(worker, grad1_density, grad2_density, grad3_density)
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad1_density
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad2_density
        type(AD_D), allocatable, dimension(:),  intent(inout)   :: grad3_density


        grad1_density = worker%get_field('Density', 'grad1', 'face interior', override_lift=.true.)
        grad2_density = worker%get_field('Density', 'grad2', 'face interior', override_lift=.true.)
        grad3_density = worker%get_field('Density', 'grad3', 'face interior', override_lift=.true.)


    end subroutine compute_density_gradient
    !*******************************************************************************

















end module bc_state_auxiliary_boundary
