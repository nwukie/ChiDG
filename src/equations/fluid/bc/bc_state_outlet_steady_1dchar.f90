module bc_state_outlet_steady_1dchar
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL, NO_ID
    use mod_fluid,              only: gam
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t, element_info
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_SUM, MPI_AllReduce, mpi_comm, &
                                      MPI_INTEGER, MPI_BCast
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
    type, public, extends(bc_state_t) :: outlet_steady_1dchar_t


    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! boundary condition function implementation
        procedure   :: compute_averages

    end type outlet_steady_1dchar_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_steady_1dchar_t),   intent(inout) :: self
        
        ! Set name, family
        call self%set_name('Outlet - Steady 1D Characteristics')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure','Required')

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
    subroutine init_bc_coupling(self,mesh,group_ID,bc_comm)
        class(outlet_steady_1dchar_t),  intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        integer(ik),                    intent(in)      :: group_ID
        type(mpi_comm),                 intent(in)      :: bc_comm

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !********************************************************************************






    !>  Compute averaged quantities over the face. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, v1_avg, v2_avg, v3_avg, vn_avg, density_avg, p_avg)
        class(outlet_steady_1dchar_t),   intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        type(mpi_comm),                     intent(in)      :: bc_COMM
        type(AD_D),                         intent(inout)   :: v1_avg
        type(AD_D),                         intent(inout)   :: v2_avg
        type(AD_D),                         intent(inout)   :: v3_avg
        type(AD_D),                         intent(inout)   :: vn_avg
        type(AD_D),                         intent(inout)   :: density_avg
        type(AD_D),                         intent(inout)   :: p_avg

        type(AD_D)          :: p_integral, v1_integral, v2_integral, v3_integral, &
                               vn_integral, density_integral, &
                               face_density, face_v1, face_v2, face_v3, face_vn, face_p
        type(element_info_t)    :: coupled_element

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3, energy, p,       &
            v1, v2, v3, vn

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r

        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
                       ielement_l_coupled, iface_coupled, dof_start_coupled, coupled_iface
        real(rk)    :: face_area, total_area

        ! Zero integrated quantities
        total_area = ZERO

        ! Get location on domain
        idomain_l  = worker%element_info%idomain_l
        ielement_l = worker%element_info%ielement_l
        iface      = worker%iface

        ! Get location on bc_patch_group
        group_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
        patch_ID = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
        face_ID  = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID

        ! Loop through coupled faces and compute their contribution to the average pressure
        do icoupled = 1,worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

            ! Get solution
            idensity = 1
            imom1    = 2
            imom2    = 3
            imom3    = 4
            ienergy  = 5
            itime    = 1


            ! Get face info from coupled element we want to interpolate from
            idomain_g_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g( icoupled)
            idomain_l_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l( icoupled)
            ielement_g_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(icoupled)
            ielement_l_coupled = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(icoupled)
            iface_coupled      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%iface(     icoupled)
            dof_start_coupled  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_start( icoupled)

            coupled_element = element_info(idomain_g       = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(icoupled),        &
                                           idomain_l       = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(icoupled),        &
                                           ielement_g      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(icoupled),       &
                                           ielement_l      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(icoupled),       &
                                           iproc           = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(icoupled),             &
                                           pelem_ID        = NO_ID,                                                                                             &
                                           eqn_ID          = NO_ID,                                                                                             &
                                           nfields         = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nfields(icoupled),          &
                                           ntime           = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ntime(icoupled),            &
                                           nterms_s        = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nterms_s(icoupled),         &
                                           nterms_c        = 0,                                                                                                 &
                                           dof_start       = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_start(icoupled),        &
                                           dof_local_start = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_local_start(icoupled),  &
                                           recv_comm       = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_comm(icoupled),        &
                                           recv_domain     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_domain(icoupled),      &
                                           recv_element    = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_element(icoupled),     &
                                           recv_dof        = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_dof(icoupled))


            coupled_iface = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%iface(icoupled)


            
            ! Interpolate coupled element solution on face of coupled element
            density = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, idensity, itime, 'value', ME)
            mom1    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom1,    itime, 'value', ME)
            mom2    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom2,    itime, 'value', ME)
            mom3    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom3,    itime, 'value', ME)
            energy  = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, ienergy,  itime, 'value', ME)


            ! TODO: fix for parallel!
            if (worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%coordinate_system == CYLINDRICAL) then
                mom2 = mom2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%interp_coords_def(:,1)
            end if
            !r = worker%coordinate('1','boundary')
            !if (worker%coordinate_system() == 'Cylindrical') then
            !    mom2 = mom2 / r
            !end if

            ! Compute velocities and pressure
            v1 = mom1/density
            v2 = mom2/density
            v3 = mom3/density
            p = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density )

            ! TODO: fix for parallel!
            vn = v1*worker%mesh%domain(idomain_l_coupled)%faces(ielement_l_coupled,iface_coupled)%unorm_def(:,1) + &
                 v2*worker%mesh%domain(idomain_l_coupled)%faces(ielement_l_coupled,iface_coupled)%unorm_def(:,2) + &
                 v3*worker%mesh%domain(idomain_l_coupled)%faces(ielement_l_coupled,iface_coupled)%unorm_def(:,3)


            ! TODO: fix for parallel!
            ! Get weights + areas
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area

            ! Integrate and contribute to average
            face_v1      = sum(v1      * areas * weights)
            face_v2      = sum(v2      * areas * weights)
            face_v3      = sum(v3      * areas * weights)
            face_vn      = sum(vn      * areas * weights)
            face_density = sum(density * areas * weights)
            face_p       = sum(p       * areas * weights)


            if (allocated(v1_integral%xp_ad_)) then
                v1_integral = v1_integral + face_v1
            else
                v1_integral = face_v1
            end if

            if (allocated(v2_integral%xp_ad_)) then
                v2_integral = v2_integral + face_v2
            else
                v2_integral = face_v2
            end if

            if (allocated(v3_integral%xp_ad_)) then
                v3_integral = v3_integral + face_v3
            else
                v3_integral = face_v3
            end if

            if (allocated(vn_integral%xp_ad_)) then
                vn_integral = vn_integral + face_vn
            else
                vn_integral = face_vn
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


        ! Compute average pressure:
        !   area-weighted pressure integral over the total area
        v1_avg      = v1_integral      / total_area
        v2_avg      = v2_integral      / total_area
        v3_avg      = v3_integral      / total_area
        vn_avg      = vn_integral      / total_area
        density_avg = density_integral / total_area
        p_avg       = p_integral       / total_area

    end subroutine compute_averages
    !*******************************************************************************************







    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. steady_1dchar
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_steady_1dchar_t),    intent(inout)   :: self
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
            v1_bc, v2_bc, v3_bc, p_bc,                                                  &
            v1_m,  v2_m,  v3_m,  vn_m, p_m,                                             &
            v1_t,  v2_t,  v3_t,                                                         &
            ddensity_c, dv1_c, dv2_c, dv3_c, dvn_c, dp_c,                               &
            c1, c2, c3, c4, c5, ddensity, dp, dv1, dv2, dv3, dvn,                       &
            dv1_avg, dv2_avg, dv3_avg, vn_avg_1, vn_avg_2, vn_avg_3


        type(AD_D)  :: p_avg, v1_avg, v2_avg, v3_avg, vn_avg, density_avg, M_avg, c_avg, c4_1d,    &
                       ddensity_avg, dvn_avg, dp_avg

        real(rk),       allocatable, dimension(:)   ::  p_user, r
        integer :: i

        ! Get back pressure from function.
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())

        ! Interpolate interior solution to face quadrature nodes
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



        ! Store boundary gradient state. Grad(Q_bc). Do this here, before we
        ! compute any transformations for cylindrical.
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


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)
        end if


        ! Update average pressure
        call self%compute_averages(worker,bc_COMM,v1_avg, v2_avg, v3_avg, vn_avg, density_avg, p_avg)
        c_avg = sqrt(gam*p_avg/density_avg)


        ! Compute velocities
        v1_m = mom1_m/density_m
        v2_m = mom2_m/density_m
        v3_m = mom3_m/density_m
        vn_m = v1_m*worker%unit_normal(1) + &
               v2_m*worker%unit_normal(2) + &
               v3_m*worker%unit_normal(3)

        ! Compute tangential velocity part
        v1_t = v1_m - vn_m*worker%unit_normal(1)
        v2_t = v2_m - vn_m*worker%unit_normal(2)
        v3_t = v3_m - vn_m*worker%unit_normal(3)

        ! Compute pressure from extrapolated data
        p_m = worker%get_field('Pressure', 'value', 'face interior')
    

        ! Compute update for average quantities
        c4_1d        = TWO*(p_user(1) - p_avg)
        ddensity_avg =  c4_1d/(TWO*c_avg*c_avg)
        dvn_avg      = -c4_1d/(TWO*density_avg*c_avg)
        dp_avg       =  HALF*c4_1d

        ! Compute the average update along physical coordinates
        dv1_avg = dvn_avg*worker%unit_normal(1)
        dv2_avg = dvn_avg*worker%unit_normal(2)
        dv3_avg = dvn_avg*worker%unit_normal(3)

        ! Compute perturbation from avg
        dv1      = v1_m      - v1_avg
        dv2      = v2_m      - v2_avg
        dv3      = v3_m      - v3_avg
        dvn      = vn_m      - vn_avg
        ddensity = density_m - density_avg
        dp       = p_m       - p_avg


        ! Compute 1D characteristics 
        allocate(c1(size(dp)), c2(size(dp)), c3(size(dp)), c4(size(dp)), c5(size(dp)))
        do i = 1,size(dp)
            c1(i) = -c_avg*c_avg*ddensity(i)   +  dp(i)
            c4(i) =  density_avg*c_avg*dvn(i)  +  dp(i)
            c5(i) = -density_avg*c_avg*dvn(i)  +  dp(i)
        end do


        ! Compute update from characteristics for perturbation quantities: No contrib from c5
        allocate(ddensity_c(size(dp)), dv1_c(size(dp)), dv2_c(size(dp)), dv3_c(size(dp)), dvn_c(size(dp)), dp_c(size(dp)))
        do i = 1,size(dp)
            ddensity_c(i) = -c1(i)/(c_avg*c_avg)  +  c4(i)/(TWO*c_avg*c_avg) 
            dvn_c(i)      =  c4(i)/(TWO*density_avg*c_avg)
            dp_c(i)       =  c4(i)/TWO
        end do


        ! Compute velocity perturbation along physical coordinates
        dv1_c = dvn_c*worker%unit_normal(1)
        dv2_c = dvn_c*worker%unit_normal(2)
        dv3_c = dvn_c*worker%unit_normal(3)


        ! Compute contribution to each direction from average normal velocity
        vn_avg_1 = vn_avg*worker%unit_normal(1)
        vn_avg_2 = vn_avg*worker%unit_normal(2)
        vn_avg_3 = vn_avg*worker%unit_normal(3)


        ! Construct boundary state from average and perturbations
        density_bc = density_m
        v1_bc = density_m
        v2_bc = density_m
        v3_bc = density_m
        p_bc = density_m
        do i = 1,size(dp)
            density_bc(i) = density_avg  +  ddensity_avg  +  ddensity_c(i)
            v1_bc(i)      = vn_avg_1(i)  +  dv1_avg(i)    +  dv1_c(i)      +   v1_t(i)
            v2_bc(i)      = vn_avg_2(i)  +  dv2_avg(i)    +  dv2_c(i)      +   v2_t(i)
            v3_bc(i)      = vn_avg_3(i)  +  dv3_avg(i)    +  dv3_c(i)      +   v3_t(i)
            p_bc(i)       = p_avg        +  dp_avg        +  dp_c(i)
        end do

        ! Form conserved variables
        density_bc = density_bc
        mom1_bc    = density_bc*v1_bc
        mom2_bc    = density_bc*v2_bc
        mom3_bc    = density_bc*v3_bc
        energy_bc  = p_bc/(gam - ONE) + (density_bc*HALF)*(v1_bc*v1_bc + &
                                                           v2_bc*v2_bc + &
                                                           v3_bc*v3_bc)

        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        ! Store boundary condition state. Q_bc
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


    end subroutine compute_bc_state
    !****************************************************************************************






end module bc_state_outlet_steady_1dchar
