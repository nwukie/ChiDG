module bc_state_outlet_neumann_pressure_fd_test
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL, &
                                      XI_MAX, ETA_MAX, ZETA_MAX
    use mod_fluid,              only: gam, Rgas

    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_SUM, MPI_AllReduce, mpi_comm, &
                                      MPI_INTEGER, MPI_BCast
    use DNAD_D
    implicit none





    !>  Pressure gradient condition via extrapolation of the pressure from some
    !!  interior location and computing a boundary-global offset parameter to
    !!  achieve an average pressure.
    !!
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!
    !!  References:
    !!              
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/23/2018
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: outlet_neumann_pressure_fd_test_t


    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Implement specialized initialization procedure
        procedure   :: compute_bc_state     ! boundary condition function implementation

        procedure   :: compute_averages

    end type outlet_neumann_pressure_fd_test_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/23/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_neumann_pressure_fd_test_t),   intent(inout) :: self
        
        ! Set name, family
        call self%set_name('Outlet - Neumann Pressure Finite Difference Test')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',  'Required')
        call self%bcproperties%add('Normal Derivative', 'Required')

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
        class(outlet_neumann_pressure_fd_test_t),    intent(inout)   :: self
        type(mesh_t),                           intent(inout)   :: mesh
        integer(ik),                            intent(in)      :: group_ID
        type(mpi_comm),                         intent(in)      :: bc_COMM

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !********************************************************************************






    !>  Update averages for the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/31/2017
    !!
    !---------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM,p_avg,c_avg,density_avg)
        class(outlet_neumann_pressure_fd_test_t),  intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        type(mpi_comm),         intent(in)      :: bc_COMM
        type(AD_D),             intent(inout)   :: p_avg
        type(AD_D),             intent(inout)   :: c_avg
        type(AD_D),             intent(inout)   :: density_avg

        type(AD_D)  :: face_p, face_c, face_density, p_integral, c_integral, density_integral

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3, energy, p, c

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r

        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime, &
                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
                       ielement_l_coupled, iface_coupled
        real(rk)    :: face_area, total_area
        type(face_info_t)   :: face_info


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

        ! Loop through coupled faces and compute their contribution to the average pressure
        do icoupled = 1,worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

            ! Get solution
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

            ! We don't want this, because this won't return the correct radius for 
            ! the current face being interpolated from.
            !r = worker%coordinate('1','boundary')
            !if (worker%coordinate_system() == 'Cylindrical') then
            !    mom_2 = mom_2 / r
            !end if


            
            ! Compute velocities and pressure
            p = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )
            c = sqrt(gam*p/density)


            ! Get weights + areas
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area


            ! Integrate and contribute to average
            face_p = sum(p*areas*weights)
            face_c = sum(c*areas*weights)
            face_density = sum(density*areas*weights)

            if (allocated(p_integral%xp_ad_)) then
                p_integral = p_integral + face_p
            else
                p_integral = face_p
            end if

            if (allocated(c_integral%xp_ad_)) then
                c_integral = c_integral + face_c
            else
                c_integral = face_c
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
        p_avg = p_integral / total_area
        c_avg = c_integral / total_area
        density_avg = density_integral / total_area


    end subroutine compute_averages
    !*******************************************************************************************







    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. neumann_pressure_fd_test
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(outlet_neumann_pressure_fd_test_t),    intent(inout)   :: self
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
            v1_bc, v2_bc, v3_bc, p_bc, T_bc,                                            &
            v1_m,  v2_m,  v3_m,  p_m,  T_m,                                             &
            density_o, mom1_o, mom2_o, mom3_o, energy_o, p_o, v1_o,                     &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,                      &
            grad1_v1, grad1_p, gradn_p, gradn_v1, gradn_density, grad1_phi1, grad1_phi4


        type(AD_D)  :: p_avg, c_avg, density_avg, delta_p_avg

        real(rk),   allocatable, dimension(:)   :: p_user, r, delta_n
        real(rk),   allocatable, dimension(:,:) :: ref_coords, phys_coords, phys_coords_offset
        real(rk)                                :: offset

        integer :: i


        ! Interpolate interior solution to face quadrature nodes
        density_m = worker%get_field('Density'    , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1' , 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2' , 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3' , 'value', 'face interior')
        energy_m  = worker%get_field('Energy'     , 'value', 'face interior')
        p_m       = worker%get_field('Pressure'   , 'value', 'face interior')
        T_m       = worker%get_field('Temperature', 'value', 'face interior')

        ! Get back pressure from function.
        gradn_p = p_m ! allocate
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        gradn_p = self%bcproperties%compute('Normal Derivative',worker%time(),worker%coords())
        
        associate( idom => worker%element_info%idomain_l, ielem => worker%element_info%ielement_l, iface => worker%iface )

            ! Retrieve reference and physical coordinates
            ref_coords = worker%mesh%domain(idom)%faces(ielem,iface)%basis_s%nodes_face(iface)
            phys_coords = worker%mesh%domain(idom)%faces(ielem,iface)%interp_coords_def(:,:)

            ! Offset axial reference nodes to lie in the element interior
            offset = 0.5    ! quarter-element distance
            !offset = 0.25    ! eight-element distance
            if ( (ref_coords(1,1) - ONE) < 1.e-7_rk ) then
                ref_coords(:,1) = ref_coords(:,1) + offset
            else if ( (ref_coords(1,2) - ONE) < 1.e-7_rk ) then
                ref_coords(:,2) = ref_coords(:,2) + offset
            else if ( (ref_coords(1,3) - ONE) < 1.e-7_rk ) then
                ref_coords(:,3) = ref_coords(:,3) + offset
            else
                call chidg_signal(FATAL,"bc_state_outlet_neumann_pressure_fd_test: invalid condition.")
            end if

            ! Allocate offset physical coords and compute from element
            phys_coords_offset = phys_coords
            do i = 1,size(phys_coords,1)
                phys_coords_offset(i,:) = worker%mesh%domain(idom)%elems(ielem)%physical_point(ref_coords(i,:),'Deformed')
            end do

        end associate

        ! Compute offset in physical space
        allocate(delta_n(size(phys_coords,1))) ! Silence debug error
        delta_n = sqrt((phys_coords(:,1) - phys_coords_offset(:,1))**TWO + &
                       (phys_coords(:,2) - phys_coords_offset(:,2))**TWO + &
                       (phys_coords(:,3) - phys_coords_offset(:,3))**TWO)


        ! Compute interpolation of solution at offset quadrature node set
        density_o = worker%interpolate_field('Density',    ref_coords)
        mom1_o    = worker%interpolate_field('Momentum-1', ref_coords)
        mom2_o    = worker%interpolate_field('Momentum-2', ref_coords)
        mom3_o    = worker%interpolate_field('Momentum-3', ref_coords)
        energy_o  = worker%interpolate_field('Energy',     ref_coords)
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_o = mom2_o / worker%coordinate('1','boundary')
        end if
        p_o = (gam-ONE)*(energy_o - HALF*( mom1_o*mom1_o + mom2_o*mom2_o + mom3_o*mom3_o )/density_o )



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
        end if


        ! Extrapolate temperature and velocity
        T_bc = T_m
        v1_bc = mom1_m/density_m
        v2_bc = mom2_m/density_m
        v3_bc = mom3_m/density_m


        ! Update average pressure
        call self%compute_averages(worker,bc_COMM,p_avg,c_avg,density_avg)

        ! interior velocity gradient
        grad1_v1 = -(v1_bc/density_m)*grad1_density_m  +  (ONE/density_m)*grad1_mom1_m

        ! interior pressure gradient
        dp_ddensity =  (gam-ONE)*HALF*(mom1_m*mom1_m + mom2_m*mom2_m + mom3_m*mom3_m)/(density_m*density_m)
        dp_dmom1    = -(gam-ONE)*mom1_m/density_m
        dp_dmom2    = -(gam-ONE)*mom2_m/density_m
        dp_dmom3    = -(gam-ONE)*mom3_m/density_m
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)

        grad1_p = dp_ddensity * grad1_density_m  + &
                  dp_dmom1    * grad1_mom1_m     + &
                  dp_dmom2    * grad1_mom2_m     + &
                  dp_dmom3    * grad1_mom3_m     + &
                  dp_denergy  * grad1_energy_m


        grad1_phi1 = density_m
        grad1_phi4 = density_m
        do i = 1,size(grad1_p)
            grad1_phi1(i) = -c_avg*c_avg*grad1_density_m(i)  +  grad1_p(i)
            grad1_phi4(i) = density_avg*c_avg*grad1_v1(i)    +  grad1_p(i)
        end do

        ! BC gradients
        gradn_p = density_m
        gradn_v1 = density_m
        gradn_density = density_m
        do i = 1,size(grad1_p)
            gradn_p(i)       = HALF*grad1_phi4(i)
            gradn_v1(i)      = (ONE/(TWO*density_avg*c_avg))*grad1_phi4(i)
            gradn_density(i) = -(ONE/(c_avg*c_avg))*grad1_phi1(i)  +  (ONE/(TWO*c_avg*c_avg))*grad1_phi4(i)
        end do


        ! Extrapolate density
        v1_o = mom1_o/density_o
        !delta_p_avg = (p_user(1) - p_avg)
        delta_p_avg = 1000._rk*(p_avg - p_user(1))
        !density_bc = density_o - delta_n*gradn_density
        !v1_bc      = v1_o      - delta_n*gradn_v1
        p_bc       = p_o       + delta_n*(gradn_p - delta_p_avg)

        ! Compute density, momentum, energy
        density_bc = p_bc/(Rgas*T_bc)
        mom1_bc    = v1_bc*density_bc
        mom2_bc    = v2_bc*density_bc
        mom3_bc    = v3_bc*density_bc
        energy_bc  = p_bc/(gam - ONE) + (density_bc*HALF)*(v1_bc*v1_bc + v2_bc*v2_bc + v3_bc*v3_bc)


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
    !*****************************************************************************************


end module bc_state_outlet_neumann_pressure_fd_test
