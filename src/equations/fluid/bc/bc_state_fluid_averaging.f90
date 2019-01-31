module bc_state_fluid_averaging
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ME, ONE, ZERO, HALF, NO_ID
    use mod_fluid,          only: gam
    use mod_interpolate,    only: interpolate_face_autodiff
    use mpi_f08,            only: mpi_comm
    use type_mesh,          only: mesh_t
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_element_info,  only: element_info_t
    use DNAD_D
    implicit none


    !>  Define a boundary state function class for fluids that provides
    !!  averaging of the primitive variables.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2018
    !!
    !--------------------------------------------------------------------------------
    type, public, abstract, extends(bc_state_t) :: bc_fluid_averaging_t

    contains

        procedure   :: init_bc_coupling
        procedure   :: compute_averages

    end type bc_fluid_averaging_t
    !********************************************************************************



contains


    !>  Initialize boundary group coupling.
    !!
    !!  Call global coupling routine to initialize implicit coupling between each
    !!  element with every other element on the boundary, a result of averaging
    !!  operations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/18/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_comm)
        class(bc_fluid_averaging_t),    intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh
        integer(ik),                    intent(in)      :: group_ID
        type(mpi_comm),                 intent(in)      :: bc_comm

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !********************************************************************************



    !>  Update the area-averaged pressure for the boundary condition.
    !!
    !!  @author Nathan A. average_pressure
    !!  @date   3/31/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_averages(self,worker,bc_COMM, vel1_avg, vel2_avg, vel3_avg, density_avg, p_avg)
        class(bc_fluid_averaging_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_COMM
        type(AD_D),                     intent(inout)   :: vel1_avg
        type(AD_D),                     intent(inout)   :: vel2_avg
        type(AD_D),                     intent(inout)   :: vel3_avg
        type(AD_D),                     intent(inout)   :: density_avg
        type(AD_D),                     intent(inout)   :: p_avg

        type(element_info_t)    :: coupled_element

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3, energy, p, v1, v2, v3

        type(AD_D)  :: p_integral, v1_integral, v2_integral, v3_integral, density_integral, &
                       face_density, face_v1, face_v2, face_v3, face_p


        integer(ik) :: ipatch, iface_bc, idomain_l, ielement_l, iface, ierr, itime,         &
                       idensity, imom1, imom2, imom3, ienergy, group_ID, patch_ID, face_ID, &
                       icoupled, idomain_g_coupled, idomain_l_coupled, ielement_g_coupled,  &
                       ielement_l_coupled, iface_coupled, dof_start_coupled, coupled_iface

        real(rk),   allocatable,    dimension(:)    :: weights, areas, r
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

            !face_info = face_info_constructor(idomain_g_coupled,  &
            !                                  idomain_l_coupled,  &
            !                                  ielement_g_coupled, &
            !                                  ielement_l_coupled, &
            !                                  iface_coupled,      &
            !                                  dof_start_coupled)

            coupled_element = element_info_t(idomain_g    = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(icoupled),     &
                                             idomain_l    = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(icoupled),     &
                                             ielement_g   = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(icoupled),    &
                                             ielement_l   = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(icoupled),    &
                                             iproc        = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(icoupled),          &
                                             pelem_ID     = NO_ID,                                                                                          &
                                             eqn_ID       = NO_ID,                                                                                          &
                                             nfields      = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%neqns(icoupled),         &
                                             nterms_s     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nterms_s(icoupled),      &
                                             nterms_c     = 0,                                                                                              &
                                             dof_start    = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_start(icoupled),     &
                                             recv_comm    = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_comm(icoupled),     &
                                             recv_domain  = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_domain(icoupled),   &
                                             recv_element = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_element(icoupled))

            coupled_iface = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%iface(icoupled)




            !
            ! Interpolate coupled element solution on face of coupled element
            !
            density = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, idensity, itime, 'value', ME)
            mom1    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom1,    itime, 'value', ME)
            mom2    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom2,    itime, 'value', ME)
            mom3    = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, imom3,    itime, 'value', ME)
            energy  = interpolate_face_autodiff(worker%mesh,worker%solverdata%q,coupled_element, worker%function_info, coupled_iface, ienergy,  itime, 'value', ME)

            !r = worker%coordinate('1','boundary')
            !if (worker%coordinate_system() == 'Cylindrical') then
            !    mom2 = mom2 / r
            !end if
            if (worker%coordinate_system() == 'Cylindrical') then
                mom2 = mom2 / worker%mesh%domain(idomain_l_coupled)%elems(ielement_l_coupled)%interp_coords_def(:,1)
            end if

            ! Compute quantities for averaging
            v1 = mom1/density
            v2 = mom2/density
            v3 = mom3/density
            p = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)

            ! Get weights + areas
            weights   = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights_face(iface_coupled)
            areas     = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%areas
            face_area = worker%mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%data(icoupled)%total_area

            ! Integrate and contribute to average
            face_density = sum(density * areas * weights)
            face_v1      = sum(v1      * areas * weights)
            face_v2      = sum(v2      * areas * weights)
            face_v3      = sum(v3      * areas * weights)
            face_p       = sum(p       * areas * weights)

            ! Allocate derivatives and clear integral for first face.
            if (icoupled == 1) then
                density_integral = face_v1
                v1_integral      = face_v1
                v2_integral      = face_v1
                v3_integral      = face_v1
                p_integral       = face_v1
                density_integral = ZERO
                v1_integral      = ZERO
                v2_integral      = ZERO
                v3_integral      = ZERO
                p_integral       = ZERO
            end if

            ! Accumulate face contribution.
            density_integral = density_integral + face_density
            v1_integral      = v1_integral      + face_v1
            v2_integral      = v2_integral      + face_v2
            v3_integral      = v3_integral      + face_v3
            p_integral       = p_integral       + face_p

            ! Accumulate surface area
            total_area = total_area + face_area

        end do !icoupled

        ! Compute average pressure:
        !   area-weighted pressure integral over the total area
        vel1_avg    = v1_integral      / total_area
        vel2_avg    = v2_integral      / total_area
        vel3_avg    = v3_integral      / total_area
        density_avg = density_integral / total_area
        p_avg       = p_integral       / total_area

    end subroutine compute_averages
    !********************************************************************************



end module bc_state_fluid_averaging
