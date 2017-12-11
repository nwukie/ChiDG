module bc_state_P_neumann
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, HALF, NO_ID, ZETA_MIN
    use mod_interpolate,    only: interpolate_edge_autodiff
    use mod_fluid,          only: gam
    use mod_chidg_mpi,      only: IRANK
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_edge_info,     only: edge_info_t
    use type_seed,          only: seed_t
    use mpi_f08,            only: mpi_comm
    use DNAD_D
    use ieee_arithmetic
    implicit none



    !>
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2017
    !!
    !-----------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: P_neumann_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type P_neumann_t
    !***********************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2017
    !!
    !------------------------------------------------------------------------------------
    subroutine init(self)    
        class(P_neumann_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('P Neumann')
        call self%set_family('Scalar')

    end subroutine init
    !************************************************************************************






    !>  Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   12/5/2017
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(P_neumann_t),     intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        character(:),   allocatable             :: user_msg
        integer(ik)                             :: idomain_g_n, idomain_l_n,    &
                                                   ielement_g_n, ielement_l_n,  &
                                                   ChiID, eqn_ID_n, iedge, ifield
        type(AD_D), allocatable, dimension(:)   ::                                      &
            p_bc, grad1_p, grad2_p, grad3_p,                                            &
            grad1_p_e, grad2_p_e, grad3_p_e
                    


!        !
!        ! Get 'p' value from face interior to extrapolate
!        !
!        p_bc = worker%get_field('Pressure', 'value', 'face interior')
!
!
!        !
!        ! TODO: 
!        ! Get identifiers for interior problem from overset data
!        !
!        associate ( idomain_g  => worker%element_info%idomain_g, &
!                    idomain_l  => worker%element_info%idomain_l, &
!                    ielement_g => worker%element_info%ielement_g, &
!                    ielement_l => worker%element_info%ielement_l, &
!                    iface      => worker%iface )
!
!        ! Assuming chimera on ZETA_MIN. Assuming a single conforming donor
!        ChiID = worker%mesh%domain(idomain_l)%faces(ielement_l,ZETA_MIN)%ChiID
!        idomain_g_n  = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%idomain_g
!        idomain_l_n  = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%idomain_l
!        ielement_g_n = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%ielement_g
!        ielement_l_n = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%ielement_l
!        eqn_ID_n     = worker%mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(1)%eqn_ID
!
!
!        !
!        ! Determine correct edge to interpolate from on the element from 
!        ! the interior problem.
!        !
!        select case (iface)
!            case(1)
!                iedge = 11
!            case(2)
!                iedge = 12
!            case(3)
!                iedge = 5
!            case(4)
!                iedge = 7
!            case default
!                user_msg = "edges can only be determined for faces 1-4."
!                call chidg_signal_one(FATAL,user_msg,iface)
!        end select
!
!
!
!        edge_info     = edge_info_t(idomain_g_n, idomain_l_n, ielement_g_n, ielement_l_n, iedge)
!        seed          = seed_t(idomain_g_n, idomain_l_n, ielement_g_n, ielement_l_n, &
!                               neqns = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%neqns, &
!                               nterms_s = worker%mesh%domain(idomain_l_n)%elems(ielement_l_n)%nterms_s, &
!                               iproc = IRANK, &
!                               recv_comm = NO_ID, &
!                               recv_domain = NO_ID, &
!                               recv_element = NO_ID)
!
!        end associate
!
!        ! Interpolate "Density" along edge
!        ifield = worker%prop(eqn_ID_n)%get_primary_field_index("Density")
!        density_e       = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'value')
!        grad1_density_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad1')
!        grad2_density_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad2')
!        grad3_density_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad3')
!
!        ! Interpolate "Momentum-1" along edge
!        ifield = worker%prop(eqn_ID_n)%get_primary_field_index("Momentum-1")
!        mom1_e       = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'value')
!        grad1_mom1_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad1')
!        grad2_mom1_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad2')
!        grad3_mom1_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad3')
!
!        ! Interpolate "Momentum-2" along edge
!        ifield = worker%prop(eqn_ID_n)%get_primary_field_index("Momentum-2")
!        mom2_e       = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'value')
!        grad1_mom2_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad1')
!        grad2_mom2_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad2')
!        grad3_mom2_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad3')
!
!        ! Interpolate "Momentum-3" along edge
!        ifield = worker%prop(eqn_ID_n)%get_primary_field_index("Momentum-3")
!        mom3_e       = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'value')
!        grad1_mom3_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad1')
!        grad2_mom3_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad2')
!        grad3_mom3_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad3')
!
!        ! Interpolate "Energy" along edge
!        ifield = worker%prop(eqn_ID_n)%get_primary_field_index("Energy")
!        energy_e       = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'value')
!        grad1_energy_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad1')
!        grad2_energy_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad2')
!        grad3_energy_e = interpolate_edge_autodiff(worker%mesh,worker%solverdata%q,edge_info,seed,ifield,worker%itime,'grad3')
!
!
!        ! Compute Jacobians of pressure
!        !p = (gam - ONE)*(energy_e - HALF*(mom1_e*mom1_e + mom2_e*mom2_e + mom3_e*mom3_e)/density_e)
!        dp_ddensity =  (gam-ONE)*HALF*(mom1_e*mom1_e + mom2_e*mom2_e + mom3_e*mom3_e)/(density_e*density_e)
!        dp_dmom1    = -(gam-ONE)*mom1_e/density_e
!        dp_dmom2    = -(gam-ONE)*mom2_e/density_e
!        dp_dmom3    = -(gam-ONE)*mom3_e/density_e
!        dp_denergy  = dp_ddensity ! init storage
!        dp_denergy  =  (gam-ONE)
!
!        ! Compute pressure gradient using Chain-rule
!        grad1_p_e = dp_ddensity * grad1_density_e  + &
!                    dp_dmom1    * grad1_mom1_e     + &
!                    dp_dmom2    * grad1_mom2_e     + &
!                    dp_dmom3    * grad1_mom3_e     + &
!                    dp_denergy  * grad1_energy_e
!
!        grad2_p_e = dp_ddensity * grad2_density_e  + &
!                    dp_dmom1    * grad2_mom1_e     + &
!                    dp_dmom2    * grad2_mom2_e     + &
!                    dp_dmom3    * grad2_mom3_e     + &
!                    dp_denergy  * grad2_energy_e
!
!        grad3_p_e = dp_ddensity * grad3_density_e  + &
!                    dp_dmom1    * grad3_mom1_e     + &
!                    dp_dmom2    * grad3_mom2_e     + &
!                    dp_dmom3    * grad3_mom3_e     + &
!                    dp_denergy  * grad3_energy_e
!
!
!
!        !
!        ! grad1_p(istart:iend) = grad1_p_e
!        ! grad2_p(istart:iend) = grad2_p_e
!        ! grad3_p(istart:iend) = grad3_p_e
!        !
!
!
!
!        if ((worker%element_info%ielement_g == 1) .and. (worker%iface == 1)) then
!           p_bc(1) = 100000._rk
!        end if 
!
!
!        call worker%store_bc_state('Presure', p_bc,    'value')
!        call worker%store_bc_state('Presure', grad1_p, 'grad1')
!        call worker%store_bc_state('Presure', grad2_p, 'grad2')
!        call worker%store_bc_state('Presure', grad3_p, 'grad3')






    end subroutine compute_bc_state
    !************************************************************************************






end module bc_state_P_neumann
