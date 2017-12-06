module bc_state_P_neumann
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
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


        type(AD_D), allocatable, dimension(:)   :: u_bc, dudx_bc, dudy_bc, dudz_bc


!        !
!        ! Get 'u' value from face interior to extrapolate
!        !
!        p_bc = worker%get_field('Pressure', 'value', 'face interior')
!
!
!        !
!        ! TODO: Probably dangerous in parallel
!        !
!        idomain_g = 
!        idomain_l = 
!        ielement_g = 
!        ielement_l = 
!
!
!
!        grad1_density = worker%get_field('Density', 'grad1', 'edge exterior', edge=int(edge_selection))
!        grad2_density = worker%get_field('Density', 'grad2', 'edge exterior', edge=int(edge_selection))
!        grad3_density = worker%get_field('Density', 'grad3', 'edge exterior', edge=int(edge_selection))
!
!
!        ! Compute Jacobians of pressure
!        !p = (gam - ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)
!        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
!        dp_dmom1    = -(gam-ONE)*mom1/density
!        dp_dmom2    = -(gam-ONE)*mom2/density
!        dp_dmom3    = -(gam-ONE)*mom3/density
!        dp_denergy  =  (gam-ONE)
!
!        ! Compute pressure gradient using Chain-rule
!        grad1_p = dp_ddensity * grad1_density    + &
!                  dp_dmom1    * grad1_mom1       + &
!                  dp_dmom2    * grad1_mom2       + &
!                  dp_dmom3    * grad1_mom3       + &
!                  dp_denergy  * grad1_energy
!
!        grad2_p = dp_ddensity * grad2_density    + &
!                  dp_dmom1    * grad2_mom1       + &
!                  dp_dmom2    * grad2_mom2       + &
!                  dp_dmom3    * grad2_mom3       + &
!                  dp_denergy  * grad2_energy
!
!        grad3_p = dp_ddensity * grad3_density    + &
!                  dp_dmom1    * grad3_mom1       + &
!                  dp_dmom2    * grad3_mom2       + &
!                  dp_dmom3    * grad3_mom3       + &
!                  dp_denergy  * grad3_energy
!
!
!
!
!
!
!        if ((worker%element_info%ielement_g == 1) .and. (worker%iface == 1)) then
!        !if (worker%iface == 1) then
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
