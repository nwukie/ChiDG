module bc_state_graddemo_gradp_extrapolate
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, TWO
    use mod_fluid,              only: gam
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: graddemo_gradp_extrapolate_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type graddemo_gradp_extrapolate_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_gradp_extrapolate_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("graddemo gradp extrapolate")
        call self%set_family("Wall")

    end subroutine init
    !********************************************************************************





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
        class(graddemo_gradp_extrapolate_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            po, grad1_po, grad2_po, grad3_po,       &
            p_bc, grad1_p_bc, grad2_p_bc, grad3_p_bc

        type(AD_D), dimension(:),   allocatable ::                              &
            density, mom1, mom2, mom3, energy,                                  &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,    &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,    &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy,    &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy





        !
        ! Retrieve values from original problem
        !
        po = worker%get_field('Pressure', 'value', 'face interior')
        !grad1_po = worker%get_field('Pressure Gradient - 1', 'value', 'face interior')
        !grad2_po = worker%get_field('Pressure Gradient - 2', 'value', 'face interior')
        !grad3_po = worker%get_field('Pressure Gradient - 3', 'value', 'face interior')



        !
        ! Interpolate solution to quadrature nodes
        !
        density       = worker%get_field('Density',    'value')
        mom1          = worker%get_field('Momentum-1', 'value')
        mom2          = worker%get_field('Momentum-2', 'value')
        mom3          = worker%get_field('Momentum-3', 'value')
        energy        = worker%get_field('Energy',     'value')

        grad1_density = worker%get_field('Density',    'grad1', override_lift=.true.)
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', override_lift=.true.)
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', override_lift=.true.)
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', override_lift=.true.)
        grad1_energy  = worker%get_field('Energy',     'grad1', override_lift=.true.)


        grad2_density = worker%get_field('Density',    'grad2', override_lift=.true.)
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', override_lift=.true.)
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', override_lift=.true.)
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', override_lift=.true.)
        grad2_energy  = worker%get_field('Energy',     'grad2', override_lift=.true.)


        grad3_density = worker%get_field('Density',    'grad3', override_lift=.true.)
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', override_lift=.true.)
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', override_lift=.true.)
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', override_lift=.true.)
        grad3_energy  = worker%get_field('Energy',     'grad3', override_lift=.true.)



        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
        dp_dmom1    = -(gam-ONE)*mom1/density
        dp_dmom2    = -(gam-ONE)*mom2/density
        dp_dmom3    = -(gam-ONE)*mom3/density
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)


        ! Compute pressure gradient using Chain-rule
        grad1_po = dp_ddensity * grad1_density  + &
                   dp_dmom1    * grad1_mom1     + &
                   dp_dmom2    * grad1_mom2     + &
                   dp_dmom3    * grad1_mom3     + &
                   dp_denergy  * grad1_energy

        grad2_po = dp_ddensity * grad2_density  + &
                   dp_dmom1    * grad2_mom1     + &
                   dp_dmom2    * grad2_mom2     + &
                   dp_dmom3    * grad2_mom3     + &
                   dp_denergy  * grad2_energy

        grad3_po = dp_ddensity * grad3_density  + &
                   dp_dmom1    * grad3_mom1     + &
                   dp_dmom2    * grad3_mom2     + &
                   dp_dmom3    * grad3_mom3     + &
                   dp_denergy  * grad3_energy



        !
        ! Set boundary condition on auxiliary problem from original problem
        !
        p_bc       = po
        grad1_p_bc = grad1_po
        grad2_p_bc = grad2_po
        grad3_p_bc = grad3_po



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Pressure_TEMP', p_bc,       'value')
        call worker%store_bc_state('Pressure_TEMP', grad1_p_bc, 'grad1')
        call worker%store_bc_state('Pressure_TEMP', grad2_p_bc, 'grad2')
        call worker%store_bc_state('Pressure_TEMP', grad3_p_bc, 'grad3')




    end subroutine compute_bc_state
    !*******************************************************************************






end module bc_state_graddemo_gradp_extrapolate
