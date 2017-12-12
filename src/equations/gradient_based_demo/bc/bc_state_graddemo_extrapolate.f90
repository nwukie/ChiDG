module bc_state_graddemo_extrapolate
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: graddemo_extrapolate_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type graddemo_extrapolate_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_extrapolate_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("graddemo extrapolate")
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
        class(graddemo_extrapolate_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
        type(mpi_comm),                 intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                              &
            density, mom1, mom2, mom3, energy,                                  &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,    &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,    &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy
        

        !
        ! Interpolate interior solution to face quadrature nodes
        !
        density = worker%get_field('Density',    'value', 'face interior')
        mom1    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy  = worker%get_field('Energy',     'value', 'face interior')


        grad1_density = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3    = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3', 'face interior')

        grad1_energy  = worker%get_field('Energy'    , 'grad1', 'face interior')
        grad2_energy  = worker%get_field('Energy'    , 'grad2', 'face interior')
        grad3_energy  = worker%get_field('Energy'    , 'grad3', 'face interior')


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density',    density, 'value')
        call worker%store_bc_state('Momentum-1', mom1,    'value')
        call worker%store_bc_state('Momentum-2', mom2,    'value')
        call worker%store_bc_state('Momentum-3', mom3,    'value')
        call worker%store_bc_state('Energy',     energy,  'value')
                                                



        call worker%store_bc_state('Density'   , grad1_density, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', grad1_mom1,    'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1,    'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1,    'grad3')
                                                
        call worker%store_bc_state('Momentum-2', grad1_mom2,    'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2,    'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2,    'grad3')
                                                
        call worker%store_bc_state('Momentum-3', grad1_mom3,    'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3,    'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3,    'grad3')

        call worker%store_bc_state('Energy'    , grad1_energy,  'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy,  'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy,  'grad3')



    end subroutine compute_bc_state
    !*******************************************************************************






end module bc_state_graddemo_extrapolate
