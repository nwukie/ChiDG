module graddemo_gradp_volume_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,ZERO
    use mod_fluid,              only: gam

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: graddemo_gradp_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_gradp_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2017
    !!
    !------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_gradp_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Graddemo GradP Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")
        !call self%set_operator_type("Volume Diffusive Flux")

        ! Set operator equations
        call self%add_primary_field("Pressure_TEMP")
!        call self%add_primary_field("Pressure Gradient - 1")
!        call self%add_primary_field("Pressure Gradient - 2")
!        call self%add_primary_field("Pressure Gradient - 3")

    end subroutine init
    !******************************************************************************



    !>  Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2017
    !!
    !!-----------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(graddemo_gradp_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                                        &
            density, mom1, mom2, mom3, energy,                                          &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,                      &
            grad1_p,                grad2_p,                grad3_p,                    &
            grad1_p_computed,       grad2_p_computed,       grad3_p_computed,           &
            grad1_p_element,        grad2_p_element,        grad3_p_element,            &
            grad1_p_residual,       grad2_p_residual,       grad3_p_residual,           &
            grad1_p_residual_modes, grad2_p_residual_modes, grad3_p_residual_modes,     &
            pressure, pressure_computed, pressure_residual, pressure_residual_modes,    &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,            &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,            &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy



        !
        ! Get Primary fields: boundary
        !
        density  = worker%get_field('Density',               'value', 'element')
        mom1     = worker%get_field('Momentum-1',            'value', 'element')
        mom2     = worker%get_field('Momentum-2',            'value', 'element')
        mom3     = worker%get_field('Momentum-3',            'value', 'element')
        energy   = worker%get_field('Energy',                'value', 'element')
        pressure = worker%get_field('Pressure_TEMP',         'value', 'element')
!        grad1_p  = worker%get_field('Pressure Gradient - 1', 'value', 'element')
!        grad2_p  = worker%get_field('Pressure Gradient - 2', 'value', 'element')
!        grad3_p  = worker%get_field('Pressure Gradient - 3', 'value', 'element')


        
!        !
!        ! Get Primary fields: interior
!        !
!        grad1_density  = worker%get_field('Density',    'grad1', 'element', override_lift=.true.)
!        grad1_mom1     = worker%get_field('Momentum-1', 'grad1', 'element', override_lift=.true.)
!        grad1_mom2     = worker%get_field('Momentum-2', 'grad1', 'element', override_lift=.true.)
!        grad1_mom3     = worker%get_field('Momentum-3', 'grad1', 'element', override_lift=.true.)
!        grad1_energy   = worker%get_field('Energy',     'grad1', 'element', override_lift=.true.)
!
!        grad2_density  = worker%get_field('Density',    'grad2', 'element', override_lift=.true.)
!        grad2_mom1     = worker%get_field('Momentum-1', 'grad2', 'element', override_lift=.true.)
!        grad2_mom2     = worker%get_field('Momentum-2', 'grad2', 'element', override_lift=.true.)
!        grad2_mom3     = worker%get_field('Momentum-3', 'grad2', 'element', override_lift=.true.)
!        grad2_energy   = worker%get_field('Energy',     'grad2', 'element', override_lift=.true.)
!
!        grad3_density  = worker%get_field('Density',    'grad3', 'element', override_lift=.true.)
!        grad3_mom1     = worker%get_field('Momentum-1', 'grad3', 'element', override_lift=.true.)
!        grad3_mom2     = worker%get_field('Momentum-2', 'grad3', 'element', override_lift=.true.)
!        grad3_mom3     = worker%get_field('Momentum-3', 'grad3', 'element', override_lift=.true.)
!        grad3_energy   = worker%get_field('Energy',     'grad3', 'element', override_lift=.true.)
!
!
!        ! Compute Jacobians of pressure
!        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
!        dp_dmom1    = -(gam-ONE)*mom1/density
!        dp_dmom2    = -(gam-ONE)*mom2/density
!        dp_dmom3    = -(gam-ONE)*mom3/density
!        dp_denergy  = dp_ddensity ! init storage
!        dp_denergy  =  (gam-ONE)



        ! Compute pressure gradient using Chain-rule
        !pressure_computed = (gam - ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)
        pressure_computed = worker%get_field('Pressure', 'value', 'element')

!        grad1_p_computed = dp_ddensity * grad1_density  + &
!                           dp_dmom1    * grad1_mom1     + &
!                           dp_dmom2    * grad1_mom2     + &
!                           dp_dmom3    * grad1_mom3     + &
!                           dp_denergy  * grad1_energy
!
!        grad2_p_computed = dp_ddensity * grad2_density  + &
!                           dp_dmom1    * grad2_mom1     + &
!                           dp_dmom2    * grad2_mom2     + &
!                           dp_dmom3    * grad2_mom3     + &
!                           dp_denergy  * grad2_energy
!
!        grad3_p_computed = dp_ddensity * grad3_density  + &
!                           dp_dmom1    * grad3_mom1     + &
!                           dp_dmom2    * grad3_mom2     + &
!                           dp_dmom3    * grad3_mom3     + &
!                           dp_denergy  * grad3_energy



        !
        ! Compute residual
        !
        pressure_residual = -(pressure - pressure_computed)
        !pressure_residual = (pressure - sqrt(pressure_computed))
!        grad1_p_residual  = -(grad1_p  - grad1_p_computed)
!        grad2_p_residual  = -(grad2_p  - grad2_p_computed)
!        grad3_p_residual  = -(grad3_p  - grad3_p_computed)
        
        !
        ! Project residual
        !
        pressure_residual_modes = worker%project_from_nodes(pressure_residual)
!        grad1_p_residual_modes  = worker%project_from_nodes(grad1_p_residual)
!        grad2_p_residual_modes  = worker%project_from_nodes(grad2_p_residual)
!        grad3_p_residual_modes  = worker%project_from_nodes(grad3_p_residual)



        !=================================================
        !                      Store
        !=================================================
        call worker%integrate_volume_source('Pressure_TEMP',pressure_residual)
!        call worker%accumulate_residual('Pressure',              pressure_residual_modes)
!        call worker%accumulate_residual('Pressure Gradient - 1', grad1_p_residual_modes )
!        call worker%accumulate_residual('Pressure Gradient - 2', grad2_p_residual_modes )
!        call worker%accumulate_residual('Pressure Gradient - 3', grad3_p_residual_modes )



    end subroutine compute
    !******************************************************************************






end module graddemo_gradp_volume_operator
