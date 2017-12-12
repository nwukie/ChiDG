module graddemo_volume_operator
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
    type, extends(operator_t), public :: graddemo_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Graddemo Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Flux")

        ! Set operator equations
        call self%add_primary_field("Density")
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy")
        call self%add_primary_field("Pressure_TEMP")
        call self%add_primary_field("Pressure Gradient - 1")
        call self%add_primary_field("Pressure Gradient - 2")
        call self%add_primary_field("Pressure Gradient - 3")

    end subroutine init
    !******************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!-----------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(graddemo_volume_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                        &
            density, mom1, mom2, mom3, energy,                          &
            density_face, mom1_face, mom2_face, mom3_face, energy_face, &
            density_element, mom1_element, mom2_element, mom3_element,  &
            density_residual, mom1_residual, mom2_residual, mom3_residual, energy_residual, &
            density_residual_modes, mom1_residual_modes, mom2_residual_modes, mom3_residual_modes, energy_residual_modes,   &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,                  &
            pressure, pressure_residual, pressure_residual_modes, pressure_modes,   &
            grad1_p,                grad2_p,                grad3_p,                &
            grad1_p_face,           grad2_p_face,           grad3_p_face,           &
            grad1_p_element,        grad2_p_element,        grad3_p_element,        &
            grad1_p_residual,       grad2_p_residual,       grad3_p_residual,       &
            grad1_p_residual_modes, grad2_p_residual_modes, grad3_p_residual_modes, &
            grad1_density_face, grad1_mom1_face, grad1_mom2_face, grad1_mom3_face, grad1_energy_face,   &
            grad2_density_face, grad2_mom1_face, grad2_mom2_face, grad2_mom3_face, grad2_energy_face,   &
            grad3_density_face, grad3_mom1_face, grad3_mom2_face, grad3_mom3_face, grad3_energy_face

        type(AD_D)  :: average_pressure

        integer(ik) :: islice, istart, iend, nnodes1d


        !
        ! Get Primary fields: boundary
        !
        density  = worker%get_field('Density',               'value', 'element')
        mom1     = worker%get_field('Momentum-1',            'value', 'element')
        mom2     = worker%get_field('Momentum-2',            'value', 'element')
        mom3     = worker%get_field('Momentum-3',            'value', 'element')
        energy   = worker%get_field('Energy',                'value', 'element')
        pressure = worker%get_field('Pressure_TEMP',         'value', 'element')
        grad1_p  = worker%get_field('Pressure Gradient - 1', 'value', 'element')
        grad2_p  = worker%get_field('Pressure Gradient - 2', 'value', 'element')
        grad3_p  = worker%get_field('Pressure Gradient - 3', 'value', 'element')


        
        !
        ! Get Primary fields: interior
        !
        density_face = worker%get_field('Density',               'value', 'face exterior', iface=5)
        mom1_face    = worker%get_field('Momentum-1',            'value', 'face exterior', iface=5)
        mom2_face    = worker%get_field('Momentum-2',            'value', 'face exterior', iface=5)
        mom3_face    = worker%get_field('Momentum-3',            'value', 'face exterior', iface=5)
        grad1_p_face = worker%get_field('Pressure Gradient - 1', 'value', 'face exterior', iface=5)
        grad2_p_face = worker%get_field('Pressure Gradient - 2', 'value', 'face exterior', iface=5)
        grad3_p_face = worker%get_field('Pressure Gradient - 3', 'value', 'face exterior', iface=5)



        !
        ! Initialize element storage
        !
        density_element = density
        mom1_element    = density
        mom2_element    = density
        mom3_element    = density
        grad1_p_element = density
        grad2_p_element = density
        grad3_p_element = density

        nnodes1d = worker%nnodes1d(density)
        do islice = 1,nnodes1d
            istart = 1 + (nnodes1d*nnodes1d)*(islice-1)
            iend = istart + (nnodes1d*nnodes1d - 1)
            density_element(istart:iend)  = density_face
            mom1_element(istart:iend)     = mom1_face
            mom2_element(istart:iend)     = mom2_face
            mom3_element(istart:iend)     = mom3_face
            grad1_p_element(istart:iend)  = grad1_p_face
            grad2_p_element(istart:iend)  = grad2_p_face
            grad3_p_element(istart:iend)  = grad3_p_face
        end do



        !
        ! Compute residual
        !
        density_residual = -(density  - density_element)
        mom1_residual    = -(mom1     - mom1_element)
        mom2_residual    = -(mom2     - mom2_element)
        mom3_residual    = -(mom3     - mom3_element)
        energy_residual  = -(energy   - ( (pressure/(gam-ONE)) + HALF*(mom1_element*mom1_element + mom2_element*mom2_element + mom3_element*mom3_element)/density_element))
        !energy_residual = -(energy   - ( (100000._rk/(gam-ONE)) + HALF*(mom1_element*mom1_element + mom2_element*mom2_element + mom3_element*mom3_element)/density_element))
        grad1_p_residual = -(grad1_p  - grad1_p_element)
        grad2_p_residual = -(grad2_p  - grad2_p_element)
        grad3_p_residual = -(grad3_p  - grad3_p_element)
!        pressure_residual = -(pressure - 100000._rk)
        
        !
        ! Project residual
        !
!        density_residual_modes  = worker%project_from_nodes(density_residual)
!        mom1_residual_modes     = worker%project_from_nodes(mom1_residual)
!        mom2_residual_modes     = worker%project_from_nodes(mom2_residual)
!        mom3_residual_modes     = worker%project_from_nodes(mom3_residual)
!        energy_residual_modes   = worker%project_from_nodes(energy_residual)
!        grad1_p_residual_modes  = worker%project_from_nodes(grad1_p_residual)
!        grad2_p_residual_modes  = worker%project_from_nodes(grad2_p_residual)
!        grad3_p_residual_modes  = worker%project_from_nodes(grad3_p_residual)
        !pressure_residual_modes = worker%project_from_nodes(pressure_residual)



!        if ( (worker%element_info%idomain_g == 2) .and. (worker%element_info%ielement_g == 1) ) then
!            average_pressure = sum(pressure)/real(size(pressure),rk)
!            print*, 'average pressure: ', average_pressure%x_ad_
!
!            pressure_residual = pressure
!            pressure_residual(:) = (average_pressure - 100000._rk)
!            pressure_residual_modes = worker%project_from_nodes(pressure_residual)
!
!!            ! We only want to affect the constant mode
!!            pressure_residual = pressure_modes
!!            print*, 'first pressure mode: ', pressure_modes(1)%x_ad_
!!            pressure_residual(1)  = (pressure_modes(1) - 100000._rk)
!!            pressure_residual(2:) = ZERO
!        else
!            pressure_residual_modes = density_residual_modes
!            pressure_residual_modes = ZERO
!        end if


        !=================================================
        !                      Store
        !=================================================
        call worker%integrate_volume_source('Density',       density_residual )
        call worker%integrate_volume_source('Momentum-1',    mom1_residual    )
        call worker%integrate_volume_source('Momentum-2',    mom2_residual    )
        call worker%integrate_volume_source('Momentum-3',    mom3_residual    )
        call worker%integrate_volume_source('Energy',        energy_residual  )
        call worker%integrate_volume_source('Pressure Gradient - 1', grad1_p_residual)
        call worker%integrate_volume_source('Pressure Gradient - 2', grad2_p_residual)
        call worker%integrate_volume_source('Pressure Gradient - 3', grad3_p_residual)
!        call worker%integrate_volume_source('Pressure_TEMP', pressure_residual)
!        call worker%accumulate_residual('Density',               density_residual_modes )
!        call worker%accumulate_residual('Momentum-1',            mom1_residual_modes    )
!        call worker%accumulate_residual('Momentum-2',            mom2_residual_modes    )
!        call worker%accumulate_residual('Momentum-3',            mom3_residual_modes    )
!        call worker%accumulate_residual('Energy',                energy_residual_modes  )
!        call worker%accumulate_residual('Pressure Gradient - 1', grad1_p_residual_modes )
!        call worker%accumulate_residual('Pressure Gradient - 2', grad2_p_residual_modes )
!        call worker%accumulate_residual('Pressure Gradient - 3', grad3_p_residual_modes )

!        call worker%accumulate_residual('Pressure',              pressure_residual_modes)


    end subroutine compute
    !******************************************************************************






end module graddemo_volume_operator
