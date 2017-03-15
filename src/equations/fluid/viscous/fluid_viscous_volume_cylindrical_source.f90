module fluid_viscous_volume_cylindrical_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !>  Volume source terms from viscous fluxes in cylindrical coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/22/2017
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: fluid_viscous_volume_cylindrical_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_volume_cylindrical_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/22/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_volume_cylindrical_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Fluid Viscous Volume Cylindrical Source')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Flux')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for viscous equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_viscous_volume_cylindrical_source_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::    &
            tau_22, source 





        !=================================================
        ! mass flux
        !=================================================


        !=================================================
        ! momentum-1 flux
        !=================================================
        if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! get shear stress
            !
            tau_22 = worker%get_model_field_element('Shear-22','value')

            source = -tau_22 / worker%coordinate('1','volume')

            call worker%integrate_volume('Momentum-1',source)

        end if


        !=================================================
        ! momentum-2 flux
        !=================================================


        !=================================================
        ! momentum-3 flux
        !=================================================


        !=================================================
        ! energy flux
        !=================================================



    end subroutine compute
    !*********************************************************************************************************






end module fluid_viscous_volume_cylindrical_source
