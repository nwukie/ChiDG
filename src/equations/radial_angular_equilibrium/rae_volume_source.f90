module rae_volume_cylindrical_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for rae_ equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: rae_volume_cylindrical_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rae_volume_cylindrical_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rae_volume_cylindrical_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("RAE Volume Cylindrical Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Pressure-1")
        call self%add_primary_field("Pressure-2")

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for rae_ equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rae_volume_cylindrical_source_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            p1, p2, density, v, source 

        real(rk),   allocatable :: r(:)


        print*, 'volume source 1'

        !=================================================
        ! momentum-1 flux
        !=================================================
        if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Get solution at quadrature nodes
            !
            p1 = worker%get_field('Pressure-1','value','element')
            p2 = worker%get_field('Pressure-2','value','element')


            !
            ! Get radial coordinate
            !
            r = worker%coordinate('1','volume')


            !
            ! Get model fields
            !
            density = worker%get_field('Density',    'value', 'element')
            v       = worker%get_field('Velocity-2', 'value', 'element')



            !
            ! Source term due to transformation to cylindrical coordinates
            !   Stationary reference frame
            !       source = (density*v*v + p)/r
            !
            !   Translating reference frame
            !       source = [(density*v*v + p) - (density*v)*Ugrid_2]/r
            !
            !   Translating/Rotating/deforming reference frame
            !       source = [ale_g*ale_Dinv(2,2)*(density*v*v + p) - ale_Dinv(2,2)*Ugrid_2*(density*v)]/r
            !
            source = (density*v*v + p1*p2)/r
            print*, 'source: ', source(:)%x_ad_

            call worker%integrate_volume_source('Pressure-1',source)

        end if


        print*, 'volume source 2'

    end subroutine compute
    !*********************************************************************************************************






end module rae_volume_cylindrical_source
