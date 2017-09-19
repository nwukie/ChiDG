module tm_volume_cylindrical_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for tm_ equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: tm_volume_cylindrical_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type tm_volume_cylindrical_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(tm_volume_cylindrical_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("TM Volume Cylindrical Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Pressure")

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for tm_ equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(tm_volume_cylindrical_source_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            p, density, u, v, source 

        real(rk),   allocatable :: r(:)



        !=================================================
        ! momentum-1 flux
        !=================================================
        if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Get solution at quadrature nodes
            !
            p = worker%get_field('Pressure','value','element')


            !
            ! Get radial coordinate
            !
            r = worker%coordinate('1','volume')


            !
            ! Get model fields
            !
            density = worker%get_field('Density',    'value', 'element')
            u       = worker%get_field('Velocity-1', 'value', 'element')
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
            source = -density*u*v/r

            call worker%integrate_volume_source('Pressure',source)

        end if



    end subroutine compute
    !*********************************************************************************************************






end module tm_volume_cylindrical_source
