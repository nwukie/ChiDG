module euler_volume_cylindrical_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
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
    type, extends(operator_t), public :: euler_volume_cylindrical_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_volume_cylindrical_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_volume_cylindrical_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Euler Volume Cylindrical Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_volume_cylindrical_source_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::    &
            density, mom2, v, p, source 

        real(rk),   allocatable :: r(:), ale_g(:), grid_vel(:,:), ale_Dinv(:,:,:)




        !=================================================
        ! mass flux
        !=================================================


        !=================================================
        ! momentum-1 flux
        !=================================================
        if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Get solution at quadrature nodes
            !
            density = worker%get_field('Density'   ,'value','element')
            mom2    = worker%get_field('Momentum-2','value','element')


            !
            ! Get grid motion/deformation information
            !
            grid_vel = worker%get_grid_velocity_element()
            ale_g    = worker%get_det_jacobian_grid_element('value')
            ale_Dinv = worker%get_inv_jacobian_grid_element()


            !
            ! Compute V_theta from angular momentum, get pressure:
            !
            r = worker%coordinate('1','volume')
            v = mom2/(r*density)
            p = worker%get_field('Pressure','value','element')


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
            !source = ( ale_g*ale_Dinv(2,2,:)*(density*v*v + p)  -  ale_Dinv(2,2,:)*grid_vel(:,2)*(density*v) )/r
            source = ( ale_g*ale_Dinv(2,2,:)*(density*v*v + p) )/r

            call worker%integrate_volume_source('Momentum-1',source)

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






end module euler_volume_cylindrical_source
