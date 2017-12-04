module tm_volume_cylindrical_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF, ZERO
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
            p, density, invdensity, u, v, vmag, t,  &
            mom1,           mom2,                   &
            grad1_density,  grad2_density,          &
            grad1_mom1,     grad2_mom1,             &
            grad1_mom2,     grad2_mom2,             &
            grad1_u,        grad2_u,                &
            grad1_v,        grad2_v,                &
            du_ddensity,    dv_ddensity,            &
            du_dmom1,       dv_dmom2,               &
            source_1, source_2, source


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
            r = worker%coordinate('1','element')


            !
            ! Get model fields
            !
            density = worker%get_field('Density',    'value', 'element')
            mom1    = worker%get_field('Momentum-1', 'value', 'element')
            mom2    = worker%get_field('Momentum-2', 'value', 'element')


            grad1_density = worker%get_field('Density : Grad1',    'value', 'element')
            grad2_density = worker%get_field('Density : Grad2',    'value', 'element')

            grad1_mom1    = worker%get_field('Momentum-1 : Grad1', 'value', 'element')
            grad2_mom1    = worker%get_field('Momentum-1 : Grad2', 'value', 'element')

            grad1_mom2    = worker%get_field('Momentum-2 : Grad1', 'value', 'element')
            grad2_mom2    = worker%get_field('Momentum-2 : Grad2', 'value', 'element')


!            if (worker%coordinate_system() == 'Cylindrical') then
!                r = worker%coordinate('1')
!                mom2       = mom2 / r
!                grad1_mom2 = (grad1_mom2/r) - mom2/r
!                grad2_mom2 = (grad2_mom2/r)
!            end if


            !
            ! Compute velocities
            !
            u = mom1 / density
            v = mom2 / density


            !
            ! compute velocity jacobians
            !
            invdensity  = ONE/density
            du_ddensity = -invdensity*invdensity*mom1
            dv_ddensity = -invdensity*invdensity*mom2

            du_dmom1 = invdensity
            dv_dmom2 = invdensity


            !
            ! compute velocity gradients via chain rule:
            !
            !   u = f(rho,rhou)
            !
            !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
            !
            grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
            grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1

            grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
            grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2





            !
            ! Compute weighting coefficients
            !
            source_1 = - ( (u*grad1_mom1 + density*u*grad1_u) + &
                           (v*grad2_mom1 + density*u*grad2_v) )
            source_2 = - ( (v*grad1_mom1 + density*u*grad1_v) + &
                           (v*grad2_mom2 + density*v*grad2_v) )

            if (worker%coordinate_system() == 'Cylindrical') then
                r = worker%coordinate('1','element')
                source_1 = source_1  +  (density*v*v + p)/r  -  (density*u*u/r)
                source_2 = source_2  -  (density*u*v)/r      -  (density*u*v/r)
            end if

            t = source_1/(source_1 + source_2)





            t = ZERO ! All tangential



            ! Source term
            source = t*source_1  +  (ONE-t)*source_2

            call worker%integrate_volume_source('Pressure',source)

        end if






!        type(AD_D), allocatable, dimension(:)   ::  &
!            p, density, mom1, mom2, u, v, source 
!
!        real(rk),   allocatable :: r(:)
!
!
!
!        !=================================================
!        ! momentum-1 flux
!        !=================================================
!        if (worker%coordinate_system() == 'Cylindrical') then
!
!            !
!            ! Get solution at quadrature nodes
!            !
!            p = worker%get_field('Pressure','value','element')
!
!
!            !
!            ! Get radial coordinate
!            !
!            r = worker%coordinate('1','volume')
!
!
!            !
!            ! Get model fields
!            !
!            density = worker%get_field('Density',    'value', 'element')
!            mom1    = worker%get_field('Momentum-1', 'value', 'element')
!            mom2    = worker%get_field('Momentum-2', 'value', 'element')
!
!            ! Get tangential momentum
!            mom2 = mom2/r
!
!            u = mom1/density
!            v = mom2/density
!
!
!
!            !
!            ! Source term due to transformation to cylindrical coordinates
!            !   Stationary reference frame
!            !       source = (density*v*v + p)/r
!            !
!            !   Translating reference frame
!            !       source = [(density*v*v + p) - (density*v)*Ugrid_2]/r
!            !
!            !   Translating/Rotating/deforming reference frame
!            !       source = [ale_g*ale_Dinv(2,2)*(density*v*v + p) - ale_Dinv(2,2)*Ugrid_2*(density*v)]/r
!            !
!            source = -density*u*v/r
!
!            call worker%integrate_volume_source('Pressure',source)
!
!        end if



    end subroutine compute
    !*********************************************************************************************************






end module tm_volume_cylindrical_source
