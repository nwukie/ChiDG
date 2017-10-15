module rac_volume_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,ZERO

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
    type, extends(operator_t), public :: rac_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rac_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !------------------------------------------------------------------------------
    subroutine init(self)
        class(rac_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("RAC Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Pressure")

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
        class(rac_volume_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            p, density, invdensity, u, v, vmag, t,  &
            mom1,           mom2,                   &
            grad1_density,  grad2_density,          &
            grad1_mom1,     grad2_mom1,             &
            grad1_mom2,     grad2_mom2,             &
            grad1_u,        grad2_u,                &
            grad1_v,        grad2_v,                &
            du_ddensity,    dv_ddensity,            &
            du_dmom1,       dv_dmom2,               &
            source_1, source_2, source,             &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r

        !
        ! Interpolate solution to quadrature nodes
        !
        p = worker%get_field('Pressure', 'value', 'element')


        !
        ! Get model fields
        !
        density = worker%get_field('Density',    'value', 'element')
        mom1    = worker%get_field('Momentum-1', 'value', 'element')
        mom2    = worker%get_field('Momentum-2', 'value', 'element')


        grad1_density = worker%get_field('Density',    'grad1', 'element')
        grad2_density = worker%get_field('Density',    'grad2', 'element')

        grad1_mom1    = worker%get_field('Momentum-1', 'grad1', 'element')
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2', 'element')

        grad1_mom2    = worker%get_field('Momentum-2', 'grad1', 'element')
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2', 'element')


        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
        end if


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




        !=================================================
        !                   Momentum-1
        !=================================================
        flux_1 = t*p
        flux_2 = (ONE-t)*p
        flux_3 = density
        flux_3 = ZERO


        call worker%integrate_volume_flux('Pressure','Advection',flux_1,flux_2,flux_3)



    end subroutine compute
    !******************************************************************************






end module rac_volume_operator
