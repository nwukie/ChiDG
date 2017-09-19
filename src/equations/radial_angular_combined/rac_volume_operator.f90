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
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rac_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("RAC Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Pressure")

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
        class(rac_volume_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            p, density, u, v,                       &
            flux_1, flux_2, flux_3

        !
        ! Interpolate solution to quadrature nodes
        !
        p = worker%get_field('Pressure', 'value', 'element')


        !
        ! Get model fields
        !
        density = worker%get_field('Density',    'value', 'element')
        u       = worker%get_field('Velocity-1', 'value', 'element')
        v       = worker%get_field('Velocity-2', 'value', 'element')


        !=================================================
        !                   Momentum-1
        !=================================================
        flux_1 = density*u*(u + v)  +  p
        flux_2 = density*v*(u + v)  +  p
        flux_3 = density*u
        flux_3 = ZERO

        call worker%integrate_volume_flux('Pressure','Advection',flux_1,flux_2,flux_3)



    end subroutine compute
    !*********************************************************************************************************






end module rac_volume_operator
