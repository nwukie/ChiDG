module tm_volume_operator
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
    type, extends(operator_t), public :: tm_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type tm_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(tm_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("TM Volume Flux")

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
        class(tm_volume_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            p, density, u, v, mom1, mom2,           &
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

        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','element')
            mom2 = mom2/r
        end if


        u = mom1/density
        v = mom2/density



        !=================================================
        !                   Momentum-1
        !=================================================
        flux_1 = density
        flux_2 = density*v*v  +  p
        flux_3 = density
        flux_1 = ZERO
        flux_3 = ZERO

        call worker%integrate_volume_flux('Pressure','Advection',flux_1,flux_2,flux_3)



    end subroutine compute
    !*********************************************************************************************************






end module tm_volume_operator
