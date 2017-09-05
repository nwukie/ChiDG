module euler_volume_operator
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
    type, extends(operator_t), public :: euler_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Euler Volume Flux")

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
        class(euler_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            p, u, v, w, enthalpy, invdensity,       &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:) :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density'   , 'value', 'element')
        mom1    = worker%get_field('Momentum-1', 'value', 'element')
        mom2    = worker%get_field('Momentum-2', 'value', 'element')
        mom3    = worker%get_field('Momentum-3', 'value', 'element')
        energy  = worker%get_field('Energy'    , 'value', 'element')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','element') 
            mom2 = mom2 / r
        end if

    
        !
        ! Compute velocity
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity


        !
        ! Compute pressure and total enthalpy
        !
        p = worker%get_field('Pressure', 'value', 'element')
        enthalpy = (energy + p)*invdensity


        !===========================
        !        Mass flux
        !===========================
        flux_1 = (density * u)
        flux_2 = (density * v)
        flux_3 = (density * w)

        call worker%integrate_volume_flux('Density','Advection',flux_1,flux_2,flux_3)

        !===========================
        !       Momentum-1 flux
        !===========================
        flux_1 = (mom1 * u)  +  p
        flux_2 = (mom1 * v)
        flux_3 = (mom1 * w)
        
        call worker%integrate_volume_flux('Momentum-1','Advection',flux_1,flux_2,flux_3)

        !============================
        !       Momentum-2 flux
        !============================
        flux_1 = (mom2 * u)
        flux_2 = (mom2 * v)  +  p
        flux_3 = (mom2 * w)
        
        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if

        call worker%integrate_volume_flux('Momentum-2','Advection',flux_1,flux_2,flux_3)

        !============================
        !       Momentum-3 flux
        !============================
        flux_1 = (mom3 * u)
        flux_2 = (mom3 * v)
        flux_3 = (mom3 * w)  +  p

        call worker%integrate_volume_flux('Momentum-3','Advection',flux_1,flux_2,flux_3)

        !============================
        !       Energy flux
        !============================
        flux_1 = (density * enthalpy * u)
        flux_2 = (density * enthalpy * v)
        flux_3 = (density * enthalpy * w)

        call worker%integrate_volume_flux('Energy','Advection',flux_1,flux_2,flux_3)

    end subroutine compute
    !*********************************************************************************************************






end module euler_volume_operator
