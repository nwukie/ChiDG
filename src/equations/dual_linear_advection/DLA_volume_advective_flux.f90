module DLA_volume_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none



    !> This equation set exists really just to test equationsets with more than one equation. 
    !! The idea is just to compute the linear advecdtion solution twice at the same time. 
    !! The equations are independent of each other. So, we can verify, for example, the volume 
    !! flux jacobians for each equation. They should be the same as for the single 
    !! LinearAdvection equation set.
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: DLA_volume_advective_flux_t

    contains

        procedure   :: init
        procedure   :: compute

    end type DLA_volume_advective_flux_t
    !***********************************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(DLA_volume_advective_flux_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("DLA Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("u_a")
        call self%add_primary_field("u_b")

    end subroutine init
    !********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(DLA_volume_advective_flux_t),     intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            ua, ub, flux_1, flux_2, flux_3

        real(rk) :: c1, c2, c3


        !
        ! Get equation set properties
        !
        c1 = 1._rk
        c2 = 0._rk
        c3 = 0._rk


        !
        ! Interpolate solution to quadrature nodes
        !
        ua = worker%get_field('u_a', 'value', 'element')
        ub = worker%get_field('u_b', 'value', 'element')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = (c1 * ua)
        flux_2 = (c2 * ua)
        flux_3 = (c3 * ua)

        call worker%integrate_volume_flux('u_a','Advection',flux_1,flux_2,flux_3)



        ! Compute volume flux at quadrature nodes
        flux_1 = (c1 * ub)
        flux_2 = (c2 * ub)
        flux_3 = (c3 * ub)

        call worker%integrate_volume_flux('u_b','Advection',flux_1,flux_2,flux_3)




    end subroutine compute
    !***************************************************************************************************




end module DLA_volume_advective_flux
