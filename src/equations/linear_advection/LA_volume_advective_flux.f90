module LA_volume_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

!    use LA_properties,          only: LA_properties_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: LA_volume_advective_flux_t


    contains
    
        procedure   :: init
        procedure   :: compute

    end type LA_volume_advective_flux_t
    !*************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(LA_volume_advective_flux_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("LA Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%set_equation("u")

    end subroutine init
    !********************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(LA_volume_advective_flux_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        integer(ik)             :: iu
        real(rk)                :: cx, cy, cz

        type(AD_D), allocatable, dimension(:)   ::  &
            u, flux_x, flux_y, flux_z



        !
        ! Get variable index from equation set
        !
        iu = prop%get_equation_index('u')


        !
        ! Get equation set properties
        !
!        select type(prop)
!            type is (LA_properties_t)
!                cx = prop%c(1)
!                cy = prop%c(2)
!                cz = prop%c(3)
!        end select
        cx = 1._rk
        cy = 0._rk
        cz = 0._rk



        !
        ! Interpolate solution to quadrature nodes
        !
        u = worker%interpolate(iu, 'value')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = cx  *  u 
        flux_y = cy  *  u
        flux_z = cz  *  u


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume(iu, flux_x, flux_y, flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LA_volume_advective_flux
