module LD_volume_diffusive_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use LD_properties,          only: LD_properties_t
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: LD_volume_diffusive_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type LD_volume_diffusive_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(LD_volume_diffusive_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Linear Diffusion Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Flux")

        ! Set operator equations
        call self%set_equation("u")

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(LD_volume_diffusive_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        integer(ik)             :: iu
!        real(rk)                :: mu_x, mu_y, mu_z

        type(AD_D), allocatable, dimension(:)   ::  &
            flux_x, flux_y, flux_z, dudx, dudy, dudz


        !
        ! Get variable index from equation set
        !
        iu = prop%get_equation_index("u")


!        !
!        ! Get equation set properties
!        !
!        select type(prop)
!            type is (LD_properties_t)
!                mu_x = prop%mu(1)
!                mu_y = prop%mu(2)
!                mu_z = prop%mu(3)
!        end select


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx = worker%get_element_variable(iu, 'ddx + lift')
        dudy = worker%get_element_variable(iu, 'ddy + lift')
        dudz = worker%get_element_variable(iu, 'ddz + lift')



        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = -dudx
        flux_y = -dudy
        flux_z = -dudz


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume(iu, flux_x, flux_y, flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LD_volume_diffusive_operator
