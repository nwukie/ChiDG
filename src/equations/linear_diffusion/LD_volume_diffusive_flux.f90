module LD_volume_diffusive_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_volume_flux,       only: volume_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use LD_properties,          only: LD_properties_t
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: LD_volume_diffusive_flux_t


    contains

        procedure   :: compute

    end type LD_volume_diffusive_flux_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(LD_volume_diffusive_flux_t),  intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        integer(ik)             :: iu
        real(rk)                :: mu_x, mu_y, mu_z

        type(AD_D), allocatable, dimension(:)   ::  &
            flux_x, flux_y, flux_z, dudx, dudy, dudz


        !
        ! Get variable index from equation set
        !
        iu = prop%get_eqn_index("u")


        !
        ! Get equation set properties
        !
        select type(prop)
            type is (LD_properties_t)
                mu_x = prop%mu(1)
                mu_y = prop%mu(2)
                mu_z = prop%mu(3)
        end select


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx = worker%interpolate(iu, 'ddx')
        dudy = worker%interpolate(iu, 'ddy')
        dudz = worker%interpolate(iu, 'ddz')



        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = -mu_x*dudx
        flux_y = -mu_y*dudy
        flux_z = -mu_z*dudz


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume(iu, flux_x, flux_y, flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LD_volume_diffusive_flux
