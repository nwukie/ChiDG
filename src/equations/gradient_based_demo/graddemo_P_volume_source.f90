module graddemo_P_volume_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/18/2018
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: graddemo_P_volume_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_P_volume_source_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   12/18/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_P_volume_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Graddemo P Volume Source')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   12/18/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(graddemo_P_volume_source_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            p, grad1_pbc, grad2_pbc, grad3_pbc,     &
            grad1_p,   grad2_p,   grad3_p,          &
            source
            
        type(AD_D)  :: p_avg

!        !if ( (worker%element_info%idomain_g == 5) .and. (worker%element_info%ielement_g == 3) ) then
!        !if ( (worker%element_info%ielement_g == 160) .and. (worker%iface == 6) ) then
!        !if ( (worker%element_info%ielement_g == 1) .and. (worker%iface == 6) ) then
!        if ( (worker%element_info%ielement_g == 40) ) then
!
!            !
!            ! Interpolate solution to quadrature nodes
!            !
!            p = worker%get_field('Pressure_TEMP', 'value', 'element')
!
!            p_avg = sum(p)/real(size(p),rk)
!            print*, 'Pressures: ', p(:)%x_ad_
!            print*, 'Avg Pressure: ', p_avg%x_ad_
!
!            source = p
!            source(:) = -(100000._rk - p_avg)
!
!
!            call worker%integrate_volume_source('Pressure_TEMP', source)
!
!        end if

    end subroutine compute
    !****************************************************************************************************






end module graddemo_P_volume_source
