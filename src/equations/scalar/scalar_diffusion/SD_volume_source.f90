module SD_volume_source
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,FOUR,PI

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use SD_properties,          only: SD_properties_t
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: SD_volume_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_volume_source_t
    !*************************************************************************

contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_volume_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Scalar Diffusion Volume Source Standard")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Source")

        ! Set operator equations
        call self%add_primary_field("u")

    end subroutine init
    !********************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SD_volume_source_t),          intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        integer(ik)                             :: iu
        type(AD_D), allocatable, dimension(:)   :: source
        real(rk),   allocatable, dimension(:)   :: x, y

        !
        ! Get variable index from equation set
        !
        iu = prop%get_primary_field_index("u")


        !
        ! Interpolate solution to quadrature nodes
        !
        source = worker%get_element_variable(iu, 'ddx + lift')


        !x = worker%x('volume')
        !source = FOUR*PI*PI*dsin(TWO*PI*x)


        x = worker%x('volume')
        y = worker%y('volume')
        source = TWO*PI*PI*dsin(PI*x)*dsin(PI*y)



        !
        ! Integrate volume flux
        !
        call worker%integrate_volume(iu, source)


    end subroutine compute
    !****************************************************************************************************






end module SD_volume_source
