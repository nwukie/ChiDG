module pressure_boundary_source
#include <messenger.h>
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
    type, extends(operator_t), public :: pressure_boundary_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type pressure_boundary_source_t
    !******************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(pressure_boundary_source_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Pressure Boundary Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("u")

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
        class(pressure_boundary_source_t),  intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::    &
            u, v, p, source 

        real(rk),   allocatable :: r(:), ale_g(:), grid_vel(:,:), ale_Dinv(:,:,:)



            !
            ! Get solution at quadrature nodes
            !
            u = worker%get_field('u','value','element')


            !
            ! Compute V_theta from angular momentum, get pressure:
            !
            r = worker%coordinate('1','volume')

            source = ( ale_g*ale_Dinv(2,2,:)*(density*v*v + p)  -  ale_Dinv(2,2,:)*grid_vel(:,2)*(density*v) )/r

            call worker%integrate_volume_source('u',source)


    end subroutine compute
    !*********************************************************************************************************






end module pressure_boundary_source
