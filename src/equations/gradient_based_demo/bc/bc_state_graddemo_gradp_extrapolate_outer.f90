module bc_state_graddemo_gradp_extrapolate_outer
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: graddemo_gradp_extrapolate_outer_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type graddemo_gradp_extrapolate_outer_t
    !********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_gradp_extrapolate_outer_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("graddemo gradp extrapolate outer")
        call self%set_family("Wall")

    end subroutine init
    !********************************************************************************





    !>  Compute routine for Pressure Outlet boundary condition state function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !---------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(graddemo_gradp_extrapolate_outer_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            pressure, grad1_pressure, grad2_pressure, grad3_pressure,   &
            grad1_p, grad2_p, grad3_p,                                  &
            grad1_grad1_p, grad2_grad1_p, grad3_grad1_p,                &
            grad1_grad2_p, grad2_grad2_p, grad3_grad2_p,                &
            grad1_grad3_p, grad2_grad3_p, grad3_grad3_p


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        pressure = worker%get_field('Pressure_TEMP',         'value', 'face interior')
        grad1_p  = worker%get_field('Pressure Gradient - 1', 'value', 'face interior')
        grad2_p  = worker%get_field('Pressure Gradient - 2', 'value', 'face interior')
        grad3_p  = worker%get_field('Pressure Gradient - 3', 'value', 'face interior')


        grad1_pressure = worker%get_field('Pressure_TEMP'        , 'grad1', 'face interior')
        grad2_pressure = worker%get_field('Pressure_TEMP'        , 'grad2', 'face interior')
        grad3_pressure = worker%get_field('Pressure_TEMP'        , 'grad3', 'face interior')

        grad1_grad1_p  = worker%get_field('Pressure Gradient - 1', 'grad1', 'face interior')
        grad2_grad1_p  = worker%get_field('Pressure Gradient - 1', 'grad2', 'face interior')
        grad3_grad1_p  = worker%get_field('Pressure Gradient - 1', 'grad3', 'face interior')

        grad1_grad2_p  = worker%get_field('Pressure Gradient - 2', 'grad1', 'face interior')
        grad2_grad2_p  = worker%get_field('Pressure Gradient - 2', 'grad2', 'face interior')
        grad3_grad2_p  = worker%get_field('Pressure Gradient - 2', 'grad3', 'face interior')

        grad1_grad3_p  = worker%get_field('Pressure Gradient - 3', 'grad1', 'face interior')
        grad2_grad3_p  = worker%get_field('Pressure Gradient - 3', 'grad2', 'face interior')
        grad3_grad3_p  = worker%get_field('Pressure Gradient - 3', 'grad3', 'face interior')

        
        !
        ! We want the gradient to be equal to the gradient from the interior problem
        !
        !grad1_pressure =  -132000._rk*sin(132._rk * worker%x('boundary'))
        grad1_pressure = grad1_p
        grad2_pressure = grad2_p
        grad3_pressure = ZERO

!        grad1_pressure = ZERO
!        grad2_pressure = ZERO
!        grad3_pressure = ZERO


        !
        ! We want to set the value of pressure at a single point
        !
!        if ( (worker%element_info%idomain_g == 6) .and. (worker%element_info%ielement_g == 1) .and. (worker%iface == 1) ) then
!            pressure(1) = 100000._rk
!            !pressure = 98000._rk
!        end if

        if ( ((worker%element_info%idomain_g == 6) .or. (worker%element_info%idomain_g == 7)) .and. &
             (worker%iface == 1) ) then
            grad1_pressure = ZERO
            grad2_pressure = ZERO
            grad3_pressure = ZERO
        end if

        if ( ((worker%element_info%idomain_g == 6) .or. (worker%element_info%idomain_g == 7)) .and. &
             (worker%iface == 2) ) then
            grad1_pressure = ZERO
            grad2_pressure = ZERO
            grad3_pressure = ZERO
        end if

!        if ( (worker%element_info%idomain_g == 6)  .and. &
!             (worker%element_info%ielement_g == 1) .and. &
!             (worker%iface == 3) ) then
!            pressure(1) = 100000._rk
!            !pressure = 98000._rk
!        end if


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Pressure_TEMP', pressure,       'value')
        call worker%store_bc_state('Pressure_TEMP', grad1_pressure, 'grad1')
        call worker%store_bc_state('Pressure_TEMP', grad2_pressure, 'grad2')
        call worker%store_bc_state('Pressure_TEMP', grad3_pressure, 'grad3')


        call worker%store_bc_state('Pressure Gradient - 1', grad1_p,        'value')
        call worker%store_bc_state('Pressure Gradient - 2', grad2_p,        'value')
        call worker%store_bc_state('Pressure Gradient - 3', grad3_p,        'value')

        call worker%store_bc_state('Pressure Gradient - 1', grad1_grad1_p,  'grad1')
        call worker%store_bc_state('Pressure Gradient - 1', grad2_grad1_p,  'grad2')
        call worker%store_bc_state('Pressure Gradient - 1', grad3_grad1_p,  'grad3')

        call worker%store_bc_state('Pressure Gradient - 2', grad1_grad2_p,  'grad1')
        call worker%store_bc_state('Pressure Gradient - 2', grad2_grad2_p,  'grad2')
        call worker%store_bc_state('Pressure Gradient - 2', grad3_grad2_p,  'grad3')

        call worker%store_bc_state('Pressure Gradient - 3', grad1_grad3_p,  'grad1')
        call worker%store_bc_state('Pressure Gradient - 3', grad2_grad3_p,  'grad2')
        call worker%store_bc_state('Pressure Gradient - 3', grad3_grad3_p,  'grad3')


    end subroutine compute_bc_state
    !*******************************************************************************






end module bc_state_graddemo_gradp_extrapolate_outer
