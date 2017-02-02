module bc_state_artificial_viscosity_outlet
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, TEN, ZERO
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: artificial_viscosity_outlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type artificial_viscosity_outlet_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_outlet_t),   intent(inout) :: self
        

        !
        ! Set name, family
        !
        call self%set_name('Artificial Viscosity Outlet')
        call self%set_family('Outlet')


    end subroutine init
    !********************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(artificial_viscosity_outlet_t),   intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            eps_m,  eps_bc,                         &
            deps_dx_m,  deps_dy_m,  deps_dz_m,      &
            deps_dx_bc, deps_dy_bc, deps_dz_bc,     &
            deriv_mag


        real(rk),   allocatable, dimension(:)   ::  &
            unormx, unormy, unormz, h



        !
        ! Interpolate interior solution to quadrature nodes
        !
        eps_m      = worker%get_primary_field_face('Artificial Viscosity' , 'value', 'face interior')
        deps_dx_m  = worker%get_primary_field_face('Artificial Viscosity' , 'ddx',   'face interior')
        deps_dy_m  = worker%get_primary_field_face('Artificial Viscosity' , 'ddy',   'face interior')
        deps_dz_m  = worker%get_primary_field_face('Artificial Viscosity' , 'ddz',   'face interior')



        !
        ! Extrapolate value
        !
        eps_bc = eps_m


        !
        ! Get unit normal vector
        !
        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)

        
        !
        ! Get approximate element bounding box size
        !
        h = worker%element_size('interior')


        deriv_mag = -eps_m/(TEN*sqrt((h(1)*unormx)**TWO + (h(2)*unormy)**TWO + (h(3)*unormz)**TWO) )
        deps_dx_m = deriv_mag * unormx
        deps_dy_m = deriv_mag * unormy
        deps_dz_m = deriv_mag * unormz


        !
        ! Store computed boundary state
        !
        call worker%store_bc_state('Artificial Viscosity' , eps_bc,    'value')

        deps_dx_m = ZERO
        deps_dy_m = ZERO
        deps_dz_m = ZERO
        call worker%store_bc_state('Artificial Viscosity' , deps_dx_m, 'ddx'  )
        call worker%store_bc_state('Artificial Viscosity' , deps_dy_m, 'ddy'  )
        call worker%store_bc_state('Artificial Viscosity' , deps_dz_m, 'ddz'  )
                                                

    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_artificial_viscosity_outlet
