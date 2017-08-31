module bc_state_artificial_viscosity_symmetry
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, TEN
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
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: artificial_viscosity_symmetry_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type artificial_viscosity_symmetry_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(artificial_viscosity_symmetry_t),   intent(inout) :: self
        

        !
        ! Set name, family
        !
        call self%set_name('Artificial Viscosity Symmetry')
        call self%set_family('Symmetry')


    end subroutine init
    !********************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(artificial_viscosity_symmetry_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop
        type(mpi_comm),                         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::      &
            eps_m,  deps_dx_m,  deps_dy_m,  deps_dz_m



        !
        ! Interpolate interior solution to quadrature nodes
        !
        eps_m      = worker%get_field('Artificial Viscosity' , 'value', 'face interior')
        deps_dx_m  = worker%get_field('Artificial Viscosity' , 'grad1', 'face interior')
        deps_dy_m  = worker%get_field('Artificial Viscosity' , 'grad2', 'face interior')
        deps_dz_m  = worker%get_field('Artificial Viscosity' , 'grad3', 'face interior')


        !
        ! Store computed boundary state
        !
        call worker%store_bc_state('Artificial Viscosity' , eps_m,     'value')
        call worker%store_bc_state('Artificial Viscosity' , deps_dx_m, 'grad1')
        call worker%store_bc_state('Artificial Viscosity' , deps_dy_m, 'grad2')
        call worker%store_bc_state('Artificial Viscosity' , deps_dz_m, 'grad3')
                                                

    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_artificial_viscosity_symmetry
