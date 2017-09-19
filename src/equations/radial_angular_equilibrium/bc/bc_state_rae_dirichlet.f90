module bc_state_rae_dirichlet
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
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: rae_dirichlet_t

    contains

        procedure   :: init                 !< Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     !< boundary condition function implementation

    end type rae_dirichlet_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rae_dirichlet_t),   intent(inout) :: self
        
        !
        ! Set name, family
        !
        call self%set_name("RAE Dirichlet")
        call self%set_family("Outlet")


        !
        ! Add functions
        !
        call self%bcproperties%add('Pressure','Required')


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
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(rae_dirichlet_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop
        type(mpi_comm),             intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            p1_bc,      p2_bc,                      &
            grad1_p1_m, grad1_p2_m,                 & 
            grad2_p1_m, grad2_p2_m,                 & 
            grad3_p1_m, grad3_p2_m
            

        real(rk),   allocatable, dimension(:)   ::  p_bc

        print*, 'dirichlet - 1'

        !
        ! Get back pressure from function.
        !
        p_bc = self%bcproperties%compute('Pressure',worker%time(),worker%coords())


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        p1_bc = worker%get_field('Pressure-1', 'value', 'face interior')
        p2_bc = worker%get_field('Pressure-2', 'value', 'face interior')

        p1_bc = sqrt(p_bc)
        p2_bc = sqrt(p_bc)



        grad1_p1_m = worker%get_field('Pressure-1', 'grad1', 'face interior')
        grad2_p1_m = worker%get_field('Pressure-1', 'grad2', 'face interior')
        grad3_p1_m = worker%get_field('Pressure-1', 'grad3', 'face interior')

        grad1_p2_m = worker%get_field('Pressure-2', 'grad1', 'face interior')
        grad2_p2_m = worker%get_field('Pressure-2', 'grad2', 'face interior')
        grad3_p2_m = worker%get_field('Pressure-2', 'grad3', 'face interior')




        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Pressure-1', p1_bc, 'value')
        call worker%store_bc_state('Pressure-2', p2_bc, 'value')





        call worker%store_bc_state('Pressure-1', grad1_p1_m, 'grad1')
        call worker%store_bc_state('Pressure-1', grad2_p1_m, 'grad2')
        call worker%store_bc_state('Pressure-1', grad3_p1_m, 'grad3')
                                                
        call worker%store_bc_state('Pressure-2', grad1_p2_m, 'grad1')
        call worker%store_bc_state('Pressure-2', grad2_p2_m, 'grad2')
        call worker%store_bc_state('Pressure-2', grad3_p2_m, 'grad3')
                                                


        print*, 'dirichlet - 2'




    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_rae_dirichlet
