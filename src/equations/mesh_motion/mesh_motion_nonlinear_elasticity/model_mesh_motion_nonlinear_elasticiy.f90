module model_mesh_motion_diffusion
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO, RKTOL
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
!    use eqn_wall_distance,  only: get_p_poisson_parameter
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Wall Distance
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: mesh_motion_diffusion_m

    contains

        procedure   :: init
        procedure   :: compute

    end type mesh_motion_diffusion_m
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(mesh_motion_diffusion_m), intent(inout)   :: self

        call self%set_name('Mesh Motion Diffusion')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Mesh Motion Grid Displacement 1')
        call self%add_model_field('Mesh Motion Grid Displacement 2')
        call self%add_model_field('Mesh Motion Grid Displacement 3')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Wall Distance using a normalization of the output
    !!  from a p-Poisson equation based method.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(mesh_motion_diffusion_m),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            d1, d2, d3, rho 


        ! Get primary field to initialize derivatives
        rho     = worker%get_field('Density', 'value')
        d1       = rho
        d2       = rho
        d3       = rho

        !
        ! Interpolate solution to quadrature nodes
        !
        ! Is this the right way to handle multiple components?
        d1       = worker%get_auxiliary_field_general('grid_displacement1', 'value')
        d2       = worker%get_auxiliary_field_general('grid_displacement2', 'value')
        d3       = worker%get_auxiliary_field_general('grid_displacement3', 'value')
        



        !
        ! Store grid displacements to model field
        !
        call worker%store_model_field('Mesh Motion Grid Displacement 1','value', d1)
        call worker%store_model_field('Mesh Motion Grid Displacement 2','value', d2)
        call worker%store_model_field('Mesh Motion Grid Displacement 3','value', d3)


    end subroutine compute
    !***************************************************************************************




end module model_mesh_motion_diffusion
