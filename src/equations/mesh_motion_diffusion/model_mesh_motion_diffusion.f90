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

        call self%set_name('Mesh Motion : Diffusion')
!        call self%set_dependency('f(Q-)')

        call self%add_model_field('Scalar Diffusion Coefficient')

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
            d, grad1_d, grad2_d, grad3_d, d_normalization, sumsqr, rho

        real(rk) :: p


        ! Get primary field to initialize derivatives
        rho     = worker%get_primary_field_general('Density', 'value')
        d       = rho
        grad1_d = rho
        grad2_d = rho
        grad3_d = rho


        !
        ! Interpolate solution to quadrature nodes
        !
        d       = worker%get_auxiliary_field_general('Mesh Motion : Diffusion', 'value')
        grad1_d = worker%get_auxiliary_field_general('Mesh Motion : Diffusion', 'grad1')
        grad2_d = worker%get_auxiliary_field_general('Mesh Motion : Diffusion', 'grad2')
        grad3_d = worker%get_auxiliary_field_general('Mesh Motion : Diffusion', 'grad3')




        d_normalization = d

        !
        ! Store Wall Distance to model field
        !
        call worker%store_model_field('Scalar Diffusion Coefficient','value', d_normalization)


    end subroutine compute
    !***************************************************************************************




end module model_mesh_motion_diffusion
