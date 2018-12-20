module fluid_laplacian_av_bc_operator
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: fluid_laplacian_av_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_laplacian_av_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_laplacian_av_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Fluid Laplacian AV BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_laplacian_av_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            grad1, grad2, grad3, pgrad1, pgrad2, pgrad3, av, &
            flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   ::  r


        
        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
        end if

        av = worker%get_field('Smoothed Artificial Viscosity', 'value', 'boundary')

        !=================================================
        ! Mass flux
        !=================================================
        grad1 = worker%get_field('Density'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density'   , 'grad3', 'boundary')
        flux_1 = -av*grad1
        flux_2 = -av*grad2
        flux_3 = -av*grad3
        

        call worker%integrate_boundary_condition('Density','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-1 flux
        !=================================================
        grad1 = worker%get_field('Momentum-1'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Momentum-1'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Momentum-1'   , 'grad3', 'boundary')
        flux_1 = -av*grad1
        flux_2 = -av*grad2
        flux_3 = -av*grad3
 
        call worker%integrate_boundary_condition('Momentum-1','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-2 flux
        !=================================================
        grad1 = worker%get_field('Momentum-2'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Momentum-2'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Momentum-2'   , 'grad3', 'boundary')
        flux_1 = -av*grad1
        flux_2 = -av*grad2
        flux_3 = -av*grad3
 

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * r
        end if

        call worker%integrate_boundary_condition('Momentum-2','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-3 flux
        !=================================================
        grad1 = worker%get_field('Momentum-3'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Momentum-3'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Momentum-3'   , 'grad3', 'boundary')
        flux_1 = -av*grad1
        flux_2 = -av*grad2
        flux_3 = -av*grad3
 

        call worker%integrate_boundary_condition('Momentum-3','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! Energy flux
        !=================================================
        grad1 = worker%get_field('Energy'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Energy'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Energy'   , 'grad3', 'boundary')
        pgrad1 = worker%get_field('Pressure Gradient - 1'   , 'value', 'boundary')
        pgrad2 = worker%get_field('Pressure Gradient - 2'   , 'value', 'boundary')
        pgrad3 = worker%get_field('Pressure Gradient - 3'   , 'value', 'boundary')
        flux_1 = -av*(grad1+pgrad1)
        flux_2 = -av*(grad2+pgrad2)
        flux_3 = -av*(grad3+pgrad3)
 
        if (any(ieee_is_nan(flux_1(:)%x_ad_))) print *, 'energy flux 1 bc is nan'
        if (any(ieee_is_nan(flux_2(:)%x_ad_))) print *, 'energy flux 2 bc is nan'
        if (any(ieee_is_nan(flux_3(:)%x_ad_))) print *, 'energy flux 3 bc is nan'

        call worker%integrate_boundary_condition('Energy','Diffusion',flux_1,flux_2,flux_3)

    end subroutine compute
    !**********************************************************************************************









end module fluid_laplacian_av_bc_operator
