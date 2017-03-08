module bc_state_spalart_allmaras_farfield
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, RKTOL
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none
    


    !>  Symmetry boundary condition for Spalart-Allmaras turbulent working variable.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: spalart_allmaras_farfield_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type spalart_allmaras_farfield_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_farfield_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Spalart Allmaras Farfield')
        call self%set_family('Farfield')


        !
        ! Add turbulent inlet parameter and default value.
        !
        call self%bcproperties%add('Turbulent Viscosity Ratio', 'Required')
        call self%set_fcn_option('Turbulent Viscosity Ratio', 'val', 3._rk)

    end subroutine init
    !********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(spalart_allmaras_farfield_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   :: &
            density_m, mom1_m, mom2_m, mom3_m, normal_momentum, &
            density_nutilde_m, density_nutilde_bc,              &
            grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m

        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3, nutilde_nu
        logical,    allocatable, dimension(:)   :: inflow, outflow


        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density',    'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')

        density_nutilde_m       = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        grad1_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad1', 'face interior')
        grad2_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad2', 'face interior')
        grad3_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'grad3', 'face interior')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
        end if


        !
        ! Get User incoming boundary condition viscosity ratio
        !
        nutilde_nu = self%bcproperties%compute('Turbulent Viscosity Ratio', worker%time(), worker%coords() )


        !
        ! Get unit normal vector
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)


        !
        ! Determine which modes are inflow/outflow
        !
        normal_momentum = mom1_m*unorm_1 + mom2_m*unorm_2 + mom3_m*unorm_3

        inflow  = ( normal_momentum <= RKTOL )
        outflow = ( normal_momentum >  RKTOL )




        !
        ! Set boundary values for density * nutilde.
        !   1: Extrapolate all values (assuming outflow)
        !   2: Where inflow is detected, set from user specified parameter
        !
        density_nutilde_bc = density_nutilde_m
        where(inflow)
            density_nutilde_bc = density_m * (nutilde_nu * 0.00018_rk)
        end where



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density * NuTilde', density_nutilde_bc,'value')


        !
        ! Store boundary condition gradient - Zero Gradient
        !
        grad1_density_nutilde_m = ZERO
        grad2_density_nutilde_m = ZERO
        grad3_density_nutilde_m = ZERO
        call worker%store_bc_state('Density * NuTilde', grad1_density_nutilde_m, 'grad1')
        call worker%store_bc_state('Density * NuTilde', grad2_density_nutilde_m, 'grad2')
        call worker%store_bc_state('Density * NuTilde', grad3_density_nutilde_m, 'grad3')
                                                
    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_spalart_allmaras_farfield
