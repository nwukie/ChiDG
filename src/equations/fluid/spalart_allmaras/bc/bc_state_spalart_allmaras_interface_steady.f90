module bc_state_spalart_allmaras_interface_steady
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, THREE, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: gam, Rgas, cp

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    !use bc_giles_HB_base,       only: giles_HB_base_t
    use bc_nonlocal_nrbc_lindblad_base,       only: nonlocal_nrbc_lindblad_base_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use ieee_arithmetic,        only: ieee_is_nan
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none



    !>  Name: inlet - 3D Giles
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!  Behavior:
    !!      
    !!  References:
    !!              
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, extends(nonlocal_nrbc_lindblad_base_t) :: spalart_allmaras_interface_steady_t

    contains

        procedure   :: init                         ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state             ! boundary condition function implementation

    end type spalart_allmaras_interface_steady_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_interface_steady_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Steady Spalart Allmaras Interface')
        call self%set_family('Inlet')

        ! Add functions
        call self%bcproperties%add('Pitch A', 'Required')
        call self%bcproperties%add('Pitch B', 'Required')

    end subroutine init
    !********************************************************************************





    !>  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_comm)
        class(spalart_allmaras_interface_steady_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_comm


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                                      &
            density_nutilde_m, grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m, density_m, mu_m, nu_m, density_nutilde_bc

        character(1)    :: side



        if (self%get_face_side(worker) == 'A') then
            ! Interpolate interior solution to face quadrature nodes
            density_nutilde_m       = worker%get_field('Density * NuTilde', 'value', 'face interior')
            grad1_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad1', 'face interior')
            grad2_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad2', 'face interior')
            grad3_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad3', 'face interior')


            ! Outlet: extrapolate value, zero gradient
            grad1_density_nutilde_m = ZERO
            grad2_density_nutilde_m = ZERO
            grad3_density_nutilde_m = ZERO
            call worker%store_bc_state('Density * NuTilde', density_nutilde_m,       'value')
            call worker%store_bc_state('Density * NuTilde', grad1_density_nutilde_m, 'grad1')
            call worker%store_bc_state('Density * NuTilde', grad2_density_nutilde_m, 'grad2')
            call worker%store_bc_state('Density * NuTilde', grad3_density_nutilde_m, 'grad3')



        else if (self%get_face_side(worker) == 'B') then

            grad1_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad1', 'face interior')
            grad2_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad2', 'face interior')
            grad3_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad3', 'face interior')

            ! Inlet: extrapolate gradient
            call worker%store_bc_state('Density * NuTilde', grad1_density_nutilde_m, 'grad1')
            call worker%store_bc_state('Density * NuTilde', grad2_density_nutilde_m, 'grad2')
            call worker%store_bc_state('Density * NuTilde', grad3_density_nutilde_m, 'grad3')

            density_m = worker%get_field('Density',           'value', 'face interior')
            mu_m      = worker%get_field('Laminar Viscosity', 'value', 'face interior')
            nu_m = mu_m/density_m
            !density_nutilde_bc = density_m * (nutilde_nu * nu_m)
            density_nutilde_bc = density_m * (THREE * nu_m)

            call worker%store_bc_state('Density * NuTilde', density_nutilde_bc,'value')

        end if


    end subroutine compute_bc_state
    !*********************************************************************************




end module bc_state_spalart_allmaras_interface_steady
