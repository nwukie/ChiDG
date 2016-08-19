module LINEULER_volume_advective_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,ZERO

    use type_mesh,              only: mesh_t
    use type_volume_flux,       only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t
    
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none

    private


    !>  Volume advective flux for Linearized Euler equations - imaginary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: LINEULER_volume_advective_flux_imag_t


    contains

        procedure  :: compute

    end type LINEULER_volume_advective_flux_imag_t
    !***********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,elem_info,function_info)
        class(LINEULER_volume_advective_flux_imag_t),   intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        type(element_info_t),                           intent(in)      :: elem_info
        type(function_info_t),                          intent(in)      :: function_info


        ! Equation indices
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe



        integer(ik)    :: idom, ielem


        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE,                        &
                    flux_x, flux_y, flux_z



        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l



        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho_i")
        irhou = prop%get_eqn_index("rhou_i")
        irhov = prop%get_eqn_index("rhov_i")
        irhow = prop%get_eqn_index("rhow_i")
        irhoE = prop%get_eqn_index("rhoE_i")




        !
        ! Interpolate solution to quadrature nodes
        !
        rho  = interpolate(mesh,sdata,elem_info,function_info,irho, 'value')
        rhou = interpolate(mesh,sdata,elem_info,function_info,irhou,'value')
        rhov = interpolate(mesh,sdata,elem_info,function_info,irhov,'value')
        rhow = interpolate(mesh,sdata,elem_info,function_info,irhow,'value')
        rhoE = interpolate(mesh,sdata,elem_info,function_info,irhoE,'value')



        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rho_x_rho  * rho  + &
                 rho_x_rhou * rhou + &
                 rho_x_rhov * rhov + &
                 rho_x_rhow * rhow + &
                 rho_x_rhoE * rhoE
        flux_y = rho_y_rho  * rho  + &
                 rho_y_rhou * rhou + &
                 rho_y_rhov * rhov + &
                 rho_y_rhow * rhow + &
                 rho_y_rhoE * rhoE
        flux_z = rho_z_rho  * rho  + &
                 rho_z_rhou * rhou + &
                 rho_z_rhov * rhov + &
                 rho_z_rhow * rhow + &
                 rho_z_rhoE * rhoE

        call integrate_volume_flux(mesh,sdata,elem_info,function_info,irho,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = rhou_x_rho  * rho  + &
                 rhou_x_rhou * rhou + &
                 rhou_x_rhov * rhov + &
                 rhou_x_rhow * rhow + &
                 rhou_x_rhoE * rhoE
        flux_y = rhou_y_rho  * rho  + &
                 rhou_y_rhou * rhou + &
                 rhou_y_rhov * rhov + &
                 rhou_y_rhow * rhow + &
                 rhou_y_rhoE * rhoE
        flux_z = rhou_z_rho  * rho  + &
                 rhou_z_rhou * rhou + &
                 rhou_z_rhov * rhov + &
                 rhou_z_rhow * rhow + &
                 rhou_z_rhoE * rhoE

        call integrate_volume_flux(mesh,sdata,elem_info,function_info,irhou,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = rhov_x_rho  * rho  + &
                 rhov_x_rhou * rhou + &
                 rhov_x_rhov * rhov + &
                 rhov_x_rhow * rhow + &
                 rhov_x_rhoE * rhoE
        flux_y = rhov_y_rho  * rho  + &
                 rhov_y_rhou * rhou + &
                 rhov_y_rhov * rhov + &
                 rhov_y_rhow * rhow + &
                 rhov_y_rhoE * rhoE
        flux_z = rhov_z_rho  * rho  + &
                 rhov_z_rhou * rhou + &
                 rhov_z_rhov * rhov + &
                 rhov_z_rhow * rhow + &
                 rhov_z_rhoE * rhoE

        call integrate_volume_flux(mesh,sdata,elem_info,function_info,irhov,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = rhow_x_rho  * rho  + &
                 rhow_x_rhou * rhou + &
                 rhow_x_rhov * rhov + &
                 rhow_x_rhow * rhow + &
                 rhow_x_rhoE * rhoE
        flux_y = rhow_y_rho  * rho  + &
                 rhow_y_rhou * rhou + &
                 rhow_y_rhov * rhov + &
                 rhow_y_rhow * rhow + &
                 rhow_y_rhoE * rhoE
        flux_z = rhow_z_rho  * rho  + &
                 rhow_z_rhou * rhou + &
                 rhow_z_rhov * rhov + &
                 rhow_z_rhow * rhow + &
                 rhow_z_rhoE * rhoE


        call integrate_volume_flux(mesh,sdata,elem_info,function_info,irhow,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhoE_x_rho  * rho  + &
                 rhoE_x_rhou * rhou + &
                 rhoE_x_rhov * rhov + &
                 rhoE_x_rhow * rhow + &
                 rhoE_x_rhoE * rhoE
        flux_y = rhoE_y_rho  * rho  + &
                 rhoE_y_rhou * rhou + &
                 rhoE_y_rhov * rhov + &
                 rhoE_y_rhow * rhow + &
                 rhoE_y_rhoE * rhoE
        flux_z = rhoE_z_rho  * rho  + &
                 rhoE_z_rhou * rhou + &
                 rhoE_z_rhov * rhov + &
                 rhoE_z_rhow * rhow + &
                 rhoE_z_rhoE * rhoE

        call integrate_volume_flux(mesh,sdata,elem_info,function_info,irhoE,flux_x,flux_y,flux_z)

    end subroutine compute
    !******************************************************************************************************






end module LINEULER_volume_advective_flux_imag
