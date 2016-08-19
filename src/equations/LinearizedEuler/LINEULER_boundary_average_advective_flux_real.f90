module LINEULER_boundary_average_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none



    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LINEULER_boundary_average_advective_flux_real_t

    contains
        procedure  :: compute

    end type LINEULER_boundary_average_advective_flux_real_t
    !*******************************************************************************************










contains



    !>   Boundary Flux routine for Euler
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(LINEULER_boundary_average_advective_flux_real_t), intent(in)      :: self
        type(mesh_t),                                           intent(in)      :: mesh(:)
        type(solverdata_t),                                     intent(inout)   :: sdata
        class(properties_t),                                    intent(inout)   :: prop
        type(face_info_t),                                      intent(in)      :: face_info
        type(function_info_t),                                  intent(in)      :: function_info

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe


        integer(ik) :: idom, ielem, iface
        integer(ik) :: ifcn, idonor, iblk, igq

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                  &
                        rhou_m,     rhou_p,                                 &
                        rhov_m,     rhov_p,                                 &
                        rhow_m,     rhow_p,                                 &
                        rhoe_m,     rhoe_p,                                 &
                        flux_x_m,   flux_y_m,   flux_z_m,                   &
                        flux_x_p,   flux_y_p,   flux_z_p,                   &
                        flux_x,     flux_y,     flux_z,                     &
                        integrand


        irho  = prop%get_eqn_index("rho_r")
        irhou = prop%get_eqn_index("rhou_r")
        irhov = prop%get_eqn_index("rhov_r")
        irhow = prop%get_eqn_index("rhow_r")
        irhoE = prop%get_eqn_index("rhoE_r")

        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, faces => mesh(idom)%faces, q => sdata%q)


            !
            ! Interpolate solution to quadrature nodes
            !
            rho_m  = interpolate(mesh,sdata,face_info,function_info, irho,  'value', ME)
            rho_p  = interpolate(mesh,sdata,face_info,function_info, irho,  'value', NEIGHBOR)

            rhou_m = interpolate(mesh,sdata,face_info,function_info, irhou, 'value', ME)
            rhou_p = interpolate(mesh,sdata,face_info,function_info, irhou, 'value', NEIGHBOR)

            rhov_m = interpolate(mesh,sdata,face_info,function_info, irhov, 'value', ME)
            rhov_p = interpolate(mesh,sdata,face_info,function_info, irhov, 'value', NEIGHBOR)

            rhow_m = interpolate(mesh,sdata,face_info,function_info, irhow, 'value', ME)
            rhow_p = interpolate(mesh,sdata,face_info,function_info, irhow, 'value', NEIGHBOR)

            rhoE_m = interpolate(mesh,sdata,face_info,function_info, irhoE, 'value', ME)
            rhoE_p = interpolate(mesh,sdata,face_info,function_info, irhoE, 'value', NEIGHBOR)




            !================================
            !       MASS FLUX
            !================================
            flux_x_m = rho_x_rho  * rho_m  + &
                       rho_x_rhou * rhou_m + &
                       rho_x_rhov * rhov_m + &
                       rho_x_rhow * rhow_m + &
                       rho_x_rhoE * rhoE_m
            flux_y_m = rho_y_rho  * rho_m  + &
                       rho_y_rhou * rhou_m + &
                       rho_y_rhov * rhov_m + &
                       rho_y_rhow * rhow_m + &
                       rho_y_rhoE * rhoE_m 
            flux_z_m = rho_z_rho  * rho_m  + &
                       rho_z_rhou * rhou_m + &
                       rho_z_rhov * rhov_m + &
                       rho_z_rhow * rhow_m + &
                       rho_z_rhoE * rhoE_m 

            flux_x_p = rho_x_rho  * rho_p  + &
                       rho_x_rhou * rhou_p + &
                       rho_x_rhov * rhov_p + &
                       rho_x_rhow * rhow_p + &
                       rho_x_rhoE * rhoE_p
            flux_y_p = rho_y_rho  * rho_p  + &
                       rho_y_rhou * rhou_p + &
                       rho_y_rhov * rhov_p + &
                       rho_y_rhow * rhow_p + &
                       rho_y_rhoE * rhoE_p 
            flux_z_p = rho_z_rho  * rho_p  + &
                       rho_z_rhou * rhou_p + &
                       rho_z_rhov * rhov_p + &
                       rho_z_rhow * rhow_p + &
                       rho_z_rhoE * rhoE_p 


            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)




            !================================
            !       X-MOMENTUM FLUX
            !================================
            flux_x_m = rhou_x_rho  * rho_m  + &
                       rhou_x_rhou * rhou_m + &
                       rhou_x_rhov * rhov_m + &
                       rhou_x_rhow * rhow_m + &
                       rhou_x_rhoE * rhoE_m
            flux_y_m = rhou_y_rho  * rho_m  + &
                       rhou_y_rhou * rhou_m + &
                       rhou_y_rhov * rhov_m + &
                       rhou_y_rhow * rhow_m + &
                       rhou_y_rhoE * rhoE_m 
            flux_z_m = rhou_z_rho  * rho_m  + &
                       rhou_z_rhou * rhou_m + &
                       rhou_z_rhov * rhov_m + &
                       rhou_z_rhow * rhow_m + &
                       rhou_z_rhoE * rhoE_m 

            flux_x_p = rhou_x_rho  * rho_p  + &
                       rhou_x_rhou * rhou_p + &
                       rhou_x_rhov * rhov_p + &
                       rhou_x_rhow * rhow_p + &
                       rhou_x_rhoE * rhoE_p
            flux_y_p = rhou_y_rho  * rho_p  + &
                       rhou_y_rhou * rhou_p + &
                       rhou_y_rhov * rhov_p + &
                       rhou_y_rhow * rhow_p + &
                       rhou_y_rhoE * rhoE_p 
            flux_z_p = rhou_z_rho  * rho_p  + &
                       rhou_z_rhou * rhou_p + &
                       rhou_z_rhov * rhov_p + &
                       rhou_z_rhow * rhow_p + &
                       rhou_z_rhoE * rhoE_p 


            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)



            !================================
            !       Y-MOMENTUM FLUX
            !================================
            flux_x_m = rhov_x_rho  * rho_m  + &
                       rhov_x_rhou * rhou_m + &
                       rhov_x_rhov * rhov_m + &
                       rhov_x_rhow * rhow_m + &
                       rhov_x_rhoE * rhoE_m
            flux_y_m = rhov_y_rho  * rho_m  + &
                       rhov_y_rhou * rhou_m + &
                       rhov_y_rhov * rhov_m + &
                       rhov_y_rhow * rhow_m + &
                       rhov_y_rhoE * rhoE_m 
            flux_z_m = rhov_z_rho  * rho_m  + &
                       rhov_z_rhou * rhou_m + &
                       rhov_z_rhov * rhov_m + &
                       rhov_z_rhow * rhow_m + &
                       rhov_z_rhoE * rhoE_m 


            flux_x_p = rhov_x_rho  * rho_p  + &
                       rhov_x_rhou * rhou_p + &
                       rhov_x_rhov * rhov_p + &
                       rhov_x_rhow * rhow_p + &
                       rhov_x_rhoE * rhoE_p
            flux_y_p = rhov_y_rho  * rho_p  + &
                       rhov_y_rhou * rhou_p + &
                       rhov_y_rhov * rhov_p + &
                       rhov_y_rhow * rhow_p + &
                       rhov_y_rhoE * rhoE_p 
            flux_z_p = rhov_z_rho  * rho_p  + &
                       rhov_z_rhou * rhou_p + &
                       rhov_z_rhov * rhov_p + &
                       rhov_z_rhow * rhow_p + &
                       rhov_z_rhoE * rhoE_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)




            !================================
            !       Z-MOMENTUM FLUX
            !================================
            flux_x_m = rhow_x_rho  * rho_m  + &
                       rhow_x_rhou * rhou_m + &
                       rhow_x_rhov * rhov_m + &
                       rhow_x_rhow * rhow_m + &
                       rhow_x_rhoE * rhoE_m
            flux_y_m = rhow_y_rho  * rho_m  + &
                       rhow_y_rhou * rhou_m + &
                       rhow_y_rhov * rhov_m + &
                       rhow_y_rhow * rhow_m + &
                       rhow_y_rhoE * rhoE_m 
            flux_z_m = rhow_z_rho  * rho_m  + &
                       rhow_z_rhou * rhou_m + &
                       rhow_z_rhov * rhov_m + &
                       rhow_z_rhow * rhow_m + &
                       rhow_z_rhoE * rhoE_m 


            flux_x_p = rhow_x_rho  * rho_p  + &
                       rhow_x_rhou * rhou_p + &
                       rhow_x_rhov * rhov_p + &
                       rhow_x_rhow * rhow_p + &
                       rhow_x_rhoE * rhoE_p
            flux_y_p = rhow_y_rho  * rho_p  + &
                       rhow_y_rhou * rhou_p + &
                       rhow_y_rhov * rhov_p + &
                       rhow_y_rhow * rhow_p + &
                       rhow_y_rhoE * rhoE_p 
            flux_z_p = rhow_z_rho  * rho_p  + &
                       rhow_z_rhou * rhou_p + &
                       rhow_z_rhov * rhov_p + &
                       rhow_z_rhow * rhow_p + &
                       rhow_z_rhoE * rhoE_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)


            !================================
            !          ENERGY FLUX
            !================================
            flux_x_m = rhoE_x_rho  * rho_m  + &
                       rhoE_x_rhou * rhou_m + &
                       rhoE_x_rhov * rhov_m + &
                       rhoE_x_rhow * rhow_m + &
                       rhoE_x_rhoE * rhoE_m
            flux_y_m = rhoE_y_rho  * rho_m  + &
                       rhoE_y_rhou * rhou_m + &
                       rhoE_y_rhov * rhov_m + &
                       rhoE_y_rhow * rhow_m + &
                       rhoE_y_rhoE * rhoE_m 
            flux_z_m = rhoE_z_rho  * rho_m  + &
                       rhoE_z_rhou * rhou_m + &
                       rhoE_z_rhov * rhov_m + &
                       rhoE_z_rhow * rhow_m + &
                       rhoE_z_rhoE * rhoE_m 

            flux_x_p = rhoE_x_rho  * rho_p  + &
                       rhoE_x_rhou * rhou_p + &
                       rhoE_x_rhov * rhov_p + &
                       rhoE_x_rhow * rhow_p + &
                       rhoE_x_rhoE * rhoE_p
            flux_y_p = rhoE_y_rho  * rho_p  + &
                       rhoE_y_rhou * rhou_p + &
                       rhoE_y_rhov * rhov_p + &
                       rhoE_y_rhow * rhow_p + &
                       rhoE_y_rhoE * rhoE_p 
            flux_z_p = rhoE_z_rho  * rho_p  + &
                       rhoE_z_rhou * rhou_p + &
                       rhoE_z_rhov * rhov_p + &
                       rhoE_z_rhow * rhow_p + &
                       rhoE_z_rhoE * rhoE_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)




        end associate

    end subroutine compute
    !********************************************************************************************************












end module LINEULER_boundary_average_advective_flux_real
