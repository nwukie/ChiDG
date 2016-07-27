module EULER_volume_advective_flux_cacheoptimized
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_timer,             only: timer_t
    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools
    use DNAD_D

    use EULER_properties,       only: EULER_properties_t
    implicit none


    type(timer_t), public   :: EULER_volume_total
    type(timer_t), public   :: EULER_volume_interpolate

    private


    
    !> Volume advective flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: EULER_volume_advective_flux_cacheoptimized_t


    contains

        procedure  :: compute

    end type EULER_volume_advective_flux_cacheoptimized_t
    !******************************************************************************










contains



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(EULER_volume_advective_flux_cacheoptimized_t),   intent(in)      :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        integer(ik),                            intent(in)      :: idom, ielem, iblk

        ! Equation indices
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe

        integer(ik)    :: ngq, block, iloop, gq_max, gq_min
        integer(ik)    :: iseed, i, idonor
        type(seed_t)   :: seed

        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                        &
                    flux_x, flux_y, flux_z, invrho


        idonor = 0

        call EULER_volume_total%start()

        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")


        !
        ! Get neighbor face and seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iblk,idonor,iblk)


        !
        ! Interpolate solution to quadrature nodes
        !
        call EULER_volume_interpolate%start()
        call interpolate_element(mesh,sdata%q,idom,ielem,irho, rho, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou,rhou,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov,rhov,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow,rhow,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE,rhoE,seed)
        call EULER_volume_interpolate%stop()


        ngq = size(rho)
        block = 10


        iloop  = 1
        gq_max = 0
        do while (gq_max < ngq)
            gq_min = 1  +  block*(iloop-1)
            gq_max = min(ngq,block*(iloop))



        invrho(gq_min:gq_max) = ONE/rho(gq_min:gq_max)



        !
        ! Compute pressure and total enthalpy
        !
        call prop%fluid%compute_pressure(rho(gq_min:gq_max),rhou(gq_min:gq_max),rhov(gq_min:gq_max),rhow(gq_min:gq_max),rhoE(gq_min:gq_max),p(gq_min:gq_max))

        H(gq_min:gq_max) = (rhoE(gq_min:gq_max) + p(gq_min:gq_max))*invrho(gq_min:gq_max)

        !===========================
        !        MASS FLUX
        !===========================
        flux_x(gq_min:gq_max) = rhou(gq_min:gq_max)
        flux_y(gq_min:gq_max) = rhov(gq_min:gq_max)
        flux_z(gq_min:gq_max) = rhow(gq_min:gq_max)


            iloop = iloop + 1
        end do


        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        iloop  = 1
        gq_max = 0
        do while (gq_max < ngq)
            gq_min = 1  +  block*(iloop-1)
            gq_max = min(ngq,block*(iloop))
        flux_x(gq_min:gq_max) = (rhou(gq_min:gq_max)*rhou(gq_min:gq_max))*invrho(gq_min:gq_max)  +  p(gq_min:gq_max)
        flux_y(gq_min:gq_max) = (rhou(gq_min:gq_max)*rhov(gq_min:gq_max))*invrho(gq_min:gq_max)
        flux_z(gq_min:gq_max) = (rhou(gq_min:gq_max)*rhow(gq_min:gq_max))*invrho(gq_min:gq_max)

            iloop = iloop + 1
        end do

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        iloop  = 1
        gq_max = 0
        do while (gq_max < ngq)
            gq_min = 1  +  block*(iloop-1)
            gq_max = min(ngq,block*(iloop))
        flux_x(gq_min:gq_max) = (rhov(gq_min:gq_max)*rhou(gq_min:gq_max))*invrho(gq_min:gq_max)
        flux_y(gq_min:gq_max) = (rhov(gq_min:gq_max)*rhov(gq_min:gq_max))*invrho(gq_min:gq_max)  +  p(gq_min:gq_max)
        flux_z(gq_min:gq_max) = (rhov(gq_min:gq_max)*rhow(gq_min:gq_max))*invrho(gq_min:gq_max)

            iloop = iloop + 1
        end do

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhov,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        iloop  = 1
        gq_max = 0
        do while (gq_max < ngq)
            gq_min = 1  +  block*(iloop-1)
            gq_max = min(ngq,block*(iloop))
        flux_x(gq_min:gq_max) = (rhow(gq_min:gq_max)*rhou(gq_min:gq_max))*invrho(gq_min:gq_max)
        flux_y(gq_min:gq_max) = (rhow(gq_min:gq_max)*rhov(gq_min:gq_max))*invrho(gq_min:gq_max)
        flux_z(gq_min:gq_max) = (rhow(gq_min:gq_max)*rhow(gq_min:gq_max))*invrho(gq_min:gq_max)  +  p(gq_min:gq_max)

            iloop = iloop + 1
        end do

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhow,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        iloop  = 1
        gq_max = 0
        do while (gq_max < ngq)
            gq_min = 1  +  block*(iloop-1)
            gq_max = min(ngq,block*(iloop))
        flux_x(gq_min:gq_max) = rhou(gq_min:gq_max)*H(gq_min:gq_max)
        flux_y(gq_min:gq_max) = rhov(gq_min:gq_max)*H(gq_min:gq_max)
        flux_z(gq_min:gq_max) = rhow(gq_min:gq_max)*H(gq_min:gq_max)

            iloop = iloop + 1
        end do

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhoE,iblk,flux_x,flux_y,flux_z)




        call EULER_volume_total%stop()

    end subroutine compute
    !*********************************************************************************************************






end module EULER_volume_advective_flux_cacheoptimized
