module type_eqnset_euler
    use mod_types,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use mod_compute_lift,       only: compute_lift_interior_face

    use atype_equationset,      only: equationset_t
    use type_equation,          only: equation_t
    use type_mesh,              only: mesh_t
    use type_element,           only: element_t
    use type_face,              only: face_t
    use type_mask,              only: mask_t


    implicit none

    private

    type, extends(equationset_t), public :: eqnset_euler_t

        ! Equation-set specific data
        type(fluid_t)   :: fluid

    contains
        ! Must define these procedures in the extended, concrete type
        procedure  :: init
        procedure  :: compute_boundary_flux_face
        procedure  :: compute_volume_flux
        procedure  :: compute_volume_source




    end type eqnset_euler_t


contains
    !==========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(eqnset_euler_t), intent(inout) :: self

        self%neqns   = 5

        ! Allocate equations
        allocate(self%eqns(self%neqns))

        ! Initialize equation parameters
        self%eqns(1)%name = "rho"
        self%eqns(1)%ind  = 1

        self%eqns(2)%name = "rhou"
        self%eqns(2)%ind  = 2

        self%eqns(3)%name = "rhov"
        self%eqns(3)%ind  = 3

        self%eqns(4)%name = "rhow"
        self%eqns(4)%ind  = 4

        self%eqns(5)%name = "rhoE"
        self%eqns(5)%ind  = 5

        ! Initialize equation set parameters
        self%rgas = 287.058_rk        ! J/(kg*K)

    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,face_m)
        class(eqnset_euler_t), intent(in)   :: self
        class(face_t),         intent(in)   :: face_m

        ! Equation indices
        integer(ik)    :: rho_ind
        integer(ik)    :: rhou_ind
        integer(ik)    :: rhov_ind
        integer(ik)    :: rhow_ind
        integer(ik)    :: rhoe_ind

        integer(ik) :: ixi_m,ieta_m,izeta_m
        integer(ik) :: ixi_p,ieta_p,izeta_p,iface_p
        integer(ik) :: inode

        ! Storage at quadrature nodes
        real(rk), dimension(face_mgq%face%nnodes)    ::     &
                        nxi_m,     neta_m,     nzeta_m,     &
                        nxi_p,     neta_p,     nzeta_p

        type(AD_D), dimension(face_m%gq%face%nnodes)    ::  &
                        rho_m,      rho_p,                  &
                        rhou_m,     rhou_p,                 &
                        rhov_m,     rhov_p,                 &
                        rhow_m,     rhow_p,                 &
                        rhoe_m,     rhoe_p,                 &
                        p_m,        p_p,                    &
                        H_m,        H_p,                    &
                        flux_x_m,   flux_y_m,   flux_z_m,   &
                        flux_x_p,   flux_y_p,   flux_z_p,   &
                        flux_x,     flux_y,     flux_z,     &
                        flux


        ! Create pointers to faces to reduce typing
        type(face_t),   pointer :: p_face_m,p_face_p
        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        rho_ind  = self%get_var("rho")
        rhou_ind = self%get_var("rhou")
        rhov_ind = self%get_var("rhov")
        rhow_ind = self%get_var("rhow")
        rhoE_ind = self%get_var("rhoE")



        ! Skip if face is a boundary face
        if (mesh%faces(ixi_m,ieta_m,izeta_m,iface_m)%ftype /= 0) cycle

        ! Get neighbor index, based on current element and face index
        call get_neighbor(iface_m,iface_p,ixi_m,ieta_m,izeta_m,ixi_p,ieta_p,izeta_p,mesh)


        ! Associate pointers to m,p faces to reduce typing
        face_m => mesh%faces(ixi_m,ieta_m,izeta_m,iface_m)
        face_p => mesh%faces(ixi_p,ieta_p,izeta_p,iface_p)


        ! Get normal vectors
        nxi_m   = face_m%norm(:,1)
        neta_m  = face_m%norm(:,2)
        nzeta_m = face_m%norm(:,3)

        nxi_p   = face_p%norm(:,1)
        neta_p  = face_p%norm(:,2)
        nzeta_p = face_p%norm(:,3)


        ! Compute GQ nodes
        call face_m%compute_var(rho_ind,rho_m)
        call face_p%compute_var(rho_ind,rho_p)

        call face_m%compute_var(rhou_ind,rhou_m)
        call face_p%compute_var(rhou_ind,rhou_p)

        call face_m%compute_var(rhov_ind,rhov_m)
        call face_p%compute_var(rhov_ind,rhov_p)

        call face_m%compute_var(rhow_ind,rhow_m)
        call face_p%compute_var(rhow_ind,rhow_p)

        call face_m%compute_var(rhoE_ind,rhoE_m)
        call face_p%compute_var(rhoE_ind,rhoE_p)


        ! Compute pressure and total enthalpy
        call self%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
        call self%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)

        H_m = (rhoE_m + p_m)/rho_m
        H_p = (rhoE_p + p_p)/rho_p

        !================================
        !       MASS FLUX
        !================================
        flux_x_m = rhou_m
        flux_y_m = rhov_m
        flux_z_m = rhow_m

        flux_x_p = rhou_p
        flux_y_p = rhov_p
        flux_z_p = rhow_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        ! Dot product with normal vector
        flux = HALF*(nxi_m*flux_x + neta_m*flux_y + nzeta_m*flux_z)

        call face_m%integrate_flux(flux)

        !================================
        !       X-MOMENTUM FLUX
        !================================
        flux_x_m = (rhou_m*rhou_m)/rho_m + p_m
        flux_y_m = (rhou_m*rhov_m)/rho_m
        flux_z_m = (rhou_m*rhow_m)/rho_m

        flux_x_p = (rhou_p*rhou_p)/rho_p + p_p
        flux_y_p = (rhou_p*rhov_p)/rho_p
        flux_z_p = (rhou_p*rhow_p)/rho_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        ! Dot product with normal vector
        flux = HALF*(nxi_m*flux_x + neta_m*flux_y + nzeta_m*flux_z)


        call face_m%integrate_flux(flux)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        flux_x_m = (rhov_m*rhou_m)/rho_m
        flux_y_m = (rhov_m*rhov_m)/rho_m + p_m
        flux_z_m = (rhov_m*rhow_m)/rho_m

        flux_x_p = (rhov_p*rhou_p)/rho_p
        flux_y_p = (rhov_p*rhov_p)/rho_p + p_p
        flux_z_p = (rhov_p*rhow_p)/rho_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        ! Dot product with normal vector
        flux = HALF*(nxi_m*flux_x + neta_m*flux_y + nzeta_m*flux_z)

        call face_m%integrate_flux(flux)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        flux_x_m = (rhow_m*rhou_m)/rho_m
        flux_y_m = (rhow_m*rhov_m)/rho_m
        flux_z_m = (rhow_m*rhow_m)/rho_m + p_m

        flux_x_p = (rhow_p*rhou_p)/rho_p
        flux_y_p = (rhow_p*rhov_p)/rho_p
        flux_z_p = (rhow_p*rhow_p)/rho_p + p_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        ! Dot product with normal vector
        flux = HALF*(nxi_m*flux_x + neta_m*flux_y + nzeta_m*flux_z)

        call face_m%integrate_flux(flux)

        !================================
        !          ENERGY FLUX
        !================================
        flux_x_m = rhou_m*H_m
        flux_y_m = rhov_m*H_m
        flux_z_m = rhow_m*H_m

        flux_x_p = rhou_p*H_p
        flux_y_p = rhov_p*H_p
        flux_z_p = rhow_p*H_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        ! Dot product with normal vector
        flux = HALF*(nxi_m*flux_x + neta_m*flux_y + nzeta_m*flux_z)

        call face_m%integrate_flux(flux)


    end subroutine



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute_volume_flux(self,elem)
        class(eqnset_euler_t),  intent(in)      :: self
        class(element_t),       intent(inout)   :: elem

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: rho_ind
        integer(ik)    :: rhou_ind
        integer(ik)    :: rhov_ind
        integer(ik)    :: rhow_ind
        integer(ik)    :: rhoe_ind


        integer(ik)    :: ixi,ieta,izeta,iterm,inode

        type(AD_D), dimension(elem%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,     &
                    flux_x, flux_y, flux_z

        !-------------------------------------------------------------
        rho_ind  = self%get_var("rho")
        rhou_ind = self%get_var("rhou")
        rhov_ind = self%get_var("rhov")
        rhow_ind = self%get_var("rhow")
        rhoE_ind = self%get_var("rhoE")

        ! Compute variables at volume GQ nodes
        call elem%compute_var(rho_ind,rho)
        call elem%compute_var(rhou_ind,rhou)
        call elem%compute_var(rhov_ind,rhov)
        call elem%compute_var(rhow_ind,rhow)
        call elem%compute_var(rhoE_ind,rhoE)

        ! Compute pressure and total enthalpy
        call self%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)/rho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow

        call elem%integrate_volume_flux(flux_x,  flux_y,  flux_z,  rho_ind)

        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)/rho  +  p
        flux_y = (rhou*rhov)/rho
        flux_z = (rhou*rhow)/rho

        call elem%integrate_volume_flux(flux_x,  flux_y,  flux_z,  rho_ind)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)/rho
        flux_y = (rhov*rhov)/rho  +  p
        flux_z = (rhov*rhow)/rho

        call elem%integrate_volume_flux(flux_x,  flux_y,  flux_z,  rho_ind)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)/rho
        flux_y = (rhow*rhov)/rho
        flux_z = (rhow*rhow)/rho  +  p

        call elem%integrate_volume_flux(flux_x,  flux_y,  flux_z,  rho_ind)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call elem%integrate_volume_flux(flux_x,  flux_y,  flux_z,  rho_ind)

    end subroutine


    !==========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,elem)
        class(eqnset_euler_t), intent(in)       :: self
        class(elements_t),     intent(inout)    :: elem

        ! Equation indices
        !---------------------------------------------------------
        integer(ik)    :: rho_ind
        integer(ik)    :: rhou_ind
        integer(ik)    :: rhov_ind
        integer(ik)    :: rhow_ind
        integer(ik)    :: rhoe_ind

        integer(ik)    :: ixi,ieta,izeta

        type(AD_D), dimension(elem%gq%vol%nnodes)  :: &
                            source, u,                &
                            x,      y,      z,        &
                            xval,   yval,   zval

        !----------------------------------------------------------
        rho_ind  = self%get_var("rho")
        rhou_ind = self%get_var("rhou")
        rhov_ind = self%get_var("rhov")
        rhow_ind = self%get_var("rhow")
        rhoE_ind = self%get_var("rhoE")


        source = 0._rk


        ! Integrate volume source function
        call elem%integrate_volume_source(source,rho_ind)
        call elem%integrate_volume_source(source,rhou_ind)
        call elem%integrate_volume_source(source,rhov_ind)
        call elem%integrate_volume_source(source,rhow_ind)
        call elem%integrate_volume_source(source,rhoE_ind)

    end subroutine



    !================================================================
    !
    !   Implements the equation of state for computing pressure
    !
    !================================================================
    subroutine compute_pressure(self,rho,rhou,rhov,rhow,rhoE,p)
        class(eqnset_euler_t),  intent(in)            :: self
        real(kind=rk),  intent(in)     :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(kind=rk),  intent(inout)  :: p(:)

        real(kind=rk), dimension(size(p))   :: gam


        call self%compute_gamma(rho,rhou,rhov,rhow,rhoE,gam)

        p = (gam-ONE)*(rhoE - HALF*rho*((rhou/rho)*(rhou/rho) + (rhov/rho)*(rhov/rho) + (rhow/rho)*(rhow/rho)))

    end subroutine


    !================================================================
    !
    !   Implements a routine for computing temperature
    !
    !================================================================
    subroutine compute_temperature(self,rho,rhou,rhov,rhow,rhoE,T)
        class(eqnset_euler_t),  intent(in)            :: self
        real(kind=rk),  intent(in)     :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(kind=rk),  intent(inout)  :: T(:)

        real(kind=rk), dimension(size(T))   :: gam,p


        call self%compute_gamma(rho,rhou,rhov,rhow,rhoE,gam)
        call self%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        T = p/(rho*self%rgas)

    end subroutine



    !================================================================
    !
    !   Implements the computing gamma for a fluid
    !
    !================================================================
    subroutine compute_gamma(self,rho,rhou,rhov,rhow,rhoE,gam)
        class(eqnset_euler_t),  intent(in)            :: self
        real(kind=rk),  intent(in)     :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(kind=rk),  intent(inout)  :: gam(:)

        gam = 1.4_rk

    end subroutine



    !================================================================
    !
    !   Implements the computing cp (specific heat at constant pressure) for a fluid
    !
    !================================================================
    subroutine compute_cp(self,rho,rhou,rhov,rhow,rhoE,cp)
        class(eqnset_euler_t),  intent(in)  :: self
        real(kind=rk),  intent(in)          :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(kind=rk),  intent(inout)       :: cp(:)

        cp = 1003.5_rk  ! J/(kg*K)

    end subroutine












    subroutine get_neighbor(iface_m,iface_p,ixi_m,ieta_m,izeta_m,ixi_p,ieta_p,izeta_p,mesh)
        integer(kind=ik),   intent(in)      :: iface_m,ixi_m,ieta_m,izeta_m
        integer(kind=ik),   intent(inout)   :: iface_p,ixi_p,ieta_p,izeta_p
        type(mesh_t),       intent(in)      :: mesh

        if (ixi_m == 1 .and. iface_m == XI_MIN) then
            ixi_p = mesh%nelem_xi
            ieta_p = ieta_m
            izeta_p = izeta_m
            iface_p = XI_MAX
        elseif (ixi_m == mesh%nelem_xi .and. iface_m == XI_MAX) then
            ixi_p = 1
            ieta_p = ieta_m
            izeta_p = izeta_m
            iface_p = XI_MIN
        elseif (ieta_m == 1 .and. iface_m == ETA_MIN) then
            ixi_p = ixi_m
            ieta_p = mesh%nelem_eta
            izeta_p = izeta_p
            iface_p = ETA_MAX
        elseif (ieta_m == mesh%nelem_eta .and. iface_m == ETA_MAX) then
            ixi_p = ixi_m
            ieta_p = 1
            izeta_p = izeta_m
            iface_p = ETA_MIN
        elseif (izeta_m == 1 .and. iface_m == ZETA_MIN) then
            ixi_p = ixi_m
            ieta_p = ieta_m
            izeta_p = mesh%nelem_zeta
            iface_p = ZETA_MAX
        elseif (izeta_m == mesh%nelem_zeta .and. iface_m == ZETA_MAX) then
            ixi_p = ixi_m
            ieta_p = ieta_m
            izeta_p = 1
            iface_p = ZETA_MIN
        else

            ! Normal Interior Faces
            select case (iface_m)
                case (XI_MIN)
                    ixi_p   = ixi_m - 1
                    ieta_p  = ieta_m
                    izeta_p = izeta_m
                    iface_p = XI_MAX
                case (XI_MAX)
                    ixi_p   = ixi_m + 1
                    ieta_p  = ieta_m
                    izeta_p = izeta_m
                    iface_p = XI_MIN
                case (ETA_MIN)
                    ixi_p   = ixi_m
                    ieta_p  = ieta_m - 1
                    izeta_p = izeta_m
                    iface_p = ETA_MAX
                case (ETA_MAX)
                    ixi_p   = ixi_m
                    ieta_p  = ieta_m + 1
                    izeta_p = izeta_m
                    iface_p = ETA_MIN
                case (ZETA_MIN)
                    ixi_p   = ixi_m
                    ieta_p  = ieta_m
                    izeta_p = izeta_m - 1
                    iface_p = ZETA_MAX
                case (ZETA_MAX)
                    ixi_p   = ixi_m
                    ieta_p  = ieta_m
                    izeta_p = izeta_m + 1
                    iface_p = ZETA_MIN
                case default
                    stop "Error: get_neighbor - invalid face index"
            end select

        end if

    end subroutine
end module type_eqnset_euler
