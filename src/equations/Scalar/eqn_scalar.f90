module eqn_scalar
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux, integrate_boundary_flux
    use DNAD_D


    implicit none

    private

    type, extends(equationset_t), public :: scalar_e

        real(rk)    :: c(3)     !> Advection velocities (cx, cy, cz)

    contains
        ! Must define these procedures in the extended, concrete type
        procedure   :: init
        procedure   :: compute_boundary_average_flux
        procedure   :: compute_boundary_upwind_flux
        procedure   :: compute_volume_flux
        procedure   :: compute_volume_source




    end type scalar_e


contains
    !===========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(scalar_e), intent(inout) :: self

        self%neqns   = 1

        ! Allocate equations
        allocate(self%eqns(self%neqns))

        ! Initialize equation parameters
        self%eqns(1)%name = "u"
        self%eqns(1)%ind  = 1


        self%c(1) = ZERO
        self%c(2) = ZERO
        self%c(3) = ONE


    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,sdata,ielem,iface,iblk)
        class(scalar_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk

        integer(ik)              :: u_i, iseed, ierr, nnodes, ineighbor, iface_p
        type(AD_D), allocatable  :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        u_i       = self%get_var('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor



        !> Allocate arrays for data at quadrature points
        allocate(u_l(nnodes),       &
                 u_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !> Get element for seeding derivatives
        if ( iblk == DIAG ) then
            iseed = ielem
        else
            iseed = mesh%faces(ielem,iblk)%ineighbor
        end if


        if (iface == XI_MIN) then
            iface_p = XI_MAX
        else if (iface == XI_MAX) then
            iface_p = XI_MIN
        else if (iface == ETA_MIN) then
            iface_p = ETA_MAX
        else if (iface == ETA_MAX) then
            iface_p = ETA_MIN
        else if (iface == ZETA_MIN) then
            iface_p = ZETA_MAX
        else if (iface == ZETA_MAX) then
            iface_p = ZETA_MIN
        end if



        !> Interpolate to quadrature points
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  u_i,u_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,u_i,u_l,iseed)



        !> Compute boundary average flux
        flux_x = self%c(1)*u_l
!        flux_x = ((self%c(1)*u_l + self%c(1)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
!        flux_y = ((self%c(2)*u_l + self%c(2)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
!        flux_z = ((self%c(3)*u_l + self%c(3)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)
!
!
!        !> Integrate flux
!        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, u_i, iblk, flux_x, flux_y, flux_z)

    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,sdata,ielem,iface,iblk)
        class(scalar_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk



        integer(ik)              :: u_i, iseed, ierr, nnodes, ineighbor, iface_p
        type(AD_D), allocatable  :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        u_i       = self%get_var('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor


        !> Allocate arrays for data at quadrature points
        allocate(u_l(nnodes),       &
                 u_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        if (iface == XI_MIN) then
            iface_p = XI_MAX
        else if (iface == XI_MAX) then
            iface_p = XI_MIN
        else if (iface == ETA_MIN) then
            iface_p = ETA_MAX
        else if (iface == ETA_MAX) then
            iface_p = ETA_MIN
        else if (iface == ZETA_MIN) then
            iface_p = ZETA_MAX
        else if (iface == ZETA_MAX) then
            iface_p = ZETA_MIN
        end if





        !> Get element for seeding derivatives
        if ( iblk == DIAG ) then
            iseed = ielem
        else
            iseed = mesh%faces(ielem,iblk)%ineighbor
        end if

        !> Interpolate to quadrature points
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  u_i,u_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,u_i,u_l,iseed)


        !> Compute boundary average flux
        flux_x = (self%c(1) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (self%c(2) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (self%c(3) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)


!        if (ielem == 22) then
!            print*, 'u_modes'
!            print*, sdata%q(ielem)%mat
!
!            print*, 'flux_x'
!            print*, flux_x(:)%x_ad_
!
!            print*, 'flux_y'
!            print*, flux_y(:)%x_ad_
!
!            print*, 'flux_z'
!            print*, flux_z(:)%x_ad_
!        end if



        !> Integrate flux
        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, u_i, iblk, flux_x, flux_y, flux_z)





    end subroutine



    !===========================================================
    !
    !   Volume Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,sdata,ielem,iblk)
        class(scalar_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk

        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        integer(ik)             :: nnodes, ierr, iseed
        integer(ik)             :: ivar_u


        associate (elem => mesh%elems(ielem), q => sdata%q)


            ! Get variable index from equation set
            ivar_u = self%get_var('u')


            ! Allocate storage for variable values at quadrature points
            nnodes = elem%gq%nnodes_v
            allocate(u(nnodes),         &
                     flux_x(nnodes),    &
                     flux_y(nnodes),    &
                     flux_z(nnodes),    stat = ierr)
            if (ierr /= 0) call AllocationError




            ! Interpolate modal varaiables to quadrature points
            !> Get element for seeding derivatives
            if ( iblk == DIAG ) then
                iseed = ielem
            else
                iseed = mesh%faces(ielem,iblk)%ineighbor
            end if

!            iseed = 0   !> no derivative tracking
            call interpolate(mesh%elems,q,ielem,ivar_u,u,iseed)


            ! Compute volume flux at quadrature nodes
            flux_x = self%c(1)  *  u
            flux_y = self%c(2)  *  u
            flux_z = self%c(3)  *  u




!            if (ielem == 22) then
!                print*, 'u_modes'
!                print*, sdata%q(ielem)%mat
!
!                print*, 'u'
!                print*, u(:)%x_ad_
!
!                print*, 'flux_x'
!                print*, flux_x(:)%x_ad_
!
!                print*, 'flux_y'
!                print*, flux_y(:)%x_ad_
!
!                print*, 'flux_z'
!                print*, flux_z(:)%x_ad_
!            end if



            ! Integrate volume flux
            call integrate_volume_flux(elem,sdata,ivar_u,iblk,flux_x,flux_y,flux_z)

        end associate

    end subroutine


    !===========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,sdata,ielem,iblk)
        class(scalar_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk


    end subroutine



end module eqn_scalar
