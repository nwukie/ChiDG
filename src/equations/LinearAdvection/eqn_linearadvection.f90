module eqn_linearadvection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux, integrate_boundary_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    implicit none

    private

    type, extends(equationset_t), public :: linearadvection_e

        real(rk)    :: c(3)     !> Advection velocities (cx, cy, cz)

    contains
        procedure   :: init
        procedure   :: compute_boundary_average_flux
        procedure   :: compute_boundary_upwind_flux
        procedure   :: compute_volume_flux
        procedure   :: compute_volume_source

    end type linearadvection_e

contains


    !===========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(linearadvection_e), intent(inout) :: self

        self%neqns   = 1

        ! Allocate equations
        allocate(self%eqns(self%neqns))

        ! Initialize equation parameters
        self%eqns(1)%name = "u"
        self%eqns(1)%ind  = 1


        self%c(1) = ONE
        self%c(2) = ZERO
        self%c(3) = ZERO


    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,sdata,ielem,iface,iblk)
        class(linearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk

        integer(ik)              :: iu, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu       = self%get_var('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor



        !> Allocate arrays for data at quadrature points
        allocate(u_l(nnodes),       &
                 u_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !> Get neighbor face and seed element for derivatives
        iface_p = compute_neighbor_face(iface)
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu,u_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu,u_l,iseed)


        !> Compute boundary average flux
        flux_x = ((self%c(1)*u_l + self%c(1)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = ((self%c(2)*u_l + self%c(2)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = ((self%c(3)*u_l + self%c(3)*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)


        !> Integrate flux
        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,sdata,ielem,iface,iblk)
        class(linearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk



        integer(ik)              :: iu, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu        = self%get_var('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor


        !> Allocate arrays for data at quadrature points
        allocate(u_l(nnodes),       &
                 u_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !> Get neighbor face and seed element for derivatives
        iface_p = compute_neighbor_face(iface)
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu,u_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu,u_l,iseed)


        !> Compute boundary upwind flux
        flux_x = (self%c(1) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (self%c(2) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (self%c(3) * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)


        !> Integrate flux
        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine



    !===========================================================
    !
    !   Volume Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,sdata,ielem,iblk)
        class(linearadvection_e),       intent(in)      :: self
        class(mesh_t),                  intent(in)      :: mesh
        class(solverdata_t),            intent(inout)   :: sdata
        integer(ik),                    intent(in)      :: ielem, iblk



        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        integer(ik)             :: nnodes, ierr, iseed
        integer(ik)             :: ivar_u, i


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


            !> Get seed element for derivatives
            iseed   = compute_seed_element(mesh,ielem,iblk)


            !> Interpolate solution to quadrature nodes
            call interpolate(mesh%elems,q,ielem,ivar_u,u,iseed)


            !> Compute volume flux at quadrature nodes
            flux_x = self%c(1)  *  u
            flux_y = self%c(2)  *  u
            flux_z = self%c(3)  *  u


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
        class(linearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk


    end subroutine









end module eqn_linearadvection
