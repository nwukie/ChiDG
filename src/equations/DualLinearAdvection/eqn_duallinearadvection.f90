module eqn_duallinearadvection
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



    !> This equation set exists really just to test equationsets with 
    !! more than one equation. The idea is just to compute the linear
    !! advecdtion solution twice at the same time. The equations are 
    !! independent of each other. So, we can verify, for example,
    !! the volume flux jacobians for each equation. They should be the
    !! same as for the single LinearAdvection equation set
    !!
    !!
    !-------------------------------------------------------------
    type, extends(equationset_t), public :: duallinearadvection_e

        real(rk)    :: c(3)     !> Advection velocities (cx, cy, cz)

    contains
        procedure   :: init
        procedure   :: compute_boundary_average_flux
        procedure   :: compute_boundary_upwind_flux
        procedure   :: compute_volume_flux
        procedure   :: compute_volume_source

    end type duallinearadvection_e

contains


    !===========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(duallinearadvection_e), intent(inout) :: self

        self%neqns   = 2

        ! Allocate equations
        allocate(self%eqns(self%neqns))

        ! Initialize equation parameters
        self%eqns(1)%name = "u_a"
        self%eqns(1)%ind  = 1

        self%eqns(2)%name = "u_b"
        self%eqns(2)%ind  = 2

        !call self%add('u_a',1)
        !call self%add('u_b',2)



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
        class(duallinearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk

        integer(ik)              :: iu_a, iu_b, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: ua_l(:), ua_r(:), ub_l(:), ub_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu_a      = self%get_var('u_a')
        iu_b      = self%get_var('u_b')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor



        !> Allocate arrays for data at quadrature points
        allocate(ua_l(nnodes),       &
                 ua_r(nnodes),       &
                 ub_l(nnodes),       &
                 ub_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !> Get neighbor face and seed element for derivatives
        iface_p = compute_neighbor_face(iface)
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_a,ua_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_a,ua_l,iseed)

        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_b,ub_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_b,ub_l,iseed)



        !> Compute boundary average flux for u_a
        flux_x = ((self%c(1)*ua_l + self%c(1)*ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = ((self%c(2)*ua_l + self%c(2)*ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = ((self%c(3)*ua_l + self%c(3)*ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_a, iblk, flux_x, flux_y, flux_z)



        !> Compute boundary average flux for u_b
        flux_x = ((self%c(1)*ub_l + self%c(1)*ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = ((self%c(2)*ub_l + self%c(2)*ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = ((self%c(3)*ub_l + self%c(3)*ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_b, iblk, flux_x, flux_y, flux_z)




    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,sdata,ielem,iface,iblk)
        class(duallinearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk



        integer(ik)              :: iu_a, iu_b, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: ua_l(:), ua_r(:), ub_l(:), ub_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu_a        = self%get_var('u_a')
        iu_b        = self%get_var('u_b')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor


        !> Allocate arrays for data at quadrature points
        allocate(ua_l(nnodes),       &
                 ua_r(nnodes),       &
                 ub_l(nnodes),       &
                 ub_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !> Get neighbor face and seed element for derivatives
        iface_p = compute_neighbor_face(iface)
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_a,ua_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_a,ua_l,iseed)

        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_b,ub_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_b,ub_l,iseed)





        !> Compute boundary upwind flux
        flux_x = (self%c(1) * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (self%c(2) * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (self%c(3) * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_a, iblk, flux_x, flux_y, flux_z)



        !> Compute boundary upwind flux
        flux_x = (self%c(1) * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (self%c(2) * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (self%c(3) * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_b, iblk, flux_x, flux_y, flux_z)






    end subroutine



    !===========================================================
    !
    !   Volume Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,sdata,ielem,iblk)
        class(duallinearadvection_e),       intent(in)      :: self
        class(mesh_t),                  intent(in)      :: mesh
        class(solverdata_t),            intent(inout)   :: sdata
        integer(ik),                    intent(in)      :: ielem, iblk



        type(AD_D), allocatable :: ua(:), ub(:), flux_x(:), flux_y(:), flux_z(:)
        integer(ik)             :: nnodes, ierr, iseed
        integer(ik)             :: iu_a, iu_b, i


        associate (elem => mesh%elems(ielem), q => sdata%q)


            ! Get variable index from equation set
            iu_a = self%get_var('u_a')
            iu_b = self%get_var('u_b')


            ! Allocate storage for variable values at quadrature points
            nnodes = elem%gq%nnodes_v
            allocate(ua(nnodes),         &
                     ub(nnodes),         &
                     flux_x(nnodes),    &
                     flux_y(nnodes),    &
                     flux_z(nnodes),    stat = ierr)
            if (ierr /= 0) call AllocationError


            !> Get seed element for derivatives
            iseed   = compute_seed_element(mesh,ielem,iblk)



            !> Interpolate solution to quadrature nodes
            call interpolate(mesh%elems,q,ielem,iu_a,ua,iseed)

            call interpolate(mesh%elems,q,ielem,iu_b,ub,iseed)

            

            !> Compute volume flux at quadrature nodes
            flux_x = self%c(1)  *  ua
            flux_y = self%c(2)  *  ua
            flux_z = self%c(3)  *  ua

            call integrate_volume_flux(elem,sdata,iu_a,iblk,flux_x,flux_y,flux_z)



            !> Compute volume flux at quadrature nodes
            flux_x = self%c(1)  *  ub
            flux_y = self%c(2)  *  ub
            flux_z = self%c(3)  *  ub

            call integrate_volume_flux(elem,sdata,iu_b,iblk,flux_x,flux_y,flux_z)



        end associate

    end subroutine


    !===========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,sdata,ielem,iblk)
        class(duallinearadvection_e),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk


    end subroutine









end module eqn_duallinearadvection
