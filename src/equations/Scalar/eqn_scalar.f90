module eqn_scalar
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solver,           only: solver_t
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


        self%c(1) = ONE
        self%c(2) = ZERO
        self%c(3) = ZERO


    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,solver,ielem,iface,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(solver_t),    intent(inout)   :: solver
        integer(ik),        intent(in)      :: ielem, iface, iblk

        type(AD_D), allocatable  :: u(:)






    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,solver,ielem,iface,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(solver_t),    intent(inout)   :: solver
        integer(ik),        intent(in)      :: ielem, iface, iblk


    end subroutine



    !==========================================================
    !
    !   Volume Flux routine for Scalar
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,solver,ielem,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(solver_t),    intent(inout)   :: solver
        integer(ik),        intent(in)      :: ielem, iblk

        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        integer(ik)             :: nnodes, ierr, iseed
        integer(ik)             :: ivar_u


        associate (elem => mesh%elems(ielem), q => solver%q)


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
            iseed = 0   !> no derivative tracking
            call interpolate(mesh%elems,q,ielem,ivar_u,u,iseed)


            ! Compute volume flux at quadrature nodes
            flux_x = self%c(1)  *  u
            flux_y = self%c(2)  *  u
            flux_z = self%c(3)  *  u


            ! Integrate volume flux
            call integrate_volume_flux(elem,solver,ivar_u,iblk,flux_x,flux_y,flux_z)

        end associate

    end subroutine


    !==========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,solver,ielem,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(solver_t),    intent(inout)   :: solver
        integer(ik),        intent(in)      :: ielem, iblk


    end subroutine



end module eqn_scalar
