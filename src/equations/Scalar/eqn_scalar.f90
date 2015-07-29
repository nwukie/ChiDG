module eqn_scalar
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use type_expansion,         only: expansion_t
    use DNAD_D,                 only: AD_D


    implicit none

    private

    type, extends(equationset_t), public :: scalar_e


    contains
        ! Must define these procedures in the extended, concrete type
        procedure   :: init
        procedure   :: compute_boundary_average_flux
        procedure   :: compute_boundary_upwind_flux
        procedure   :: compute_volume_flux
        procedure   :: compute_volume_source




    end type scalar_e


contains
    !==========================================================
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

    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,q,rhs,ielem,iface,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        class(expansion_t), intent(inout)   :: rhs(:)

        integer(ik),        intent(in)      :: ielem, iface, iblk

        type(AD_D), allocatable  :: u(:)










    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,q,rhs,ielem,iface,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        class(expansion_t), intent(inout)   :: rhs(:)
        integer(ik),        intent(in)      :: ielem, iface, iblk


    end subroutine



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,q,rhs,ielem,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        class(expansion_t), intent(inout)   :: rhs(:)
        integer(ik),        intent(in)      :: ielem, iblk


    end subroutine


    !==========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,q,rhs,ielem,iblk)
        class(scalar_e),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        class(expansion_t), intent(inout)   :: rhs(:)
        integer(ik),        intent(in)      :: ielem, iblk


    end subroutine



end module eqn_scalar
