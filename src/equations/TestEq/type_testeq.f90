module type_testeq
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use type_expansion,         only: expansion_t


    implicit none

    private

    type, extends(equationset_t), public :: testeq_t

        ! Equation-set specific data
!        type(fluid_t)   :: fluid

    contains
        ! Must define these procedures in the extended, concrete type
        procedure   :: init
        procedure   :: compute_boundary_average_flux
        procedure   :: compute_boundary_upwind_flux
        procedure   :: compute_volume_flux
        procedure   :: compute_volume_source




    end type testeq_t


contains
    !==========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(testeq_t), intent(inout) :: self
!
!        self%neqns   = 5
!
!        ! Allocate equations
!        allocate(self%eqns(self%neqns))
!
!        ! Initialize equation parameters
!        self%eqns(1)%name = "rho"
!        self%eqns(1)%ind  = 1
!
!        self%eqns(2)%name = "rhou"
!        self%eqns(2)%ind  = 2
!
!        self%eqns(3)%name = "rhov"
!        self%eqns(3)%ind  = 3
!
!        self%eqns(4)%name = "rhow"
!        self%eqns(4)%ind  = 4
!
!        self%eqns(5)%name = "rhoE"
!        self%eqns(5)%ind  = 5
!
!        ! Initialize equation set parameters
!        self%rgas = 287.058_rk        ! J/(kg*K)

    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,q,ielem,iface,iblk)
        class(testeq_t),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        integer(ik),        intent(in)      :: ielem, iface, iblk

    end subroutine


    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,q,ielem,iface,iblk)
        class(testeq_t),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        integer(ik),        intent(in)      :: ielem, iface, iblk


    end subroutine



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,q,ielem,iblk)
        class(testeq_t),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        integer(ik),        intent(in)      :: ielem, iblk

    end subroutine


    !==========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,q,ielem,iblk)
        class(testeq_t),    intent(in)      :: self
        class(mesh_t),      intent(in)      :: mesh
        class(expansion_t), intent(inout)   :: q(:)
        integer(ik),        intent(in)      :: ielem, iblk

    end subroutine



end module type_testeq
