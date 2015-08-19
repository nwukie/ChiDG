module atype_bc
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_mesh,          only: mesh_t
    use type_element,       only: element_t
    use atype_equationset,  only: equationset_t
    use atype_solverdata,   only: solverdata_t
    implicit none
    private


    !> Abstract base-type for boundary conditions
    !!
    !!
    !!
    !-------------------------------------------------
    type, public, abstract :: bc_t
        private
        integer(ik), allocatable :: ielems(:)    !> Indices of elements associated with boundary condition
        integer(ik), allocatable :: ifaces(:)    !> Indices of the boundary face for elements elems(ielems)


        logical, public :: isInitialized = .false.

    contains
        procedure :: init               !> Boundary condition initialization
        procedure :: init_spec          !> Call specialized initialization routine
        procedure :: apply              !> Spatial application of the boundary condition

        procedure(compute_interface), deferred :: compute  !> Implements boundary condition calculation
    end type bc_t



    abstract interface
        subroutine compute_interface(self,mesh,iface)
            use mod_kinds,  only: ik
            import bc_t
            import mesh_t
            class(bc_t),    intent(inout)   :: self
            type(mesh_t),   intent(in)      :: mesh
            integer(ik),    intent(in)      :: iface
        end subroutine
    end interface



contains

    !> Initialize boundary condition routine
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  mesh    mesh_t object containing elements and faces
    !!  @param[in]  iface   block face index to which the boundary condition is being applied
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,iface)
        class(bc_t),    intent(inout)       :: self
        type(mesh_t),   intent(in), target  :: mesh
        integer(ik),    intent(in)          :: iface

        type(element_t), pointer    :: elems_p(:,:)   !> Pointer to remap face elements to a plane
        integer(ik)                 :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc, &
                                       xi_begin,  eta_begin, zeta_begin, xi_end, eta_end, zeta_end, &
                                       ixi, ieta, izeta, ierr


        nelem_xi   = mesh%nelem_xi
        nelem_eta  = mesh%nelem_eta
        nelem_zeta = mesh%nelem_zeta

        xi_begin   = 1
        eta_begin  = 1
        zeta_begin = 1

        xi_end   = nelem_xi
        eta_end  = nelem_eta
        zeta_end = nelem_zeta


        !> Compute number of elements associated with the boundary condition
        !! Constrain index ranges for a particular face on the block
        select case (iface)
            case (XI_MIN)                           !> XI_MIN constant
                nelem_bc = nelem_eta * nelem_zeta
                xi_end = 1
            case (XI_MAX)                           !> XI_MAX constant
                nelem_bc = nelem_eta * nelem_zeta
                xi_begin = nelem_xi
            case (ETA_MIN)                          !> ETA_MIN constant
                nelem_bc = nelem_xi * nelem_zeta
                eta_end = 1
            case (ETA_MAX)                          !> ETA_MAX constant
                nelem_bc = nelem_xi * nelem_zeta
                eta_begin = nelem_eta
            case (ZETA_MIN)                         !> ZETA_MIN constant
                nelem_bc = nelem_xi * nelem_eta
                zeta_end = 1
            case (ZETA_MAX)                         !> ZETA_MAX constant
                nelem_bc = nelem_xi * nelem_eta
                zeta_begin = nelem_zeta
            case default
                call signal(FATAL,"bc%init: Invalid block face 'iface'. Valid face indices are iface = [1-6]")
        end select

        !> Allocate storage for element and face indices
        allocate(self%ielems(nelem_bc), self%ifaces(nelem_bc), stat=ierr)
        if (ierr /= 0) call AllocationError

        ielem_bc = 1
        !> Loop over a face of the block and store element indices
        do izeta = zeta_begin,zeta_end
            do ieta = eta_begin,eta_end
                do ixi = xi_begin,xi_end
                    self%ielems(ielem_bc) = mesh%elems_m(ixi,ieta,izeta)%ielem
                    self%ifaces(ielem_bc) = iface
                    ielem_bc = ielem_bc + 1
                end do ! ixi
            end do ! ieta
        end do ! izeta

        call self%init_spec()

        self%isInitialized = .true. !> Set initialization confirmation
    end subroutine
    



    !>
    subroutine apply(self,eqnset,mesh,sdata,iblk)
        class(bc_t),            intent(inout)   :: self
        class(equationset_t),   intent(in)      :: eqnset
        type(mesh_t),           intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: iblk


    end subroutine


    !> Default specialized initialization procedure. This is called from the base bc%init procedure
    !! and can be overwritten by derived types to implement specialized initiailization details.
    subroutine init_spec(self)
        class(bc_t),    intent(inout)   :: self

    end subroutine



end module atype_bc
