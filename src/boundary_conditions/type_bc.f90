module type_bc
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, BOUNDARY
    use type_mesh,          only: mesh_t
    use type_element,       only: element_t
    use type_equationset,   only: equationset_t
    use type_solverdata,    only: solverdata_t
    use type_dict,          only: dict_t
    use type_properties,    only: properties_t
    implicit none
    private


    !> Abstract base-type for boundary conditions
    !!  - contains a list of associated element indices
    !!  - contains a list of face indices
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------
    type, public, abstract :: bc_t
        integer(ik), allocatable    :: ielems(:)                !< Indices of elements associated with boundary condition
        integer(ik), allocatable    :: ifaces(:)                !< Indices of the boundary face for elements elems(ielems)
        logical, public             :: isInitialized = .false.  !< Logical switch for indicating the boundary condition initializaiton status

    contains
        procedure :: init                                       !< Boundary condition initialization
        procedure :: init_spec                                  !< Call specialized initialization routine
        procedure :: apply                                      !< Spatial application of the boundary condition
        procedure(compute_interface), deferred :: compute       !< Implements boundary condition calculation

    end type bc_t



    abstract interface
        subroutine compute_interface(self,mesh,sdata,prop,idom,ielem,iface,iblk)
            use mod_kinds,  only: ik
            import bc_t
            import mesh_t
            import solverdata_t
            import properties_t

            class(bc_t),            intent(inout)   :: self
            type(mesh_t),           intent(in)      :: mesh(:)
            type(solverdata_t),     intent(inout)   :: sdata
            class(properties_t),    intent(inout)   :: prop
            integer(ik),            intent(in)      :: idom
            integer(ik),            intent(in)      :: ielem
            integer(ik),            intent(in)      :: iface
            integer(ik),            intent(in)      :: iblk
        end subroutine
    end interface



contains

    !> Initialize boundary condition routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    mesh_t object containing elements and faces
    !!  @param[in]  iface   block face index to which the boundary condition is being applied
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,iface,options)
        class(bc_t),            intent(inout)       :: self
        type(mesh_t),           intent(inout)       :: mesh
        integer(ik),            intent(in)          :: iface 
        type(dict_t), optional, intent(in)          :: options

        
        integer(ik)                 :: nelem_xi, nelem_eta, nelem_zeta, nelem_bc, ielem_bc, & 
                                       xi_begin, eta_begin, zeta_begin, xi_end, eta_end, zeta_end, & 
                                       ixi, ieta, izeta, ierr, ielem, ielem_test
        
        nelem_xi   = mesh%nelem_xi
        nelem_eta  = mesh%nelem_eta
        nelem_zeta = mesh%nelem_zeta

        xi_begin   = 1
        eta_begin  = 1
        zeta_begin = 1

        xi_end   = nelem_xi
        eta_end  = nelem_eta
        zeta_end = nelem_zeta


        !
        ! Compute number of elements associated with the boundary condition
        ! Constrain index ranges for a particular face on the block
        !
        select case (iface)
            case (XI_MIN)                           ! XI_MIN constant
                nelem_bc = nelem_eta * nelem_zeta
                xi_end = 1
            case (XI_MAX)                           ! XI_MAX constant
                nelem_bc = nelem_eta * nelem_zeta
                xi_begin = nelem_xi
            case (ETA_MIN)                          ! ETA_MIN constant
                nelem_bc = nelem_xi * nelem_zeta
                eta_end = 1
            case (ETA_MAX)                          ! ETA_MAX constant
                nelem_bc = nelem_xi * nelem_zeta
                eta_begin = nelem_eta
            case (ZETA_MIN)                         ! ZETA_MIN constant
                nelem_bc = nelem_xi * nelem_eta
                zeta_end = 1
            case (ZETA_MAX)                         ! ZETA_MAX constant
                nelem_bc = nelem_xi * nelem_eta
                zeta_begin = nelem_zeta
            case default
                call chidg_signal(FATAL,"bc%init: Invalid block face 'iface'. Valid face indices are iface = [1-6]")
        end select


        !
        ! Allocate storage for element and face indices
        !
        allocate(self%ielems(nelem_bc), self%ifaces(nelem_bc), stat=ierr)
        if (ierr /= 0) call AllocationError


        ielem_bc = 1
        !
        ! Loop over a face of the block and store element indices
        !
        do izeta = zeta_begin,zeta_end
            do ieta = eta_begin,eta_end
                do ixi = xi_begin,xi_end
                    ielem = ixi + nelem_xi*(ieta-1) + nelem_xi*nelem_eta*(izeta-1)

                    self%ielems(ielem_bc) = ielem
                    self%ifaces(ielem_bc) = iface
                    ielem_bc = ielem_bc + 1


                    !
                    ! Set face to boundary condition face
                    !
                    mesh%faces(ielem,iface)%ftype = BOUNDARY
                end do ! ixi
            end do ! ieta
        end do ! izeta


        !
        ! Call user-specialized boundary condition initializatio        
        !
        call self%init_spec(mesh,iface,options)

        self%isInitialized = .true. ! Set initialization confirmation
    end subroutine init
    






    !>  Apply boundary condition to the mesh and solution
    !!      - Loops through the associated elements(faces) and calls the specialized bc_t%compute
    !!        procedure for computing the rhs and linearization.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    mesh_t defining elements and faces
    !!  @param[inout]   sdata   solverdata_t containing solution, rhs, and linearization(lin) data
    !!  @param[in]      iblk    Block of the linearization for the current element that is being computed (XI_MIN, XI_MAX, eta.)
    !!  @param[inout]   prop    properties_t object containing equationset properties and material_t objects
    !!
    !--------------------------------------------------------------------
    subroutine apply(self,mesh,sdata,prop,idom,iblk)
        class(bc_t),            intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        class(solverdata_t),    intent(inout)   :: sdata
        class(properties_t),    intent(inout)   :: prop
        integer(ik),            intent(in)      :: idom
        integer(ik),            intent(in)      :: iblk

        integer(ik) :: ielem_bc, ielem, iface

        !
        ! Loop through associated boundary condition elements and call compute routine for the boundary flux calculation
        !
        do ielem_bc = 1,size(self%ielems)
            ielem = self%ielems(ielem_bc)   ! Get index of the element being operated on
            iface = self%ifaces(ielem_bc)   ! Get face index of element 'ielem' that is being operated on

            !
            ! For the current boundary element(face), call specialized compute procedure
            !
            call self%compute(mesh,sdata,prop,idom,ielem,iface,iblk)

        end do


    end subroutine apply







    !> Default specialized initialization procedure. This is called from the base bc%init procedure
    !! and can be overwritten by derived types to implement specialized initiailization details.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   mesh        mesh_t object containing elements and faces
    !!  @param[in]      iface       block face index to which the boundary condition is being applied
    !!  @param[in]      options     dictionary object containing boundary condition options
    !!
    !--------------------------------------------------------------------------------
    subroutine init_spec(self,mesh,iface,options)
        class(bc_t),            intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: iface
        type(dict_t), optional, intent(in)      :: options




    end subroutine init_spec






end module type_bc
