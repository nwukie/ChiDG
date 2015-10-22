!> This module contains procedures for initializing and maintaining the Chimera
!! interfaces.
!!
!!  @author Nathan A. Wukie
!!
!!
!---------------------------------------------------------------
module mod_chimera
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: NFACES, ORPHAN, CHIMERA

    use type_mesh,      only: mesh_t









contains


    !> Routine for detecting Chimera faces. 
    !!
    !! Routine flags face as a Chimera face if it has an ftype==ORPHAN, 
    !! indicating it is not an interior face and it has not been assigned a 
    !! boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   mesh    Array of mesh types. One for each domain.
    !!
    !------------------------------------------------------------------------
    subroutine detect_chimera_interfaces(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: idom, ndom, ielem, iface
        logical     :: orphan_face = .false.

        !
        ! Get number of domains
        !
        ndom = size(mesh)

        
        !
        ! Loop through each domain
        !
        do idom = 1,ndom

            !
            ! Loop through each element
            !
            do ielem = 1,mesh(idom)%nelem


                !
                ! Loop through each face
                !
                do iface = 1,NFACES

                    !
                    ! Test if the current face is unattached
                    !
                    orphan_face = ( mesh(idom)%faces(ielem,iface)%ftype == ORPHAN ) 

                    !
                    ! If orphan_face, set as Chimera face so it can search for donors in other domains
                    !
                    if (orphan_face) then

                        mesh(idom)%faces(ielem,iface)%ftype = CHIMERA

                    end if


                end do ! iface

            end do ! ielem

        end do ! idom


    end subroutine










end module mod_chimera
