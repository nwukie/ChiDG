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

        integer(ik) :: idom, ndom, ielem, iface, nchimera_faces
        logical     :: orphan_face = .false.

        !
        ! Get number of domains
        !
        ndom = size(mesh)

        
        !
        ! Loop through each domain
        !
        do idom = 1,ndom


            nchimera_faces = 0
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



                        ! Increment domain-local chimera face count
                        nchimera_faces = nchimera_faces + 1

                        ! Set face-type to CHIMERA
                        mesh(idom)%faces(ielem,iface)%ftype = CHIMERA

                        ! Set domain-local Chimera identifier. Really, just the index order which they were detected in, starting from 1.
                        mesh(idom)%faces(ielem,iface)%ChiID = nchimera_faces



                    end if


                end do ! iface

            end do ! ielem

            ! Set total number of Chimera faces detected for domain - idom
            mesh(idom)%chimera%recv%nfaces = nchimera_faces

        end do ! idom


    end subroutine













    subroutine accumulate_chimera_interfaces(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)






    end subroutine












































end module mod_chimera
