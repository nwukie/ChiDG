module mod_spatial
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: NFACES
    use type_domain,    only: domain_t

    implicit none


contains

    subroutine update_space(domain)
        type(domain_t), intent(inout)   :: domain

        integer(ik) :: iblk, ielem, iface, nelem


        associate ( mesh => domain%mesh, sdata => domain%sdata)

            nelem = mesh%nelem
            !------------------------------------------------------------------------------------------
            !                                      Interior Scheme
            !------------------------------------------------------------------------------------------

            !> Loop through given element and neighbors and compute the corresponding linearization
            !> XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG
            do iblk = 1,7




                !> Loop through elements in the domain
                do ielem = 1,nelem



                    ! For the current element, compute the contributions from boundary integrals
                    do iface = 1,NFACES
                        call domain%eqnset%compute_boundary_average_flux(mesh,sdata,ielem,iface,iblk)
                        call domain%eqnset%compute_boundary_upwind_flux( mesh,sdata,ielem,iface,iblk)
                    end do !face



                    ! For the current element, compute the contributions from volume integrals
                    call domain%eqnset%compute_volume_flux(  mesh,sdata,ielem,iblk)
                    call domain%eqnset%compute_volume_source(mesh,sdata,ielem,iblk)



                end do !elem




            end do !block


            !-------------------------------------------------------------------------------------------
            !                                     Boundary Scheme
            !-------------------------------------------------------------------------------------------


        end associate

    end subroutine



end module mod_spatial
