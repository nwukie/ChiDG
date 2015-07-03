module mod_prototype
    implicit none


    do blocks (1,7)

        do elements

            ! Seed derivative arrays
            call seed_derivatives(block)


            ! Compute boundary flux
            do faces (1,6)
                call eqnset%compute_boundary_average_flux(mesh,iface)
                call eqnset%compute_boundary_upwind_flux(mesh,iface)
            end do faces


            call eqnset%compute_volume_flux(mesh,ielem)
            call eqnset%compute_volume_source(mesh,ielem)


        end do ! elements

    end do ! lin blocks




end module mod_prototype
