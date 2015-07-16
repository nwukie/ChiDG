!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie


program chidg
    use mod_kinds,              only: ik,rk
    use mod_constants,          only: ZERO, ONE, FIVE
    use mod_tecio,              only: write_tecio_variable
    use type_point,             only: point_t
    use type_domain,            only: domain_t
    use mod_element_mapping,    only: compute_element_mappings
    use messenger

    !-----------------------------------------------------------
    !-----------------------------------------------------------
    ! Variable declarations
    implicit none
    type(domain_t)  :: dom(1)
    integer(ik)     :: ipt, ixi, ieta, izeta
    real(rk)        :: x(8), y(8), z(8)
    type(point_t)   :: pts(2,2,2)


    call compute_element_mappings

    !> Initialize coredg static data
    x = [ZERO, FIVE, ZERO, FIVE, ZERO, FIVE, ZERO, FIVE]
    y = [ZERO, ZERO, ONE, ONE, ZERO, ZERO, ONE, ONE]
    z = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE]

    ipt = 1
    do izeta = 1,2
        do ieta = 1,2
            do ixi = 1,2
                call pts(ixi,ieta,izeta)%set(x(ipt),y(ipt),z(ipt))
                ipt = ipt + 1
            end do
        end do
    end do

    call dom(1)%init_geom(8,pts)
    call dom(1)%init_sol('testeq',27)


    call write_tecio_variable(dom(1),'test.plt',1)




    !======================================================
    !
    !   Solve
    !
    !======================================================







    !=======================================================
    !
    !   Post-execution tasks
    !
    !=======================================================

!    print*, "Completed in: ", elapsed, " Seconds"
!    print*, "======================="
!    print*, "Done!"
!    print*, "======================="

end program chidg
