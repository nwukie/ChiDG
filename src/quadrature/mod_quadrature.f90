module mod_quadrature
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI

    use type_quadrature,    only: quadrature_t

    implicit none


    ! Arrays of quadrature types, so an element can choose the order
    ! of the quadrature, based on it's polynomial approximation
    integer(ik), parameter :: ngq = 5

    type(quadrature_t), target, save :: GQ(ngq)
    type(quadrature_t), target, save :: GQMESH(ngq)


contains

!    subroutine initialize_quadrature()
!        use mod_io,     only: nterms_sol1d, nterms_sol2d, nterms_sol3d, &
!                              nterms_mesh1d,nterms_mesh2d,nterms_mesh3d
!        integer(ik)    :: igq
!        integer(ik)    :: nnodes_face, nnodes_vol
!        integer(ik)    :: ierr
!
!        call compute_nnodes_integration(1,nterms_sol3d,nnodes_face,nnodes_vol)
!
!        do igq = 1,ngq
!            call GQ(igq)%init(nnodes_face,nnodes_vol,nterms_sol3d)
!            call GQMESH(igq)%init(nnodes_face,nnodes_vol,nterms_mesh3d)
!        end do
!
!    end subroutine initialize_quadrature







    !>  Compute the number of quadrature nodes to use in order to accurately compute integrals
    !!  in the DG discretization. Based on solution and coordinate polynomial orders.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in] nterms_s     Number of terms in the solution expansion
    !!  @param[in] nterms_c     Number of terms in the coordinate expansion
    !!  @param[out] nnodes_face Number of quadrature nodes defined for a face
    !!  @param[out] nnodes_vol  Number of quadrature nodes defined for a volume
    !------------------------------------------------------------------------------------
    subroutine compute_nnodes_gq(nterms_s,nterms_c,nnodes_face,nnodes_vol)
        use mod_io,                     only: gq_rule

        integer(ik), intent(in)        :: nterms_s, nterms_c
        integer(ik), intent(out)       :: nnodes_face, nnodes_vol
        integer(ik)                    :: nterms1d,nnodes1d,nnodes2d,nnodes3d

        ! Find number of terms in the 1d expansion
        nterms1d = 0
        do while (nterms1d*nterms1d*nterms1d /= nterms_s)
            nterms1d = nterms1d + 1
        end do
        if (nterms1d*nterms1d*nterms1d > nterms_s) stop "Error: compute_nnodes_integration in term count"


        ! Compute number of 1D nodes, based on integration rule
        select case (gq_rule)
            case(1)
                ! Collocation quadrature
                nnodes1d = nterms1d
            case(2)
!                nnodes1d = 3*nterms1d/2 + 1
                nnodes1d = 3*nterms1d+1
            case(3)
                nnodes1d = 2*nterms1d + 1
            case default
                print*, "Error: compute_nnodes_integration - valid rules are (1,2,3)"
                stop
        end select

        ! Compute face and volume nodes
        nnodes2d = nnodes1d*nnodes1d
        nnodes3d = nnodes1d*nnodes1d*nnodes1d

        nnodes_face = nnodes2d
        nnodes_vol  = nnodes3d


    end subroutine

end module mod_quadrature
