module mod_quadrature
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI

    use type_quadrature,    only: quadrature_t

    implicit none


    ! Arrays of quadrature types, so an element can choose the order
    ! of the quadrature, based on it's polynomial approximation
    integer(ik), parameter :: ngq = 100

    type(quadrature_t), target, save :: GQ(ngq)


contains

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






    !>  Routine to find a quadrature instance or initialize a new one and return
    !!  its location in the global quadrature array 'igq' in GQ(igq)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]      nterms  Number of terms in a polynomial expansion
    !!  @param[in]      nn_v    Number of volume quadrature nodes required
    !!  @param[in]      nn_f    Number of face quadrature nodes required
    !!  @param[inout]   gqout   Integer index of the selected quadrature instance in
    !!                          the global quadrature instance array, GQ
    !-------------------------------------------------------------------------------------
    subroutine get_quadrature(nterms,nn_v,nn_f,gqout)
        integer(ik),    intent(in)       :: nterms, nn_v, nn_f
        integer(ik),    intent(inout)    :: gqout

        integer(ik) :: igq
        logical     :: has_correct_nodes_terms

        do igq = 1,size(GQ)

            if (GQ(igq)%isInitialized) then
                has_correct_nodes_terms = (GQ(igq)%nterms == nterms) .and. (GQ(igq)%nnodes_v == nn_v)

                if (has_correct_nodes_terms) then
                    gqout = igq
                    exit
                end if
            else
                ! If we are here, then no initialized GQ instance was found that met the requirements,
                ! so, we initialize a new one.
                call GQ(igq)%init(nn_f,nn_v,nterms)
                gqout = igq
                exit
            end if

        end do

    end subroutine












end module mod_quadrature
