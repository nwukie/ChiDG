module mod_project
    use mod_kinds,      only: rk,ik
    use mod_quadrature, only: GQ, get_quadrature, compute_nnodes_gq
    use mod_grid_tools, only: compute_discrete_coordinates
    use type_point,     only: point_t
    use type_expansion, only: expansion_t
    use atype_function, only: function_t

    implicit none



contains

    !>  Project values from a function evaluation to the polynomial basis
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  fcn     Incoming function to evaluate. Arguments should be in cartesian coordinates(x, y, z).
    !!  @param[in]  nterms  Number of terms in the basis being projected to
    !!  @param[in]  cmodes  Modal coefficients of coordinate expansion to evaluate cartesian coordinates
    !!  @param[out] fmodes  Modal coefficients of the projected function
    !------------------------------------------------------------------------------
    subroutine project_function_xyz(fcn,nterms,cmodes,fmodes)
        class(function_t),  intent(in)  :: fcn
        integer(ik),        intent(in)  :: nterms
        type(expansion_t),  intent(in)  :: cmodes           ! Expansion contains x-modes, y-modes, and z-modes
        real(rk),           intent(out) :: fmodes(nterms)

        type(point_t),  allocatable     :: pts(:)
        real(rk),       allocatable     :: fvals(:)
        integer(ik)                     :: nterms_1d, nn_face, nn_vol, igq, ierr
        integer(ik)                     :: gq_p, gq_c
        logical                         :: has_correct_nodes_terms


        ! Compute number of face and volume quadrature nodes
        call compute_nnodes_gq(nterms,cmodes%nterms,nn_face,nn_vol)


        ! Find the correct quadrature instance for projecting the function values
        ! and for evaluating the coordinates at discrete points.
        call get_quadrature(nterms,nn_vol,nn_face,gq_p)
        call get_quadrature(cmodes%nterms,nn_vol,nn_face,gq_c)


        ! Compute discrete cartesian coordinates
        allocate(pts(nn_vol), fvals(nn_vol), stat=ierr)
        if (ierr /= 0) stop "Error: project_function_xyz -- allocation error"
        call compute_discrete_coordinates(cmodes,gq_c,pts)


        ! Call function for evaluation and multiply by quadrature weights
        fvals = fcn%calc(pts)  *  GQ(gq_p)%vol%weights

        ! Project
        fmodes = matmul(transpose(GQ(gq_p)%vol%val),fvals) / GQ(gq_p)%vol%dmass

    end subroutine








    subroutine project_volume()



    end subroutine

















end module mod_project
