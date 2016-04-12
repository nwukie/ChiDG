module mod_grid_tools_two
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, TWO, X_DIR, Y_DIR, Z_DIR, XI_DIR, ETA_DIR, ZETA_DIR, TWO_DIM, THREE_DIM
    use mod_inv

    use type_mesh,          only: mesh_t
    use type_point,         only: point_t
    use type_ivector,       only: ivector_t
    implicit none




contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine compute_element_donor(mesh,point_phys,idom_d,ielem_d,point_comp)
        type(mesh_t),   intent(in)      :: mesh(:)
        type(point_t),  intent(in)      :: point_phys
        integer(ik),    intent(inout)   :: idom_d
        integer(ik),    intent(inout)   :: ielem_d
        type(point_t),  intent(inout)   :: point_comp

        real(rk)    :: res, tol, xgq, ygq, zgq,                 &
                       xmin, xmax, ymin, ymax, zmin, zmax,      &
                       dx, dy, dz, xn, yn, zn, xi, eta, zeta

        integer(ik) :: ncandidates, idom, ielem, inewton, icandidate, spacedim

        type(ivector_t)         :: candidate_domains
        type(ivector_t)         :: candidate_elements

        logical                 :: contained = .false.
        logical                 :: donor_found = .false.

        real(rk)    :: mat(3,3), minv(3,3)
        real(rk)    :: R(3)
        real(rk)    :: dcoord(3)


        tol = 1.e-12_rk

        xgq = point_phys%c1_
        ygq = point_phys%c2_
        zgq = point_phys%c3_

        !
        ! Loop through domains and search for potential donor candidates
        !
        ncandidates = 0
        do idom = 1,size(mesh)


            !
            ! Loop through elements in the current domain
            !
            do ielem = 1,mesh(idom)%nelem


                !
                ! Get bounding coordinates for the current element
                !
                xmin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c1_)
                xmax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c1_)

                ymin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c2_)
                ymax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c2_)

                zmin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c3_)
                zmax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c3_)


                !
                ! Grow bounding box by 10%. Use delta x,y,z instead of scaling xmin etc. in case xmin is 0
                !
                dx = abs(xmax - xmin)  
                dy = abs(ymax - ymin)
                dz = abs(zmax - zmin)


                xmin = xmin - 0.1*dx
                xmax = xmax + 0.1*dx
                ymin = ymin - 0.1*dy
                ymax = ymax + 0.1*dy
                zmin = (zmin-0.001) - 0.1*dz    ! This is to help 2D
                zmax = (zmax+0.001) + 0.1*dz    ! This is to help 2D

                !
                ! Test if gq_node is contained within the bounding coordinates
                !
                contained = ( (xmin < xgq) .and. (xgq < xmax ) .and. &
                              (ymin < ygq) .and. (ygq < ymax ) .and. &
                              (zmin < zgq) .and. (zgq < zmax ) )


                !
                ! If the node was within the bounding coordinates, flag the element as a potential donor
                !
                if ( contained ) then
                   call candidate_domains%push_back(idom) 
                   call candidate_elements%push_back(ielem)
                   ncandidates = ncandidates + 1
                end if


            end do ! ielem

        end do ! idom








        !
        ! Test gq_node on candidate element volume using Newton's method to map to donor local coordinates
        !
        donor_found = .false.
        do icandidate = 1,ncandidates

            idom  = candidate_domains%at(icandidate)
            ielem = candidate_elements%at(icandidate)
            spacedim = mesh(idom)%spacedim

            !
            ! Newton iteration to find the donor local coordinates
            !
            xi   = 0._rk
            eta  = 0._rk
            zeta = 0._rk
            do inewton = 1,20

                !
                ! Compute local cartesian coordinates as a function of xi,eta,zeta
                !
!                xn = mesh_point(mesh(idom)%elems(ielem),X_DIR,xi,eta,zeta)
!                yn = mesh_point(mesh(idom)%elems(ielem),Y_DIR,xi,eta,zeta)
!                zn = mesh_point(mesh(idom)%elems(ielem),Z_DIR,xi,eta,zeta)
                xn = mesh(idom)%elems(ielem)%x(xi,eta,zeta)
                yn = mesh(idom)%elems(ielem)%y(xi,eta,zeta)
                zn = mesh(idom)%elems(ielem)%z(xi,eta,zeta)



                !
                ! Assemble residual vector
                !
                R(1) = -(xn - xgq)
                R(2) = -(yn - ygq)
                R(3) = -(zn - zgq)


                !
                ! Assemble coordinate jacobian matrix
                !
!                mat(1,1) = metric_point(mesh(idom)%elems(ielem),X_DIR,XI_DIR,  xi,eta,zeta)
!                mat(2,1) = metric_point(mesh(idom)%elems(ielem),Y_DIR,XI_DIR,  xi,eta,zeta)
!                mat(3,1) = metric_point(mesh(idom)%elems(ielem),Z_DIR,XI_DIR,  xi,eta,zeta)
!                mat(1,2) = metric_point(mesh(idom)%elems(ielem),X_DIR,ETA_DIR, xi,eta,zeta)
!                mat(2,2) = metric_point(mesh(idom)%elems(ielem),Y_DIR,ETA_DIR, xi,eta,zeta)
!                mat(3,2) = metric_point(mesh(idom)%elems(ielem),Z_DIR,ETA_DIR, xi,eta,zeta)
!                mat(1,3) = metric_point(mesh(idom)%elems(ielem),X_DIR,ZETA_DIR,xi,eta,zeta)
!                mat(2,3) = metric_point(mesh(idom)%elems(ielem),Y_DIR,ZETA_DIR,xi,eta,zeta)
!                mat(3,3) = metric_point(mesh(idom)%elems(ielem),Z_DIR,ZETA_DIR,xi,eta,zeta)

                if ( spacedim == THREE_DIM ) then
                    mat(1,1) = mesh(idom)%elems(ielem)%compute_metric(X_DIR,XI_DIR,  xi,eta,zeta)
                    mat(2,1) = mesh(idom)%elems(ielem)%compute_metric(Y_DIR,XI_DIR,  xi,eta,zeta)
                    mat(3,1) = mesh(idom)%elems(ielem)%compute_metric(Z_DIR,XI_DIR,  xi,eta,zeta)
                    mat(1,2) = mesh(idom)%elems(ielem)%compute_metric(X_DIR,ETA_DIR, xi,eta,zeta)
                    mat(2,2) = mesh(idom)%elems(ielem)%compute_metric(Y_DIR,ETA_DIR, xi,eta,zeta)
                    mat(3,2) = mesh(idom)%elems(ielem)%compute_metric(Z_DIR,ETA_DIR, xi,eta,zeta)
                    mat(1,3) = mesh(idom)%elems(ielem)%compute_metric(X_DIR,ZETA_DIR,xi,eta,zeta)
                    mat(2,3) = mesh(idom)%elems(ielem)%compute_metric(Y_DIR,ZETA_DIR,xi,eta,zeta)
                    mat(3,3) = mesh(idom)%elems(ielem)%compute_metric(Z_DIR,ZETA_DIR,xi,eta,zeta)

                else if ( spacedim == TWO_DIM ) then
                    mat(1,1) = mesh(idom)%elems(ielem)%compute_metric(X_DIR,XI_DIR,  xi,eta,zeta)
                    mat(2,1) = mesh(idom)%elems(ielem)%compute_metric(Y_DIR,XI_DIR,  xi,eta,zeta)
                    mat(3,1) = mesh(idom)%elems(ielem)%compute_metric(Z_DIR,XI_DIR,  xi,eta,zeta)
                    mat(1,2) = mesh(idom)%elems(ielem)%compute_metric(X_DIR,ETA_DIR, xi,eta,zeta)
                    mat(2,2) = mesh(idom)%elems(ielem)%compute_metric(Y_DIR,ETA_DIR, xi,eta,zeta)
                    mat(3,2) = mesh(idom)%elems(ielem)%compute_metric(Z_DIR,ETA_DIR, xi,eta,zeta)
                    mat(1,3) = ZERO
                    mat(2,3) = ZERO
                    mat(3,3) = ONE

                end if



                !
                ! Invert jacobian matrix
                !
                minv = inv(mat)


                !
                ! Compute coordinate update
                !
                dcoord = matmul(minv,R)


                !
                ! Update coordinates
                !
                xi   = xi   + dcoord(1)
                eta  = eta  + dcoord(2)
                zeta = zeta + dcoord(3)


                !
                ! Compute residual coordinate norm
                !
                res = norm2(R)


                !
                ! Exit if converged
                !
                if ( res < tol ) then
                    donor_found = .true.

                    idom_d  = idom
                    ielem_d = ielem
                    call point_comp%set(xi,eta,zeta)

                    exit
                end if


                !
                ! Limit computational coordinates, in case they go out of bounds.
                !
                if ( xi   >  ONE ) xi   =  ONE
                if ( xi   < -ONE ) xi   = -ONE
                if ( eta  >  ONE ) eta  =  ONE
                if ( eta  < -ONE ) eta  = -ONE
                if ( zeta >  ONE ) zeta =  ONE
                if ( zeta < -ONE ) zeta = -ONE

            end do ! inewton

            if (donor_found) then
                exit
            end if

        end do ! icandidate




        !
        ! Sanity check on donors and set donor_element location
        !
        if ( .not. donor_found ) call chidg_signal(FATAL,"compute_element_donor: No donor found for gq_node")




    end subroutine compute_element_donor
    !******************************************************************************************







end module mod_grid_tools_two
