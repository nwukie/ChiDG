module mod_kirchoffs
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: NFACES, ZERO, ONE, TWO, PI
    use type_chidg,                 only: chidg_t
    use type_mesh,                  only: mesh_t
    use type_point,                 only: point_t
    use type_bcset,                 only: bcset_t
    use type_equationset_wrapper,   only: equationset_wrapper_t
    use type_solverdata,            only: solverdata_t
    use type_point,                 only: point_t

    use mod_interpolate,            only: interpolate_face
    use mod_primitive_linearized_euler,     only: cbar, omega
    use mod_io,                             only: nterms_s
    implicit none




contains


    subroutine kirchoff(chidg_file)
        character(*),   intent(in)  :: chidg_file

        type(chidg_t)   :: chidg
        integer(ik)     :: res, ierr, itheta, fileunit
        real(rk)        :: theta_min, theta_max, r, dtheta

        real(rk),       dimension(:),   allocatable :: theta
        complex(rk),    dimension(:),   allocatable :: pressures
        type(point_t),  dimension(:),   allocatable :: points

        nterms_s = 5*5*5

        !
        ! Initialize ChiDG environment
        !
        call chidg%init('env')


        !
        ! Read grid data from file
        !
        call chidg%read_grid(chidg_file,3)


        !
        ! Read boundary conditions
        !
        call chidg%read_boundaryconditions(chidg_file)


        !
        ! Initialize solution data storage
        !
        !call chidg%data%init_sdata()
        call chidg%initialize_solution_domains(nterms_s)
        call chidg%initialize_solution_solver()

        print*, chidg%data%eqnset(1)%item%name

        !
        ! Initialize solution
        !
        call chidg%read_solution(chidg_file)

        

        print*, chidg%data%eqnset(1)%item%name
        !
        ! Wrap-up initialization activities
        !
        !call chidg%init('finalize')

        print*, chidg%data%eqnset(1)%item%name

        res       = 4000

        allocate(theta(res), points(res), stat=ierr)
        if ( ierr /= 0 ) call AllocationError

        theta_min = ZERO
        !theta_max = 120._rk*PI/180._rk
        theta_max = PI/TWO

        r        = 46._rk
        theta    = ZERO
        dtheta = (theta_max - theta_min)/real(res,kind=rk)

        do itheta = 2,res
            theta(itheta) = theta(itheta-1) + dtheta
        end do

        points(:)%c1_ = r*cos(theta)
        points(:)%c2_ = r*sin(theta)
        points(:)%c3_ = ZERO

        print*, chidg%data%eqnset(1)%item%name
        pressures = compute_kirchoffs_integral(chidg%data%mesh,chidg%data%bcset,chidg%data%eqnset,chidg%data%sdata,points) 


        print*, pressures


        !
        ! Close ChiDG
        !
        call chidg%close()



        open(newunit=fileunit, file='p.out')

        do itheta = 1,res
            write(fileunit,*) theta(itheta), sqrt( real(pressures(itheta))**TWO + imag(pressures(itheta))**TWO ) / sqrt(TWO)
        end do


    end subroutine kirchoff
    !******************************************************************************************************


















    !>  Routine for computing a Kirchoffs integral for propagating acoustic disturbances
    !!  to far-field observation points.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/9/2016
    !!
    !!  @param[in]  mesh    Array of mesh data structures
    !!  @param[in]  eqnset  Array of equationsets, one for each mesh
    !!  @param[in]  sdata   Solver data for all domains
    !!  @param[in]  points  Observer points that pressure is being propagated to
    !!
    !!  @return     pressures   Array of complex pressure associated with observer points
    !!
    !---------------------------------------------------------------------------------
    function compute_kirchoffs_integral(mesh,bcset,eqnset,sdata,points) result(pressures)
        type(mesh_t),                   intent(in)  :: mesh(:)
        type(bcset_t),                  intent(in)  :: bcset(:)
        type(equationset_wrapper_t),    intent(in)  :: eqnset(:)
        type(solverdata_t),             intent(in)  :: sdata
        type(point_t),                  intent(in)  :: points(:)

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: irho_r, irho_i
        integer(ik)    :: iu_r,   iu_i
        integer(ik)    :: iv_r,   iv_i
        integer(ik)    :: iw_r,   iw_i
        integer(ik)    :: ip_r,   ip_i

        integer(ik) :: ndomains, idom, ielem, iface, ibc, nnodes, iobs, ierr
        integer(ik) :: idom_l, ielem_l, iface_l
        logical     :: primitive_linearized_euler = .false.
        logical     :: kirchoff_surface = .false.
        complex(rk) :: pressures(size(points)), imag, integral
        real(rk)    :: xo, yo, zo

        real(rk), dimension(:), allocatable     ::  &
                    rho_r, u_r, v_r, w_r, p_r,      &
                    rho_i, u_i, v_i, w_i, p_i,      &
                    xs, ys, zs, r, rx, ry, rz,      &
                    dpi_dx, dpi_dy, dpi_dz,                   &
                    dpr_dx, dpr_dy, dpr_dz,                   &
                    nx, ny, nz, nxi, neta, nzeta, rdotn

        complex(rk), dimension(:), allocatable  ::                  &
                    rho, u, v, w, p,                                &
                    dpdx, dpdy, dpdz, dpdn, integrand,    &
                    leading_term, one_over_r_term, one_over_r2_term







        imag      = cmplx(ZERO, ONE )
        pressures = cmplx(ZERO, ZERO)


        ndomains = size(mesh)



        !
        ! Loop domains
        !
        do idom = 1,ndomains


            irho_r = eqnset(idom)%item%prop%get_eqn_index("rho_r")
            iu_r   = eqnset(idom)%item%prop%get_eqn_index("u_r")
            iv_r   = eqnset(idom)%item%prop%get_eqn_index("v_r")
            iw_r   = eqnset(idom)%item%prop%get_eqn_index("w_r")
            ip_r   = eqnset(idom)%item%prop%get_eqn_index("p_r")

            irho_i = eqnset(idom)%item%prop%get_eqn_index("rho_i")
            iu_i   = eqnset(idom)%item%prop%get_eqn_index("u_i")
            iv_i   = eqnset(idom)%item%prop%get_eqn_index("v_i")
            iw_i   = eqnset(idom)%item%prop%get_eqn_index("w_i")
            ip_i   = eqnset(idom)%item%prop%get_eqn_index("p_i")





            !
            ! Check domain has the correct equation set for the integra
            !
            primitive_linearized_euler = ( eqnset(idom)%item%get_name() == 'PrimitiveLinearizedEuler' )
            if ( .not. primitive_linearized_euler ) call chidg_signal(FATAL,"compute_kirchoffs_integral: domain does not have correct equation set")



            ! 
            ! Loop boundary conditions for domain
            !
            do ibc = 1,size(bcset(idom)%bcs)

            
                ! Check if boundary condition is Kirchoff boundary 
                kirchoff_surface = ( bcset(idom)%bcs(ibc)%bc%get_name() == "Kirchoff" )



                if ( kirchoff_surface ) then




                    ! For each face, compute Kirchoff integral
                    do iface = 1,size(bcset(idom)%bcs(ibc)%bc%faces)

                        ! Get face-local indices
                        idom_l  = idom
                        ielem_l = bcset(idom)%bcs(ibc)%bc%elems(iface)
                        iface_l = bcset(idom)%bcs(ibc)%bc%faces(iface)



                        ! Get gq nodes and allocate storage
                        nnodes = mesh(idom_l)%faces(ielem_l,iface_l)%gq%face%nnodes
                        allocate(rho_r(nnodes), u_r(nnodes), v_r(nnodes), w_r(nnodes), p_r(nnodes), &
                                 rho_i(nnodes), u_i(nnodes), v_i(nnodes), w_i(nnodes), p_i(nnodes), &
                                 dpdx(nnodes), dpdy(nnodes), dpdz(nnodes),                          &
                                 r(nnodes), &
                                 leading_term(nnodes), one_over_r_term(nnodes), one_over_r2_term(nnodes),  &
                                 rho(nnodes),   u(nnodes),   v(nnodes),   w(nnodes),   p(nnodes), stat=ierr)
                        if ( ierr /= 0 ) call AllocationError






                        ! Get variables at quadrature nodes
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, irho_r,  rho_r )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iu_r,    u_r   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iv_r,    v_r   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iw_r,    w_r   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, ip_r,    p_r   )

                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, irho_i,  rho_i )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iu_i,    u_i   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iv_i,    v_i   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, iw_i,    w_i   )
                        call interpolate_face(mesh,sdata%q, idom_l, ielem_l, iface_l, ip_i,    p_i   )



                        ! Assemble complex variables
                        rho = cmplx(rho_r, rho_i)
                        u   = cmplx(u_r, u_i)
                        v   = cmplx(v_r, v_i)
                        w   = cmplx(w_r, w_i)
                        p   = cmplx(p_r, p_i)

                        ! Get surface(source) coordinates
                        xs = mesh(idom_l)%faces(ielem_l,iface_l)%quad_pts(:)%c1_
                        ys = mesh(idom_l)%faces(ielem_l,iface_l)%quad_pts(:)%c2_
                        zs = mesh(idom_l)%faces(ielem_l,iface_l)%quad_pts(:)%c3_



                        ! Get unit normal vectors
                        nx = mesh(idom_l)%faces(ielem_l,iface_l)%unorm(:,1)
                        ny = mesh(idom_l)%faces(ielem_l,iface_l)%unorm(:,2)
                        nz = mesh(idom_l)%faces(ielem_l,iface_l)%unorm(:,3)



                        ! Compute derivatives of pressure
                        dpr_dx = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdx, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_r) )
                        dpi_dx = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdx, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_i) )
                        dpr_dy = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdy, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_r) )
                        dpi_dy = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdy, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_i) )
                        dpr_dz = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdz, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_r) )
                        dpi_dz = matmul(mesh(idom_l)%faces(ielem_l,iface_l)%dtdz, sdata%q%dom(idom_l)%lvecs(ielem_l)%getvar(ip_i) )



        

                        dpdx = cmplx(dpr_dx, dpi_dx)
                        dpdy = cmplx(dpr_dy, dpi_dy)
                        dpdz = cmplx(dpr_dz, dpi_dz)


                        ! Compute derivative of pressure in the direction of the face normal
                        dpdn = dpdx*nx + dpdy*ny + dpdz*nz




                        ! Compute integral contribution to each observer point
                        do iobs = 1,size(points)


                            ! Get coordinates for current observer
                            xo = points(iobs)%c1_
                            yo = points(iobs)%c2_
                            zo = points(iobs)%c3_

                            ! Compute distance from source to current observer for each quadrature point
                            r = sqrt( (xo-xs)**TWO  +  (yo-ys)**TWO  +  (zo-zs)**TWO )
    
                            ! Compute radius vector at each quadrature point
                            rx = xo-xs
                            ry = yo-ys
                            rz = zo-zs


                            ! Compute cos(theta), r dot n
                            !
                            !   TODO: FIX THIS FOR CARTESIAN NORMALS
                            !
                            rdotn = rx*nx + ry*ny + rz*nz

                            
                            leading_term     = exp(imag*omega*r/cbar)
                            one_over_r_term  = (ONE/r)*(-(imag*omega/cbar)*rdotn*p - dpdn) 
                            one_over_r2_term = p*rdotn/(r**TWO)


                            ! Compute Kirchhoff Integrand
                            integrand = leading_term*( one_over_r_term  +  one_over_r2_term )

                            ! Compute integral contribution from current face
                            integral = sum( integrand * mesh(idom_l)%faces(ielem_l,iface_l)%jinv * mesh(idom_l)%faces(ielem_l,iface_l)%gq%face%weights(:,iface_l) )

                            ! Contribute to current observer
                            pressures(iobs) = pressures(iobs) + integral


                        end do !iobs






                        ! Get gq nodes and allocate storage
                        deallocate(rho_r, u_r, v_r, w_r, p_r, &
                                   rho_i, u_i, v_i, w_i, p_i, &
                                   dpdx,  dpdy, dpdz,         &
                                   r, &
                                   leading_term,    one_over_r_term, one_over_r2_term, &
                                   rho,   u,   v,   w,   p)




                    end do ! iface








                end if  ! if kirchoff_surface



            end do ! ibc




        end do ! idom




    end function compute_kirchoffs_integral
    !**********************************************************************************






end module mod_kirchoffs
