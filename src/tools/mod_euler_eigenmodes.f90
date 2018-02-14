module mod_euler_eigenmodes
#include <messenger.h>
    use mod_constants,              only: ZERO, TWO, ONE
    use mod_function,               only: create_function
    use mod_spatial,                only: update_space
    use type_chidg,                 only: chidg_t
    use type_function,              only: function_t
    use mod_file_utilities,         only: delete_file
    use mod_inv,                    only: zinv
    use mod_gauss_legendre,         only: gl_nodes, gl_weights
    use mod_legendre,               only: legendre_val1D, dlegendre_val1D
    use mod_gridspace,              only: linspace
    use mod_fluid,                  only: gam
!    use mod_primitive_linearized_euler
    use mod_io
    implicit none

    external ZGGEV
    external ZGEEV


    !
    ! Azimuthal orther
    !
    integer(ik) :: mod_m = 2   ! Could get set by boundary conditions
    !integer(ik) :: mod_m = 0    ! Could get set by boundary conditions
    integer(ik) :: mod_n = 1    ! Could get set by boundary conditions
    real(rk), parameter :: romega = 3441.9548_rk
    !real(rk), parameter :: romega = 0._rk
    real(rk), parameter :: iomega = -10.e-5_rk*romega


    !
    ! Dimensional mean flow
    !
    real(rk), parameter :: rho_d = 1.2212179_rk
    real(rk), parameter :: u_d   = 0._rk
    real(rk), parameter :: v_d   = 0._rk
    real(rk), parameter :: w_d   = 103.2586448_rk
    real(rk), parameter :: p_d   = 103341.6659_rk


    !
    ! Mean primitives
    !
    real(rk), parameter :: rhobar = rho_d
    real(rk), parameter :: ubar   = u_d
    real(rk), parameter :: vbar   = v_d
    real(rk), parameter :: wbar   = w_d
    real(rk), parameter :: pbar   = p_d
    real(rk), parameter :: cbar   = sqrt(gam * pbar / rhobar)
    real(rk), parameter :: Mbar   = ubar/cbar






contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !--------------------------------------------------------------------
    subroutine compute_euler_eigenmodes()

        integer,  parameter :: res = 50
        real(rk), parameter :: ri = 0.25_rk
        real(rk), parameter :: ro = 1.0_rk
        real(rk), parameter :: dr = (ro-ri)/real(res-1,rk)
        integer,  parameter :: nfields = 5
        integer             :: matrix_size, work_size, handle, info, dof, nterms, iterm, inode, idof
        integer :: iia, iib, iic, iid, iie, iif, iiia, iiib, iiic, iiid, iiie, iiif

        real(rk),         allocatable,   dimension(:)   :: nodes, r
        real(rk),         allocatable,   dimension(:,:) :: mass, rmass, stiff, bound1, bound2, zero_matrix, eigenvector_mag, val, dat, stencil, identity, ridentity, bc
        complex(kind=8),  allocatable,   dimension(:,:) :: vl, vr, A, B, Binv, system, M1, M2, M3, M4, M5, A_reduced, B_reduced, A_rr, B_rr
        complex(kind=8),  allocatable,   dimension(:)   :: work, lambda, alpha, beta
        double precision, allocatable,   dimension(:)   :: rwork
        integer(ik),      allocatable,   dimension(:)   :: ind
        integer(ik) :: ifield, ideriv, ierr, nterms_1d, ind_col, ind_row, irow_small_start, icol_small_start, irow, i, ivec, istart_big, istart_small, icol, icol_start, irow_start
        logical :: file_exists


        !
        ! Compute grid
        !
        r = linspace(ri,ro,res)


        !
        ! Compute degrees of freedom
        !
        dof = res*nfields


        !
        ! Allocate global storage
        !
        allocate(                   &
                 A(   dof,dof),     &
                 B(   dof,dof),     &
                 Binv(dof,dof),     &
                 M1(  dof,dof),     &
                 M2(  dof,dof),     &
                 M3(  dof,dof),     &
                 M4(  dof,dof),     &
                 M5(  dof,dof),     &
                 system(dof,dof),   &
                 vl(    dof,dof),   &
                 vr(    dof,dof),   &
                 lambda(dof),       &
                 alpha(dof),        &
                 beta(dof),         &
                 work(10*dof),      &
                 rwork(10*dof),     &
                 stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Allocate Finite Element sub-matrices
        !
        allocate(stencil(    res,res),  &
                 identity(   res,res),  &
                 ridentity(  res,res),  &
                 zero_matrix(res,res),  &
                 stat=ierr)
        if (ierr /= 0) call AllocationError

        M1 = ZERO
        M2 = ZERO
        M3 = ZERO
        M4 = ZERO
        M5 = ZERO
        A  = ZERO
        B  = ZERO
        zero_matrix = ZERO



        !
        ! Construct Fourth-order stencil matrix
        !
        stencil = 0._rk
        stencil(1,  1:5) = [-25._rk, 48._rk, -36._rk, 16._rk, -3._rk]
        stencil(2,  1:5) = [-3._rk, -10._rk,  18._rk, -6._rk,  1._rk]
        stencil(res-1,res-4:res) = [-1._rk,  6._rk, -18._rk,  10._rk,  3._rk]
        stencil(res,  res-4:res) = [3._rk, -16._rk,  36._rk, -48._rk, 25._rk]
        do i = 3,res-2
            stencil(i,i-2) =  1._rk
            stencil(i,i-1) = -8._rk
            stencil(i,i+1) =  8._rk
            stencil(i,i+2) = -1._rk
        end do
        
        ! Scale by dr
        stencil = (1._rk/(12._rk*dr))*stencil



        !
        ! Construct identity matrix for source terms
        !
        identity = 0._rk
        do i = 1,res
            identity(i,i) = 1._rk
        end do



        !
        ! Construct identity scaled by 1/r
        !
        ridentity = 0._rk
        do i = 1,res
            ridentity(i,i) = 1._rk/r(i)
        end do



        !
        ! Circumferential derivative source
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*vbar*real(mod_m,rk))

        irow = 2
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,bc*vbar*real(mod_m,rk))

        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*vbar*real(mod_m,rk))

        irow = 4
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*vbar*real(mod_m,rk))

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*vbar*real(mod_m,rk))

        irow = 1
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*rhobar*real(mod_m,rk))

        irow = 3
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*(real(mod_m,rk)/rhobar))

        irow = 5
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M1(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M1(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,ridentity*gam*pbar*real(mod_m,rk))





        !
        ! Temporal derivative source(real, gets added as imaginary contribution)
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*romega)

        irow = 2
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = identity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,bc*romega)

        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*romega)

        irow = 4
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*romega)

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*romega)



        !
        ! Temporal derivative source(imag, just a small perturbation to facilitate determining 
        ! propagation direction for modes later. Get added as real source term.
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(-identity*iomega,zero_matrix)
                                                                                                                                                                                           
        irow = 2                                                                                                                                                                           
        icol = 2                                                                                                                                                                           
        irow_start = 1 + res*(irow-1)                                                                                                                                                      
        icol_start = 1 + res*(icol-1)                                                                                                                                                      
        bc = identity                                                                                                                                                                      
        bc(:,1)   = ZERO                                                                                                                                                                  
        bc(:,res) = ZERO                                                                                                                                                                  
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(-bc*iomega,zero_matrix)
                                                                                                                                                                                           
        irow = 3                                                                                                                                                                           
        icol = 3                                                                                                                                                                           
        irow_start = 1 + res*(irow-1)                                                                                                                                                      
        icol_start = 1 + res*(icol-1)                                                                                                                                                      
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(-identity*iomega,zero_matrix)
                                                                                                                                                                                           
        irow = 4                                                                                                                                                                           
        icol = 4                                                                                                                                                                           
        irow_start = 1 + res*(irow-1)                                                                                                                                                      
        icol_start = 1 + res*(icol-1)                                                                                                                                                      
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(-identity*iomega,zero_matrix)
                                                                                                                                                                                           
        irow = 5                                                                                                                                                                           
        icol = 5                                                                                                                                                                           
        irow_start = 1 + res*(irow-1)                                                                                                                                                      
        icol_start = 1 + res*(icol-1)                                                                                                                                                      
        M2(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M2(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(-identity*iomega,zero_matrix)






        !
        ! Equation/coordinate system source
        !
        irow = 2
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(ridentity*vbar*vbar/rhobar ,zero_matrix)

        irow = 2
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(ridentity*TWO*vbar ,zero_matrix)

        irow = 2
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(ridentity*(ONE/rhobar) ,zero_matrix)

        irow = 3
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(-ridentity*ubar*vbar/rhobar ,zero_matrix)

        irow = 3
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(-bc*vbar ,zero_matrix)


        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(-ridentity*ubar ,zero_matrix)

        irow = 5
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(bc*(gam-ONE)*rhobar*vbar*vbar ,zero_matrix)

        irow = 5
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(-ridentity*(gam-ONE)*rhobar*ubar*vbar ,zero_matrix)

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M3(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M3(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) - cmplx(-ridentity*(gam-ONE)*ubar ,zero_matrix)





        !
        ! Axial derivative source
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*wbar)

        irow = 2
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = identity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,bc*wbar)

        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*wbar)

        irow = 4
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*wbar)

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*wbar)

        irow = 1
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*rhobar)

        irow = 4
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity/rhobar)

        irow = 5
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        B(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = B(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(zero_matrix,identity*gam*pbar)






        !
        ! Divergence derivative
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(stencil*ubar,zero_matrix)

        irow = 2
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = stencil
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*ubar,zero_matrix)

        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(stencil*ubar,zero_matrix)

        irow = 4
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(stencil*ubar,zero_matrix)

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(stencil*ubar,zero_matrix)

        irow = 1
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = stencil
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*rhobar,zero_matrix)

        irow = 2
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(stencil*(ONE/rhobar),zero_matrix)

        irow = 5
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = stencil
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M4(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M4(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*gam*pbar,zero_matrix)



        !
        ! Divergence source
        !
        irow = 1
        icol = 1
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(ridentity*ubar,zero_matrix)

        irow = 2
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*ubar,zero_matrix)

        irow = 3
        icol = 3
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(ridentity*ubar,zero_matrix)

        irow = 4
        icol = 4
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(ridentity*ubar,zero_matrix)

        irow = 5
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(ridentity*ubar,zero_matrix)

        irow = 1
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*rhobar,zero_matrix)

        irow = 2
        icol = 5
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(ridentity*(ONE/rhobar),zero_matrix)

        irow = 5
        icol = 2
        irow_start = 1 + res*(irow-1)
        icol_start = 1 + res*(icol-1)
        bc = ridentity
        bc(:,1)   = ZERO
        bc(:,res) = ZERO
        M5(irow_start:irow_start + (res-1),icol_start:icol_start + (res-1)) = M5(irow_start:irow_start+(res-1),icol_start:icol_start+(res-1)) + cmplx(bc*gam*pbar,zero_matrix)









        !
        ! Compose A as contributions from separate derivative and source terms
        !
        A = M1 + M2 + M3 + M4 + M5


        !
        ! Compute access/storage indices for removing rows/columns
        !
        iia = 1
        iib = res
        iic = res+2
        iid = 2*res-1
        iie = 2*res+1
        iif = 5*res

        iiia = 1
        iiib = res
        iiic = res+1
        iiid = 2*res-2
        iiie = 2*res-2+1
        iiif = 5*res-2

        ! Remove rows
        allocate(A_reduced(size(A,1)-2,size(A,2)), B_reduced(size(B,1)-2,size(A,2)))
        A_reduced(iiia:iiib,:) = A(iia:iib,:)
        A_reduced(iiic:iiid,:) = A(iic:iid,:)
        A_reduced(iiie:iiif,:) = A(iie:iif,:)

        B_reduced(iiia:iiib,:) = B(iia:iib,:)
        B_reduced(iiic:iiid,:) = B(iic:iid,:)
        B_reduced(iiie:iiif,:) = B(iie:iif,:)



!        iia = 1
!        iib = res
!        iic = res+2
!        iid = 2*res-1
!        iie = 2*res+1
!        iif = 5*res
!
!        iiia = 1
!        iiib = res
!        iiic = res+1
!        iiid = 2*res-2
!        iiie = 2*res-2+1
!        iiif = 5*res-2

        ! Remove columns
        allocate(A_rr(size(A_reduced,1),size(A_reduced,2)-2), B_rr(size(B_reduced,1),size(B_reduced,2)-2))
        A_rr(:,iiia:iiib) = A_reduced(:,iia:iib)
        A_rr(:,iiic:iiid) = A_reduced(:,iic:iid)
        A_rr(:,iiie:iiif) = A_reduced(:,iie:iif)

        B_rr(:,iiia:iiib) = B_reduced(:,iia:iib)
        B_rr(:,iiic:iiid) = B_reduced(:,iic:iid)
        B_rr(:,iiie:iiif) = B_reduced(:,iie:iif)







!        !
!        ! Generalized Eigenvalue Problem
!        !
!        B = -B
!
!        matrix_size = size(A,1)
!        work_size = size(work,1)
!        print*, 'Calling eigenvalue problem: '
!        !call ZGGEV( JOBVL, JOBVR,      N,        A,       LDA,    B,     LDB,     ALPHA, BETA, VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!        call ZGGEV( 'N'   , 'V'  , matrix_size,   A,  matrix_size, B, matrix_size, alpha, beta, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!        print*, 'ZGGEV Info: ', info
!        lambda = alpha/beta
!        do i = 1,size(lambda)
!            print*, alpha(i), beta(i), lambda(i)
!        end do

        !
        ! Generalized Eigenvalue Problem: REDUCED
        !
        B_rr = -B_rr

        deallocate(alpha,beta,vl,vr,work,rwork)
        allocate(alpha(size(A_rr,1)),               &
                 beta(size(A_rr,1)),                &
                 vl(size(A_rr,1),size(A_rr,1)),     &
                 vr(size(A_rr,1),size(A_rr,1)),     &
                 work(10*size(A_rr,1)),             &
                 rwork(10*size(A_rr,1)) )
        matrix_size = size(A_rr,1)
        work_size = size(work,1)
        print*, 'Calling eigenvalue problem: '
        !call ZGGEV( JOBVL, JOBVR,      N,        A,         LDA,      B,       LDB,     ALPHA, BETA, VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
        call ZGGEV( 'N'   , 'V'  , matrix_size,   A_rr,  matrix_size, B_rr, matrix_size, alpha, beta, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
        print*, 'ZGGEV Info: ', info
        lambda = alpha/beta
        do i = 1,size(lambda)
            print*, alpha(i), beta(i), lambda(i)
        end do






!        !
!        ! Standard Eigenvalue Problem
!        !
!        B = -B
!        Binv = zinv(B)
!        system = matmul(Binv,A)
!
!        matrix_size = size(system,1)
!        work_size = size(work,1)
!        print*, 'Calling eigenvalue problem: '
!        !call ZGEEV( JOBVL, JOBVR,      N,      A,          LDA,       W,    VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!        call ZGEEV( 'N'   , 'V'  , matrix_size, system, matrix_size, lambda, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!        print*, 'ZGEEV Info: ', info
!        do i = 1,size(lambda)
!            print*, lambda(i)
!        end do




        !
        ! Write eigenvalues to file for plotting
        !
        print*, 'Writing eigenvalues...'
        inquire(file='eigenvalues.dat', exist=file_exists)
        if (file_exists) call delete_file('eigenvalues.dat')
        open(newunit=handle, file='eigenvalues.dat', iostat=ierr)
        do i = 1,size(lambda)
            !print*, lambda(i)
            write(handle,*) realpart(lambda(i)), ',', imagpart(lambda(i))
        end do
        close(handle)


        !
        ! Compute magnitude of eigenvector
        !
        eigenvector_mag = sqrt(realpart(vr)**TWO  +  imagpart(vr)**TWO)
        !eigenvector_mag = realpart(vr)


        ! Write eigenvectors to file
        print*, 'Writing eigenvectors...'
        inquire(file='eigenvectors.dat', exist=file_exists)
        if (file_exists) call delete_file('eigenvectors.dat')
        open(newunit=handle, file='eigenvectors.dat', iostat=ierr)
        do i = 1,size(r,1)
            write(handle,*) r(i), realpart(vr(res*4 + i - 2,:))     ! pressure(real)    reduced
            !write(handle,*) r(i), realpart(vr(res*4 + i,:))        ! pressure(real)
            !write(handle,*) r(i), eigenvector_mag(res*4 + i,:)     ! pressure(mag)
            !write(handle,*) r(i), eigenvector_mag(res*1 + i,:)     ! radial velocity
        end do
        close(handle)



    end subroutine compute_euler_eigenmodes
    !********************************************************************




















!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   1/23/2018
!    !!
!    !--------------------------------------------------------------------
!    subroutine compute_euler_eigenmodes_old()
!
!        real(rk), parameter :: ri = 1._rk
!        real(rk), parameter :: ro = 2.5_rk
!        integer,  parameter :: nfields = 5
!        integer             :: matrix_size, work_size, handle, info, dof, nterms, res, iterm, inode, idof
!
!        real(rk),         allocatable,   dimension(:)   :: nodes
!        real(rk),         allocatable,   dimension(:,:) :: mass, rmass, stiff, bound1, bound2, zero_matrix, eigenvector_mag, val, dat
!        complex(kind=8),  allocatable,   dimension(:,:) :: vl, vr, A, B, Binv, system
!        complex(kind=8),  allocatable,   dimension(:)   :: work, lambda, alpha, beta
!        double precision, allocatable,   dimension(:)   :: rwork
!        integer(ik),      allocatable,   dimension(:)   :: ind
!        integer(ik) :: ifield, ideriv, ierr, nterms_1d, ind_col, ind_row, irow_small_start, icol_small_start, irow, i, ivec, istart_big, istart_small, icol, icol_start, irow_start
!        logical :: file_exists
!
!
!        !
!        ! Read solution order
!        !
!        print*, 'Enter the solution order of accuracy: '
!        ierr = 1
!        do while (ierr /= 0)
!            read(*,'(I8)', iostat=ierr) nterms
!            if ( (ierr/=0) ) print*, "Invalid input: expecting an integer for the solution order."
!        end do
!
!
!        !
!        ! Compute degrees of freedom
!        !
!        dof = nterms*nfields
!
!
!        !
!        ! Allocate global storage
!        !
!        allocate(                   &
!                 A(     dof,dof),   &
!                 B(     dof,dof),   &
!                 Binv(  dof,dof),   &
!                 system(dof,dof),   &
!                 vl(    dof,dof),   &
!                 vr(    dof,dof),   &
!                 lambda(dof),       &
!                 alpha(dof),        &
!                 beta(dof),         &
!                 work(10*dof),      &
!                 rwork(10*dof),     &
!                 stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        !
!        ! Allocate Finite Element sub-matrices
!        !
!        allocate(mass(  nterms,nterms), &
!                 stiff( nterms,nterms), &
!                 bound1(nterms,nterms), &
!                 bound2(nterms,nterms), &
!                 stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!
!
!
!        !
!        ! Compute Finite Element sub-matrices
!        !
!        mass   = compute_mass_matrix(ri,ro,nterms)
!        rmass  = compute_rmass_matrix(ri,ro,nterms)
!        stiff  = compute_stiffness_matrix(ri,ro,nterms)
!        bound1 = compute_boundary_matrix(nterms,'min')
!        bound2 = compute_boundary_matrix(nterms,'max')
!        zero_matrix = mass
!        zero_matrix = ZERO
!
!        
!        !
!        ! Testing transpose
!        !
!        !stiff = transpose(stiff)
!
!
!
!        print*, 'Mass: '
!        do irow = 1,nterms
!            print*, mass(irow,:)
!        end do
!
!        print*, 'RMass: '
!        do irow = 1,nterms
!            print*, rmass(irow,:)
!        end do
!
!        print*, 'Stiffness: '
!        do irow = 1,nterms
!            print*, stiff(irow,:)
!        end do
!
!
!        print*, 'Boundary 1: '
!        do irow = 1,nterms
!            print*, bound1(irow,:)
!        end do
!
!        print*, 'Boundary 2: '
!        do irow = 1,nterms
!            print*, bound2(irow,:)
!        end do
!
!
!        !
!        ! Zero out storage
!        !
!        A = ZERO
!        B = ZERO
!
!
!
!        !
!        ! Circumferential derivative source
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*vbar*real(mod_m,rk))
!
!        irow = 2
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*vbar*real(mod_m,rk))
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*vbar*real(mod_m,rk))
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*vbar*real(mod_m,rk))
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*vbar*real(mod_m,rk))
!
!        irow = 1
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*rhobar*real(mod_m,rk))
!
!        irow = 3
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*(real(mod_m,rk)/rhobar))
!
!        irow = 5
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,rmass*gam*pbar*real(mod_m,rk))
!
!
!
!
!
!
!        !
!        ! Temporal derivative source
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,-mass*omega)
!
!        irow = 2
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,-mass*omega)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,-mass*omega)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,-mass*omega)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,-mass*omega)
!
!
!
!        !
!        ! Equation/coordinate system source
!        !
!        irow = 2
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(rmass*vbar*vbar/rhobar ,zero_matrix)
!
!        irow = 2
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(rmass*TWO*vbar ,zero_matrix)
!
!        irow = 2
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(rmass*(ONE/rhobar) ,zero_matrix)
!
!        irow = 3
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(-rmass*ubar*vbar/rhobar ,zero_matrix)
!
!        irow = 3
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(-rmass*vbar ,zero_matrix)
!
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(-rmass*ubar ,zero_matrix)
!
!        irow = 5
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(rmass*(gam-ONE)*rhobar*vbar*vbar ,zero_matrix)
!
!        irow = 5
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(-rmass*(gam-ONE)*rhobar*ubar*vbar ,zero_matrix)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(-rmass*(gam-ONE)*ubar ,zero_matrix)
!
!
!
!
!
!        !
!        ! Axial derivative source
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*wbar)
!
!        irow = 2
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*wbar)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*wbar)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*wbar)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*wbar)
!
!        irow = 1
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*rhobar)
!
!        irow = 4
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass/rhobar)
!
!        irow = 5
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        B(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = B(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(zero_matrix,mass*gam*pbar)
!
!
!
!
!
!
!        !
!        ! Divergence volume derivative
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*ubar,zero_matrix)
!
!        irow = 2
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*ubar,zero_matrix)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*ubar,zero_matrix)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*ubar,zero_matrix)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*ubar,zero_matrix)
!
!        irow = 1
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*rhobar,zero_matrix)
!
!        irow = 2
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*(ONE/rhobar),zero_matrix)
!
!        irow = 5
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) - cmplx(stiff*gam*pbar,zero_matrix)
!
!
!
!        !
!        ! Divergence source
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*ubar,zero_matrix)
!
!        irow = 2
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*ubar,zero_matrix)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*ubar,zero_matrix)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*ubar,zero_matrix)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*ubar,zero_matrix)
!
!        irow = 1
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*rhobar,zero_matrix)
!
!        irow = 2
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*(ONE/rhobar),zero_matrix)
!
!        irow = 5
!        icol = 2
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(rmass*gam*pbar,zero_matrix)
!
!
!
!
!        !
!        ! Divergence boundary(face 1)
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*ubar,zero_matrix)
!
!!        irow = 2
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*ubar,zero_matrix)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*ubar,zero_matrix)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*ubar,zero_matrix)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*ubar,zero_matrix)
!
!!        irow = 1
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*rhobar,zero_matrix)
!
!        irow = 2
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*(ONE/rhobar),zero_matrix)
!
!!        irow = 5
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound1*gam*pbar,zero_matrix)
!
!
!
!        !
!        ! Divergence boundary(face 2)
!        !
!        irow = 1
!        icol = 1
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*ubar,zero_matrix)
!
!!        irow = 2
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*ubar,zero_matrix)
!
!        irow = 3
!        icol = 3
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*ubar,zero_matrix)
!
!        irow = 4
!        icol = 4
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*ubar,zero_matrix)
!
!        irow = 5
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*ubar,zero_matrix)
!
!!        irow = 1
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*rhobar,zero_matrix)
!
!        irow = 2
!        icol = 5
!        irow_start = 1 + nterms*(irow-1)
!        icol_start = 1 + nterms*(icol-1)
!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*(ONE/rhobar),zero_matrix)
!
!!        irow = 5
!!        icol = 2
!!        irow_start = 1 + nterms*(irow-1)
!!        icol_start = 1 + nterms*(icol-1)
!!        A(irow_start:irow_start + (nterms-1),icol_start:icol_start + (nterms-1)) = A(irow_start:irow_start+(nterms-1),icol_start:icol_start+(nterms-1)) + cmplx(bound2*gam*pbar,zero_matrix)
!
!
!
!
!
!
!
!
!
!        print*, 'A(real):'
!        do irow = 1,size(A,1)
!            print*, realpart(A(irow,:))
!        end do
!
!        print*, 'A(imag):'
!        do irow = 1,size(A,1)
!            print*, imagpart(A(irow,:))
!        end do
!
!
!
!
!
!        !
!        ! Generalized Eigenvalue Problem
!        !
!        B = -B
!
!        matrix_size = size(A,1)
!        work_size = size(work,1)
!        print*, 'Calling eigenvalue problem: '
!        !call ZGGEV( JOBVL, JOBVR,      N,        A,       LDA,    B,     LDB,     ALPHA, BETA, VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!        call ZGGEV( 'N'   , 'V'  , matrix_size,   A,  matrix_size, B, matrix_size, alpha, beta, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!        print*, 'ZGGEV Info: ', info
!        lambda = alpha/beta
!        do i = 1,size(lambda)
!            print*, alpha(i), beta(i), lambda(i)
!        end do
!
!
!
!
!
!
!
!!        !
!!        ! Standard Eigenvalue Problem
!!        !
!!        B = -B
!!        Binv = zinv(B)
!!        system = matmul(Binv,A)
!!
!!        matrix_size = size(system,1)
!!        work_size = size(work,1)
!!        print*, 'Calling eigenvalue problem: '
!!        !call ZGEEV( JOBVL, JOBVR,      N,      A,          LDA,       W,    VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!!        call ZGEEV( 'N'   , 'V'  , matrix_size, system, matrix_size, lambda, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!!        print*, 'ZGEEV Info: ', info
!!        do i = 1,size(lambda)
!!            print*, lambda(i)
!!        end do
!
!
!
!
!        !
!        ! Write eigenvalues to file for plotting
!        !
!        inquire(file='eigenvalues.dat', exist=file_exists)
!        if (file_exists) call delete_file('eigenvalues.dat')
!        open(newunit=handle, file='eigenvalues.dat', iostat=ierr)
!        do i = 1,size(lambda)
!            !print*, lambda(i)
!            write(handle,*) realpart(lambda(i)), ',', imagpart(lambda(i))
!        end do
!        close(handle)
!
!
!        !
!        ! Compute magnitude of eigenvector
!        !
!        eigenvector_mag = sqrt(realpart(vr)**TWO  +  imagpart(vr)**TWO)
!        !eigenvector_mag = realpart(vr)
!
!
!
!        !
!        ! Write pressure eigenvectors to file for plotting
!        !
!        res = 50
!        nodes = linspace(-ONE,ONE,res)
!        allocate(val(res,nterms), &
!                 dat(res,dof), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        ! Populate interpolation matrix
!        do iterm = 1,nterms
!            do inode = 1,size(nodes)
!                val(inode,iterm) = legendre_val1D(iterm,nodes(inode))
!            end do
!        end do
!
!        print*, 'val: '
!        do irow = 1,size(val,1)
!            print*, val(irow,:)
!        end do
!
!        
!        ! Interpolate pressure eigenvectors
!        do idof = 1,dof
!            dat(:,idof) = matmul(val,eigenvector_mag((dof-(nterms-1)):dof,idof))    ! pressure
!            !dat(:,idof) = matmul(val,eigenvector_mag((nterms+1):(nterms*2),idof))    ! u
!        end do
!        
!        ! Write eigenvectors to file
!        inquire(file='eigenvectors.dat', exist=file_exists)
!        if (file_exists) call delete_file('eigenvectors.dat')
!        open(newunit=handle, file='eigenvectors.dat', iostat=ierr)
!        do i = 1,size(dat,1)
!            write(handle,*) nodes(i), dat(i,:)
!        end do
!        close(handle)
!
!
!
!
!
!        do i = 1,size(vr,2)
!            print*, 'ith eigenvector: ', i
!            print*, vr(:,i)
!        end do
!
!
!!        !
!!        ! Scatter eigenvector to chidg_vector for output
!!        !
!!
!!        !
!!        ! Read solution order
!!        !
!!        print*, 'Enter the eigenvector to store: '
!!        ierr = 1
!!        do while (ierr /= 0)
!!            read(*,'(I8)', iostat=ierr) ivec
!!            if ( (ierr/=0) ) print*, "Invalid input: expecting an integer for the eigenvector index."
!!        end do
!!
!!        do ifield = 1,5
!!            istart_big   = 1 + (ifield-1)*nterms_3d
!!            istart_small = 1 + (ifield-1)*nterms_1d
!!            do i = 1,nterms_1d
!!                print*, 'Storing dof at: ', istart_big + ind(i)-1
!!                chidg%data%sdata%q%dom(1)%vecs(1)%vec(istart_big + ind(i)-1) = sqrt(realpart(vr(istart_small+i-1,ivec))**TWO  +  imagpart(vr(istart_small+i-1,ivec))**TWO)
!!            end do
!!        end do
!
!
!
!
!    end subroutine compute_euler_eigenmodes_old
!    !********************************************************************
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   1/23/2018
!    !!
!    !--------------------------------------------------------------------
!    subroutine compute_euler_eigenmodes_old_old(grid_file)
!        character(*),   intent(in)  :: grid_file
!
!        type(chidg_t)                   :: chidg
!        class(function_t),  allocatable :: fcn
!
!        real(rk),    allocatable,   dimension(:,:)  :: A1, A2a, A2b, A3a, A3b, B, C, D, E, F, zero_big, zero_small
!        !complex(rk), allocatable,   dimension(:,:)  :: Mc, Nc
!        !complex*16,  allocatable,   dimension(:,:)  :: vl, vr, Mbig, Nbig, Msmall, Nsmall
!        !complex*16,  allocatable,   dimension(:)    :: alpha, beta, work, lambda
!        complex(kind=8),  allocatable,   dimension(:,:)  :: vl, vr, Mbig, Nbig, Msmall, Nsmall, Ninv, A
!        complex(kind=8),  allocatable,   dimension(:)    :: alpha, beta, work, lambda
!        double precision, allocatable, dimension(:) :: rwork
!        integer(ik), allocatable,   dimension(:)    :: ind
!        integer(ik) :: ifield, ideriv, ierr, nterms_1d, nterms_3d, ind_col, ind_row, irow_small_start, irow_big_start, icol_small_start, icol_big_start, irow, info, i, ivec, istart_big, istart_small
!        integer :: matrix_size, work_size, handle
!        logical :: file_exists
!
!
!        !
!        ! Read solution order
!        !
!        print*, 'Enter the solution order of accuracy: '
!        ierr = 1
!        do while (ierr /= 0)
!            read(*,'(I8)', iostat=ierr) solution_order
!            if ( (ierr/=0) ) print*, "Invalid input: expecting an integer for the solution order."
!        end do
!
!
!
!
!
!        !
!        ! Initialize ChiDG environment
!        !
!        call chidg%start_up('mpi')
!!        call chidg%start_up('namelist')
!        call chidg%start_up('core')
!
!
!        !
!        ! Set ChiDG Algorithms, Accuracy
!        !
!        call chidg%set('Time Integrator' , algorithm=time_integrator                   )
!!        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions)
!!        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions)
!!        call chidg%set('Preconditioner'  , algorithm=preconditioner                    )
!        call chidg%set('Solution Order'  , integer_input=solution_order)
!
!
!
!        !
!        ! Read grid and boundary condition data
!        !
!        call chidg%read_mesh(grid_file)
!
!
!
!        !
!        ! Initialize solution
!        !   1: 'none', init fields with values from mod_io module variable initial_fields(:)
!        !   2: read initial solution from ChiDG hdf5 file
!        !
!        initial_fields = ZERO
!        call create_function(fcn,'constant')
!        do ifield = 1,chidg%data%mesh%domain(1)%neqns
!            call fcn%set_option('val',initial_fields(ifield))
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,fcn,ifield)
!        end do
!
!
!
!        !
!        ! Initialize algorithms and state
!        !
!        call chidg%time_integrator%initialize_state(chidg%data)
!
!
!
!        !
!        ! Update spatial operators
!        !
!        chidg%data%time_manager%itime = 1
!        chidg%data%time_manager%t     = ZERO
!        call update_space(chidg%data,differentiate=.true.)
!        
!
!
!
!        !
!        ! Get matrix contributions
!        !
!        call chidg%data%sdata%A1%to_real(A1)        ! Volume advection
!        call chidg%data%sdata%A2a%to_real(A2a)      ! Boundary Average
!        call chidg%data%sdata%A2b%to_real(A2b)      ! Boundary Dissipation
!        call chidg%data%sdata%A3a%to_real(A3a)      ! Boundary condition(face1)
!        call chidg%data%sdata%A3b%to_real(A3b)      ! Boundary condition(face2)
!        call chidg%data%sdata%B%to_real(B)          ! Circumferential derivative
!        call chidg%data%sdata%C%to_real(C)          ! Axial derivative
!        call chidg%data%sdata%D%to_real(D)          ! Equation/coordinate source terms
!        call chidg%data%sdata%E%to_real(E)          ! Temporal derivative
!        call chidg%data%sdata%F%to_real(F)          ! Divergence source
!
!
!
!        allocate(zero_big(size(A1,1),size(A1,2)), stat=ierr)
!        if (ierr /= 0) call AllocationError
!        zero_big = ZERO
!
!
!
!!        print*, 'A3a:'
!!        do irow = 1,size(A3a,1)
!!            print*, A3a(irow,:)
!!        end do
!!
!!        print*, 'A3b:'
!!        do irow = 1,size(A3b,1)
!!            print*, A3b(irow,:)
!!        end do
!
!
!
!
!
!
!        !
!        ! Construct left and right matrices for Generalized Eigenvalue Problem
!        !
!        Mbig = cmplx(A1 + A3a + A3b + D - F,-B+E)
!        Nbig = cmplx(zero_big,C)
!
!        !
!        ! Allocate reduced matrices
!        !
!        allocate(Msmall(solution_order*5,solution_order*5), &
!                 Nsmall(solution_order*5,solution_order*5), &
!                 vl(solution_order*5,solution_order*5),     &
!                 vr(solution_order*5,solution_order*5),     &
!                 work(10*solution_order*5),                  &
!                 rwork(10*solution_order*5),                 &
!                 stat=ierr)
!        if (ierr /= 0) call AllocationError
!        
!
!        !
!        ! Compute number of terms in expansions
!        !
!        nterms_1d = solution_order
!        nterms_3d = solution_order*solution_order*solution_order
!
!
!
!        !
!        ! Reorder
!        !
!        ind = [1,4,11,28,67,127,217,344]
!        do ideriv = 1,5
!            do ifield = 1,5
!
!                irow_small_start = 1 + nterms_1d*(ifield-1)
!                icol_small_start = 1 + nterms_1d*(ideriv-1)
!
!                irow_big_start = 1 + nterms_3d*(ifield-1)
!                icol_big_start = 1 + nterms_3d*(ideriv-1)
!                
!                ! Assign entries from big to small
!                do ind_col = 1,nterms_1d
!                    do ind_row = 1,nterms_1d
!                        Msmall(irow_small_start+ind_row-1, icol_small_start+ind_col-1) = Mbig(irow_big_start+ind(ind_row)-1, icol_big_start+ind(ind_col)-1)
!                        Nsmall(irow_small_start+ind_row-1, icol_small_start+ind_col-1) = Nbig(irow_big_start+ind(ind_row)-1, icol_big_start+ind(ind_col)-1)
!                    end do
!                end do
!
!                
!            end do
!        end do
!
!    
!
!!        !
!!        ! Print
!!        !
!!        do ideriv = 1,5
!!            do ifield = 1,5
!!
!!                irow_big_start = 1 + nterms_3d*(ifield-1)
!!                icol_big_start = 1 + nterms_3d*(ideriv-1)
!!                
!!                print*, ifield, ideriv
!!                do ind_row = 1,nterms_3d
!!                    print*, Mbig(irow_big_start+(ind_row-1),icol_big_start:icol_big_start+(nterms_3d-1))
!!                end do
!!                
!!            end do
!!        end do
!
!        print*, 'real(Msmall):'
!        do irow = 1,size(Msmall,1)
!            print*, realpart(Msmall(irow,:))
!        end do
!
!        print*, 'imag(Msmall):'
!        do irow = 1,size(Msmall,1)
!            print*, imagpart(Msmall(irow,:))
!        end do
!
!        print*, 'real(Nsmall):'
!        do irow = 1,size(Nsmall,1)
!            print*, realpart(Nsmall(irow,:))
!        end do
!
!        print*, 'imag(Nsmall):'
!        do irow = 1,size(Nsmall,1)
!            print*, imagpart(Nsmall(irow,:))
!        end do
!
!
!
!        allocate(alpha(size(Msmall,1)), beta(size(Nsmall,1)), lambda(size(Msmall,1)), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!!        !
!!        ! Generalized Eigenvalue Problem
!!        !
!!        matrix_size = size(Msmall,1)
!!        work_size = size(work,1)
!!        print*, 'Calling eigenvalue problem: '
!!        !call ZGGEV( JOBVL, JOBVR,      N,        A,        LDA,       B,        LDB,     ALPHA, BETA, VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!!        call ZGGEV( 'N'   , 'V'  , matrix_size, Msmall, matrix_size, Nsmall, matrix_size, alpha, beta, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!!        print*, 'ZGGEV Info: ', info
!!        lambda = alpha/beta
!!        do i = 1,size(lambda)
!!            print*, alpha(i), beta(i), lambda(i)
!!        end do
!
!
!        !
!        ! Standard Eigenvalue Problem
!        !
!        Ninv = zinv(Nsmall)
!        A = matmul(Ninv,Msmall)
!
!        matrix_size = size(Msmall,1)
!        work_size = size(work,1)
!        print*, 'Calling eigenvalue problem: '
!        !call ZGEEV( JOBVL, JOBVR,      N,      A,        LDA,       W,    VL,     LDVL,    VR,     LDVR,    WORK,   LWORK,   RWORK, INFO )
!        call ZGEEV( 'N'   , 'V'  , matrix_size, A, matrix_size, lambda, vl, matrix_size, vr, matrix_size, work, work_size, rwork, info )
!        print*, 'ZGEEV Info: ', info
!        do i = 1,size(lambda)
!            print*, lambda(i)
!        end do
!
!
!
!
!        inquire(file='eigenvalues.dat', exist=file_exists)
!        if (file_exists) call delete_file('eigenvalues.dat')
!        open(newunit=handle, file='eigenvalues.dat', iostat=ierr)
!        do i = 1,size(lambda)
!            !print*, lambda(i)
!            write(handle,*) realpart(lambda(i)), ',', imagpart(lambda(i))
!        end do
!        close(handle)
!
!
!        do i = 1,size(vr,2)
!            print*, 'ith eigenvector: ', i
!            print*, vr(:,i)
!        end do
!
!
!        !
!        ! Scatter eigenvector to chidg_vector for output
!        !
!
!        !
!        ! Read solution order
!        !
!        print*, 'Enter the eigenvector to store: '
!        ierr = 1
!        do while (ierr /= 0)
!            read(*,'(I8)', iostat=ierr) ivec
!            if ( (ierr/=0) ) print*, "Invalid input: expecting an integer for the eigenvector index."
!        end do
!
!        do ifield = 1,5
!            istart_big   = 1 + (ifield-1)*nterms_3d
!            istart_small = 1 + (ifield-1)*nterms_1d
!            do i = 1,nterms_1d
!                print*, 'Storing dof at: ', istart_big + ind(i)-1
!                chidg%data%sdata%q%dom(1)%vecs(1)%vec(istart_big + ind(i)-1) = sqrt(realpart(vr(istart_small+i-1,ivec))**TWO  +  imagpart(vr(istart_small+i-1,ivec))**TWO)
!            end do
!        end do
!
!
!        call chidg%write_mesh('eigenvector.h5')
!        call chidg%write_fields('eigenvector.h5')
!
!
!
!        !
!        ! Close ChiDG
!        !
!        call chidg%shut_down('core')
!        call chidg%shut_down('mpi')
!
!
!
!!      integer, parameter :: N=4, nb=64, Nmax=10
!!      integer :: lda,ldb,ldvr,lwork
!!      parameter (lda=Nmax, ldb=Nmax, ldvr=Nmax,lwork=Nmax+Nmax*nb)
!!      integer :: i,j,info
!!      !complex(kind=8) :: A(lda,Nmax), alpha(Nmax), B(ldb,Nmax), beta(Nmax), dummy(1,1), vr(ldvr,Nmax), work(lwork), eig(Nmax)
!!      !complex*16:: A(lda,Nmax), alpha(Nmax), B(ldb,Nmax), beta(Nmax), dummy(1,1), vr(ldvr,Nmax), work(lwork), eig(Nmax)
!!      complex(kind=rk) :: A(lda,Nmax), alpha(Nmax), B(ldb,Nmax), beta(Nmax), dummy(1,1), vr(ldvr,Nmax), work(lwork), eig(Nmax)
!!      double precision :: rwork(8*Nmax)
!!
!!
!!      A(1,1)=(-21.10,-22.50);A(1,2)=(53.50,-50.50)
!!      A(1,3)=(-34.50,127.50);A(1,4)=(7.50,0.50)
!!      A(2,1)=(-0.46,-7.78);A(2,2)=(-3.5,-37.5)
!!      A(2,3)=(-15.5,58.5);A(2,4)=(-10.5,-1.5)
!!      A(3,1)=(4.3,-5.5);A(3,2)=(39.7,-17.1)
!!      A(3,3)=(-68.5,12.5);A(3,4)=(-7.5,-3.5)
!!      A(4,1)=(5.5,4.4);A(4,2)=(14.4,43.3)
!!      A(4,3)=(-32.5,-46);A(4,4)=(-19,-32.5)
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      B(1,1)=(1,-5);B(1,2)=(1.6,1.2)
!!      B(1,3)=(-3,0);B(1,4)=(0,-1)
!!      B(2,1)=(0.8,-0.6);B(2,2)=(3,-5)
!!      B(2,3)=(-4,3);B(2,4)=(-2.4,-3.2)
!!      B(3,1)=(1,0);B(3,2)=(2.4,1.8)
!!      B(3,3)=(-4,-5);B(3,4)=(0,-3)
!!      B(4,1)=(0,1);B(4,2)=(-1.8,2.4)
!!      B(4,3)=(0,-4);B(4,4)=(4,-5)
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!      call zggev('n','v',N,A,lda,B,ldb,alpha,beta,dummy,1,vr,ldvr,work,lwork,rwork,info)
!!
!!      eig=alpha/beta
!!
!!      !here skip the heading to the file
!!
!!
!!
!!      do j=1,N   
!!         print*, eig(j)
!!      end do 
!
!
!
!
!    end subroutine compute_euler_eigenmodes_old_old
!    !********************************************************************








    !>  Compute the mass matrix for the 1D DG discretization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/26/2018
    !!
    !------------------------------------------------------
    function compute_mass_matrix(ri,ro,nterms) result(mass)
        real(rk),       intent(in)  :: ri
        real(rk),       intent(in)  :: ro
        integer(ik),    intent(in)  :: nterms

        real(rk) :: jinv, nodes(nterms), weights(nterms), mass(nterms,nterms), temp(nterms,nterms), val(nterms,nterms)
        integer(ik) :: iterm, inode

        !
        ! Compute Gauss-Legendre Quadrature nodes/weights
        !
        nodes   = gl_nodes(nterms)
        weights = gl_weights(nterms)

        !
        ! Compute metric terms
        !
        jinv = (ro-ri)/TWO


        !
        ! Compute polynomial evaluations at nodes
        !
        do iterm = 1,nterms
            do inode = 1,nterms
                val(inode,iterm) = legendre_val1D(iterm,nodes(inode))
            end do
        end do
        

        !
        ! Pre-multiply quadrature weights and jacobian scaling
        !
        do iterm = 1,nterms
            temp(:,iterm) = val(:,iterm)*weights*jinv
        end do

        
        !
        ! Compute mass matrix as matrix-vector multiplication
        !
        mass = matmul(transpose(val),temp)

    end function compute_mass_matrix
    !*********************************************************************





    !>  Compute the mass matrix with 1/r scaling for the 1D DG discretization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/26/2018
    !!
    !-------------------------------------------------------------------------
    function compute_rmass_matrix(ri,ro,nterms) result(rmass)
        real(rk),       intent(in)  :: ri
        real(rk),       intent(in)  :: ro
        integer(ik),    intent(in)  :: nterms

        real(rk),   allocatable :: temp(:,:), val(:,:), r(:), nodes(:), weights(:)
        real(rk)                :: jinv, rmass(nterms,nterms)
        integer(ik)             :: iterm, inode, nnodes, ierr

        !
        ! Compute Gauss-Legendre Quadrature nodes/weights
        !
        !nnodes = nterms*5   ! this function is nonlinear, so we need more points.
        nnodes = nterms*2   ! this function is nonlinear, so we need more points.
        nodes   = gl_nodes(nnodes)
        weights = gl_weights(nnodes)
        r       = ((ro-ri)/TWO)*nodes + (ro+ri)/TWO


        !
        ! Compute metric terms
        !
        jinv = (ro-ri)/TWO


        !
        ! Compute polynomial evaluations at nodes
        !
        allocate(val(nnodes,nterms), stat=ierr)
        do iterm = 1,nterms
            do inode = 1,nnodes
                val(inode,iterm) = legendre_val1D(iterm,nodes(inode))
            end do
        end do
        

        !
        ! Pre-multiply quadrature weights and jacobian scaling
        !
        allocate(temp(nnodes,nterms), stat=ierr)
        do iterm = 1,nterms
            temp(:,iterm) = val(:,iterm)*weights*(ONE/r)*jinv
        end do

        
        !
        ! Compute mass matrix as with 1/r scaling as matrix-vector multiplication
        !
        rmass = matmul(transpose(val),temp)

    end function compute_rmass_matrix
    !*********************************************************************






    !>  Compute the mass matrix for the 1D DG discretization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------
    function compute_boundary_matrix(nterms,side) result(bnd)
        integer(ik),    intent(in)  :: nterms
        character(*),   intent(in)  :: side

        real(rk) :: node, bnd(nterms,nterms), temp(nterms,nterms), val(nterms,1)
        integer(ik) :: iterm, inode

        !
        ! Compute Gauss-Legendre Quadrature nodes/weights
        !
        if (trim(side) == 'min') then
            node = -ONE
        else if (trim(side) == 'max') then
            node = ONE
        else
            stop
        end if


        !
        ! Compute polynomial evaluations at boundary
        !
        do iterm = 1,nterms
            val(iterm,1) = legendre_val1D(iterm,node)
        end do
        
        
        !
        ! Compute boundary matrix as matrix-vector multiplication
        !
        bnd = matmul(val,transpose(val))


        !
        ! Negate if normal vector is -1
        !
        if (trim(side) == 'min') then
            bnd = -bnd
        end if

    end function compute_boundary_matrix
    !************************************************************************

















    !>  Compute the stiffness matrix for the 1D DG discretization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/26/2018
    !!
    !------------------------------------------------------------------------
    function compute_stiffness_matrix(ri,ro,nterms) result(stiff)
        real(rk),       intent(in)  :: ri
        real(rk),       intent(in)  :: ro
        integer(ik),    intent(in)  :: nterms

        real(rk)    :: nodes(nterms), weights(nterms), stiff(nterms,nterms), temp(nterms,nterms), val(nterms,nterms), dval(nterms,nterms)
        integer(ik) :: iterm, inode

        !
        ! Compute Gauss-Legendre Quadrature nodes/weights
        !
        nodes   = gl_nodes(nterms)
        weights = gl_weights(nterms)


        !
        ! Compute polynomial evaluations at nodes
        !
        do iterm = 1,nterms
            do inode = 1,nterms
                val(inode,iterm) = legendre_val1D(iterm,nodes(inode))
            end do
        end do


        !
        ! Compute polynomial derivatives at nodes
        !
        do iterm = 1,nterms
            do inode = 1,nterms
                dval(inode,iterm) = dlegendre_val1D(iterm,nodes(inode))
            end do
        end do


        !
        ! Pre-multiply quadrature weights and scaling
        !
        do iterm = 1,nterms
            temp(:,iterm) = val(:,iterm)*weights
        end do

        
        !
        ! Compute stiffness matrix as matrix-vector multiplication
        !
        stiff = matmul(transpose(dval),temp)

    end function compute_stiffness_matrix
    !**************************************************************************











!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   1/29/2017
!    !!
!    !!
!    !-------------------------------------------------------------------------
!    function filter(r,eigenvalues, eigenvectors)
!        real(rk),           intent(in)  :: r(:)
!        complex(kind=8),    intent(in)  :: eigenvalues(:)
!        complex(kind=8),    intent(in)  :: eigenvectors(:,:)
!
!
!        do ivec = 1,size(eigenvectors,2)
!
!            mag = realpart(eigenvectors(:,ivec))**TWO  +  imagpart(eigenvectors(:,ivec)
!
!        end do
!
!
!
!
!    end function filter
!    !*************************************************************************











end module mod_euler_eigenmodes
