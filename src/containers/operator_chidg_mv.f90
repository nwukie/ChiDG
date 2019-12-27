module operator_chidg_mv
#include <messenger.h>
#include "petsc/finclude/petscmat.h"
    use petscmat,           only: MatMult, MatMultTranspose, MatNorm, NORM_INFINITY

    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO, ONE
    use mod_chidg_mpi,      only: IRANK, ChiDG_COMM
    use type_chidg_matrix,  only: chidg_matrix_t
    use mod_time,           only: time_manager_global
    use mpi_f08,            only: MPI_BCast, MPI_Send, MPI_Recv, MPI_INTEGER4,MPI_REAL8, &
                                  MPI_LOGICAL, MPI_STATUS_IGNORE
    use type_chidg_vector
    use type_timer,         only: timer_t
    implicit none

    external DGEMV

    type(timer_t)   :: timer_comm, timer_blas


contains




    !>  Select whether the matrix-vector or the transposed-matrix-vector product needs to be 
    !!  carried out. 
    !!
    !!  The transposed matrix is needed for adjoint
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/6/2018
    !!
    !!  @param[in]     A                Matrix
    !!  @param[in]     x                Vector
    !!  @param[in]     resvec_init      Vector for initializing the result vector.
    !!                                  This is necessary if the result vector
    !!                                  has a structure different from vector x.
    !!                                  This is needed for grid-node sensitivities.
    !!
    !!  @param[out]    res_vec          Resultant chidg vector with the shape of 
    !!                                  resvec_init (if present) or x
    !!
    !------------------------------------------------------------------------------------
    function chidg_mv(A,x,resvec_init,pattern) result(res_vec)
        type(chidg_matrix_t),   intent(inout)               :: A
        type(chidg_vector_t),   intent(inout)               :: x
        type(chidg_vector_t),   intent(inout),  optional    :: resvec_init
        character(*),           intent(in),     optional    :: pattern
        
        type(chidg_vector_t)        :: res_vec

        ! Select the specialized mv function
        if (.not. A%transposed) then
            res_vec = matrix_vector_product(A,x,pattern=pattern)
        else if (A%transposed) then
            res_vec = transposed_matrix_vector_product(A,x,resvec_init)
        end if

    end function chidg_mv
    !****************************************************************************************








    !>  This function implements the important matrix-vector multiplication 
    !!  operation : A*x :  for chidg_matrix * chidg_vector operations.
    !!
    !!  The structure of the operation is as follows:
    !!      A (in-out)
    !!      x (in-out)
    !!
    !!      Communication:
    !!      ----------------------------------
    !!          : x (initiate non-blocking send)
    !!
    !!      Processor local portion:
    !!      ----------------------------------
    !!          : A*x (INTERIOR - lblks   )
    !!          : A*x (CHIMERA  - chi_blks)
    !!          : A*x (BOUNDARY - bc_blks )
    !!
    !!      Matrix-free contributions:
    !!      ----------------------------------
    !!          : A*x (Harmonic Balance)
    !!
    !!      Communication:
    !!      ----------------------------------
    !!          : x (initiate blocking recv)
    !!
    !!      Processor global(parallel) portion
    !!      ----------------------------------
    !!          : A*x (INTERIOR - lblks   )
    !!          : A*x (CHIMERA  - chi_blks)
    !!          : A*x (BOUNDARY - bc_blks )
    !!
    !!      Communication:
    !!      ----------------------------------
    !!          : x (wait on non-blocking send requests)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!  @date   6/6/2016 (AFRL) - parallelization
    !!  @date   4/13/2017 - boundary coupling
    !!
    !!  @author Mayank Sharma
    !!  @date   2/23/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    function matrix_vector_product(A,x,pattern) result(res)
        type(chidg_matrix_t),   intent(inout)           :: A
        type(chidg_vector_t),   intent(inout)           :: x
        character(*),           intent(in), optional    :: pattern

        type(chidg_vector_t)        :: res
        integer(ik)                 :: idom, ielem, iblk, imat, itime, ivar,  &
                                       itime_i, recv_comm, recv_domain, recv_element, idiag
        integer(ik)                 :: dparent_g, dparent_l, eparent_g, eparent_l
        integer(ik)                 :: matrix_proc, vector_proc, nrows, ncols, ierr
        integer(ik)                 :: res_istart, res_iend, x_istart, x_iend, xstart, xend, tparent
        logical                     :: local_multiply, parallel_multiply
        logical                     :: nonconforming = .false.
        logical                     :: HB_flag_status = .false.
        logical                     :: pattern_multiply
        character(:), allocatable   :: HB_flag, multiply_pattern
        real(rk),     allocatable   :: D(:,:)
        real(rk),     allocatable   :: temp_1(:), temp_2(:)

        PetscErrorCode :: perr
        PetscReal :: petsc_norm


        if (allocated(x%wrapped_petsc_vector)) then

            call timer_comm%start()
            call timer_comm%stop()
            call timer_blas%start()

            res = x

            call MatMult(A%wrapped_petsc_matrix%petsc_matrix,x%wrapped_petsc_vector%petsc_vector,res%wrapped_petsc_vector%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_mv: error calling petsc MatMult.')

            call timer_blas%stop()

        else



            if (present(pattern)) then
                multiply_pattern = pattern
            else
                multiply_pattern = 'all'
            end if


            !
            ! Allocate result and clear
            !
            res = x
            call res%clear


            !
            ! Check to see if matrix has been initialized with information about where to locate 
            ! vector information being received from other processors.
            !
            if ( .not. A%recv_initialized ) then
                call A%init_recv(x)
            end if


            !
            ! Begin non-blocking send of parallel vector information
            !
            call timer_comm%start()
            call x%comm_send()
            call timer_comm%stop()


            do itime = 1,time_manager_global%ntime


                !
                ! Compute A*x for local matrix-vector product
                !
                do idom = 1,size(A%dom)


                    !
                    ! Routine for proc-local INTERIOR coupling (lblks)
                    !
                    do ielem = 1,size(A%dom(idom)%lblks,1)
                        idiag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                        do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()

                        
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                            local_multiply   = ( matrix_proc == vector_proc )
                            pattern_multiply = (multiply_pattern == 'all') .or. (multiply_pattern == 'off-diagonal' .and. (imat /= idiag))
            
                            !if ( local_multiply .and. pattern_multiply) then
                            if ( local_multiply ) then
                                dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                                eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)

                                res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                x_istart = x%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                x_iend   = x%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),   &
                                            xvec   => x%dom(dparent_l)%vecs(eparent_l)%vec(x_istart:x_iend),     &
                                            Amat   => A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat )

                                    nonconforming = ( size(Amat,2) /= size(xvec) )
                                    if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Local m-v operation")

                                    call timer_blas%start()
                                    resvec = resvec + matmul(Amat,xvec)
                                    call timer_blas%stop()

                                end associate

                            end if

                        end do !imat
                    end do !ielem







                    !
                    ! Routine for proc-local CHIMERA coupling (chi_blks)
                    !
                    if (allocated(A%dom(idom)%chi_blks)) then
                        do ielem = 1,size(A%dom(idom)%chi_blks,1)
                            do imat = 1,A%dom(idom)%chi_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)

                                local_multiply    = ( matrix_proc == vector_proc )

                                if ( local_multiply ) then
                                    dparent_l = A%dom(idom)%chi_blks(ielem,itime)%dparent_l(imat)
                                    eparent_l = A%dom(idom)%chi_blks(ielem,itime)%eparent_l(imat)

                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                    x_iend     = x%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),    &
                                                xvec   => x%dom(dparent_l)%vecs(eparent_l)%vec(x_istart:x_iend), &
                                                Amat   => A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%mat  ) 

                                        !
                                        ! Test matrix vector sizes
                                        !
                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if






                            end do !imat
                        end do ! ielem
                    end if  ! allocated



                    !
                    ! Routine for proc-local BOUNDARY coupling (bc_blks)
                    !
                    if (allocated(A%dom(idom)%bc_blks)) then
                        do ielem = 1,size(A%dom(idom)%bc_blks,1)
                            do imat = 1,A%dom(idom)%bc_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)

                                local_multiply    = ( matrix_proc == vector_proc )

                                if ( local_multiply ) then
                                    dparent_l = A%dom(idom)%bc_blks(ielem,itime)%dparent_l(imat)
                                    eparent_l = A%dom(idom)%bc_blks(ielem,itime)%eparent_l(imat)

                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                    x_iend     = x%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),    &
                                                xvec   => x%dom(dparent_l)%vecs(eparent_l)%vec(x_istart:x_iend), &
                                                Amat   => A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%mat ) 

                                        !
                                        ! Test matrix vector sizes
                                        !
                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Boundary m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated


                    !
                    ! Routine for proc-local HARMONIC BALANCE coupling (hb_blks)
                    !
                    if (allocated(A%dom(idom)%hb_blks)) then
                        do ielem = 1,size(A%dom(idom)%hb_blks,1)
                            do imat = 1,A%dom(idom)%hb_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%hb_blks(ielem,itime)%parent_proc(imat)

                                local_multiply    = ( matrix_proc == vector_proc )

                                if ( local_multiply ) then
                                    dparent_l = A%dom(idom)%hb_blks(ielem,itime)%dparent_l(imat)
                                    eparent_l = A%dom(idom)%hb_blks(ielem,itime)%eparent_l(imat)
                                    tparent   = A%dom(idom)%hb_blks(ielem,itime)%tparent(imat)

                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%dom(dparent_l)%vecs(eparent_l)%get_time_start(tparent)
                                    x_iend     = x%dom(dparent_l)%vecs(eparent_l)%get_time_end(tparent)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),    &
                                                xvec   => x%dom(dparent_l)%vecs(eparent_l)%vec(x_istart:x_iend), &
                                                Amat   => A%dom(idom)%hb_blks(ielem,itime)%data_(imat)%mat ) 

                                        ! Test matrix vector sizes
                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Harmonic Balance m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated






                    
    !                !
    !                ! Routine for harmonic balance computations
    !                ! Used only when harmonic balance is specified
    !                !
    !                HB_flag = time_manager_global%get_name()
    !
    !                if (HB_flag == 'Harmonic Balance' .or. HB_flag == 'Harmonic_Balance' .or. HB_flag == 'harmonic balance' &
    !                    .or. HB_flag == 'harmonic_balance' .or. HB_flag == 'HB') then
    !
    !                    D = time_manager_global%D
    !
    !                    do ielem = 1,size(A%dom(idom)%lblks,1) 
    !
    !                        imat = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
    !
    !                        matrix_proc = IRANK
    !                        vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)
    !
    !                        local_multiply    = (matrix_proc == vector_proc)
    !                        parallel_multiply = (matrix_proc /= vector_proc)
    !
    !                        if (local_multiply) then
    !                            do itime_i = 1,size(A%dom(idom)%lblks,2)
    !                                if (itime_i /= itime) then
    !
    !                                    associate ( nvars  => x%dom(idom)%vecs(ielem)%nvars(), &
    !                                                nterms => x%dom(idom)%vecs(ielem)%nterms(), &
    !                                                mass   => A%dom(idom)%lblks(ielem,itime_i)%mass )
    !
    !                                    if (allocated(temp_1)) deallocate(temp_1)
    !                                    if (allocated(temp_2)) deallocate(temp_2)
    !                                    allocate(temp_1(nterms),temp_2(nterms), stat=ierr)
    !                                    if (ierr /= 0) call AllocationError
    !
    !                                    call timer_blas%start()
    !                                    do ivar = 1,nvars
    !
    !                                        temp_1 = D(itime,itime_i)*matmul(mass,x%dom(idom)%vecs(ielem)%getvar(ivar,itime_i))
    !                                        temp_2 = res%dom(idom)%vecs(ielem)%getvar(ivar,itime) + temp_1
    !                                        call res%dom(idom)%vecs(ielem)%setvar(ivar,itime,temp_2)
    !
    !                                    end do  ! ivar
    !                                    call timer_blas%stop()
    !
    !                                    end associate
    !
    !                                end if 
    !                            end do  ! itime_i
    !                        end if  ! local_multiply
    !
    !                    end do  ! ielem
    !
    !                end if  ! HB_flag



                end do ! idom

            end do ! itime


        


            !
            ! Begin blocking recv of parallel vector information
            !
            call timer_comm%start()
            call x%comm_recv()
            call timer_comm%stop()




            do itime = 1,time_manager_global%ntime

                !
                ! Compute A*x for parallel matrix-vector product
                !
                do idom = 1,size(A%dom)

                    !
                    ! Routine for global(parallel) INTERIOR coupling (lblks)
                    !
                    do ielem = 1,size(A%dom(idom)%lblks,1)
                        do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()
                        
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                            parallel_multiply = ( matrix_proc /= vector_proc )
            
                            if ( parallel_multiply ) then
                                recv_comm    = A%dom(idom)%lblks(ielem,itime)%get_recv_comm(imat)
                                recv_domain  = A%dom(idom)%lblks(ielem,itime)%get_recv_domain(imat)
                                recv_element = A%dom(idom)%lblks(ielem,itime)%get_recv_element(imat)

                                res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(itime)
                                x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(itime)
                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),                               &
                                            xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec(x_istart:x_iend),  &
                                            Amat   => A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat )

                                    call timer_blas%start()
                                    resvec = resvec + matmul(Amat,xvec)
                                    call timer_blas%stop()

                                end associate
                            end if

                        end do !imat
                    end do !ielem



                    !
                    ! Routine for global(parallel) CHIMERA coupling (chi_blks)
                    !
                    if (allocated(A%dom(idom)%chi_blks)) then
                        do ielem = 1,size(A%dom(idom)%chi_blks,1)
                            do imat = 1,A%dom(idom)%chi_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)
                                parallel_multiply = ( matrix_proc /= vector_proc )

                                if ( parallel_multiply ) then
                                    recv_comm    = A%dom(idom)%chi_blks(ielem,itime)%get_recv_comm(imat)
                                    recv_domain  = A%dom(idom)%chi_blks(ielem,itime)%get_recv_domain(imat)
                                    recv_element = A%dom(idom)%chi_blks(ielem,itime)%get_recv_element(imat)


                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(itime)
                                    x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(itime)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),                                &
                                                xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec(x_istart:x_iend),   &
                                                Amat   => A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%mat )


                                        !
                                        ! Test matrix vector sizes
                                        !
                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated




                    !
                    ! Routine for global(parallel) BOUNDARY coupling (bc_blks)
                    !
                    if (allocated(A%dom(idom)%bc_blks)) then
                        do ielem = 1,size(A%dom(idom)%bc_blks,1)
                            do imat = 1,A%dom(idom)%bc_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)

                                parallel_multiply = ( matrix_proc /= vector_proc )

                                if ( parallel_multiply ) then

                                    recv_comm    = A%dom(idom)%bc_blks(ielem,itime)%get_recv_comm(imat)
                                    recv_domain  = A%dom(idom)%bc_blks(ielem,itime)%get_recv_domain(imat)
                                    recv_element = A%dom(idom)%bc_blks(ielem,itime)%get_recv_element(imat)


                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(itime)
                                    x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(itime)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),                               &
                                                xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec(x_istart:x_iend),  &
                                                Amat   => A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%mat )

                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Boundary m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated


                    !
                    ! Routine for global(parallel) HARMONIC BALANCE coupling (hb_blks)
                    !
                    if (allocated(A%dom(idom)%hb_blks)) then
                        do ielem = 1,size(A%dom(idom)%hb_blks,1)
                            do imat = 1,A%dom(idom)%hb_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%hb_blks(ielem,itime)%parent_proc(imat)
                                parallel_multiply = ( matrix_proc /= vector_proc )

                                if ( parallel_multiply ) then
                                    tparent      = A%dom(idom)%hb_blks(ielem,itime)%tparent(imat)
                                    recv_comm    = A%dom(idom)%hb_blks(ielem,itime)%get_recv_comm(imat)
                                    recv_domain  = A%dom(idom)%hb_blks(ielem,itime)%get_recv_domain(imat)
                                    recv_element = A%dom(idom)%hb_blks(ielem,itime)%get_recv_element(imat)

                                    res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                    x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(tparent)
                                    x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(tparent)
                                    associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),    &
                                                xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec(x_istart:x_iend),  &
                                                Amat   => A%dom(idom)%hb_blks(ielem,itime)%data_(imat)%mat ) 

                                        ! Test matrix vector sizes
                                        nonconforming = ( size(Amat,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Harmonic Balance m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated


                end do ! idom
            
            end do ! itime


            !
            ! Wait until all sends have been recieved
            !
            call timer_comm%start()
            call x%comm_wait()
            call timer_comm%stop()

        end if


    end function matrix_vector_product
    !****************************************************************************************












    !>  This function implements the important matrix-vector multiplication 
    !!  operation : [(A)^T]*x :  for chidg_matrix * chidg_vector operations.
    !!
    !!  The structure of the operation is as follows:
    !!      A (in-out)
    !!      x (in-out)
    !!
    !!
    !!      Processor local portion:
    !!      ----------------------------------
    !!          : A*x (INTERIOR - lblks   )
    !!          : A*x (CHIMERA  - chi_blks)
    !!          : A*x (BOUNDARY - bc_blks )
    !!
    !!      Matrix-free contributions:
    !!      ----------------------------------
    !!          : A*x (Harmonic Balance)
    !!
    !!      Processor global(parallel) portion
    !!      ----------------------------------
    !!          : A*x (INTERIOR - lblks   )
    !!          : A*x (CHIMERA  - chi_blks)
    !!          : A*x (BOUNDARY - bc_blks )
    !!
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   04/30/2017
    !!
    !!  TODO: find a way to use non-blocking MPI procedures to speed up the multiplication  
    !!
    !------------------------------------------------------------------------------------
    function transposed_matrix_vector_product(A,x,res_init) result(res)
        type(chidg_matrix_t),    intent(inout)              :: A
        type(chidg_vector_t),    intent(inout)              :: x
        type(chidg_vector_t),    intent(inout),  optional   :: res_init

        type(chidg_vector_t)        :: res
        integer(ik)                 :: idom, ielem, iblk, imat, itime, ivar,  &
                                       itime_i, recv_comm, recv_domain, recv_element
        integer(ik)                 :: dparent_g, dparent_l, eparent_g, eparent_l
        integer(ik)                 :: matrix_proc, vector_proc, nrows, ncols, ierr
        integer(ik)                 :: res_istart, res_iend, x_istart, x_iend, res_istart_l, res_iend_l
        integer(ik)                 :: Amat_rows, Amat_cols
        integer(ik)                 :: iproc, dindex_l, eindex_l, res_size, resvec_p_size
        logical                     :: local_multiply, parallel_multiply
        logical                     :: nonconforming = .false.
        logical                     :: HB_flag_status = .false.
        logical                     :: searching_res_vector = .false.
        character(:), allocatable   :: HB_flag
        real(rk),     allocatable   :: D(:,:)
        real(rk),     allocatable   :: temp_1(:), temp_2(:), resvec_add(:), resvec_p(:)
        real(rk),     allocatable   :: Amat_t(:,:), mass_t(:,:)
        integer(ik)                 :: parent_index(2), local_index(2)



        PetscErrorCode :: perr
        PetscReal :: petsc_norm


        if (allocated(x%wrapped_petsc_vector)) then

            call timer_comm%start()
            call timer_comm%stop()
            call timer_blas%start()

            res = x

            call MatMultTranspose(A%wrapped_petsc_matrix%petsc_matrix,x%wrapped_petsc_vector%petsc_vector,res%wrapped_petsc_vector%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'chidg_mv: error calling petsc MatMultTranspose.')

            call timer_blas%stop()

        else


            ! Allocate result and clear
            if (present(res_init)) then
                ! Res has the shape of a given chidg vector res_init
                res = res_init
            else
                ! By default, res is initialized as the input argument vector x
                res = x
            end if
            call res%clear


            ! Check to see if matrix has been initialized with information about where to locate 
            ! vector information being received from other processors.
            if ( .not. A%recv_initialized ) then
                call A%init_recv(x)
            end if


            !
            ! Begin non-blocking send of parallel vector information
            !
            !call timer_comm%start()
            !call x%comm_send()
            !call timer_comm%stop()


            do itime = 1,time_manager_global%ntime

                !
                ! Compute A*x for local matrix-vector product
                !
                do idom = 1,size(A%dom)

                    !
                    ! Routine for proc-local INTERIOR coupling (lblks)
                    !
                    ! Modified the standard MV multiplication by switching idom and ielem with dparent_l and eparent_l in:
                    ! res_istart, res_iend, x_istart, x_iend, res_vec, x_vec
                    ! Transposed Amat -> Amat_t
                    !
                    ! Watch out here! each dense matrix is multiplied by the dense vector of the element itslef and the product is 
                    ! stored in resultant chidg vector at location dparent_l and eparent_l (indexes of the element wrt the 
                    ! linearization is carried out).
                    !
                    do ielem = 1,size(A%dom(idom)%lblks,1)

                        do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()
                        
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )

            
                            if ( local_multiply ) then
                                dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                                eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)


                                res_istart  = res%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                res_iend    = res%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                x_istart    = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                x_iend      = x%dom(idom)%vecs(ielem)%get_time_end(itime)
                                
                                ! Get the number of rows and columns of the densematrix
                                Amat_rows   = A%dom(idom)%lblks(ielem,itime)%data_(imat)%nrows_
                                Amat_cols   = A%dom(idom)%lblks(ielem,itime)%data_(imat)%ncols_   
                                

                                ! Reallocate Amat_t
                                if ( allocated(Amat_t) ) deallocate (Amat_t)
                                allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                if (ierr /= 0 ) call AllocationError 
                                
                                ! Transpose densematrix
                                Amat_t      = A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat_transpose()
                               
                                associate ( resvec => res%dom(dparent_l)%vecs(eparent_l)%vec(res_istart:res_iend),   &
                                            xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) )
                                            !Amat   => A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat )

                                    call timer_blas%start()
                                    resvec = resvec + matmul(Amat_t,xvec)
                                    call timer_blas%stop()


                                end associate

                            end if

                        end do !imat
                    end do !ielem



                    !
                    ! Routine for proc-local CHIMERA coupling (chi_blks)
                    !
                    ! Modified the standard MV multiplication by switching idom and ielem with dparent_l and eparent_l in:
                    ! res_istart, res_iend, x_istart, x_iend, res_vec, x_vec
                    ! Transposed Amat -> Amat_t
                    !
                    ! Watch out here! each dense matrix is multiplied by the dense vector of the element itslef and the product is 
                    ! stored in resultant chidg vector at location dparent_l and eparent_l (indexes of the element wrt the 
                    ! linearization is carried out).
                    !
                    if (allocated(A%dom(idom)%chi_blks)) then
                        do ielem = 1,size(A%dom(idom)%chi_blks,1)
                            do imat = 1,A%dom(idom)%chi_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)

                                local_multiply    = ( matrix_proc == vector_proc )
                                parallel_multiply = ( matrix_proc /= vector_proc )


                                if ( local_multiply ) then
                                    dparent_l = A%dom(idom)%chi_blks(ielem,itime)%dparent_l(imat)
                                    eparent_l = A%dom(idom)%chi_blks(ielem,itime)%eparent_l(imat)


                                    res_istart = res%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                    res_iend   = res%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                    x_istart   = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    x_iend     = x%dom(idom)%vecs(ielem)%get_time_end(itime)

                                    ! Get the number of rows and columns of the densematrix
                                    Amat_rows   = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%nrows_
                                    Amat_cols   = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%ncols_   
                                    

                                    ! Reallocate Amat_t
                                    if ( allocated(Amat_t) ) deallocate (Amat_t)
                                    allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                    if (ierr /= 0 ) call AllocationError 
                                    
                                    ! Transpose densematrix
                                    Amat_t      = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%mat_transpose()
                                    
                                    associate ( resvec => res%dom(dparent_l)%vecs(eparent_l)%vec(res_istart:res_iend),    &
                                                xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) )!, &
                                                !Amat   => A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%mat  ) 

                                        !
                                        ! Test matrix vector sizes
                                        !
                                        nonconforming = ( size(Amat_t,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat_t,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated



                    !
                    ! Routine for proc-local BOUNDARY coupling (bc_blks)
                    !
                    ! Modified the standard MV multiplication by switching idom and ielem with dparent_l and eparent_l in:
                    ! res_istart, res_iend, x_istart, x_iend, res_vec, x_vec
                    ! Transposed Amat -> Amat_t
                    !
                    ! Watch out here! each dense matrix is multiplied by the dense vector of the element itslef and the product is 
                    ! stored in resultant chidg vector at location dparent_l and eparent_l (indexes of the element wrt the 
                    ! linearization is carried out).
                    !
                    if (allocated(A%dom(idom)%bc_blks)) then
                        do ielem = 1,size(A%dom(idom)%bc_blks,1)
                            do imat = 1,A%dom(idom)%bc_blks(ielem,itime)%size()

                                matrix_proc = IRANK
                                vector_proc = A%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)

                                local_multiply    = ( matrix_proc == vector_proc )
                                parallel_multiply = ( matrix_proc /= vector_proc )


                                if ( local_multiply ) then
                                    dparent_l = A%dom(idom)%bc_blks(ielem,itime)%dparent_l(imat)
                                    eparent_l = A%dom(idom)%bc_blks(ielem,itime)%eparent_l(imat)


                                    res_istart = res%dom(dparent_l)%vecs(eparent_l)%get_time_start(itime)
                                    res_iend   = res%dom(dparent_l)%vecs(eparent_l)%get_time_end(itime)
                                    x_istart   = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                    x_iend     = x%dom(idom)%vecs(ielem)%get_time_end(itime)

                                    ! Get the number of rows and columns of the densematrix
                                    Amat_rows   = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%nrows_
                                    Amat_cols   = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%ncols_   
                                    

                                    ! Reallocate Amat_t
                                    if ( allocated(Amat_t) ) deallocate (Amat_t)
                                    allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                    if (ierr /= 0 ) call AllocationError 
                                    
                                    ! Transpose densematrix
                                    Amat_t      = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%mat_transpose()


                                    associate ( resvec => res%dom(dparent_l)%vecs(eparent_l)%vec(res_istart:res_iend),    &
                                                xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) ) !, &
                                                !Amat   => A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%mat ) 

                                        !
                                        ! Test matrix vector sizes
                                        !
                                        nonconforming = ( size(Amat_t,2) /= size(xvec) )
                                        if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                        call timer_blas%start()
                                        resvec = resvec + matmul(Amat_t,xvec)
                                        call timer_blas%stop()

                                    end associate
                                end if


                            end do !imat
                        end do ! ielem
                    end if  ! allocated





                    
                    !
                    ! Routine for harmonic balance computations
                    ! Used only when harmonic balance is specified
                    !
                    ! Modified the standard HB MV multiplication by switching itime with itime_i in D and mass matrix, 
                    ! and transposing the mass matrix itslef.
                    !
                    ! Watch out here! each dense matrix is multiplied by the dense vector of the element itslef and the product is 
                    ! stored in resultant chidg vector at location dparent_l and eparent_l (indexes of the element wrt the 
                    ! linearization is carried out).
                    !
                    ! TODO: make sure it is correct when it comes to adjoint HB
                    !
                    HB_flag = time_manager_global%get_name()

                    if (HB_flag == 'Harmonic Balance' .or. HB_flag == 'Harmonic_Balance' .or. HB_flag == 'harmonic balance' &
                        .or. HB_flag == 'harmonic_balance' .or. HB_flag == 'HB') then

                        D = time_manager_global%D

                        do ielem = 1,size(A%dom(idom)%lblks,1) 

                            imat = A%dom(idom)%lblks(ielem,itime)%get_diagonal()

                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                            local_multiply    = (matrix_proc == vector_proc)
                            parallel_multiply = (matrix_proc /= vector_proc)

                            if (local_multiply) then
                                do itime_i = 1,size(A%dom(idom)%lblks,2)
                                    if (itime_i /= itime) then
                                        

                                        associate ( nvars  => x%dom(idom)%vecs(ielem)%nvars(), &
                                                    nterms => x%dom(idom)%vecs(ielem)%nterms() )!, &
                                                    !mass   => A%dom(idom)%lblks(ielem,itime_i)%mass )
                                        
                                        ! Reallocate mass matrix transpose
                                        !if ( allocated(mass_t) ) deallocate (mass_t)
                                        !allocate(mass_t(nterms,nterms), stat=ierr)
                                        !if (ierr /= 0) call AllocationError
                                        
                                        ! Evaluate mass matrix transpose
                                        mass_t = transpose(A%dom(idom)%lblks(ielem,itime)%mass)

                                        if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                                        allocate(temp_1(nterms),temp_2(nterms), stat=ierr)
                                        if (ierr /= 0) call AllocationError

                                        call timer_blas%start()
                                        do ivar = 1,nvars

                                            temp_1 = D(itime_i,itime)*matmul(mass_t,x%dom(idom)%vecs(ielem)% &
                                                                                  getvar(ivar,itime_i))
                                            temp_2 = res%dom(idom)%vecs(ielem)%getvar(ivar,itime) + temp_1
                                            call res%dom(idom)%vecs(ielem)%setvar(ivar,itime,temp_2)

                                        end do  ! ivar
                                        call timer_blas%stop()

                                        end associate

                                    end if 
                                end do  ! itime_i
                            end if  ! locla_multiply

                        end do  ! ielem

                    end if  ! HB_flag


                end do ! idom

            end do ! itime


        

            !
            ! Begin blocking recv of parallel vector information
            !
            call timer_comm%start()
            !call x%comm_recv()
            !call timer_comm%stop()




            do itime = 1,time_manager_global%ntime

                do iproc = 0,NRANK-1

                    if (iproc == IRANK ) then
                        

                        !
                        ! Compute A*x for parallel matrix-vector product
                        !
                        do idom = 1,size(A%dom)

                            !!
                            !! Routine for global(parallel) INTERIOR coupling (lblks)
                            !!
                            do ielem = 1,size(A%dom(idom)%lblks,1)
                                do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()
                                
                                    matrix_proc = IRANK
                                    vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                                    local_multiply    = ( matrix_proc == vector_proc )
                                    parallel_multiply = ( matrix_proc /= vector_proc )

                    
                                    if ( parallel_multiply ) then
                                        dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                                        eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)
                                            
                                        parent_index(1) = dparent_l 
                                        parent_index(2) = eparent_l
                        
                                        !
                                        ! searching solution dense_vector on parallel processor
                                        !
                                        searching_res_vector = .true.
                                        call MPI_BCast(searching_res_vector,1,MPI_LOGICAL, iproc , ChiDG_COMM, ierr)
                            
                                        ! Broadcast the IRANK of the processor I need to talk to
                                        call MPI_BCast(vector_proc,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)
                                        
                                        ! Send to the parent's processor the parent's domain and element indices
                                        call MPI_Send(parent_index,2,MPI_INTEGER4,vector_proc,0,ChiDG_COMM,ierr)

                                            
                                        ! Get the number of rows and columns of the densematrix
                                        Amat_rows   = A%dom(idom)%lblks(ielem,itime)%data_(imat)%nrows_
                                        Amat_cols   = A%dom(idom)%lblks(ielem,itime)%data_(imat)%ncols_   
                                        

                                        ! Reallocate Amat_t
                                        if ( allocated(Amat_t) ) deallocate (Amat_t)
                                        allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                        if (ierr /= 0 ) call AllocationError 
                                        
                                        ! Transpose densematrix
                                        Amat_t      = A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat_transpose()


                                        x_istart   = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                        x_iend     = x%dom(idom)%vecs(ielem)%get_time_end(itime)
                                        associate ( xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) )

                                            if ( allocated (resvec_p) ) deallocate (resvec_p)
                                            allocate (resvec_p(Amat_cols), stat=ierr)
                                            if (ierr/=0) call AllocationError
                                            
                                            ! Get the size of resvec_p 
                                            resvec_p_size = size(resvec_p)
                                            
                                            ! Carry out MV multiplication
                                            call timer_blas%start()
                                            resvec_p = matmul(Amat_t,xvec)
                                            call timer_blas%stop()

                                            ! Send size of resvec_p that need to be sent to the parent
                                            call MPI_Send(resvec_p_size,1,MPI_INTEGER4,vector_proc,1,ChiDG_COMM,ierr)

                                            ! Send the resvec_p to the parent
                                            call MPI_Send(resvec_p,resvec_p_size,MPI_REAL8,vector_proc,2,ChiDG_COMM,ierr)                                        

                                        end associate
                                        

                                    end if

                                end do !imat
                            end do !ielem



                            !!
                            !! Routine for global(parallel) CHIMERA coupling (chi_blks)
                            !!
                            if (allocated(A%dom(idom)%chi_blks)) then
                                do ielem = 1,size(A%dom(idom)%chi_blks,1)
                                    do imat = 1,A%dom(idom)%chi_blks(ielem,itime)%size()


                                        matrix_proc = IRANK
                                        vector_proc = A%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)

                                        local_multiply    = ( matrix_proc == vector_proc )
                                        parallel_multiply = ( matrix_proc /= vector_proc )


                                        if ( parallel_multiply ) then
                                            dparent_l = A%dom(idom)%chi_blks(ielem,itime)%dparent_l(imat)
                                            eparent_l = A%dom(idom)%chi_blks(ielem,itime)%eparent_l(imat)
                                                
                                            parent_index(1) = dparent_l 
                                            parent_index(2) = eparent_l
                            
                                            !
                                            ! searching solution dense_vector on parallel processor
                                            !
                                            searching_res_vector = .true.
                                            call MPI_BCast(searching_res_vector,1,MPI_LOGICAL, iproc , ChiDG_COMM, ierr)
                                
                                            ! Broadcast the IRANK of the processor I need to talk to
                                            call MPI_BCast(vector_proc,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)
                                            
                                            ! Send to the parent's processor the parent's domain and element indices
                                            call MPI_Send(parent_index,2,MPI_INTEGER4,vector_proc,0,ChiDG_COMM,ierr)

                                            
                                                
                                            ! Get the number of rows and columns of the densematrix
                                            Amat_rows   = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%nrows_
                                            Amat_cols   = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%ncols_   
                                            

                                            ! Reallocate Amat_t
                                            if ( allocated(Amat_t) ) deallocate (Amat_t)
                                            allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                            if (ierr /= 0 ) call AllocationError 
                                            
                                            ! Transpose densematrix
                                            Amat_t      = A%dom(idom)%chi_blks(ielem,itime)%data_(imat)%mat_transpose()


                                            x_istart   = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                            x_iend     = x%dom(idom)%vecs(ielem)%get_time_end(itime)
                                            associate ( xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) )

                                                if ( allocated (resvec_p) ) deallocate (resvec_p)
                                                allocate (resvec_p(Amat_cols), stat=ierr)
                                                if (ierr/=0) call AllocationError
                                                
                                                ! Get the size of resvec_p 
                                                resvec_p_size = size(resvec_p)

                                                ! Carry out MV multiplication
                                                call timer_blas%start()
                                                resvec_p = matmul(Amat_t,xvec)
                                                call timer_blas%stop()

                                                ! Send size of resvec_p that need to be sent to the parent
                                                call MPI_send(resvec_p_size,1,MPI_INTEGER4,vector_proc,1,ChiDG_COMM,ierr)

                                                ! Send the resvec_p to the parent
                                                call MPI_Send(resvec_p,resvec_p_size,MPI_REAL8,vector_proc,2,ChiDG_COMM,ierr)                                        

                                            end associate
                                        
                                        end if


                                    end do !imat
                                end do ! ielem
                            end if  ! allocated



                            !
                            ! Routine for global(parallel) BOUNDARY coupling (bc_blks)
                            !
                            if (allocated(A%dom(idom)%bc_blks)) then
                                do ielem = 1,size(A%dom(idom)%bc_blks,1)
                                    do imat = 1,A%dom(idom)%bc_blks(ielem,itime)%size()


                                        matrix_proc = IRANK
                                        vector_proc = A%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)

                                        local_multiply    = ( matrix_proc == vector_proc )
                                        parallel_multiply = ( matrix_proc /= vector_proc )


                                        if ( parallel_multiply ) then
                                            dparent_l = A%dom(idom)%bc_blks(ielem,itime)%dparent_l(imat)
                                            eparent_l = A%dom(idom)%bc_blks(ielem,itime)%eparent_l(imat)

                                            parent_index(1) = dparent_l 
                                            parent_index(2) = eparent_l
                            
                                            !
                                            ! searching solution dense_vector on parallel processor
                                            !
                                            searching_res_vector = .true.
                                            call MPI_BCast(searching_res_vector,1,MPI_LOGICAL, iproc , ChiDG_COMM, ierr)
                                
                                            ! Broadcast the IRANK of the processor I need to talk to
                                            call MPI_BCast(vector_proc,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)
                                            
                                            ! Send to the parent's processor the parent's domain and element indices
                                            call MPI_Send(parent_index,2,MPI_INTEGER4,vector_proc,0,ChiDG_COMM,ierr)

                                            
                                                
                                            ! Get the number of rows and columns of the densematrix
                                            Amat_rows   = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%nrows_
                                            Amat_cols   = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%ncols_   
                                            

                                            ! Reallocate Amat_t
                                            if ( allocated(Amat_t) ) deallocate (Amat_t)
                                            allocate ( Amat_t(Amat_cols,Amat_rows) , stat=ierr )
                                            if (ierr /= 0 ) call AllocationError 
                                            
                                            ! Transpose densematrix
                                            Amat_t      = A%dom(idom)%bc_blks(ielem,itime)%data_(imat)%mat_transpose()


                                            x_istart   = x%dom(idom)%vecs(ielem)%get_time_start(itime)
                                            x_iend     = x%dom(idom)%vecs(ielem)%get_time_end(itime)
                                            associate ( xvec   => x%dom(idom)%vecs(ielem)%vec(x_istart:x_iend) )

                                                if ( allocated (resvec_p) ) deallocate (resvec_p)
                                                allocate (resvec_p(Amat_cols), stat=ierr)
                                                if (ierr/=0) call AllocationError
                                                
                                                ! Get the size of resvec_p 
                                                resvec_p_size = size(resvec_p)

                                                ! Carry out MV multiplication
                                                call timer_blas%start()
                                                resvec_p = matmul(Amat_t,xvec)
                                                call timer_blas%stop()

                                                ! Send size of resvec_p that need to be sent to the parent
                                                call MPI_send(resvec_p_size,1,MPI_INTEGER4,vector_proc,1,ChiDG_COMM,ierr)

                                                ! Send the resvec_p to the parent
                                                call MPI_Send(resvec_p,resvec_p_size,MPI_REAL8,vector_proc,2,ChiDG_COMM,ierr)                                        

                                            end associate


                                        end if


                                    end do !imat
                                end do ! ielem
                            end if  ! allocated


                        end do ! idom
                    
                        ! Stop looking for chunk of res_vec
                        searching_res_vector = .false.
                        call MPI_BCast(searching_res_vector,1,MPI_LOGICAL, iproc , ChiDG_COMM, ierr)

                    end if !iproc=IRANK

                    
                    
                    ! Receive chunks of res vector from iproc and add them to current res
                    if (iproc /= IRANK) then

                        ! Check if iproc is searching for res_vector
                        call MPI_BCast(searching_res_vector,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                        do while (searching_res_vector)
                            
                            ! Processor that needs to get data from iproc
                            call MPI_BCast(vector_proc,1,MPI_INTEGER4,iproc,ChiDG_COMM,ierr)

                            if ( vector_proc == IRANK ) then
                                                    
                                ! Receive local domain and element indices to find the chunk of res vector required by iproc
                                call MPI_Recv(local_index,2,MPI_INTEGER4,iproc,0,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                
                                dindex_l = local_index(1)
                                eindex_l = local_index(2)
                                
                                
                                res_istart_l = res%dom(dindex_l)%vecs(eindex_l)%get_time_start(itime)
                                res_iend_l   = res%dom(dindex_l)%vecs(eindex_l)%get_time_end(itime)
                                

                                associate ( resvec_l => res%dom(dindex_l)%vecs(eindex_l)%vec(res_istart_l:res_iend_l) )
                                    
                                    ! Receive size of the chunk of the res vector to add
                                    call MPI_Recv(res_size,1,MPI_INTEGER4,iproc,1,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)
                                    
                                    ! Allocate resvec_add
                                    if ( allocated(resvec_add) ) deallocate (resvec_add)
                                    allocate (resvec_add(res_size), stat=ierr)
                                    if (ierr/=0) call AllocationError
                                    
                                    ! Receive chunk of res vector to add to the current res vector 
                                    call MPI_Recv(resvec_add,res_size,MPI_REAL8,iproc,2,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)

                                    ! Test vector chunk sizes
                                    nonconforming = ( res_size /= size(resvec_l) )
                                    if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mtv: nonconforming Chimera m-v operation")

                                    ! Add quantity calculated by iproc to the chunk of res vector
                                    resvec_l = resvec_l + resvec_add                            

                                end associate
                            
                            
                            end if !proc_needed

                            call MPI_BCast(searching_res_vector,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)


                        end do ! while searching_res_vector


                    end if !iproc/=IRANK
                
                end do !iproc

            end do ! itime


            ! Wait until all sends have been recieved
            call timer_comm%stop()


        end if ! petsc_wrapper

    end function transposed_matrix_vector_product
    !*************************************************************************************************






end module operator_chidg_mv
