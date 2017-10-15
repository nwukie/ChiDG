module operator_chidg_mv
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO, ONE
    use mod_chidg_mpi,      only: IRANK, ChiDG_COMM
    use type_chidg_matrix,  only: chidg_matrix_t
    use mod_time,           only: time_manager_global
    use type_chidg_vector

    use type_timer,         only: timer_t

    implicit none

    external DGEMV

    type(timer_t)   :: timer_comm, timer_blas


!    public operator(*)
!    interface operator(*)
!        module procedure MULT_chidgMatrix_chidgVector
!    end interface

contains


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
    function chidg_mv(A,x) result(res)
        type(chidg_matrix_t),    intent(inout)   :: A
        type(chidg_vector_t),    intent(inout)   :: x

        type(chidg_vector_t)        :: res
        integer(ik)                 :: idom, ielem, iblk, imat, itime, ivar,  &
                                       itime_i, recv_comm, recv_domain, recv_element
        integer(ik)                 :: dparent_g, dparent_l, eparent_g, eparent_l
        integer(ik)                 :: matrix_proc, vector_proc, nrows, ncols, ierr
        integer(ik)                 :: res_istart, res_iend, x_istart, x_iend
        logical                     :: local_multiply, parallel_multiply
        logical                     :: nonconforming = .false.
        logical                     :: HB_flag_status = .false.
        character(:), allocatable   :: HB_flag
        real(rk),     allocatable   :: D(:,:)
        real(rk),     allocatable   :: temp_1(:), temp_2(:)




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
                    do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()
                    
                        matrix_proc = IRANK
                        vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )

        
                        if ( local_multiply ) then
                            dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                            eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)


                            res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                            res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                            x_istart = x%dom(idom)%vecs(eparent_l)%get_time_start(itime)
                            x_iend   = x%dom(idom)%vecs(eparent_l)%get_time_end(itime)
                            associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),   &
                                        xvec   => x%dom(idom)%vecs(eparent_l)%vec(x_istart:x_iend),     &
                                        Amat   => A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat )

                               
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
                            parallel_multiply = ( matrix_proc /= vector_proc )


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
                            parallel_multiply = ( matrix_proc /= vector_proc )


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
                ! Routine for harmonic balance computations
                ! Used only when harmonic balance is specified
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
                                                nterms => x%dom(idom)%vecs(ielem)%nterms(), &
                                                mass   => A%dom(idom)%lblks(ielem,itime_i)%mass )

                                    if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                                    allocate(temp_1(nterms),temp_2(nterms), stat=ierr)
                                    if (ierr /= 0) call AllocationError

                                    call timer_blas%start()
                                    do ivar = 1,nvars

                                        temp_1 = D(itime,itime_i)*matmul(mass,x%dom(idom)%vecs(ielem)%getvar(ivar,itime_i))
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

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )

        
                        if ( parallel_multiply ) then
                            dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                            eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)

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

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                dparent_l = A%dom(idom)%chi_blks(ielem,itime)%dparent_l(imat)
                                eparent_l = A%dom(idom)%chi_blks(ielem,itime)%eparent_l(imat)

                                recv_comm    = A%dom(idom)%chi_blks(ielem,itime)%get_recv_comm(imat)
                                recv_domain  = A%dom(idom)%chi_blks(ielem,itime)%get_recv_domain(imat)
                                recv_element = A%dom(idom)%chi_blks(ielem,itime)%get_recv_element(imat)


                                res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(itime)
                                x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(itime)
                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec,                                    &
                                            xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec,   &
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

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                dparent_l = A%dom(idom)%bc_blks(ielem,itime)%dparent_l(imat)
                                eparent_l = A%dom(idom)%bc_blks(ielem,itime)%eparent_l(imat)

                                recv_comm    = A%dom(idom)%bc_blks(ielem,itime)%get_recv_comm(imat)
                                recv_domain  = A%dom(idom)%bc_blks(ielem,itime)%get_recv_domain(imat)
                                recv_element = A%dom(idom)%bc_blks(ielem,itime)%get_recv_element(imat)


                                res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                x_istart   = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_start(itime)
                                x_iend     = x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%get_time_end(itime)
                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec,                                    &
                                            xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec,   &
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


            end do ! idom
        
        end do ! itime


        !
        ! Wait until all sends have been recieved
        !
        call timer_comm%start()
        call x%comm_wait()
        call timer_comm%stop()



    end function chidg_mv
    !****************************************************************************************


end module operator_chidg_mv
