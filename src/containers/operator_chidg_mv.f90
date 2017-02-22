module operator_chidg_mv
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO, ONE
    use mod_chidg_mpi,      only: IRANK, ChiDG_COMM
    use type_chidg_matrix,  only: chidg_matrix_t
    use type_time_manager,  only: time_manager_t
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


    !> This function implements the important matrix-vector multiplication 
    !! operation : A*x : for multi-domain configurations, which use the chidg'Container' 
    !! type containers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/6/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    function chidg_mv(A,x) result(res)
        type(chidg_matrix_t),    intent(inout)   :: A
        type(chidg_vector_t),    intent(inout)   :: x

        type(chidg_vector_t)        :: res
        type(time_manager_t)        :: time_manager
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


        HB_flag = time_manager%get_name()

        if ( HB_flag == 'Harmonic Balance' .or. HB_flag == 'Harmonic_Balance' .or. HB_flag == 'harmonic balance' &
             .or. HB_flag == 'harmonic_balance' .or. HB_flag == 'HB') then

            HB_flag_status = .true.
            D = time_manager%D

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


        do itime = 1,time_manager%ntime

            !
            ! Compute A*x for local matrix-vector product
            !
            do idom = 1,size(A%dom)

                !
                ! Routine for neighbor/diag blocks (lblks)
                !
                do ielem = 1,size(A%dom(idom)%lblks,1)
!                    do itime_i = 1,size(A%dom(idom)%lblks,2)
                        do imat = 1,A%dom(idom)%lblks(ielem,itime)%size()
                        
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )

            
                            if ( local_multiply ) then
                                dparent_l = A%dom(idom)%lblks(ielem,itime)%dparent_l(imat)
                                eparent_l = A%dom(idom)%lblks(ielem,itime)%eparent_l(imat)

    !                            associate ( resvec => res%dom(idom)%vecs(ielem)%vec,    &
    !                                        xvec   => x%dom(idom)%vecs(eparent_l)%vec,  &
    !                                        Amat   => A%dom(idom)%lblks(ielem,itime)%mat )
    !
    !                                resvec = resvec + matmul(Amat,xvec)
    !
    !                            end associate

                                res_istart = res%dom(idom)%vecs(ielem)%get_time_start(itime)
                                res_iend   = res%dom(idom)%vecs(ielem)%get_time_end(itime)
                                x_istart = x%dom(idom)%vecs(eparent_l)%get_time_start(itime)
                                x_iend   = x%dom(idom)%vecs(eparent_l)%get_time_end(itime)
                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec(res_istart:res_iend),    &
                                            xvec   => x%dom(idom)%vecs(eparent_l)%vec(x_istart:x_iend),  &
                                            Amat   => A%dom(idom)%lblks(ielem,itime)%data_(imat)%mat )

                                   
                                    !
                                    ! TODO: Where does timer move?
                                    !
                                    call timer_blas%start()
 
                                    resvec = resvec + matmul(Amat,xvec)
                                    
                                    if (HB_flag_status .and. &
                                        imat == A%dom(idom)%lblks(ielem,itime)%get_diagonal()) then

                                        do itime_i = 1,size(A%dom(idom)%lblks,2)
                                          
                                            if (itime_i /= itime) then
                                         
                                                !
                                                ! TODO: Need another associate statement?
                                                ! TODO: ielem = eparent_l when imat = DIAG?
                                                !
                                                associate ( nvars  => x%dom(idom)%vecs(ielem)%nvars(), &
                                                            nterms => x%dom(idom)%vecs(ielem)%nterms(), &
                                                            mass   => A%dom(idom)%lblks(ielem,itime_i)%mass )

                                                if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                                                allocate(temp_1(nterms), temp_2(nterms), stat=ierr)
                                                if (ierr /= 0) call AllocationError

                                                do ivar = 1,nvars

                                                    temp_1 = D(itime,itime_i)*matmul(mass,x%dom(idom)%vecs(ielem)% &
                                                                                          getvar(ivar,itime_i))
                                                    temp_2 = res%dom(idom)%vecs(ielem)%getvar(ivar,itime_i) + temp_1
                                                    call res%dom(idom)%vecs(ielem)%setvar(ivar,itime_i,temp_2)

                                                end do  ! ivar

                                                end associate

                                            end if

                                        end do  ! itime_i

                                    end if

                                    call timer_blas%stop()


                                end associate

                            end if

                        end do !imat
!                    end do !itime
                end do !ielem



                !
                ! Routine for off-diagonal, chimera blocks
                !
                if (allocated(A%dom(idom)%chi_blks)) then
                    do ielem = 1,size(A%dom(idom)%chi_blks,1)
                        !do itime = 1,size(A%dom(idom)%chi_blks,2)
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
                        !end do ! itime
                    end do ! ielem
                end if  ! allocated



!            !
!            ! Routine for boundary condition blocks
!            !
!            if ( allocated(A%dom(idom)%bc_blks) ) then
!                do ielem = 1,size(A%dom(idom)%bc_blks,1)
!                    do iblk = 1,size(A%dom(idom)%bc_blks,2)
!
!
!                        if ( allocated(A%dom(idom)%bc_blks(ielem,iblk)%mat) ) then
!                             dparent = A%dom(idom)%bc_blks(ielem,iblk)%dparent()
!                             eparent = A%dom(idom)%bc_blks(ielem,iblk)%eparent()
!
!                            associate ( resvec => res%dom(idom)%vecs(ielem)%vec,        &
!                                        xvec   => x%dom(dparent)%vecs(eparent)%vec,     &
!                                        Amat   => A%dom(idom)%bc_blks(ielem,iblk)%mat   ) 
!
!                                !
!                                ! Test matrix vector sizes
!                                !
!                                nonconforming = ( size(Amat,2) /= size(xvec) )
!                                if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")
!
!
!                                !
!                                ! Do MV multiply and add to vector
!                                !
!                                resvec = resvec + matmul(Amat,xvec)
!
!                                ! Test without global coupling
!                                !if (ielem == eparent) then
!                                !    resvec = resvec + matmul(Amat,xvec)
!                                !end if
!
!
!                            end associate
!                        end if
!
!                    end do ! iblk
!                end do ! ielem
!            end if  ! allocated

            end do ! idom

        end do ! itime


    


        !
        ! Begin blocking recv of parallel vector information
        !
        call timer_comm%start()
        call x%comm_recv()
        call timer_comm%stop()




        do itime = 1,time_manager%ntime

            !
            ! Compute A*x for parallel matrix-vector product
            !
            do idom = 1,size(A%dom)

                !
                ! Routine for neighbor/diag blocks (lblks)
                !
                do ielem = 1,size(A%dom(idom)%lblks,1)
                    !do itime = 1,size(A%dom(idom)%lblks,2)
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


                                    if (HB_flag_status .and. &
                                        imat == A%dom(idom)%lblks(ielem,itime)%get_diagonal()) then

                                        do itime_i = 1,size(A%dom(idom)%lblks,2)
                                          
                                            if (itime_i /= itime) then
                                           
                                                !
                                                ! TODO: Need another associate statement?
                                                ! TODO: recv_element = ielem when imat = DIAG?
                                                !

                                                associate ( nvars  => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(ielem)%nvars(), &
                                                            nterms => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(ielem)%nterms(), &
                                                            mass   => A%dom(idom)%lblks(ielem,itime_i)%mass )

                                                if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                                                allocate(temp_1(nterms), temp_2(nterms), stat=ierr)
                                                if (ierr /= 0) call AllocationError

                                                do ivar = 1,nvars

                                                    temp_1 = D(itime,itime_i)*matmul(mass,x%recv%comm(recv_comm)% &
                                                                                    dom(recv_domain)%vecs(ielem)%getvar(ivar,itime_i))
                                                    temp_2 = res%dom(idom)%vecs(ielem)%getvar(ivar,itime_i) + temp_1
                                                    call res%dom(idom)%vecs(ielem)%setvar(ivar,itime_i,temp_2)

                                                end do  ! ivar

                                                end associate

                                            end if

                                        end do  ! itime_i

                                    end if

                                end associate
                            end if

                        end do !imat
                    !end do !itime
                end do !ielem



                !
                ! Routine for off-diagonal, chimera blocks
                !
                if (allocated(A%dom(idom)%chi_blks)) then
                    do ielem = 1,size(A%dom(idom)%chi_blks,1)
                        !do itime = 1,size(A%dom(idom)%chi_blks,2)
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
                        !end do ! itime
                    end do ! ielem
                end if  ! allocated




                !
                ! TODO: ADD COUPLED BC DATA
                !

            end do ! idom
        
        end do ! itime


        !
        ! Wait until all sends have been recieved
        !
        call timer_comm%start()
        call x%comm_wait()
        call timer_comm%stop()



    !end function MULT_chidgMatrix_chidgVector
    end function chidg_mv
    !****************************************************************************************


end module operator_chidg_mv
