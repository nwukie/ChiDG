module operator_chidg_mv
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ONE
    use mod_chidg_mpi,      only: IRANK, ChiDG_COMM
    use type_chidgMatrix,   only: chidgMatrix_t
    use type_chidgVector

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
    !function MULT_chidgMatrix_chidgVector(A,x) result(res)
    function chidg_mv(A,x) result(res)
        type(chidgMatrix_t),    intent(inout)   :: A
        type(chidgVector_t),    intent(inout)   :: x

        type(chidgVector_t)     :: res
        integer(ik)             :: idom, ielem, iblk, recv_comm, recv_domain, recv_element
        integer(ik)             :: dparent_g, dparent_l, eparent_g, eparent_l
        integer(ik)             :: matrix_proc, vector_proc, nrows, ncols, ierr
        logical                 :: local_multiply, parallel_multiply
        logical                 :: nonconforming = .false.



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


        !
        ! Compute A*x for local matrix-vector product
        !
        do idom = 1,size(A%dom)

            !
            ! Routine for neighbor/diag blocks (lblks)
            !
            !$OMP PARALLEL DO PRIVATE(matrix_proc, vector_proc, local_multiply, parallel_multiply, dparent_l, eparent_l, iblk)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do iblk = 1,size(A%dom(idom)%lblks,2)
                    
                    if (allocated(A%dom(idom)%lblks(ielem,iblk)%mat)) then
                        matrix_proc = IRANK
                        vector_proc = A%dom(idom)%lblks(ielem,iblk)%parent_proc()

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )

        
                        if ( local_multiply ) then
                            dparent_l = A%dom(idom)%lblks(ielem,iblk)%dparent_l()
                            eparent_l = A%dom(idom)%lblks(ielem,iblk)%eparent_l()

                            associate ( resvec => res%dom(idom)%vecs(ielem)%vec,    &
                                        xvec   => x%dom(idom)%vecs(eparent_l)%vec,  &
                                        Amat   => A%dom(idom)%lblks(ielem,iblk)%mat )

                                !resvec = resvec + matmul(Amat,xvec)
                                !res%dom(idom)%vecs(ielem)%vec = res%dom(idom)%vecs(ielem)%vec + matmul(A%dom(idom)%lblks(ielem,iblk)%mat, x%dom(idom)%vecs(eparent_l)%vec)
                                
                                !nrows = size(Amat,1)
                                !ncols = size(Amat,2)
                                !call DGEMV('N', nrows, ncols, ONE, Amat, nrows, xvec, 1, ONE, resvec, 1)
                                nrows = size(A%dom(idom)%lblks(ielem,iblk)%mat,1)
                                ncols = size(A%dom(idom)%lblks(ielem,iblk)%mat,2)
                                call timer_blas%start()
                                call DGEMV('N', nrows, ncols, ONE, A%dom(idom)%lblks(ielem,iblk)%mat, nrows, x%dom(idom)%vecs(eparent_l)%vec, 1, ONE, res%dom(idom)%vecs(ielem)%vec, 1)
                                !call chidg_matmul(A%dom(idom)%lblks(ielem,iblk)%mat,x%dom(idom)%vecs(eparent_l)%vec, res%dom(idom)%vecs(ielem)%vec)
                                call timer_blas%stop()

                            end associate
                        end if

                    end if

                end do
            end do
            !$OMP END PARALLEL DO



            !
            ! Routine for off-diagonal, chimera blocks
            !
            if (allocated(A%dom(idom)%chi_blks)) then
                !$OMP PARALLEL DO PRIVATE(matrix_proc, vector_proc, local_multiply, parallel_multiply, dparent_l, eparent_l, iblk, nonconforming)
                do ielem = 1,size(A%dom(idom)%chi_blks,1)
                    do iblk = 1,size(A%dom(idom)%chi_blks,2)


                        if (allocated(A%dom(idom)%chi_blks(ielem,iblk)%mat)) then
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%chi_blks(ielem,iblk)%parent_proc()

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( local_multiply ) then
                                dparent_l = A%dom(idom)%chi_blks(ielem,iblk)%dparent_l()
                                eparent_l = A%dom(idom)%chi_blks(ielem,iblk)%eparent_l()

                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec,        &
                                            xvec   => x%dom(dparent_l)%vecs(eparent_l)%vec, &
                                            Amat   => A%dom(idom)%chi_blks(ielem,iblk)%mat  ) 

                                    !
                                    ! Test matrix vector sizes
                                    !
                                    nonconforming = ( size(Amat,2) /= size(xvec) )
                                    if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                    resvec = resvec + matmul(Amat,xvec)
                                    !res%dom(idom)%vecs(ielem)%vec = res%dom(idom)%vecs(ielem)%vec + matmul(A%dom(idom)%chi_blks(ielem,iblk)%mat, x%dom(dparent_l)%vecs(eparent_l)%vec)


                                    !nrows = size(A%dom(idom)%chi_blks(ielem,iblk)%mat,1)
                                    !ncols = size(A%dom(idom)%chi_blks(ielem,iblk)%mat,2)
                                    call timer_blas%start()
                                    !call DGEMV('N', nrows, ncols, ONE, A%dom(idom)%chi_blks(ielem,iblk)%mat, nrows, x%dom(dparent_l)%vecs(eparent_l)%vec, 1, ONE, res%dom(idom)%vecs(ielem)%vec, 1)
                                    call timer_blas%stop()

                                end associate
                            end if

                        end if

                    end do ! iblk
                end do ! ielem
                !OMP END PARALLEL DO
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




    


        !
        ! Begin blocking recv of parallel vector information
        !
        call timer_comm%start()
        call x%comm_recv()
        call timer_comm%stop()






        !
        ! Compute A*x for parallel matrix-vector product
        !
        do idom = 1,size(A%dom)

            !
            ! Routine for neighbor/diag blocks (lblks)
            !
            !$OMP PARALLEL DO PRIVATE(matrix_proc, vector_proc, local_multiply, parallel_multiply, dparent_l, eparent_l, iblk, recv_comm, recv_domain, recv_element)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do iblk = 1,size(A%dom(idom)%lblks,2)
                    
                    if (allocated(A%dom(idom)%lblks(ielem,iblk)%mat)) then
                        matrix_proc = IRANK
                        vector_proc = A%dom(idom)%lblks(ielem,iblk)%parent_proc()

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )

        
                        if ( parallel_multiply ) then
                            dparent_l = A%dom(idom)%lblks(ielem,iblk)%dparent_l()
                            eparent_l = A%dom(idom)%lblks(ielem,iblk)%eparent_l()

                            recv_comm    = A%dom(idom)%lblks(ielem,iblk)%recv_comm
                            recv_domain  = A%dom(idom)%lblks(ielem,iblk)%recv_domain
                            recv_element = A%dom(idom)%lblks(ielem,iblk)%recv_element

                            associate ( resvec => res%dom(idom)%vecs(ielem)%vec,                                    &
                                        xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec,   &
                                        Amat   => A%dom(idom)%lblks(ielem,iblk)%mat )

                                !resvec = resvec + matmul(Amat,xvec)
                                !res%dom(idom)%vecs(ielem)%vec = res%dom(idom)%vecs(ielem)%vec + matmul(A%dom(idom)%lblks(ielem,iblk)%mat, x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec)

                                nrows = size(A%dom(idom)%lblks(ielem,iblk)%mat,1)
                                ncols = size(A%dom(idom)%lblks(ielem,iblk)%mat,2)
                                call timer_blas%start()
                                call DGEMV('N', nrows, ncols, ONE, A%dom(idom)%lblks(ielem,iblk)%mat, nrows, x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec, 1, ONE, res%dom(idom)%vecs(ielem)%vec, 1)
                                !call chidg_matmul(A%dom(idom)%lblks(ielem,iblk)%mat,x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec, res%dom(idom)%vecs(ielem)%vec)
                                call timer_blas%stop()

                            end associate
                        end if

                    end if

                end do
            end do
            !$OMP END PARALLEL DO



            !
            ! Routine for off-diagonal, chimera blocks
            !
            if (allocated(A%dom(idom)%chi_blks)) then
                !$OMP PARALLEL DO PRIVATE(matrix_proc, vector_proc, local_multiply, parallel_multiply, dparent_l, eparent_l, iblk, recv_comm, recv_domain, recv_element, nonconforming)
                do ielem = 1,size(A%dom(idom)%chi_blks,1)
                    do iblk = 1,size(A%dom(idom)%chi_blks,2)


                        if (allocated(A%dom(idom)%chi_blks(ielem,iblk)%mat)) then
                            matrix_proc = IRANK
                            vector_proc = A%dom(idom)%chi_blks(ielem,iblk)%parent_proc()

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                dparent_l = A%dom(idom)%chi_blks(ielem,iblk)%dparent_l()
                                eparent_l = A%dom(idom)%chi_blks(ielem,iblk)%eparent_l()

                                recv_comm    = A%dom(idom)%chi_blks(ielem,iblk)%recv_comm
                                recv_domain  = A%dom(idom)%chi_blks(ielem,iblk)%recv_domain
                                recv_element = A%dom(idom)%chi_blks(ielem,iblk)%recv_element


                                associate ( resvec => res%dom(idom)%vecs(ielem)%vec,                                    &
                                            xvec   => x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec,   &
                                            Amat   => A%dom(idom)%chi_blks(ielem,iblk)%mat )


                                    !
                                    ! Test matrix vector sizes
                                    !
                                    nonconforming = ( size(Amat,2) /= size(xvec) )
                                    if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                    resvec = resvec + matmul(Amat,xvec)
                                    !res%dom(idom)%vecs(ielem)%vec = res%dom(idom)%vecs(ielem)%vec + matmul(A%dom(idom)%chi_blks(ielem,iblk)%mat, x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec)


                                    !nrows = size(A%dom(idom)%chi_blks(ielem,iblk)%mat,1)
                                    !ncols = size(A%dom(idom)%chi_blks(ielem,iblk)%mat,2)
                                    !call DGEMV('N', nrows, ncols, ONE, A%dom(idom)%chi_blks(ielem,iblk)%mat, nrows, x%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec, 1, ONE, res%dom(idom)%vecs(ielem)%vec, 1)


                                end associate
                            end if

                        end if

                    end do ! iblk
                end do ! ielem
                !OMP END PARALLEL DO
            end if  ! allocated




            !
            ! TODO: ADD COUPLED BC DATA
            !

        end do ! idom



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
