module precon_ILU0_overset
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, ONE
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK
    use mod_io,                 only: verbosity

    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_chidg_vector
    implicit none






    !>  ILU0_overset preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_ILU0_overset_t

        type(chidg_matrix_t)     :: LD

    contains
    
        procedure   :: init
        procedure   :: update
        procedure   :: apply

!        procedure   :: restrict

    end type precon_ILU0_overset_t
    !*******************************************************************************************




contains



    !>  Initialize the ILU0_overset preconditioner. This is for allocating storage. In this case, 
    !!  we allocate a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  @param[inout]   domain  domain_t instance containing a mesh component used to 
    !!                          initialize the block matrix
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_ILU0_overset_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data

        call self%LD%init(mesh=data%mesh, mtype='LowerDiagonal')
        call self%LD%clear()

        self%initialized = .true.

    end subroutine init
    !*******************************************************************************************








    !>  Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  Updated for spectral time integrators
    !!
    !!  @author Mayank Sharma + Nathan Wukie
    !!  @date   10/09/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_ILU0_overset_t),   intent(inout)   :: self
        type(chidg_matrix_t),           intent(in)      :: A
        type(chidg_vector_t),           intent(in)      :: b


        integer(ik) :: idom, ielem, itime, idiagA, idiagLD, irow, icol, &
                       eparent_l, ilowerA, ilowerLD, itranspose, dparent_g_lower, &
                       eparent_g_lower


        call write_line(' Computing ILU0_overset factorization', io_proc=GLOBAL_MASTER, silence=(verbosity<5))


        !
        ! Test preconditioner initialization
        !
        if ( .not. self%initialized ) call chidg_signal(FATAL,'ILU0_overset%update: preconditioner has not yet been initialized.')

       
        do itime = 1,size(A%dom(1)%lblks,2)
        
            !
            ! For each domain
            !
            do idom = 1,size(A%dom)



                !
                ! Store diagonal blocks of A
                !
                do ielem = 1,size(A%dom(idom)%lblks,1)

                        idiagA = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                        idiagLD = self%LD%dom(idom)%lblks(ielem,itime)%get_diagonal()

                        self%LD%dom(idom)%lblks(ielem,itime)%data_(idiagLD)%mat = A%dom(idom)%lblks(ielem,itime)%data_(idiagA)%mat

                end do !ielem


                !
                ! Invert first diagonal block
                !
                idiagLD = self%LD%dom(idom)%lblks(1,itime)%get_diagonal()
                self%LD%dom(idom)%lblks(1,itime)%data_(idiagLD)%mat = inv(self%LD%dom(idom)%lblks(1,itime)%data_(idiagLD)%mat)


                !
                ! Loop through all Proc-Local rows
                !
                do irow = 2,size(A%dom(idom)%lblks,1)


                    !
                    ! Operate on all the L blocks for the current row
                    !
                    do icol = 1,A%dom(idom)%local_lower_blocks(irow,itime)%size()

                        ilowerA = A%dom(idom)%local_lower_blocks(irow,itime)%at(icol)

                        dparent_g_lower = A%dom(idom)%lblks(irow,itime)%dparent_g(ilowerA)
                        eparent_g_lower = A%dom(idom)%lblks(irow,itime)%eparent_g(ilowerA)

                        ilowerLD = self%LD%dom(idom)%lblks(irow,itime)%loc(dparent_g_lower,eparent_g_lower,itime)

                        if (A%dom(idom)%lblks(irow,itime)%parent_proc(ilowerA) == IRANK) then

                            ! Get parent index
                            eparent_l = A%dom(idom)%lblks(irow,itime)%eparent_l(ilowerA)

                            ! Get diagonal entry
                            idiagLD = self%LD%dom(idom)%lblks(eparent_l,itime)%get_diagonal()

                            ! Compute and store the contribution to the lower-triangular part of LD
                            self%LD%dom(idom)%lblks(irow,itime)%data_(ilowerLD)%mat = matmul(A%dom(idom)%lblks(irow,itime)%data_(ilowerA)%mat,self%LD%dom(idom)%lblks(eparent_l,itime)%data_(idiagLD)%mat)

                            ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                            itranspose = A%dom(idom)%lblks(irow,itime)%itranspose(ilowerA)
                            idiagLD = self%LD%dom(idom)%lblks(irow,itime)%get_diagonal()

                            ! Compute and store the contribution to the lower-triangular part of LD
                            self%LD%dom(idom)%lblks(irow,itime)%data_(idiagLD)%mat = self%LD%dom(idom)%lblks(irow,itime)%data_(idiagLD)%mat  -  &
                                     matmul(self%LD%dom(idom)%lblks(irow,itime)%data_(ilowerLD)%mat,  A%dom(idom)%lblks(eparent_l,itime)%data_(itranspose)%mat)

                        end if

                    end do ! icol


                    !
                    ! Pre-Invert current diagonal block and store
                    !
                    idiagLD = self%LD%dom(idom)%lblks(irow,itime)%get_diagonal()
                    self%LD%dom(idom)%lblks(irow,itime)%data_(idiagLD)%mat = inv(self%LD%dom(idom)%lblks(irow,itime)%data_(idiagLD)%mat)


                end do !irow

            end do ! idom

        end do  ! itime

        call write_line(' Done Computing ILU0_overset factorization', io_proc=GLOBAL_MASTER, silence=(verbosity<5))


    end subroutine update
    !*******************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  Updated for spectral time integrators
    !!
    !!  @author Mayank Sharma + Nathan Wukie
    !!  @date   10/09/2017
    !!
    !-------------------------------------------------------------------------------------------
    function apply(self,A,v,z_old) result(z)
        class(precon_ILU0_overset_t),   intent(inout)           :: self
        type(chidg_matrix_t),           intent(in)              :: A
        type(chidg_vector_t),           intent(in)              :: v
        type(chidg_vector_t),           intent(in), optional    :: z_old

        type(chidg_vector_t)         :: z, overset

        integer(ik)             :: ielem, itime, idiag, dparent_l, eparent_l, idom, irow, icol, &
                                   ilowerA, ilowerLD, iupper, dparent_g_lower, eparent_g_lower, &
                                   inner_time, precon_ntime, matrix_proc, vector_proc, recv_comm, recv_domain, recv_element
        real(rk),   allocatable :: temp(:)

        logical :: local_multiply


        call self%timer%start()


        !
        ! Initialize z for preconditioning
        !
        z = v


        !
        ! Set ntime for preconditioner computations
        ! TODO: Can also be set using time manager data
        !
        precon_ntime = z%get_ntime()


        do itime = 1,precon_ntime
            do idom = 1,size(A%dom)


                ! Forward Solve - Local
                do irow = 1,size(self%LD%dom(idom)%lblks,1)

                    ! Lower-Triangular blocks
                    do icol = 1,A%dom(idom)%local_lower_blocks(irow,itime)%size()

                        ilowerA = A%dom(idom)%local_lower_blocks(irow,itime)%at(icol)
                        dparent_g_lower = A%dom(idom)%lblks(irow,itime)%dparent_g(ilowerA)
                        eparent_g_lower = A%dom(idom)%lblks(irow,itime)%eparent_g(ilowerA)
                        ilowerLD = self%LD%dom(idom)%lblks(irow,itime)%loc(dparent_g_lower,eparent_g_lower,itime)

                        if ( A%dom(idom)%lblks(irow,itime)%parent_proc(ilowerA) == IRANK ) then
                                
                            eparent_l = self%LD%dom(idom)%lblks(irow,itime)%eparent_l(ilowerLD)
                            if (allocated(temp)) deallocate(temp)
                            allocate(temp(size(z%dom(idom)%vecs(irow)%gettime(itime))))
                            temp = z%dom(idom)%vecs(irow)%gettime(itime) - matmul(self%LD%dom(idom)%lblks(irow,itime)%data_(ilowerLD)%mat, z%dom(idom)%vecs(eparent_l)%gettime(itime))
                            call z%dom(idom)%vecs(irow)%settime(itime,temp)

                        end if

                    end do

                end do ! irow



                ! Backward Solve
                do irow = size(A%dom(idom)%lblks,1),1,-1

                    ! Upper-Triangular blocks
                    do icol = 1,A%dom(idom)%local_upper_blocks(irow,itime)%size()

                        iupper = A%dom(idom)%local_upper_blocks(irow,itime)%at(icol)

                        if (A%dom(idom)%lblks(irow,itime)%parent_proc(iupper) == IRANK) then

                                eparent_l = A%dom(idom)%lblks(irow,itime)%eparent_l(iupper)
                                if (allocated(temp)) deallocate(temp)
                                allocate(temp(size(z%dom(idom)%vecs(irow)%gettime(itime))))
                                temp = z%dom(idom)%vecs(irow)%gettime(itime) - matmul(A%dom(idom)%lblks(irow,itime)%data_(iupper)%mat, z%dom(idom)%vecs(eparent_l)%gettime(itime))
                                call z%dom(idom)%vecs(irow)%settime(itime,temp)

                        end if


                    end do


                    ! Diagonal block
                    idiag = self%LD%dom(idom)%lblks(irow,itime)%get_diagonal()
                    if (allocated(temp)) deallocate(temp)
                    allocate(temp(size(z%dom(idom)%vecs(irow)%gettime(itime))))
                    temp = matmul(self%LD%dom(idom)%lblks(irow,itime)%data_(idiag)%mat, z%dom(idom)%vecs(irow)%gettime(itime))
                    call z%dom(idom)%vecs(irow)%settime(itime,temp)


                end do ! irow


            end do ! idom

        end do ! itime



        !
        ! Communicate z, so we have parallel local and overset parts of the vector
        !
        call z%comm_send()
        call z%comm_recv()
        call z%comm_wait()


        ! Copy result from block-local solves
        overset = v
        call overset%clear()



        ! Handle Block-local PARALLEL coupling in an explicit manner
        do itime = 1,precon_ntime
            do idom = 1,size(A%dom)
                do ielem = 1,size(A%dom(idom)%lblks,1)
                    do icol = 1,A%dom(idom)%lblks(ielem,itime)%size()

                        matrix_proc = IRANK
                        vector_proc = A%dom(idom)%lblks(ielem,itime)%parent_proc(icol)
                        local_multiply = (matrix_proc == vector_proc)

                        if ( local_multiply ) then
                            ! local coupling was already handled in ILU decomposition
                        else
                            ! handling parallel block-local coupling in explicit manner
                            recv_comm    = A%dom(idom)%lblks(ielem,itime)%get_recv_comm(icol)
                            recv_domain  = A%dom(idom)%lblks(ielem,itime)%get_recv_domain(icol)
                            recv_element = A%dom(idom)%lblks(ielem,itime)%get_recv_element(icol)
                            associate ( ovec => overset%dom(idom)%vecs(ielem)%vec, &
                                        Amat => A%dom(idom)%lblks(ielem,itime)%data_(icol)%mat, &
                                        zvec => z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec )
                                ovec = ovec + matmul(Amat,zvec)
                            end associate
                        end if 

                    end do !icol
                end do !irow
            end do !idom
        end do !itime


        
        ! Handle ALL(local and parallel) OVERSET coupling in an explicit manner
        do itime = 1,precon_ntime
            do idom = 1,size(A%dom)
                do ielem = 1,size(A%dom(idom)%chi_blks,1)
                    do icol = 1,A%dom(idom)%chi_blks(ielem,itime)%size()

                        matrix_proc = IRANK
                        vector_proc = A%dom(idom)%chi_blks(ielem,itime)%parent_proc(icol)
                        local_multiply = (matrix_proc == vector_proc)

                        ! Handle BOTH local and parallel overset coupling in an explicit manner
                        if ( local_multiply ) then
                            dparent_l = A%dom(idom)%chi_blks(ielem,itime)%data_(icol)%dparent_l()
                            eparent_l = A%dom(idom)%chi_blks(ielem,itime)%data_(icol)%eparent_l()
                            associate ( ovec => overset%dom(idom)%vecs(ielem)%vec, &
                                        Amat => A%dom(idom)%chi_blks(ielem,itime)%data_(icol)%mat, &
                                        zvec => z%dom(dparent_l)%vecs(eparent_l)%vec )
                                ovec = ovec + matmul(Amat,zvec)
                            end associate
                        else
                            recv_comm    = A%dom(idom)%chi_blks(ielem,itime)%get_recv_comm(icol)
                            recv_domain  = A%dom(idom)%chi_blks(ielem,itime)%get_recv_domain(icol)
                            recv_element = A%dom(idom)%chi_blks(ielem,itime)%get_recv_element(icol)
                            associate ( ovec => overset%dom(idom)%vecs(ielem)%vec, &
                                        Amat => A%dom(idom)%chi_blks(ielem,itime)%data_(icol)%mat, &
                                        zvec => z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec )
                                ovec = ovec + matmul(Amat,zvec)
                            end associate
                        end if 

                    end do !icol
                end do !irow
            end do !idom
        end do !itime






        !
        ! OVERSET:
        !
        do itime = 1,precon_ntime
            do idom = 1,size(A%dom)

                ! Forward Solve - Local
                do irow = 1,size(self%LD%dom(idom)%lblks,1)

                    ! Lower-Triangular blocks
                    do icol = 1,A%dom(idom)%local_lower_blocks(irow,itime)%size()

                        ilowerA = A%dom(idom)%local_lower_blocks(irow,itime)%at(icol)
                        dparent_g_lower = A%dom(idom)%lblks(irow,itime)%dparent_g(ilowerA)
                        eparent_g_lower = A%dom(idom)%lblks(irow,itime)%eparent_g(ilowerA)
                        ilowerLD = self%LD%dom(idom)%lblks(irow,itime)%loc(dparent_g_lower,eparent_g_lower,itime)

                        if ( A%dom(idom)%lblks(irow,itime)%parent_proc(ilowerA) == IRANK ) then
                            eparent_l = self%LD%dom(idom)%lblks(irow,itime)%eparent_l(ilowerLD)
                            associate ( ovec => overset%dom(idom)%vecs(irow)%vec, &
                                        mat  => self%LD%dom(idom)%lblks(irow,itime)%data_(ilowerLD)%mat, &
                                        vec  => overset%dom(idom)%vecs(eparent_l)%vec )
                                ovec = ovec - matmul(mat,vec)
                            end associate
                        end if

                    end do

                end do ! irow


                ! Backward Solve
                do irow = size(A%dom(idom)%lblks,1),1,-1
                    ! Upper-Triangular blocks
                    do icol = 1,A%dom(idom)%local_upper_blocks(irow,itime)%size()
                        iupper = A%dom(idom)%local_upper_blocks(irow,itime)%at(icol)

                        if (A%dom(idom)%lblks(irow,itime)%parent_proc(iupper) == IRANK) then
                            eparent_l = A%dom(idom)%lblks(irow,itime)%eparent_l(iupper)
                            associate ( ovec => overset%dom(idom)%vecs(irow)%vec, &
                                        mat  => A%dom(idom)%lblks(irow,itime)%data_(iupper)%mat, &
                                        vec  => overset%dom(idom)%vecs(eparent_l)%vec )
                                ovec = ovec - matmul(mat,vec)
                            end associate
                        end if
                    end do

                    ! Diagonal block
                    idiag = self%LD%dom(idom)%lblks(irow,itime)%get_diagonal()
                    overset%dom(idom)%vecs(irow)%vec = matmul(self%LD%dom(idom)%lblks(irow,itime)%data_(idiag)%mat, overset%dom(idom)%vecs(irow)%vec)

                end do ! irow

            end do ! idom
        end do ! itime


        ! Compute: U^-1 L^-1 y  -  U^-1 L^-1 overset
        z = z - overset


        call self%timer%stop()

    end function apply
    !-----------------------------------------------------------------------------------------








!    !>  Produce a restricted version of the current preconditioner.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   7/24/2017
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    function restrict(self,nterms_r) result(restricted)
!        class(precon_ILU0_overset_t),   intent(in)  :: self
!        integer(ik),            intent(in)  :: nterms_r
!
!        type(precon_ILU0_overset_t) :: restricted
!
!        restricted%LD = self%LD%restrict(nterms_r)
!        restricted%initialized = .true.
!
!    end function restrict
!    !****************************************************************************************









end module precon_ILU0_overset
