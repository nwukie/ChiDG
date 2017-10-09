!>  Restricted Additive Schwarz(RAS) preconditioner. ILU0 local solve.
!!
!!
!!  Legend:
!!  =================================================================
!!  local = The current processor's local elements and their coupling
!!  proc# = Another processor's elements
!!  cpl   = Coupling between the elements of neighboring processors
!!
!!
!!  In a domain decomposition block-Jacobi preconditioner, the processor-local problem is solved 
!!  without taking the inter-processor coupling of the problem into account.
!!
!!  |--------|---------|--------|
!!  |        |         |        |
!!  | proc 0 |   cpl   |  cpl   |
!!  |        |         |        |
!!  |--------|---------|--------|                                      |---------|  |-|     |-|
!!  |        |         |        |           block-Jacobi               |         |  | |     | |
!!  |  cpl   |  local  |  cpl   |           ----------->               |  local  |  |x|  =  |b|
!!  |        |         |        |                                      |         |  | |     | |
!!  |--------|---------|--------|                                      |---------|  |-|     |-|
!!  |        |         |        |
!!  |  cpl   |   cpl   | proc 2 |
!!  |        |         |        |
!!  |--------|---------|--------|
!!
!!
!!
!!  The RAS preconditioner considers the local problem and also includes a 1-element overlap 
!!  with the neighboring processors and their coupling. The implementation here rearranges the 
!!  problem so that all overlap data is put at the end of the preconditioning matrix. The 
!!  restricted part of the preconditioner means that the solution of the preconditioning matrix 
!!  is only applied to the processor-local portion of the global vector. This is why only a 
!!  portion of the x-vector is shown here in the diagram.
!!
!!  |--------|---------|--------|
!!  |        |         |        |
!!  | proc 0 |   cpl   |  cpl   |
!!  |        |         |        |
!!  |--------|---------|--------|                               |---------|--------| |-|    |-|
!!  |        |         |        |  Restricted Additive Schwarz  |         |        | | |    | |
!!  |  cpl   |  local  |  cpl   |          ----------->         |  local  |  cpl   | |x|    | |
!!  |        |         |        |                               |         |        | | |    | |
!!  |--------|---------|--------|                               |---------|--------| |-|  = |b|
!!  |        |         |        |                               |         | proc#  |        | |
!!  |  cpl   |   cpl   | proc 2 |                               |  cpl    | overlap|        | |
!!  |        |         |        |                               |---------|--------|        |-|
!!  |--------|---------|--------|
!!
!!
!!
!!
!-----------------------------------------------------------------------------------------------
module precon_RASILU0
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, &
                                          ETA_MAX, ZETA_MAX, ONE
    use mod_inv,                    only: inv
    use mod_chidg_mpi,              only: IRANK, NRANK, ChiDG_COMM
    use mod_io,                     only: verbosity

    use type_RASILU0_send,          only: RASILU0_send_t
    use type_RASILU0_recv,          only: RASILU0_recv_t
    use type_preconditioner,        only: preconditioner_t
    use type_chidg_data,            only: chidg_data_t
    use type_chidg_matrix,          only: chidg_matrix_t
    use type_chidg_vector,          only: chidg_vector_t

    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08,                    only: MPI_ISend, MPI_Recv, MPI_REAL8, MPI_REQUEST, &
                                          MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, MPI_Barrier
    implicit none

    external DGEMV

    !>  Restricted Additive Schwarz(RAS) preconditioner using ILU0 to solve the local problems 
    !!  includeing the RAS overlap data.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_RASILU0_t

        type(chidg_matrix_t)        :: LD       !< Lower-Diagonal matrix for local problems

        type(RASILU0_send_t)        :: send     !< Overlapping data to send to other processors
        type(RASILU0_recv_t)        :: recv     !< Overlapping data to receive from other processors

        type(mpi_request_vector_t)  :: mpi_requests

    contains
    
        procedure   :: init
        procedure   :: update
        procedure   :: apply

        procedure   :: comm_send
        procedure   :: comm_recv
        procedure   :: comm_wait

    end type precon_RASILU0_t
    !******************************************************************************************




contains



    !>  Initialize the ILU0 preconditioner. This is for allocating storage. In this case, we 
    !!  allocate a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/10/2016
    !!
    !!  @param[inout]   data      chidg data container, with mesh, solution, etc.
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_RASILU0_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data

        integer :: ierr

        integer(ik) :: iread, ielem, iblk, diag

        call write_line('   Restricted Additive Schwarz(RAS) preconditioner: ', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        !
        ! Initialize Lower-Diagonal matrix for processor-local data
        !
        call self%LD%init(mesh=data%mesh, mtype='LowerDiagonal')
        self%initialized = .true.

 
        !
        ! Initialize the overlap data
        !
        call write_line('       RAS: initializing send pattern...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call self%send%init(data%mesh, data%sdata%lhs)
        call write_line('       RAS: initializing receive pattern...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call self%recv%init(data%mesh, data%sdata%lhs, data%sdata%rhs)


        !
        ! Release nonblocking send buffers
        !
        call write_line('       RAS: waiting on remaining communication buffers ...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call self%send%init_wait()
        call write_line('       RAS: initialization complete!', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))


    end subroutine init
    !******************************************************************************************








    !>  Update the preconditioner.
    !!
    !!  For the Restricted Additive Schwarz algorithm, this exchanges the processor overlap data 
    !!  and computes information that is reused every iteration. Here, it is premultiplying some
    !!  of the blocks and computing some matrix inversions that do not depend on the incoming 
    !!  vector.
    !!
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_RASILU0_t),    intent(inout)   :: self
        type(chidg_matrix_t),       intent(in)      :: A
        type(chidg_vector_t),       intent(in)      :: b


        character(:),   allocatable :: user_msg
        integer(ik)                 :: ielem, irow, icol, eparent_l, idom, ndom, ilower, ilowerA, ilowerLD,    &
                                       trans_elem, trans_blk, itranspose, nrowsA, ncolsA, ncolsB,      &
                                       iblk_diag_parent, iblk_diag, iblk, icomm,            &
                                       parent_proc, ierr, iproc, idiagLD, idiagA, dparent_g_lower, eparent_g_lower, itime


        call write_line('   RAS: Computing ILU0 factorization', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))



        !
        ! Communicate matrix overlapping components
        !
        call self%mpi_requests%clear()

        call self%comm_send(A)
        call self%comm_recv()
        call self%comm_wait()




        !
        ! Test preconditioner initialization
        !
        user_msg = 'RAS-ILU0%update: preconditioner has not yet been initialized.'
        if ( .not. self%initialized ) call chidg_signal(FATAL,user_msg)


        !
        ! For each domain
        !
        ndom = size(A%dom)
        do idom = 1,ndom


            !
            ! Store diagonal blocks of A
            !
            do ielem = 1,size(A%dom(idom)%lblks,1)
                idiagA  =       A%dom(idom)%lblks(ielem,1)%get_diagonal()
                idiagLD = self%LD%dom(idom)%lblks(ielem,1)%get_diagonal()
                self%LD%dom(idom)%lblks(ielem,1)%data_(idiagLD)%mat = A%dom(idom)%lblks(ielem,1)%data_(idiagA)%mat
            end do


            !
            ! Invert first diagonal block
            !
            idiagLD = self%LD%dom(idom)%lblks(1,1)%get_diagonal()
            self%LD%dom(idom)%lblks(1,1)%data_(idiagLD)%mat = inv(self%LD%dom(idom)%lblks(1,1)%data_(idiagLD)%mat)


            !
            ! Loop through all Proc-Local rows
            !
            itime = 1
            do irow = 2,size(A%dom(idom)%lblks,1)


                ! Operate on all the L blocks for the current row
                do icol = 1,A%dom(idom)%local_lower_blocks(irow,itime)%size()

                    ilowerA = A%dom(idom)%local_lower_blocks(irow,itime)%at(icol)

                    dparent_g_lower = A%dom(idom)%lblks(irow,1)%dparent_g(ilowerA)
                    eparent_g_lower = A%dom(idom)%lblks(irow,1)%eparent_g(ilowerA)

                    ilowerLD = self%LD%dom(idom)%lblks(irow,1)%loc(dparent_g_lower,eparent_g_lower)

                    if (A%dom(idom)%lblks(irow,1)%parent_proc(ilowerA) == IRANK) then

                        ! Get parent index and transpose block
                        eparent_l   = A%dom(idom)%lblks(irow,1)%eparent_l(ilowerA)

                        ! Compute and store the contribution to the lower-triangular part of LD
                        idiagLD = self%LD%dom(idom)%lblks(eparent_l,1)%get_diagonal()
                        self%LD%dom(idom)%lblks(irow,1)%data_(ilowerLD)%mat = matmul(A%dom(idom)%lblks(irow,1)%data_(ilowerA)%mat,self%LD%dom(idom)%lblks(eparent_l,1)%data_(idiagLD)%mat)

                        ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. (The component in the transposed position)
                        itranspose  = A%dom(idom)%lblks(irow,1)%itranspose(ilowerA)
                        idiagLD = self%LD%dom(idom)%lblks(irow,1)%get_diagonal()
                        self%LD%dom(idom)%lblks(irow,1)%data_(idiagLD)%mat = self%LD%dom(idom)%lblks(irow,1)%data_(idiagLD)%mat  -  &
                                    matmul(self%LD%dom(idom)%lblks(irow,1)%data_(ilowerLD)%mat,  A%dom(idom)%lblks(eparent_l,1)%data_(itranspose)%mat)


                    end if

                end do ! icol


                ! Pre-Invert current diagonal block and store
                idiagLD = self%LD%dom(idom)%lblks(irow,1)%get_diagonal()
                self%LD%dom(idom)%lblks(irow,1)%data_(idiagLD)%mat = inv(self%LD%dom(idom)%lblks(irow,1)%data_(idiagLD)%mat)


            end do !irow

        end do ! idom




        !
        ! Loop through all Proc-overlap lower blocks
        !
        do idom = 1,size(self%recv%dom)
            do icomm = 1,size(self%recv%dom(idom)%comm)
                do ielem = 1,size(self%recv%dom(idom)%comm(icomm)%elem)
                    do iblk = 1,self%recv%dom(idom)%comm(icomm)%elem(ielem)%lower%size()

                        iblk_diag = self%recv%dom(idom)%comm(icomm)%elem(ielem)%diag%at(1)
                        ilower    = self%recv%dom(idom)%comm(icomm)%elem(ielem)%lower%at(iblk)


                        parent_proc = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%parent_proc()
                        eparent_l   = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%eparent_l()
                        trans_elem  = self%recv%dom(idom)%comm(icomm)%elem(ielem)%trans_elem(ilower)
                        trans_blk   = self%recv%dom(idom)%comm(icomm)%elem(ielem)%trans_blk(ilower)



                        ! Compute and store the contribution to the lower-triangular part of A comm, since A comm shouldn't 
                        ! get used anywhere else
                        if ( parent_proc /= IRANK ) then
                            ! If lower block is coupled with another block in the overlap, get DIAGONAL component 
                            ! sent from A since overlap data is stored there
                            iblk_diag_parent = self%recv%dom(idom)%comm(icomm)%elem(trans_elem)%diag%at(1)

                            associate ( lower      => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%mat,  &
                                        upper_diag => self%recv%dom(idom)%comm(icomm)%elem(trans_elem)%blks(iblk_diag_parent)%mat )
                                lower = matmul(lower,upper_diag)
                            end associate

                        else
                            iblk_diag_parent = self%LD%dom(idom)%lblks(eparent_l,1)%get_diagonal()
                            ! If lower block is coupled with an interior block, get DIAGONAL component from LD, since we 
                            ! don't want to overwrite these entries in A as that would blow away data we need for the MV product
                            associate( lower      => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%mat,  &
                                       upper_diag => self%LD%dom(idom)%lblks(eparent_l,1)%data_(iblk_diag_parent)%mat )
                                lower = matmul(lower,upper_diag)
                            end associate
                        end if




                        ! Modify the current diagonal by this lower-triangular part multiplied by opposite upper-triangular part. 
                        ! (The component in the transposed position)
                        if ( parent_proc /= IRANK ) then

                            associate ( diag  => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%mat,   &
                                        lower => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%mat,      &
                                        trans => self%recv%dom(idom)%comm(icomm)%elem(trans_elem)%blks(trans_blk)%mat )
                                diag = diag - matmul(lower,trans)
                            end associate

                        else
                            associate ( diag  => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%mat,   &
                                        lower => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%mat,      &
                                        trans => A%dom(idom)%lblks(eparent_l,1)%data_(trans_blk)%mat )
                                diag = diag - matmul(lower,trans)
                            end associate
                        end if



                    end do !iblk


                    ! Pre-Invert current diagonal block and store
                    iblk_diag = self%recv%dom(idom)%comm(icomm)%elem(ielem)%diag%at(1)
                    self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%mat = inv(self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%mat)


                end do !ielem
            end do !icomm
        end do !idom



        call write_line(' Done Computing RAS-ILU0 factorization', io_proc=GLOBAL_MASTER, silence=(verbosity<5))


    end subroutine update
    !*****************************************************************************************








    !>  Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'.
    !!
    !!  We are using ILU0 for the local solve, so we perform the exchange of boundary vector 
    !!  data, compute the forward solve of the incomplete factorization and then compute the 
    !!  backward solve of the incomplete factorization.
    !!  
    !!  Assumption: All overlap data is placed at the end of the local system of equations.
    !!          
    !!          |         |        |
    !!          |  local  |coupling|
    !!          |         |        |
    !!          |---------|--------|
    !!          |         |        |
    !!          |coupling |overlap |
    !!          |         |        |
    !!
    !!          So, the algorithm is structured as:
    !!              - forward solve  (local)
    !!              - forward solve  (overlap)
    !!              - backward solve (overlap)
    !!              - backward solve (local)
    !!
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_RASILU0_t),   intent(inout)   :: self
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: v

        type(chidg_vector_t) :: z

        integer(ik)         :: ielem, eparent_l, irow, icol, idom, ndom, ilower, iupper, &
                               iblk, iblk_diag, icomm, parent_proc, recv_comm, recv_domain, &
                               recv_element, recv_comm_diag, recv_domain_diag, &
                               recv_element_diag, ncols, nrows, diag_irow, dparent_g_lower, &
                               eparent_g_lower, ilowerA, ilowerLD, itime
        logical             :: interior_block


        call self%timer%start()



        !
        ! Initialize z for preconditioning
        !
        z = v



        !
        ! Exchange boundary vector data
        !
        call z%comm_send()
        call z%comm_recv()
        call z%comm_wait()




        !
        ! For each domain
        !
        ndom = size(A%dom)
        itime = 1
        do idom = 1,ndom


            !
            ! Forward Solve - Local
            !
            do irow = 1,size(self%LD%dom(idom)%lblks,1)


                !
                ! Lower-Triangular blocks
                !
                do icol = 1,A%dom(idom)%local_lower_blocks(irow,itime)%size()

                    ilowerA = A%dom(idom)%local_lower_blocks(irow,itime)%at(icol)

                    dparent_g_lower = A%dom(idom)%lblks(irow,1)%dparent_g(ilowerA)
                    eparent_g_lower = A%dom(idom)%lblks(irow,1)%eparent_g(ilowerA)

                    ilowerLD = self%LD%dom(idom)%lblks(irow,1)%loc(dparent_g_lower,eparent_g_lower)

                    if ( A%dom(idom)%lblks(irow,1)%parent_proc(ilowerA) == IRANK) then
                            ! Get associated parent block index
                            eparent_l = self%LD%dom(idom)%lblks(irow,1)%eparent_l(ilowerLD)
                            !z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(self%LD%dom(idom)%lblks(irow,ilower)%mat, z%dom(idom)%vecs(eparent_l)%vec)

                            nrows = size(self%LD%dom(idom)%lblks(irow,1)%data_(ilowerLD)%mat,1)
                            ncols = size(self%LD%dom(idom)%lblks(irow,1)%data_(ilowerLD)%mat,2)
                            call DGEMV('N', nrows, ncols, -ONE, self%LD%dom(idom)%lblks(irow,1)%data_(ilowerLD)%mat, nrows, z%dom(idom)%vecs(eparent_l)%vec, 1, ONE, z%dom(idom)%vecs(irow)%vec, 1)



                    end if
                end do

            end do ! irow


            !
            ! Forward Solve - Overlap
            !
            do icomm = 1,size(self%recv%dom(idom)%comm)

                do ielem = 1,size(self%recv%dom(idom)%comm(icomm)%elem)
                    iblk_diag = self%recv%dom(idom)%comm(icomm)%elem(ielem)%diag%at(1)

                    ! Location in comm vector
                    recv_comm_diag    = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_comm
                    recv_domain_diag  = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_domain
                    recv_element_diag = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_element


                    do iblk = 1,self%recv%dom(idom)%comm(icomm)%elem(ielem)%lower%size()
                        ilower = self%recv%dom(idom)%comm(icomm)%elem(ielem)%lower%at(iblk)

                        parent_proc = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%parent_proc()
                        eparent_l   = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%eparent_l()
                    
                        associate ( zvec_diag => z%recv%comm(recv_comm_diag)%dom(recv_domain_diag)%vecs(recv_element_diag)%vec, &
                                    lower     => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%mat )

                        if ( parent_proc == IRANK ) then

                            !zvec_diag = zvec_diag - matmul(lower,z%dom(idom)%vecs(eparent_l)%vec)
                            nrows = size(lower,1)
                            ncols = size(lower,2)
                            call DGEMV('N', nrows, ncols, -ONE, lower, nrows, z%dom(idom)%vecs(eparent_l)%vec, 1, ONE, zvec_diag, 1)

                        else
                            recv_comm    = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%recv_comm
                            recv_domain  = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%recv_domain
                            recv_element = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(ilower)%recv_element

                            !zvec_diag = zvec_diag - matmul(lower,z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec)
                            nrows = size(lower,1)
                            ncols = size(lower,2)
                            call DGEMV('N', nrows, ncols, -ONE, lower, nrows, z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec, 1, ONE, zvec_diag, 1)
                        end if

                        end associate


                    end do !icol
                end do !irow

            end do !icomm






            !
            ! Backward solve - Overlap
            !
            do icomm = size(self%recv%dom(idom)%comm),1,-1
                do ielem = size(self%recv%dom(idom)%comm(icomm)%elem),1,-1

                    iblk_diag         = self%recv%dom(idom)%comm(icomm)%elem(ielem)%diag%at(1)
                    recv_comm_diag    = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_comm
                    recv_domain_diag  = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_domain
                    recv_element_diag = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%recv_element


                    do iblk = 1,self%recv%dom(idom)%comm(icomm)%elem(ielem)%upper%size()
                        iupper = self%recv%dom(idom)%comm(icomm)%elem(ielem)%upper%at(iblk)

                        recv_comm    = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iupper)%recv_comm
                        recv_domain  = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iupper)%recv_domain
                        recv_element = self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iupper)%recv_element

                        associate ( zvec_diag  => z%recv%comm(recv_comm_diag)%dom(recv_domain_diag)%vecs(recv_element_diag)%vec, &
                                    upper      => self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iupper)%mat,                  &
                                    zvec_upper => z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec )


                            !zvec_diag = zvec_diag - matmul(upper, zvec_upper)
                            nrows = size(upper,1)
                            ncols = size(upper,2)
                            call DGEMV('N', nrows, ncols, -ONE, upper, nrows, zvec_upper, 1, ONE, zvec_diag, 1)

                        end associate

                    end do !iblk


                    !
                    ! Diagonal block
                    !
                    associate ( zvec_diag  => z%recv%comm(recv_comm_diag)%dom(recv_domain_diag)%vecs(recv_element_diag)%vec )
                        zvec_diag = matmul(self%recv%dom(idom)%comm(icomm)%elem(ielem)%blks(iblk_diag)%mat, zvec_diag)
                    end associate


                end do !ielem
            end do !icomm




            !
            ! Backward Solve - Local
            !
            itime = 1
            do irow = size(A%dom(idom)%lblks,1),1,-1

                !
                ! Upper-Triangular blocks
                !
                do icol = 1,A%dom(idom)%local_upper_blocks(irow,itime)%size()

                    iupper = A%dom(idom)%local_upper_blocks(irow,itime)%at(icol)


                    interior_block = (A%dom(idom)%lblks(irow,1)%parent_proc(iupper) == IRANK)

                    if (interior_block) then

                        ! Get associated parent block index
                        eparent_l = A%dom(idom)%lblks(irow,1)%eparent_l(iupper)
                        !z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(A%dom(idom)%lblks(irow,iupper)%mat, z%dom(idom)%vecs(eparent_l)%vec)
                        nrows = size(A%dom(idom)%lblks(irow,1)%data_(iupper)%mat,1)
                        ncols = size(A%dom(idom)%lblks(irow,1)%data_(iupper)%mat,2)
                        call DGEMV('N', nrows, ncols, -ONE, A%dom(idom)%lblks(irow,1)%data_(iupper)%mat, nrows, z%dom(idom)%vecs(eparent_l)%vec, 1, ONE, z%dom(idom)%vecs(irow)%vec, 1)

                    else

                        ! Get associated parent block index
                        recv_comm    = A%dom(idom)%lblks(irow,1)%data_(iupper)%recv_comm
                        recv_domain  = A%dom(idom)%lblks(irow,1)%data_(iupper)%recv_domain
                        recv_element = A%dom(idom)%lblks(irow,1)%data_(iupper)%recv_element
                        !z%dom(idom)%vecs(irow)%vec = z%dom(idom)%vecs(irow)%vec - matmul(A%dom(idom)%lblks(irow,iupper)%mat, z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec)

                        nrows = size(A%dom(idom)%lblks(irow,1)%data_(iupper)%mat,1)
                        ncols = size(A%dom(idom)%lblks(irow,1)%data_(iupper)%mat,2)
                        call DGEMV('N', nrows, ncols, -ONE, A%dom(idom)%lblks(irow,1)%data_(iupper)%mat, nrows, z%recv%comm(recv_comm)%dom(recv_domain)%vecs(recv_element)%vec, 1, ONE, z%dom(idom)%vecs(irow)%vec, 1)

                    end if


                end do


                !
                ! Diagonal block
                !
                diag_irow = self%LD%dom(idom)%lblks(irow,1)%get_diagonal()
                z%dom(idom)%vecs(irow)%vec = matmul(self%LD%dom(idom)%lblks(irow,1)%data_(diag_irow)%mat, z%dom(idom)%vecs(irow)%vec)

            end do ! irow




        end do ! idom



        call self%timer%stop()

    end function apply
    !------------------------------------------------------------------------------------------












    !>  Send the matrix blocks from the current processor that couple with neighboring 
    !!  processors to those processors so they can be included in the ILU0 solve of the neighbor 
    !!  processor.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine comm_send(self,A)
        class(precon_RASILU0_t),    intent(inout)               :: self
        type(chidg_matrix_t),       intent(in), asynchronous    :: A

        integer(ik)                 :: icomm, idom_send, ielem_send, iblk_send, idom, ielem, &
                                       iblk, proc, nrows, ncols, idomain_g, ierr
        integer(ik), allocatable    :: send_blocks(:)
        type(mpi_request)           :: request


        call write_line('       RAS: sending...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))


        do icomm = 1,size(self%send%comm)
            proc = self%send%comm(icomm)%proc

            do idom_send = 1,size(self%send%comm(icomm)%dom)
                idom      = self%send%comm(icomm)%dom(idom_send)%idomain_l
                idomain_g = self%send%comm(icomm)%dom(idom_send)%idomain_g

                !
                ! Loop through element faces and find neighbors that are off-processor on 'proc' 
                ! and send overlap element data from chidg_matrix.
                !
                do ielem_send = 1,self%send%comm(icomm)%dom(idom_send)%elem_send%size()
                    ielem = self%send%comm(icomm)%dom(idom_send)%elem_send%at(ielem_send) 

                    !
                    ! Which blocks are being sent
                    !
                    send_blocks = self%send%comm(icomm)%dom(idom_send)%blk_send(ielem_send)%data()

                    do iblk_send = 1,size(send_blocks)
                        iblk = send_blocks(iblk_send)

                        nrows = size(A%dom(idom)%lblks(ielem,1)%data_(iblk)%mat,1)
                        ncols = size(A%dom(idom)%lblks(ielem,1)%data_(iblk)%mat,2)
                        call MPI_ISend(A%dom(idom)%lblks(ielem,1)%data_(iblk)%mat, nrows*ncols, MPI_REAL8, proc, idomain_g, ChiDG_COMM, request, ierr)

                        ! Store requests to be checked by MPI_Wait
                        call self%mpi_requests%push_back(request)

                    end do !iblk

                end do !ielem



            end do !idom_send
        end do !icomm




    end subroutine comm_send
    !******************************************************************************************
       







    !>  Receive matrix blocks from neighboring processors for neighboring elements.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine comm_recv(self)
        class(precon_RASILU0_t),   intent(inout)   :: self

        integer(ik), allocatable    :: send_blocks(:)
        integer(ik)                 :: icomm, idom, ielem, iblk, proc, nrows, ncols, idomain_g, &
                                       idom_recv, ielem_recv, iblk_recv, ierr

        type(mpi_request)           :: request

        call write_line('       RAS: receiving...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        do idom_recv = 1,size(self%recv%dom)

            do icomm = 1,size(self%recv%dom(idom_recv)%comm)
                proc = self%recv%dom(idom_recv)%comm(icomm)%proc

                do ielem_recv = 1,size(self%recv%dom(idom_recv)%comm(icomm)%elem)
                
                    do iblk_recv = 1,size(self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks)
                        nrows     = self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks(iblk_recv)%nrows_
                        ncols     = self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks(iblk_recv)%ncols_
                        idomain_g = self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks(iblk_recv)%dparent_g_

                        call MPI_Recv(self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks(iblk_recv)%mat, nrows*ncols, MPI_REAL8, proc, idomain_g, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                        !call MPI_IRecv(self%recv%dom(idom_recv)%comm(icomm)%elem(ielem_recv)%blks(iblk_recv)%mat, nrows*ncols, MPI_REAL8, proc, idomain_g, ChiDG_COMM, request, ierr)
                        !call self%mpi_requests%push_back(request)

                    end do !iblk_recv

                end do !ielem_recv

            end do !icomm

        end do !idom_recv



    end subroutine comm_recv
    !******************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/10/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine comm_wait(self)
        class(precon_RASILU0_t),    intent(inout)   :: self

        integer(ik) :: nwait, iwait, ierr

        call write_line('       RAS: waiting...', ltrim=.false., io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        nwait = self%mpi_requests%size()
        if (nwait > 0) then

            call MPI_Waitall(nwait, self%mpi_requests%data(1:nwait), MPI_STATUSES_IGNORE, ierr)
            call self%mpi_requests%clear()

        end if


    end subroutine comm_wait
    !******************************************************************************************























end module precon_RASILU0
