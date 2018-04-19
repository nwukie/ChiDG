module type_chidg_matrix
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: NO_ID, ZERO
    use mod_chidg_mpi,          only: IRANK
    use type_densematrix,       only: densematrix_t
    use type_domain_matrix,     only: domain_matrix_t
    use type_mesh,              only: mesh_t
    use type_face_info,         only: face_info_t
    use type_seed,              only: seed_t
    use type_chidg_vector,      only: chidg_vector_t
    use DNAD_D
    implicit none




    !>  ChiDG matrix type. Contains an array of domain_matrix_t types, each corresponding to a 
    !!  domain.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  TODO: TEST self%stamp component and behavior
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chidg_matrix_t

        type(domain_matrix_t), allocatable    :: dom(:) ! Array of domain-matrices. One for each domain

        logical     :: local_initialized = .false.      ! Has the matrix processor-local data been initialized
        logical     :: recv_initialized  = .false.      ! Has matrix been initialized with information about chidg_vector%recv

        integer     :: stamp(8) = NO_ID                 ! Stamp from date_and_time that gets updated when store routines are called

    contains
        ! Initializers
        generic,    public  :: init => initialize
        procedure,  private :: initialize               ! chidg_matrix initialization

        procedure, public   :: init_recv                ! Initialize with information about chidg_vector%recv for mv multiply

        ! Setters
        procedure   :: store                            ! Store interior coupling
        procedure   :: store_chimera                    ! Store chimera coupling
        procedure   :: store_bc                         ! Store boundary condition coupling
        procedure   :: store_hb                         ! Store linearization data for coupling across harmonic balance levels 
        procedure   :: clear                            ! Zero matrix-values

        ! Processors
        procedure   :: restrict

        procedure   :: get_ntime
        procedure   :: ndomains
        procedure   :: to_file                          ! Write matrix to file
        procedure   :: to_real                          ! Return the matrix stored as real

        procedure   :: release
        final       :: destructor

    end type chidg_matrix_t
    !*****************************************************************************************



    private
contains




    !>  Subroutine for initializing chidg_matrix_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  domains     Array of domain_t instances
    !!  
    !!
    !------------------------------------------------------------------------------------------
    subroutine initialize(self,mesh,mtype)
        class(chidg_matrix_t),  intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh
        character(*),           intent(in)      :: mtype

        integer(ik) :: ierr, ndomains, idom

        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        if (allocated(self%dom)) deallocate(self%dom)

        ! Allocate domain_matrix_t for each domain
        ndomains = mesh%ndomains()
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Call initialization procedure for each domain_matrix_t
        do idom = 1,ndomains
             call self%dom(idom)%init(mesh,idom,mtype=mtype)
        end do

        ! Set initialization to true
        self%local_initialized = .true.

    end subroutine initialize
    !******************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_recv(self,x)
        class(chidg_matrix_t),   intent(inout)   :: self
        type(chidg_vector_t),    intent(in)      :: x

        integer(ik) :: idom, ielem, iblk, itime, matrix_proc, vector_proc, comm_proc
        integer(ik) :: dparent_g, eparent_g, parent_proc, icomm, idom_recv, ielem_recv, &
                       drecv_g, erecv_g, recv_domain, recv_elem, imat
        logical     :: local_multiply, parallel_multiply, match_found

        
        !
        ! Loop through INTERIOR coupling and look for parallel multiply
        !
        do idom = 1,size(self%dom)
            do ielem = 1,size(self%dom(idom)%lblks,1)
                do itime = 1,size(self%dom(idom)%lblks,2)
                    do imat = 1,self%dom(idom)%lblks(ielem,itime)%size()

                        matrix_proc = IRANK
                        vector_proc = self%dom(idom)%lblks(ielem,itime)%parent_proc(imat)

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )


                        if ( parallel_multiply ) then
                            !
                            ! Get information about element we need to multiply with
                            !
                            dparent_g   = self%dom(idom)%lblks(ielem,itime)%dparent_g(imat)
                            eparent_g   = self%dom(idom)%lblks(ielem,itime)%eparent_g(imat)
                            parent_proc = self%dom(idom)%lblks(ielem,itime)%parent_proc(imat)



                            !
                            ! Loop through chidg_vector%recv to find match
                            !
                            match_found = .false.
                            do icomm = 1,size(x%recv%comm)

                                comm_proc = x%recv%comm(icomm)%proc

                                if ( comm_proc == parent_proc ) then
                                    do idom_recv = 1,size(x%recv%comm(icomm)%dom)
                                        do ielem_recv = 1,size(x%recv%comm(icomm)%dom(idom_recv)%vecs)

                                            ! Get recv element indices
                                            drecv_g = x%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%dparent_g()
                                            erecv_g = x%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%eparent_g()

                                            ! If they match the domain_matrix, set the recv indices so chidg_mv knows how to compute matrix-vector product
                                            if ( (drecv_g == dparent_g) .and. (erecv_g == eparent_g) ) then
                                                call self%dom(idom)%lblks(ielem,itime)%set_recv_comm(imat,icomm)
                                                call self%dom(idom)%lblks(ielem,itime)%set_recv_domain(imat,idom_recv)
                                                call self%dom(idom)%lblks(ielem,itime)%set_recv_element(imat,ielem_recv)
                                                match_found = .true.
                                            end if

                                        end do !ielem_recv
                                    end do !idom_recv
                                end if

                            end do ! icomm

                            if (.not. match_found) call chidg_signal(FATAL,"chidg_matrix%init_recv: no matching recv element found in vector")



                        end if



                    end do !imat
                end do !itime
            end do !ielem
        end do ! idom









        !
        ! Loop through CHIMERA coupling and look for parallel multiply
        !
        do idom = 1,size(self%dom)

            if (allocated(self%dom(idom)%chi_blks)) then
                do ielem = 1,size(self%dom(idom)%chi_blks,1)
                    do itime = 1,size(self%dom(idom)%chi_blks,2)
                        do imat = 1,self%dom(idom)%chi_blks(ielem,itime)%size()
                        
                            matrix_proc = IRANK
                            vector_proc = self%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                !
                                ! Get information about element we need to multiply with
                                !
                                dparent_g   = self%dom(idom)%chi_blks(ielem,itime)%dparent_g(imat)
                                eparent_g   = self%dom(idom)%chi_blks(ielem,itime)%eparent_g(imat)
                                parent_proc = self%dom(idom)%chi_blks(ielem,itime)%parent_proc(imat)



                                !
                                ! Loop through chidg_vector%recv to find match
                                !
                                match_found = .false.
                                do icomm = 1,size(x%recv%comm)

                                    comm_proc = x%recv%comm(icomm)%proc

                                    if ( comm_proc == parent_proc ) then
                                        do idom_recv = 1,size(x%recv%comm(icomm)%dom)
                                            
                                            recv_domain = x%recv%comm(icomm)%dom(idom_recv)%vecs(1)%dparent_g()
                                            if ( recv_domain == dparent_g ) then

                                            do ielem_recv = 1,size(x%recv%comm(icomm)%dom(idom_recv)%vecs)

                                                ! Get recv element indices
                                                recv_elem = x%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%eparent_g()

                                    

                                                ! If they match the domain_matrix, set the recv indices so chidg_mv knows how to compute matrix-vector product
                                                if ( recv_elem == eparent_g )  then
                                                    call self%dom(idom)%chi_blks(ielem,itime)%set_recv_comm(imat,icomm)
                                                    call self%dom(idom)%chi_blks(ielem,itime)%set_recv_domain(imat,idom_recv)
                                                    call self%dom(idom)%chi_blks(ielem,itime)%set_recv_element(imat,ielem_recv)
                                                    match_found = .true.
                                                    exit
                                                end if

                                            end do !ielem_recv

                                            end if ! recv_domain == dparent

                                            if (match_found) exit
                                        end do !idom_recv
                                    end if

                                    if (match_found) exit

                                end do ! icomm

                                if (.not. match_found) call chidg_signal(FATAL,"chidg_matrix%init_recv: no matching recv element found in vector")

                            end if


                        end do !imat
                    end do !itime
                end do !ielem
            end if

        end do ! idom






        !
        ! Loop through BOUNDARY coupling and look for parallel multiply
        !
        do idom = 1,size(self%dom)

            if (allocated(self%dom(idom)%bc_blks)) then
                do ielem = 1,size(self%dom(idom)%bc_blks,1)
                    do itime = 1,size(self%dom(idom)%bc_blks,2)
                        do imat = 1,self%dom(idom)%bc_blks(ielem,itime)%size()
                        
                            matrix_proc = IRANK
                            vector_proc = self%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                !
                                ! Get information about element we need to multiply with
                                !
                                dparent_g   = self%dom(idom)%bc_blks(ielem,itime)%dparent_g(imat)
                                eparent_g   = self%dom(idom)%bc_blks(ielem,itime)%eparent_g(imat)
                                parent_proc = self%dom(idom)%bc_blks(ielem,itime)%parent_proc(imat)



                                !
                                ! Loop through chidg_vector%recv to find match
                                !
                                match_found = .false.
                                do icomm = 1,size(x%recv%comm)

                                    comm_proc = x%recv%comm(icomm)%proc

                                    if ( comm_proc == parent_proc ) then
                                        do idom_recv = 1,size(x%recv%comm(icomm)%dom)
                                            
                                            recv_domain = x%recv%comm(icomm)%dom(idom_recv)%vecs(1)%dparent_g()
                                            if ( recv_domain == dparent_g ) then

                                            do ielem_recv = 1,size(x%recv%comm(icomm)%dom(idom_recv)%vecs)

                                                ! Get recv element indices
                                                recv_elem = x%recv%comm(icomm)%dom(idom_recv)%vecs(ielem_recv)%eparent_g()

                                    

                                                ! If they match the domain_matrix, set the recv indices so chidg_mv can compute m-v product
                                                if ( recv_elem == eparent_g )  then
                                                    call self%dom(idom)%bc_blks(ielem,itime)%set_recv_comm(imat,icomm)
                                                    call self%dom(idom)%bc_blks(ielem,itime)%set_recv_domain(imat,idom_recv)
                                                    call self%dom(idom)%bc_blks(ielem,itime)%set_recv_element(imat,ielem_recv)
                                                    match_found = .true.
                                                    exit
                                                end if

                                            end do !ielem_recv

                                            end if ! recv_domain == dparent

                                            if (match_found) exit
                                        end do !idom_recv
                                    end if

                                    if (match_found) exit

                                end do ! icomm

                                if (.not. match_found) call chidg_signal(FATAL,"chidg_matrix%init_recv: no matching recv element found in vector")

                            end if


                        end do !imat
                    end do !itime
                end do !ielem
            end if

        end do ! idom


        ! Set recv initialization to true
        self%recv_initialized = .true.

    end subroutine init_recv
    !*******************************************************************************************








    !> Procedure for storing linearization information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial 
    !!                          derivatives for the linearization matrix.
    !!  @param[in]  face_info   Information about where the coupling was computed and whom it
    !!                          was computed with respect to.
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed.
    !!
    !------------------------------------------------------------------------------------------
    subroutine store(self,integral,face_info,seed,ivar,itime)
        class(chidg_matrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        type(face_info_t),      intent(in)      :: face_info
        type(seed_t),           intent(in)      :: seed
        integer(ik),            intent(in)      :: ivar 
        integer(ik),            intent(in)      :: itime

        integer(ik) :: idomain_l

        idomain_l = face_info%idomain_l

        ! Store linearization in associated domain domain_matrix_t
        call self%dom(idomain_l)%store(integral,face_info,seed,ivar,itime)

        ! Update stamp
        call date_and_time(values=self%stamp)

    end subroutine store
    !*******************************************************************************************






    !> Procedure for stiring linearization information for Chimera faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial 
    !!                          derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing the indices defining the Chimera face
    !!  @param[in]  seed        seed_t containing the indices defining the element against 
    !!                          which the Chimera face was linearized
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !------------------------------------------------------------------------------------------
    subroutine store_chimera(self,integral,face_info,seed,ivar,itime)
        class(chidg_matrix_t),      intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar 
        integer(ik),                intent(in)      :: itime

        integer(ik) :: idomain_l

        idomain_l = face_info%idomain_l

        ! Store linearization in associated domain domain_matrix_t
        call self%dom(idomain_l)%store_chimera(integral,face_info,seed,ivar,itime)

        ! Update stamp
        call date_and_time(values=self%stamp)

    end subroutine store_chimera
    !******************************************************************************************









    !>  Procedure for storing linearization information for boundary condition faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial
    !!                          derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing the indices defining the Chimera face
    !!  @param[in]  seed        seed_t containing the indices defining the element against 
    !!                          which the Chimera face was linearized
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !------------------------------------------------------------------------------------------
    subroutine store_bc(self,integral,face_info,seed,ivar,itime)
        class(chidg_matrix_t),      intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar 
        integer(ik),                intent(in)      :: itime

        integer(ik) :: idomain_l

        idomain_l = face_info%idomain_l

        ! Store linearization in associated domain domain_matrix_t
        call self%dom(idomain_l)%store_bc(integral,face_info,seed,ivar,itime)

        ! Update stamp
        call date_and_time(values=self%stamp)

    end subroutine store_bc
    !*******************************************************************************************





    !>  Procedure for storing linearization information for harmonic-balance cross-timelevel
    !!  coupling.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2018
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial
    !!                          derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing the indices defining the Chimera face
    !!  @param[in]  seed        seed_t containing the indices defining the element against 
    !!                          which the Chimera face was linearized
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !------------------------------------------------------------------------------------------
    subroutine store_hb(self,integral,face_info,seed,ivar,itime)
        class(chidg_matrix_t),      intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar 
        integer(ik),                intent(in)      :: itime

        integer(ik) :: idomain_l

        idomain_l = face_info%idomain_l

        ! Store linearization in associated domain domain_matrix_t
        call self%dom(idomain_l)%store_hb(integral,face_info,seed,ivar,itime)

        ! Update stamp
        call date_and_time(values=self%stamp)

    end subroutine store_hb
    !*******************************************************************************************






    !> Set all chidg_matrix matrix-values to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !! 
    !----------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidg_matrix_t),   intent(inout)   :: self

        integer(ik) :: idom
    
        ! Call domain_matrix_t%clear() on all matrices
        do idom = 1,size(self%dom)
           call self%dom(idom)%clear() 
        end do
    
    end subroutine clear
    !**********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/23/2017
    !!
    !!
    !----------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(chidg_matrix_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(chidg_matrix_t)    :: restricted
        integer(ik)             :: ierr, idom

        ! Allocate storage
        allocate(restricted%dom(size(self%dom)), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Copy restricted versions of domain_matrix_t's
        do idom = 1,size(self%dom)
            restricted%dom(idom) = self%dom(idom)%restrict(nterms_r)
        end do !idom

        restricted%local_initialized = self%local_initialized
        restricted%recv_initialized  = self%recv_initialized
        restricted%stamp             = self%stamp

    end function restrict
    !**********************************************************************************









    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   10/05/2017
    !!
    !!  TODO: Put in checks for ntime with relation to different blocks
    !!
    !----------------------------------------------------------------------------------
    function get_ntime(self) result(ntime)
        class(chidg_matrix_t),  intent(in)  :: self

        integer(ik)         :: ntime

        ! Get ntime from densematrix vector array of the 1st domain
        ntime = size(self%dom(1)%lblks,2)

    end function get_ntime
    !**********************************************************************************







    !>  Return matrix as a single real allocation.
    !!
    !!  *** WARNING ***
    !!  only use this for matrices that are known to be small. Calling this for a large system
    !!  will results in very large memory usage.
    !!  *** WARNING ***
    !!
    !!  WARNING: 
    !!      - ASSUMES 1 DOMAIN
    !!      - ASSUMES ALL ELEMENT HAVE SAME ORDER AND NUMBER OF FIELDS.
    !!      - ASSUMES 1 TIME INSTANCE
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine to_real(self,real_matrix)
        class(chidg_matrix_t),      intent(in)      :: self
        real(rk),   allocatable,    intent(inout)   :: real_matrix(:,:)

        type(densematrix_t)     :: dmat
        integer(ik)             :: idom, itime, ielem, iblk, handle, idomain_g, ielement_g, irow, icol, &
                                   nelements, row_offset, col_offset, irow_start, icol_start, block_size, nblk, ierr, ndof


        if (allocated(self%dom)) then

        ! Compute total number of dof's
        ndof = 0
        do idom = 1,self%ndomains()
            do itime = 1,self%dom(idom)%ntime()
                do ielem = 1,self%dom(idom)%nelements()
                    dmat = self%dom(idom)%lblks(ielem,itime)%at(1)
                    ndof = ndof + size(dmat%mat,1)
                end do
            end do
        end do


        ! Allocate storage for LHS matrix
        allocate(real_matrix(ndof,ndof), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Store matrix entries in LHS
        real_matrix = ZERO
        do idom = 1,self%ndomains()
            do itime = 1,self%dom(idom)%ntime()
                do ielem = 1,self%dom(idom)%nelements()
                    do iblk = 1,self%dom(idom)%lblks(ielem,itime)%size()

                        nblk = self%dom(idom)%lblks(ielem,itime)%size()
                        dmat = self%dom(idom)%lblks(ielem,itime)%at(iblk)
                        block_size = size(dmat%mat,1)


                        ! Coupled element
                        idomain_g  = dmat%dparent_g()
                        ielement_g = dmat%eparent_g()

                        
                        ! Compute global indices for the current block index (1,1)
                        irow_start = 1  +  block_size*(ielem-1)
                        icol_start = 1  +  block_size*(ielement_g-1)

                        
                        ! Compute offset indices
                        do col_offset = 1,block_size
                            do row_offset = 1,block_size

                                irow = irow_start + (row_offset-1)
                                icol = icol_start + (col_offset-1)
                                real_matrix(irow,icol) = dmat%mat(row_offset,col_offset)

                            end do
                        end do 

                    end do !iblk
                end do !ielem
            end do !itime
        end do !idom

        end if !allocated

    end subroutine to_real
    !**************************************************************************






    !>  Write matrix to file.
    !!
    !!  WARNING: 
    !!      - ASSUMES 1 DOMAIN
    !!      - ASSUMES ALL ELEMENT HAVE SAME ORDER AND NUMBER OF FIELDS.
    !!      - ASSUMES 1 TIME INSTANCE
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/27/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine to_file(self,filename)
        class(chidg_matrix_t),  intent(in)  :: self
        character(*),           intent(in)  :: filename

        type(densematrix_t)     :: dmat
        integer(ik)             :: idom, itime, ielem, iblk, handle, idomain_g, ielement_g, irow, icol, &
                                   nelements, row_offset, col_offset, irow_start, icol_start, block_size, nblk

        ! Open file
        open(newunit=handle, file=trim(filename))

        ! Write format
        call write_line('(irow, icol, nblk, value)')


        do idom = 1,self%ndomains()
            do itime = 1,self%dom(idom)%ntime()
                do ielem = 1,self%dom(idom)%nelements()
                    do iblk = 1,self%dom(idom)%lblks(ielem,itime)%size()

                        nblk = self%dom(idom)%lblks(ielem,itime)%size()
                        dmat = self%dom(idom)%lblks(ielem,itime)%at(iblk)
                        block_size = size(dmat%mat,1)


                        ! Coupled element
                        idomain_g  = dmat%dparent_g()
                        ielement_g = dmat%eparent_g()

                        
                        ! Compute global indices for the current block index (1,1)
                        irow_start = 1  +  block_size*(ielem-1)
                        icol_start = 1  +  block_size*(ielement_g-1)

                        
                        ! Compute offset indices
                        do col_offset = 1,block_size
                            do row_offset = 1,block_size

                                irow = irow_start + (row_offset-1)
                                icol = icol_start + (col_offset-1)
                                call write_line(irow,icol,nblk,dmat%mat(row_offset,col_offset), columns=.true., column_width=25, delimiter=',', ltrim=.true., handle=handle)

                            end do
                        end do 

                    end do !iblk
                end do !ielem
            end do !itime
        end do !idom



        ! Close file
        close(handle)

    end subroutine to_file
    !**************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/27/2017
    !!
    !------------------------------------------------------------------------------
    function ndomains(self) result(ndomains_)
       class(chidg_matrix_t),   intent(in)  :: self
       
       integer(ik)  :: ndomains_
       
       ndomains_ = size(self%dom) 

    end function ndomains
    !******************************************************************************








    !>  Release allocated resources.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !----------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_matrix_t),  intent(inout)   :: self

        if (allocated(self%dom)) deallocate(self%dom)

    end subroutine release
    !**********************************************************************************





    !>  chidg_matrix destructor.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidg_matrix_t),    intent(inout)   :: self

    end subroutine destructor
    !*********************************************************************************





end module type_chidg_matrix
