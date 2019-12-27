module type_chidg_adjointx
#include<messenger.h>
    use type_chidg_matrix,      only: chidg_matrix_t, chidg_matrix
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: ZERO, dX_DIFF
    use mod_io,                 only: backend
    use type_mesh,              only: mesh_t
    use type_nvector,           only: nvector_t
    use type_node,              only: node_t
    use type_storage_flags,     only: storage_flags_t
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mpi_f08,                only: MPI_BCast, MPI_ISend, MPI_Recv, MPI_INTEGER4,MPI_REAL8,   &
                                      MPI_ISend, MPI_AllReduce, MPI_statuses_ignore, MPI_Request, &
                                      MPI_Status, MPI_Status_size, MPI_Waitall, MPI_status_ignore
    use type_mpi_request_vector,only: mpi_request_vector_t
    implicit none


    !>  Storage for sensitivities of funtionals wrt grid nodes movements
    !!
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/26/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    type,   public  :: chidg_adjointx_t

        ! dR/dX matrix
        type(chidg_matrix_t)                    :: Rx
        
        ! Resultant vector from [v]*[dR/dx], v being adjoint variables realtive to each 
        ! functional at a particular istep
        type(chidg_vector_t),   allocatable     :: vRx(:)                   ! [ifunc]

        ! Functional partial/total derivatives wrt grid nodes movement at a praticular istep
        type(chidg_vector_t),   allocatable     :: Jx(:)                    ! [ifunc]

        ! Sum of functional total derivatives contribution at each istep
        type(chidg_vector_t),   allocatable     :: Jx_unsteady(:)           ! [ifunc]

        ! Sensitivities contribution from the auxiliary problem
        type(chidg_vector_t),   allocatable     :: wAx(:)                   ! [ifunc]

        ! List of nodes sensitivities for post-processing
        type(nvector_t),        allocatable     :: node_sensitivities(:,:)  ! [ifunc,idom_g]



        ! Initialization completed
        logical         ::  adjointx_initialized = .false.


    contains
        
        ! Initialization
        procedure       :: init
        procedure       :: init_containers

        ! Data process
        procedure       :: gather_all
        procedure       :: dump

        ! Clear 
        procedure       :: Jx_clear
        procedure       :: vRx_clear
        procedure       :: release

    end type chidg_adjointx_t
    !************************************************************************************


contains


    !>  Allocate storage for adjointx derivatives
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/3/2018
    !!
    !!  param[in]   nstep       Number of steps of the solution (nstep = 1 for steady and HB)
    !!  param[in]   nfunc       Number of functionals
    !!
    !------------------------------------------------------------------------------------
    subroutine init(self,nfunc,nstep,sflags)
        class(chidg_adjointx_t), intent(inout)   :: self
        integer(ik),             intent(in)      :: nstep
        integer(ik),             intent(in)      :: nfunc
        type(storage_flags_t),   intent(in)      :: sflags

        integer(ik)     :: ierr

        ! Deallocate any data previously allocated
        call self%release()

        ! Default, no error
        ierr = 0


        ! Initialize vRx array of chidg_vectors
        if (sflags%vRx) allocate( self%vRx(nfunc), stat=ierr )
        if (ierr/=0) call AllocationError

        ! Initialize Jx array of chidg_vectors
        if (sflags%Jx) allocate( self%Jx(nfunc), stat=ierr )
        if (ierr/=0) call AllocationError

        ! Initialize Junsteady array of chidg_vectors
        if (sflags%Jx_unsteady) allocate( self%Jx_unsteady(nfunc), stat=ierr )
        if (ierr/=0) call AllocationError
        
        ! Initialize array of sensitivities from auxiliary problem 
        if (sflags%wAx) allocate( self%wAx(nfunc), stat=ierr )
        if (ierr/=0) call AllocationError


        self%adjointx_initialized = .true.

    end subroutine init
    !************************************************************************************






    !>  Allocate storage necessary for adjointx chidg_vectors
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/3/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine init_containers(self,mesh,ntime,sflags)
        class(chidg_adjointx_t),    intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(inout)   :: ntime  
        type(storage_flags_t),      intent(in)      :: sflags
             
        integer(ik)                 :: ifunc,istep,nfuncs
!        integer(ik)                 :: specialization
        character(2)                :: mtype

        
!        specialization = dX_DIFF
        

        if (allocated(self%Jx)) then
            nfuncs = size(self%Jx)
        else
            nfuncs = 0
        end if

        ! Initialize vectors
        do ifunc = 1,nfuncs 
            !if (sflags%vRx)         call self%vRx(ifunc)%init(mesh,ntime,specialization)
            !if (sflags%Jx)          call self%Jx(ifunc)%init(mesh,ntime,specialization)
            !if (sflags%Jx_unsteady) call self%Jx_unsteady(ifunc)%init(mesh,ntime,specialization)

            if (sflags%vRx)         self%vRx(ifunc)         = chidg_vector(trim(backend))
            if (sflags%Jx)          self%Jx(ifunc)          = chidg_vector(trim(backend))
            if (sflags%Jx_unsteady) self%Jx_unsteady(ifunc) = chidg_vector(trim(backend))

            if (sflags%vRx)         call self%vRx(ifunc)%init(mesh,ntime,'grid differentiation')
            if (sflags%Jx)          call self%Jx(ifunc)%init(mesh,ntime,'grid differentiation')
            if (sflags%Jx_unsteady) call self%Jx_unsteady(ifunc)%init(mesh,ntime,'grid differentiation')

        end do !ifunc


        if (sflags%Rx) then
            ! Initialize matrix
            self%Rx = chidg_matrix(trim(backend))
            mtype = 'dX'
            call self%Rx%init(mesh,mtype)
            call self%Rx%init_recv(self%vRx(1))

            ! Set Rx transpose flag on for matrix-vector product
            if (sflags%Rx_trans) self%Rx%transposed = .true.
        
        end if

    end subroutine init_containers
    !************************************************************************************








    !>  Gather all node sensitivities for post-processing
    !!  The master proc receives and manage the post-processing 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/25/2018
    !!
    !!  TODO: communication can be improved
    !!
    !------------------------------------------------------------------------------------
    subroutine gather_all(self,mesh)
        class(chidg_adjointx_t),    intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
             
        integer(ik)                     :: ifunc, istep, ndomains_g, nfuncs, idom, ielem, inode,   &
                                           idomain_g, node_ID_l, coords_system, node_index,     &
                                           itag, ntags, ierr, iproc, n_handles
        real(rk)                        :: node_coords(3), node_sens(3)
        type(node_t)                    :: temp_node
        type(mpi_request_vector_t)      :: isend_requests
        type(MPI_Request)               :: irequest0, irequest1, irequest2, irequest3, irequest4, irequest5
        type(MPI_status), allocatable   :: isend_status(:)
        
        print*, 'gather_all - 1'

        ! Find the max domain global index
        !ndoms_g = size(mesh%npoints,1)
        ndomains_g = mesh%ndomains_g()


        ! Deduce number of functionals computed
        nfuncs = size(self%Jx_unsteady)

        print*, 'gather_all - 2'

        ! Allocate storage for node sensitivities
        if (allocated(self%node_sensitivities)) deallocate (self%node_sensitivities)
        allocate(self%node_sensitivities(nfuncs,ndomains_g), stat=ierr)
        if (ierr/=0) call AllocationError

        print*, 'gather_all - 3'

        ! Loop trhough the functionals and gather all sensitivities
        ! This is done by each processor individually.
        do ifunc = 1,nfuncs

            !! Loop through the local nodes and gather sensitivities locally
            !do idom = 1,size(self%Jx_unsteady(ifunc)%dom)
            !    idomain_g    = mesh%domain(idom)%idomain_g
            !    do ielem = 1,size(self%Jx_unsteady(ifunc)%dom(idom)%vecs)
            ! Loop through the local nodes and gather sensitivities locally
            do idom = 1,mesh%ndomains()
                idomain_g = mesh%domain(idom)%idomain_g
                do ielem = 1,mesh%domain(idom)%nelements()

                    ! NOTE: Jx has been initialized with specialization 'dX' this means that
                    !       nterms_ attribute of the densevector corresponds to nnodes_r 
                    !       and nvars_ attribute to 3 (3 directions). Therefore, nterms()
                    !       actually return the number of support/reference node of the 
                    !       element
                    !do inode = 1,self%Jx_unsteady(ifunc)%dom(idom)%vecs(ielem)%nterms()
                    do inode = 1,mesh%domain(idom)%elems(ielem)%nterms_c

        print*, 'gather_all - 4'
                        ! Find support node coordinates, coordinate system and connectivity (ie local node ID)
                        node_coords(1) = mesh%domain(idom)%elems(ielem)%node_coords(inode,1)
                        node_coords(2) = mesh%domain(idom)%elems(ielem)%node_coords(inode,2)
                        node_coords(3) = mesh%domain(idom)%elems(ielem)%node_coords(inode,3)
                        node_ID_l      = mesh%domain(idom)%elems(ielem)%connectivity(inode)
                        coords_system  = mesh%domain(idom)%elems(ielem)%coordinate_system
                        
                        ! Read sensitivities in the chidg vector
                        node_sens(1) = self%Jx_unsteady(ifunc)%dom(idom)%vecs(ielem)%getterm(1,inode,1)
                        node_sens(2) = self%Jx_unsteady(ifunc)%dom(idom)%vecs(ielem)%getterm(2,inode,1)
                        node_sens(3) = self%Jx_unsteady(ifunc)%dom(idom)%vecs(ielem)%getterm(3,inode,1)
                        
                        ! Find out if the same node has been already added
                        node_index = self%node_sensitivities(ifunc,idomain_g)%search_by_coords(node_coords,node_ID_l)

                        if (node_index == 0) then
                            ! The node is new. Add it
                            call temp_node%init_node(node_ID_l,idomain_g,coords_system,node_coords,node_sens)
                            call self%node_sensitivities(ifunc,idomain_g)%push_back(temp_node)
                        else
                            ! The node has been already added, add sensitivities contribution to it.
                            call self%node_sensitivities(ifunc,idomain_g)%data(node_index)%add_sensitivities(node_sens)
                        end if

        print*, 'gather_all - 5'
                    end do !inode
                end do !ielem
            end do !idom
        

            ! All the processors send the node data to the master proc
            if (IRANK /= GLOBAL_MASTER) then
                
        print*, 'gather_all - 6'
                ! Count over all information to send and allocate vector of requests
                n_handles = 1
                do idom = 1,size(self%Jx_unsteady(ifunc)%dom)
                    idomain_g = mesh%domain(idom)%idomain_g
                    n_handles = n_handles + 5*self%node_sensitivities(ifunc,idomain_g)%size()
                end do   

                ! Start sending messages
                itag = 1
                call MPI_ISend(n_handles,1,MPI_INTEGER4,GLOBAL_MASTER,itag,ChiDG_COMM,irequest0,ierr)
                call isend_requests%push_back(irequest0)

                do idom = 1,size(self%Jx_unsteady(ifunc)%dom)
                    idomain_g = mesh%domain(idom)%idomain_g
                    do inode = 1,self%node_sensitivities(ifunc,idomain_g)%size()
                        
                        itag = itag + 1
                        call MPI_ISend(self%node_sensitivities(ifunc,idomain_g)%data(inode)%node_ID_l        ,1,MPI_INTEGER4,GLOBAL_MASTER,itag,ChiDG_COMM,irequest1,ierr)
                        itag = itag + 1
                        call MPI_ISend(self%node_sensitivities(ifunc,idomain_g)%data(inode)%domain_g         ,1,MPI_INTEGER4,GLOBAL_MASTER,itag,ChiDG_COMM,irequest2,ierr)
                        itag = itag + 1
                        call MPI_ISend(self%node_sensitivities(ifunc,idomain_g)%data(inode)%coordinate_system,1,MPI_INTEGER4,GLOBAL_MASTER,itag,ChiDG_COMM,irequest3,ierr)
                        itag = itag + 1
                        call MPI_ISend(self%node_sensitivities(ifunc,idomain_g)%data(inode)%coords           ,3,MPI_REAL8   ,GLOBAL_MASTER,itag,ChiDG_COMM,irequest4,ierr)
                        itag = itag + 1
                        call MPI_ISend(self%node_sensitivities(ifunc,idomain_g)%data(inode)%sensitivities    ,3,MPI_REAL8   ,GLOBAL_MASTER,itag,ChiDG_COMM,irequest5,ierr)

                        call isend_requests%push_back(irequest1)
                        call isend_requests%push_back(irequest2)
                        call isend_requests%push_back(irequest3)
                        call isend_requests%push_back(irequest4)
                        call isend_requests%push_back(irequest5)

                    end do !inode
                end do !idom

            end if ! All procs but GLOBAL MASTER

        print*, 'gather_all - 7'

            ! Synchronize all procs, otherwise the master proc will go on and look for new messages
            call MPI_Barrier(ChiDG_COMM,ierr)



            ! Master processor receives and stores the node data from all the processors
            if (IRANK == GLOBAL_MASTER) then
                
                do iproc = 0,NRANK-1
                    if (iproc /= IRANK) then

                        call MPI_Recv(ntags,1,MPI_INTEGER4,iproc,1,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                        
                        itag = 2
                        do while (itag .le. ntags)

                            call MPI_Recv(node_ID_l    ,1,MPI_INTEGER4,iproc,itag,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                            itag = itag + 1
                            call MPI_Recv(idomain_g    ,1,MPI_INTEGER4,iproc,itag,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                            itag = itag + 1
                            call MPI_Recv(coords_system,1,MPI_INTEGER4,iproc,itag,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                            itag = itag + 1
                            call MPI_Recv(node_coords  ,3,MPI_REAL8   ,iproc,itag,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                            itag = itag + 1
                            call MPI_Recv(node_sens    ,3,MPI_REAL8   ,iproc,itag,ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                            itag = itag + 1

                            
                            ! Define coordinate system

                            ! Find if the same node as been already added
                            node_index = self%node_sensitivities(ifunc,idomain_g)%search_by_coords(node_coords,node_ID_l)

                            if (node_index == 0) then
                                ! The node is new. Add it
                                call temp_node%init_node(node_ID_l,idomain_g,coords_system,node_coords,node_sens)
                                call self%node_sensitivities(ifunc,idomain_g)%push_back(temp_node)

                            else
                                ! The node has been already added, add sensitivities contribution.
                                call self%node_sensitivities(ifunc,idomain_g)%data(node_index)%add_sensitivities(node_sens)

                            end if
                    
                        end do ! itag

                    end if !iproc /= IRANK
                end do !iproc
                
            end if ! GLOBAL MASTER

        print*, 'gather_all - 8'
            ! Synchronize all procs
            call MPI_Barrier(ChiDG_COMM,ierr)

            ! The sender processors (IRANK/=GLOBAL_MASTER) check the send request status
            ! and return an error if any of the request faild.
            if (isend_requests%size() > 0) then
                if (allocated(isend_status)) deallocate(isend_status)
                allocate(isend_status(isend_requests%size()),stat=ierr)
                call MPI_Waitall(n_handles,isend_requests%data,isend_status,ierr)
            end if

        print*, 'gather_all - 9'

            ! Master processor reorder the nodes for each each block in node_index order 
            if (IRANK == GLOBAL_MASTER) then
                do idomain_g = 1,ndomains_g
                    call self%node_sensitivities(ifunc,idomain_g)%reorder_by_index()
                end do
            end if

        print*, 'gather_all - 10'

            ! Synchronize all procs
            call MPI_Barrier(ChiDG_COMM,ierr)

            ! Clear send request buffer for next functional
            call isend_requests%clear()

        end do !ifunc
        
        !
        ! Optional, dump sensitivities for a given (domain,node,ifunc*)
        !
        !call self%dump(1,2)

        print*, 'gather_all - 11'

    end subroutine gather_all
    !************************************************************************************







    !>  Dump node info for a given id. Only global master. 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   10/14/2019
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine dump(self,idom,inode,ifunc) 
        class(chidg_adjointx_t),     intent(inout)   :: self
        integer(ik),                 intent(in)      :: idom
        integer(ik),                 intent(in)      :: inode
        integer(ik), optional,       intent(in)      :: ifunc

        integer(ik)     :: idomain_g,ndoms_g,i,ierr

        if (IRANK == GLOBAL_MASTER) then

            if (present(ifunc)) then
                
                print*, '----------------------------------------------------------'
                print*, 'Functional ID      ', ifunc
                print*, 'Node ID:           ', self%node_sensitivities(ifunc,idom)%data(inode)%node_ID_l
                print*, 'Node global domain:', self%node_sensitivities(ifunc,idom)%data(inode)%domain_g
                print*, 'Node coord system: ', self%node_sensitivities(ifunc,idom)%data(inode)%coordinate_system
                print*, 'Node coordinates:  ', self%node_sensitivities(ifunc,idom)%data(inode)%coords
                print*, 'Node sensitivities:', self%node_sensitivities(ifunc,idom)%data(inode)%sensitivities
                print*, '----------------------------------------------------------'
            
            else
                
                do i = 1,size(self%node_sensitivities(:,idom))
                    print*, '----------------------------------------------------------'
                    print*, 'Functional ID      ', i
                    print*, 'Node ID:           ', self%node_sensitivities(i,idom)%data(inode)%node_ID_l
                    print*, 'Node global domain:', self%node_sensitivities(i,idom)%data(inode)%domain_g
                    print*, 'Node coord system: ', self%node_sensitivities(i,idom)%data(inode)%coordinate_system
                    print*, 'Node coordinates:  ', self%node_sensitivities(i,idom)%data(inode)%coords
                    print*, 'Node sensitivities:', self%node_sensitivities(i,idom)%data(inode)%sensitivities
                    print*, '----------------------------------------------------------'
                end do
            end if

        end if
            

    end subroutine dump 
    !************************************************************************************





    !>  Set all Jx vectors to ZERO 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine Jx_clear(self)
        class(chidg_adjointx_t),     intent(inout)   :: self
        
        integer(ik)     :: ifunc 

        do ifunc = 1,size(self%Jx)
            call self%Jx(ifunc)%clear()
        end do

    end subroutine Jx_clear
    !************************************************************************************







    !>  Set all vRx vectors to ZERO 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine vRx_clear(self)
        class(chidg_adjointx_t),     intent(inout)   :: self
       
        integer(ik)     :: ifunc 

        do ifunc = 1,size(self%Jx)
            call self%vRx(ifunc)%clear()
        end do

    end subroutine vRx_clear
    !************************************************************************************






    !>  Release allocated data 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_adjointx_t),     intent(inout)   :: self

        integer(ik) :: i
        
        if (allocated(self%vRx))                deallocate (self%vRx)
        if (allocated(self%Jx))                 deallocate (self%Jx)
        if (allocated(self%Jx_unsteady))        deallocate (self%Jx_unsteady)
        if (allocated(self%wAx))                deallocate (self%wAx)
        if (allocated(self%node_sensitivities)) deallocate (self%node_sensitivities)

        call self%Rx%release()

    end subroutine release
    !************************************************************************************




end module type_chidg_adjointx
