module type_chidg_matrix
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_chidg_mpi,          only: IRANK
    use type_blockmatrix,       only: blockmatrix_t
    use type_mesh,              only: mesh_t
    use type_mesh_new,          only: mesh_new_t
    use type_face_info,         only: face_info_t
    use type_seed,              only: seed_t
    use type_chidg_vector,       only: chidg_vector_t
    use DNAD_D
    implicit none




    !>  ChiDG matrix type. Contains an array of blockmatrix_t types, each corresponding to a 
    !!  domain.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chidg_matrix_t

        type(blockmatrix_t), allocatable    :: dom(:)                       !< Array of block-matrices. One for each domain

        logical                             :: local_initialized = .false.  !< Has the matrix processor-local data been initialized
        logical                             :: recv_initialized  = .false.  !< Has matrix been initialized with information about chidgVector%recv

    contains
        ! Initializers
        generic,    public  :: init => initialize
        procedure,  private :: initialize                   !< ChiDGMatrix initialization

        procedure, public   :: init_recv                    !< Initialize with information about chidgVector%recv for mv multiply

        ! Setters
        procedure   :: store                                !< Store linearization data for local blocks
        procedure   :: store_chimera                        !< Store linearization data for chimera blocks
        procedure   :: store_bc                             !< Store linearization data for boundary condition blocks
        procedure   :: clear                                !< Zero matrix-values


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
    !subroutine initialize(self,mesh,bcset_coupling,mtype)
    !subroutine initialize(self,mesh,bc,mtype)
    subroutine initialize(self,mesh,mtype)
        class(chidg_matrix_t),   intent(inout)          :: self
        type(mesh_new_t),        intent(in)             :: mesh
        character(*),           intent(in)              :: mtype
        !type(mesh_t),           intent(in)              :: mesh(:)
        !type(bc_t),             intent(in), optional    :: bc(:)
        !type(bcset_coupling_t), intent(in), optional    :: bcset_coupling(:)

        integer(ik) :: ierr, ndomains, idom


        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%dom)) deallocate(self%dom)


        !
        ! Allocate blockmatrix_t for each domain
        !
        ndomains = mesh%ndomains()
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization procedure for each blockmatrix_t
        !
        do idom = 1,ndomains

! WITH BC COUPLING
!            if ( present(bcset_coupling) ) then
!                call self%dom(idom)%init(mesh(idom),bcset_coupling(idom),mtype)
!            else
!                call self%dom(idom)%init(mesh(idom),mtype=mtype)
!            end if

! WITHOUT BC COUPLING
             call self%dom(idom)%init(mesh%domain(idom),mtype=mtype)

        end do


        !
        ! Set initialization to true
        !
        self%local_initialized = .true.

    end subroutine initialize
    !******************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
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
        ! Loop through LOCAL blocks and look for parallel multiply
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
                            ! Loop through chidgVector%recv to find match
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

                                            ! If they match the blockmatrix, set the recv indices so chidg_mv knows how to compute matrix-vector product
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
        ! Loop through CHIMERA blocks and look for parallel multiply
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
                                ! Loop through chidgVector%recv to find match
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

                                    

                                                ! If they match the blockmatrix, set the recv indices so chidg_mv knows how to compute matrix-vector product
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
        ! Set recv initialization to true
        !
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
    !!                          derivatives for the linearization matrix
    !!  @param[in]  idom        Domain index for storing the linearization
    !!  @param[in]  ielem       Element index for which the linearization was computed
    !!  @param[in]  iblk        Index of the block for the linearization of the given elemen
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
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

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idomain_l)%store(integral,face_info,seed,ivar,itime)

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

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idomain_l)%store_chimera(integral,face_info,seed,ivar,itime)

    end subroutine store_chimera
    !******************************************************************************************









    !>  Procedure for stiring linearization information for boundary condition faces
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

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idomain_l)%store_bc(integral,face_info,seed,ivar,itime)

    end subroutine store_bc
    !*******************************************************************************************








    !> Set all ChiDGMatrix matrix-values to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !! 
    !! 
    !-------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidg_matrix_t),   intent(inout)   :: self

        integer(ik) :: idom
    

        !
        ! Call blockmatrix_t%clear() on all matrices
        !
        do idom = 1,size(self%dom)
           call self%dom(idom)%clear() 
        end do
    
    
    end subroutine clear
    !*******************************************************************************************







    !>  Release allocated resources.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/3/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_matrix_t),  intent(inout)   :: self

        if (allocated(self%dom)) deallocate(self%dom)

    end subroutine release
    !*******************************************************************************************





    !>  chidg_matrix destructor.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidg_matrix_t),    intent(inout)   :: self

    end subroutine destructor
    !******************************************************************************************





end module type_chidg_matrix
