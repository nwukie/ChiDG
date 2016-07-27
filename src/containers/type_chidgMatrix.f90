module type_chidgMatrix
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_chidg_mpi,          only: IRANK
    use type_blockmatrix,       only: blockmatrix_t
    use type_mesh,              only: mesh_t
    use type_face_info,         only: face_info_t
    use type_seed,              only: seed_t
    use type_bcset_coupling,    only: bcset_coupling_t
    use type_chidgVector,       only: chidgVector_t
    use DNAD_D
    implicit none




    !> ChiDG matrix type. Contains an array of blockmatrix_t types, each corresponding to a domain.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    type, public :: chidgMatrix_t

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


        final       :: destructor

    end type chidgMatrix_t
    !***********************************************************************************************************



    private
contains




    !>  Subroutine for initializing chidgMatrix_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  domains     Array of domain_t instances
    !!  
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine initialize(self,mesh,bcset_coupling,mtype)
        class(chidgMatrix_t),   intent(inout)           :: self
        type(mesh_t),           intent(in)              :: mesh(:)
        type(bcset_coupling_t), intent(in), optional    :: bcset_coupling(:)
        character(*),           intent(in)              :: mtype

        integer(ik) :: ierr, ndomains, idom


        !
        ! Allocate blockmatrix_t for each domain
        !
        ndomains = size(mesh)
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
             call self%dom(idom)%init(mesh(idom),mtype=mtype)

        end do


        !
        ! Set initialization to true
        !
        self%local_initialized = .true.

    end subroutine initialize
    !***********************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine init_recv(self,x)
        class(chidgMatrix_t),   intent(inout)   :: self
        type(chidgVector_t),    intent(in)      :: x

        integer(ik) :: idom, ielem, iblk, matrix_proc, vector_proc, comm_proc
        integer(ik) :: dparent_g, eparent_g, parent_proc, icomm, idom_recv, ielem_recv, drecv_g, erecv_g
        logical     :: local_multiply, parallel_multiply, match_found

        
        !
        ! Loop through LOCAL blocks and look for parallel multiply
        !
        do idom = 1,size(self%dom)
            do ielem = 1,size(self%dom(idom)%lblks,1)
                do iblk = 1,size(self%dom(idom)%lblks,2)
                    
                    if (allocated(self%dom(idom)%lblks(ielem,iblk)%mat)) then
                        matrix_proc = IRANK
                        vector_proc = self%dom(idom)%lblks(ielem,iblk)%parent_proc()

                        local_multiply    = ( matrix_proc == vector_proc )
                        parallel_multiply = ( matrix_proc /= vector_proc )


                        if ( parallel_multiply ) then
                            !
                            ! Get information about element we need to multiply with
                            !
                            dparent_g   = self%dom(idom)%lblks(ielem,iblk)%dparent_g()
                            eparent_g   = self%dom(idom)%lblks(ielem,iblk)%eparent_g()
                            parent_proc = self%dom(idom)%lblks(ielem,iblk)%parent_proc()



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
                                                self%dom(idom)%lblks(ielem,iblk)%recv_comm    = icomm
                                                self%dom(idom)%lblks(ielem,iblk)%recv_domain  = idom_recv
                                                self%dom(idom)%lblks(ielem,iblk)%recv_element = ielem_recv
                                                match_found = .true.
                                            end if

                                        end do !ielem_recv
                                    end do !idom_recv
                                end if

                            end do ! icomm

                            if (.not. match_found) call chidg_signal(FATAL,"chidgMatrix%init_recv: no matching recv element found in vector")



                        end if

                    end if

                end do !iblk
            end do !ielem

        end do ! idom









        !
        ! Loop through CHIMERA blocks and look for parallel multiply
        !
        do idom = 1,size(self%dom)

            if (allocated(self%dom(idom)%chi_blks)) then
                do ielem = 1,size(self%dom(idom)%chi_blks,1)
                    do iblk = 1,size(self%dom(idom)%chi_blks,2)
                        
                        if (allocated(self%dom(idom)%chi_blks(ielem,iblk)%mat)) then
                            matrix_proc = IRANK
                            vector_proc = self%dom(idom)%chi_blks(ielem,iblk)%parent_proc()

                            local_multiply    = ( matrix_proc == vector_proc )
                            parallel_multiply = ( matrix_proc /= vector_proc )


                            if ( parallel_multiply ) then
                                !
                                ! Get information about element we need to multiply with
                                !
                                dparent_g   = self%dom(idom)%chi_blks(ielem,iblk)%dparent_g()
                                eparent_g   = self%dom(idom)%chi_blks(ielem,iblk)%eparent_g()
                                parent_proc = self%dom(idom)%chi_blks(ielem,iblk)%parent_proc()



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
                                                    self%dom(idom)%chi_blks(ielem,iblk)%recv_comm    = icomm
                                                    self%dom(idom)%chi_blks(ielem,iblk)%recv_domain  = idom_recv
                                                    self%dom(idom)%chi_blks(ielem,iblk)%recv_element = ielem_recv
                                                    match_found = .true.
                                                end if

                                            end do !ielem_recv
                                        end do !idom_recv
                                    end if

                                end do ! icomm

                                if (.not. match_found) call chidg_signal(FATAL,"chidgMatrix%init_recv: no matching recv element found in vector")



                            end if

                        end if

                    end do !iblk
                end do !ielem
            end if

        end do ! idom











        !
        ! Set recv initialization to true
        !
        self%recv_initialized = .true.


    end subroutine init_recv
    !***********************************************************************************************************














    !> Procedure for storing linearization information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  idom        Domain index for storing the linearization
    !!  @param[in]  ielem       Element index for which the linearization was computed
    !!  @param[in]  iblk        Index of the block for the linearization of the given elemen
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine store(self, integral, idom, ielem, iblk, ivar)
        class(chidgMatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        integer(ik),            intent(in)      :: idom, ielem, iblk, ivar

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idom)%store(integral,ielem,iblk,ivar)

    end subroutine store
    !***********************************************************************************************************









    !> Procedure for stiring linearization information for Chimera faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing the indices defining the Chimera face
    !!  @param[in]  seed        seed_t containing the indices defining the element against which the Chimera face was linearized
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine store_chimera(self,integral,face,seed,ivar)
        class(chidgMatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar 

        integer(ik) :: idomain_l

        idomain_l = face%idomain_l

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idomain_l)%store_chimera(integral,face,seed,ivar)

    end subroutine store_chimera
    !***********************************************************************************************************









    !> Procedure for stiring linearization information for boundary condition faces
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing the indices defining the Chimera face
    !!  @param[in]  seed        seed_t containing the indices defining the element against which the Chimera face was linearized
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine store_bc(self,integral,face,seed,ivar)
        class(chidgMatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar 

        integer(ik) :: idomain_l

        idomain_l = face%idomain_l

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idomain_l)%store_bc(integral,face,seed,ivar)

    end subroutine store_bc
    !***********************************************************************************************************


















    !> Set all ChiDGMatrix matrix-values to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !! 
    !! 
    !----------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidgMatrix_t),   intent(inout)   :: self

        integer(ik) :: idom
    

        !
        ! Call blockmatrix_t%clear() on all matrices
        !
        do idom = 1,size(self%dom)
           call self%dom(idom)%clear() 
        end do
    
    
    end subroutine clear
    !***********************************************************************************************************











    !> ChiDGMatrix destructor.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidgMatrix_t),    intent(inout)   :: self

    end subroutine
    !***********************************************************************************************************



end module type_chidgMatrix
