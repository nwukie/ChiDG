module type_blockmatrix
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: DIAG, ZERO, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, &
                                      NFACES, CHIMERA, NO_INTERIOR_NEIGHBOR
    use type_mesh,              only: mesh_t
    use type_densematrix,       only: densematrix_t
    use type_densematrix_vector,only: densematrix_vector_t
    use type_face_info,         only: face_info_t
    use type_seed,              only: seed_t
    use type_bcset_coupling,    only: bcset_coupling_t
    use type_ivector,           only: ivector_t
    use mod_chidg_mpi,          only: IRANK
    use DNAD_D
    implicit none


    !> Container for storing denseblock linearizations that make up a blockmatrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !! localblocks (nelem x 7)
    !!
    !!            xi_min   xi_max   eta_min   eta_max   zeta_min    zeta_max    diag
    !!
    !!  elem #1:
    !!  elem #2:
    !!  elem #3:
    !!    .
    !!    .
    !-------------------------------------------------------------------------------------------
    type, public :: blockmatrix_t

        !
        ! Primary storage, (nelem,ntime)
        !
        type(densematrix_vector_t), allocatable :: lblks(:,:)       
        type(densematrix_vector_t), allocatable :: chi_blks(:,:)       
        type(densematrix_vector_t), allocatable :: bc_blks(:,:)       

        !
        ! Supporting data
        !
        integer(ik),            allocatable :: ldata(:,:)               !< Block-local data, (nelem,3) nvars, nterms, ntime.
        !integer(ik),            allocatable :: local_transpose(:,:)     !< Block index of the transposed location (nelem,6)
        !type(ivector_t),        allocatable :: local_transpose(:)
        type(ivector_t),        allocatable :: local_lower_blocks(:)    !< For each element, which blocks (1-6) are lower blocks
        type(ivector_t),        allocatable :: local_upper_blocks(:)    !< For each element, which blocks (1-6) are upper blocks


    contains

        ! Initializers
        generic,   public  :: init => initialize_linearization          !< Initialize full linearization matrix
        procedure, private :: initialize_linearization


        ! Setters
        procedure :: store                                              !< Store linearization data for local blocks
        procedure :: store_chimera                                      !< Store linearization data for chimera blocks
        procedure :: store_bc                                           !< Store linearization data for boundary condition blocks
        procedure :: clear                                              !< Zero all data storage

        final :: destructor

    end type blockmatrix_t
    !*******************************************************************************************



    private
contains



    !> Subroutine for initializing local linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    mesh_t containing arrays of elements and faces
    !!  @param[in]  mtype   character string indicating the type of matrix to be initialized
    !!                      (ie. Full, Lower-Diagonal, Upper-Diagonal
    !!
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date   11/10/2016
    !!
    !-------------------------------------------------------------------------------------------
    subroutine initialize_linearization(self,mesh,bcset_coupling,mtype)
        class(blockmatrix_t),   intent(inout)             :: self
        class(mesh_t),          intent(in)                :: mesh
        type(bcset_coupling_t), intent(in), optional      :: bcset_coupling
        character(*),           intent(in)                :: mtype

        character(:),   allocatable :: user_msg
        integer(ik),    allocatable :: blocks(:)
        integer(ik)                 :: nelem, ierr, ielem, iblk, size1d, parent,  &
                                       block_index, neqns, nterms_s, ntime,             &
                                       nchimera_elements, maxdonors, idonor, iface,     &
                                       itime, dparent_g, dparent_l, eparent_g,          &
                                       eparent_l, parent_proc, eparent_l_trans,         &
                                       imat_trans, ChiID, ndonors, max_coupled_elems,   &
                                       ncoupled_elems, icoupled_elem, icoupled_elem_bc, &
                                       ielem_bc, ibc, imat
        logical                     :: new_elements, chimera_face, more_donors,         &
                                       donor_already_called, contains_chimera_face,     &
                                       block_initialized, lower_block, upper_block,     &
                                       transposed_block
        logical                     :: init_chimera = .false.
        logical                     :: init_bc      = .false.

        type(densematrix_t)         :: temp_blk, temp1
        !
        ! Select matrix blocks to initialize 
        !
        select case (trim(mtype))
            case ('full','Full','FULL')
                blocks       = [XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG]
                init_chimera = .true.
                init_bc      = .true.

            case ('L','l','Lower','lower')
                blocks       = [XI_MIN,ETA_MIN,ZETA_MIN]
                init_chimera = .false.
                init_bc      = .false.

            case ('U','u','Upper','upper')
                blocks       = [XI_MAX,ETA_MAX,ZETA_MAX]
                init_chimera = .false.
                init_bc      = .false.

            case ('LD','ld','LowerDiagonal','lowerdiagonal')
                blocks       = [XI_MIN,ETA_MIN,ZETA_MIN,DIAG]
                init_chimera = .false.
                init_bc      = .false.
                
            case ('UD','ud','UpperDiagonal','upperdiagonal')
                blocks       = [XI_MAX,ETA_MAX,ZETA_MAX,DIAG]
                init_chimera = .false.
                init_bc      = .false.

            case ('D', 'd', 'Diagonal', 'diagonal')
                blocks       = [DIAG]
                init_chimera = .false.
                init_bc      = .false.

            case default
                call chidg_signal(FATAL,'blockmatrix%init: unrecognized matrix type')

        end select



        ! Check to make sure the mesh numerics were initialized
        user_msg = "blockmatrix%initialize_linearization: Incoming mesh_t was not &
                    initialized. Make sure to call mesh%init_sol"
        if (.not. mesh%solInitialized) call chidg_signal(FATAL,user_msg)



        nelem = mesh%nelem      ! Number of elements in the local block
        ntime = mesh%ntime      ! Number of time levels
        !
        ! Allocation for 'local blocks'
        !
        if (allocated(self%lblks)) deallocate(self%lblks)
        allocate(self%lblks(nelem,ntime),           &
                 self%ldata(nelem,3),               &
                 self%local_lower_blocks(nelem),    &
                 self%local_upper_blocks(nelem), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Allocation for 'chimera blocks'
        !
        if (allocated(self%chi_blks)) deallocate(self%chi_blks)
        if (init_chimera) allocate(self%chi_blks(nelem,ntime), stat=ierr)
        if (ierr  /=   0) call AllocationError





!        !------------------------------------------------------------------------------
!        !
!        !                  Allocation for 'boundary condition blocks'
!        !
!        ! TODO: If the only 'coupled' element to the local face flux is the local element,
!        !       then a block is still allocated for it here. However, the linearization of the
!        !       local element is stored in self%lblks so the block here is not used. It
!        !       would be wasted compute time in the matrix-vector product because it is 
!        !       just zeros.
!        !
!        !------------------------------------------------------------------------------
!        if ( init_bc .and. present(bcset_coupling) ) then
!            if (allocated(self%bc_blks)) deallocate(self%bc_blks)
!
!            !
!            ! Get maximum number of coupled elements across all boundary conditions.
!            !
!            max_coupled_elems = 0
!            do ibc = 1,size(bcset_coupling%bc)
!
!                !
!                ! Loop through bc elems and test number of coupled elements against current maximum
!                !
!                do ielem = 1,size(bcset_coupling%bc(ibc)%elems)
!                    ncoupled_elems = bcset_coupling%bc(ibc)%coupled_elems(ielem)%size()
!                
!                    if ( ncoupled_elems > max_coupled_elems ) then
!                        max_coupled_elems = ncoupled_elems
!                    end if
!
!                end do ! ielem
!
!            end do ! ibc
!
!
!            !
!            ! Allocate boundary condition blocks
!            !
!            allocate(self%bc_blks(nelem,max_coupled_elems), stat=ierr)
!            if (ierr /= 0) call AllocationError
!
!        end if










        !
        ! Loop through elements and call initialization for 'local', 'chimera', and 
        ! 'boundary condition' denseblock matrices.
        !
        do ielem = 1,mesh%nelem
            do itime = 1,mesh%ntime 
                imat = 1
                
                ! Set the element indices that the densematrix_vector is associated with.
                self%lblks(ielem,itime)%idomain_g  = mesh%elems(ielem)%idomain_g
                self%lblks(ielem,itime)%idomain_l  = mesh%elems(ielem)%idomain_l
                self%lblks(ielem,itime)%ielement_g = mesh%elems(ielem)%ielement_g
                self%lblks(ielem,itime)%ielement_l = mesh%elems(ielem)%ielement_l
                self%lblks(ielem,itime)%mass       = mesh%elems(ielem)%mass


                !--------------------------------------------
                !
                ! Initialization  --  'local blocks'
                !
                !--------------------------------------------
                do block_index = 1,size(blocks)
                    iblk = blocks(block_index)
                    size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s

                    !
                    ! Parent is the element with respect to which the linearization is computed
                    !
                    dparent_l = mesh%idomain_l
                    if (iblk == DIAG) then
                        dparent_g   = mesh%elems(ielem)%idomain_g
                        dparent_l   = mesh%elems(ielem)%idomain_l
                        eparent_g   = mesh%elems(ielem)%ielement_g
                        eparent_l   = mesh%elems(ielem)%ielement_l
                        parent_proc = IRANK
                    else
                        dparent_g   = mesh%faces(ielem,iblk)%ineighbor_domain_g
                        dparent_l   = mesh%faces(ielem,iblk)%ineighbor_domain_l
                        eparent_g   = mesh%faces(ielem,iblk)%ineighbor_element_g
                        eparent_l   = mesh%faces(ielem,iblk)%ineighbor_element_l
                        parent_proc = mesh%faces(ielem,iblk)%ineighbor_proc
                    end if



                    !
                    ! Call initialization procedure if parent is not 0 (0 meaning there is no parent, probably a boundary)
                    !
                    if (eparent_l /= NO_INTERIOR_NEIGHBOR) then

                        ! Initialize dense block
                        call temp_blk%init(size1d,size1d,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc)
                        call self%lblks(ielem,itime)%push_back(temp_blk)

                        ! Store data about number of equations and number of terms in solution expansion
                        self%ldata(ielem,1) = mesh%elems(ielem)%neqns
                        self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
                        self%ldata(ielem,3) = mesh%elems(ielem)%ntime


                        !
                        ! If off-diagonal, store block index as 'upper' or 'lower'
                        !
                        ! TODO: Add consideration for ntime and how that affects upper/lower status.
                        !

                        lower_block = ( (parent_proc == IRANK .and. eparent_l < ielem) .or. (parent_proc < IRANK) )
                        upper_block = ( (parent_proc == IRANK .and. eparent_l > ielem) .or. (parent_proc > IRANK) )

                        if ( lower_block ) then
                            call self%local_lower_blocks(ielem)%push_back(imat)
                        else if ( upper_block ) then
                            call self%local_upper_blocks(ielem)%push_back(imat)
                        end if

                        imat = imat + 1
                    end if

                end do ! init local
                !********************************************
                



                !--------------------------------------------
                !
                ! Initialization  --  'chimera blocks'
                !
                !--------------------------------------------
                if (init_chimera) then
                    ! Set the element indices that the densematrix_vector is associated with.
                    self%chi_blks(ielem,itime)%idomain_g  = mesh%elems(ielem)%idomain_g
                    self%chi_blks(ielem,itime)%idomain_l  = mesh%elems(ielem)%idomain_l
                    self%chi_blks(ielem,itime)%ielement_g = mesh%elems(ielem)%ielement_g
                    self%chi_blks(ielem,itime)%ielement_l = mesh%elems(ielem)%ielement_l

                    do iface = 1,NFACES

                        !
                        ! If facetype is CHIMERA
                        !
                        chimera_face = ( mesh%faces(ielem,iface)%ftype == CHIMERA )
                        if (chimera_face) then

                            !
                            ! Get ChiID and number of donor elements
                            !
                            ChiID   = mesh%faces(ielem,iface)%ChiID
                            ndonors = mesh%chimera%recv%data(ChiID)%ndonors()

                            !
                            ! Call block initialization for each Chimera donor
                            !
                            do idonor = 1,ndonors
                                neqns       = mesh%chimera%recv%data(ChiID)%donor_neqns%at(idonor)
                                nterms_s    = mesh%chimera%recv%data(ChiID)%donor_nterms_s%at(idonor)
                                dparent_g   = mesh%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)
                                dparent_l   = mesh%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)
                                eparent_g   = mesh%chimera%recv%data(ChiID)%donor_element_g%at(idonor)
                                eparent_l   = mesh%chimera%recv%data(ChiID)%donor_element_l%at(idonor)
                                parent_proc = mesh%chimera%recv%data(ChiID)%donor_proc%at(idonor)

                                size1d = neqns * nterms_s

                                !
                                ! Check if block initialization was already called for current donor
                                !
                                donor_already_called = .false.
                                do imat = 1,self%chi_blks(ielem,itime)%size()
                                    
                                    ! dummy densematrix to get a specific densematrix inside the chi_blks
                                    ! densematrix_vector temporary variable to access densematrix routine
                                    temp1 = self%chi_blks(ielem,itime)%at(imat)
                                    
                                    donor_already_called = ( dparent_g == temp1%dparent_g() .and. &
                                                             dparent_l == temp1%dparent_l() .and. &
                                                             eparent_g == temp1%eparent_g() .and. &
                                                             eparent_l == temp1%eparent_l() )
                                    if (donor_already_called) exit
                                end do

                            
                                !
                                ! If a block for the donor element hasn't yet been initialized, call initialization procedure
                                !
                                if (.not. donor_already_called) then

                                    !
                                    ! Call block initialization, store
                                    !
                                    call temp_blk%init(size1d,size1d,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc)
                                    call self%chi_blks(ielem,itime)%push_back(temp_blk)

                                end if

                            end do ! idonor

                        end if

                    end do ! iface

                end if  ! init_chimera
                !********************************************
            
            end do  ! itime

        end do ! ielem



!        !--------------------------------------------
!        !
!        ! Initialization  --  'boundary condition blocks'
!        !
!        !--------------------------------------------
!        !
!        ! Loop through boundary conditions and initialize blocks for coupling
!        !
!        if ( init_bc .and. present(bcset_coupling) ) then
!            do ibc = 1,size(bcset_coupling%bc)
!
!                !
!                ! For the current boundary condition, loop through bc elements.
!                !
!                do ielem_bc = 1,size(bcset_coupling%bc(ibc)%elems)
!
!                    ncoupled_elems = bcset_coupling%bc(ibc)%coupled_elems(ielem_bc)%size()
!                    !
!                    ! Initialize block storage for each coupled element
!                    !
!                    do icoupled_elem_bc = 1,ncoupled_elems
!
!                        !
!                        ! Get block indices
!                        !
!                        ielem         = bcset_coupling%bc(ibc)%elems(ielem_bc)
!                        icoupled_elem = bcset_coupling%bc(ibc)%coupled_elems(ielem_bc)%at(icoupled_elem_bc)
!
!
!                        !
!                        ! Check if block has already been initialized for the coupled element
!                        !
!                        block_initialized = .false.
!                        do iblk = 1,size(self%bc_blks,2)
!                            if ( self%bc_blks(ielem,iblk)%eparent() == icoupled_elem ) then
!                                block_initialized = .true.
!                                exit
!                            end if
!                        end do
!
!
!                        if ( .not. block_initialized ) then
!                            !
!                            ! Compute block size
!                            !
!                            size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s
!
!                            !
!                            ! Call boundary condition block initialization
!                            !
!                            dparent_l = mesh%idomain_l
!                            eparent_l = icoupled_elem
!                            call self%bc_blks(ielem,icoupled_elem_bc)%init(size1d,dparent_l,eparent_l)
!                        end if
!
!
!                    end do ! icoupled_elem
!
!                end do ! ielem
!
!            end do ! ibc
!
!        end if ! init_bc







        ! 
        ! Initialize transpose data
        !
        select case (trim(mtype))
            case ('full','Full','FULL')

                do ielem = 1,mesh%nelem
                    do itime = 1,mesh%ntime
                        do imat = 1,self%lblks(ielem,itime)%size()

                            if ( self%lblks(ielem,itime)%parent_proc(imat) == IRANK ) then

                                    ! Get parent element of off-diagonal block
                                    eparent_l = self%lblks(ielem,itime)%eparent_l(imat)

                                    !
                                    ! Find block index of transposed location in parent lblks
                                    !
                                    do imat_trans = 1,self%lblks(eparent_l,itime)%size()
                                        ! Make sure the block we are seaching is on-proc
                                        if ( self%lblks(eparent_l,itime)%parent_proc(imat_trans) == IRANK ) then
                                            eparent_l_trans = self%lblks(eparent_l,itime)%eparent_l(imat_trans)

                                            transposed_block = ( eparent_l_trans == ielem )
                                            if ( transposed_block ) then
                                                call self%lblks(ielem,itime)%set_itranspose(imat,imat_trans)
                                                exit
                                            end if
                                        end if

                                        if ( (imat_trans == self%lblks(eparent_l,itime)%size() ) .and. &
                                             (transposed_block .eqv. .false.) ) call chidg_signal(FATAL,"blockmatrix%init: no transposed element found")
                                    end do !imat_trans

                            end if

                        end do !imat
                    end do !itime
                end do !ielem

        end select




    end subroutine initialize_linearization
    !******************************************************************************************















    !>  Stores derivative data to the linearization matrix
    !!
    !!      -- Given the integral data computed from the spatial discretization,
    !!         store the derivative values from the AD data types
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/10/2016
    !!
    !!  @param[in]  integral    Array of modes with embedded partial derivatives for the linearization
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!  @param[in]  itime       Index of a time level for the linearization of the given element
    !!
    !-----------------------------------------------------------------------------------------
    subroutine store(self,integral,face_info,seed,ivar,itime)
        class(blockmatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        type(face_info_t),      intent(in)      :: face_info
        type(seed_t),           intent(in)      :: seed
        integer(ik),            intent(in)      :: ivar
        integer(ik),            intent(in)      :: itime

        integer(ik) :: nterms, imat
        integer(ik) :: ielement_l, idonor_domain_g, idonor_element_g
        

        ielement_l = face_info%ielement_l

        idonor_domain_g  = seed%idomain_g
        idonor_element_g = seed%ielement_g


        ! Get stored information for the block
        nterms = self%ldata(ielement_l,2)  

        ! Find donor densematrix location 
        imat = self%lblks(ielement_l,itime)%loc(idonor_domain_g,idonor_element_g)

        ! Call subroutine on densematrix 
        call self%lblks(ielement_l,itime)%store(imat,ivar,nterms,integral)

    end subroutine store
    !******************************************************************************************










    !>  Stores derivative data from Chimera faces to the linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/10/2015
    !!
    !!  @param[in]  integral    Array of modes with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!  @param[in]  itime       Index of a time level for the linearization of the given element
    !!
    !------------------------------------------------------------------------------------------
    subroutine store_chimera(self,integral,face_info,seed,ivar,itime)
        class(blockmatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar
        integer(ik),                intent(in)      :: itime

        integer(ik) :: nterms, imat
        integer(ik) :: ielement_l, idonor_domain_g, idonor_element_g
        

        ielement_l = face_info%ielement_l

        idonor_domain_g  = seed%idomain_g
        idonor_element_g = seed%ielement_g


        ! Get stored information for the block
        nterms = self%ldata(ielement_l,2)


        ! Find donor densematrix location 
        imat = self%chi_blks(ielement_l,itime)%loc(idonor_domain_g,idonor_element_g)


        ! Store derivatives
        call self%chi_blks(ielement_l,itime)%store(imat,ivar,nterms,integral)


    end subroutine store_chimera
    !******************************************************************************************










    !>  Stores derivative data from boundary condition coupling to the linearization matrix.
    !!
    !!  The linearization of a flux wrt the local element is stored in self%lblks in the DIAG location.
    !!  The linearization of a flux wrt other elements on the boundary is stored in self%bc_blks.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/14/2016
    !!
    !!  @param[in]  integral    Array of modes with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!  @param[in]  itime       Index of a time level for the linearization of the given element
    !!
    !------------------------------------------------------------------------------------------
    subroutine store_bc(self,integral,face,seed,ivar,itime)
        class(blockmatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar
        integer(ik),                intent(in)      :: itime

        integer(ik) :: idomain_l, ielement_l
        integer(ik) :: idonor_domain_l, idonor_element_l, imat
        integer(ik) :: nterms, size_integral
        logical     :: local_element_linearization = .false.


        idomain_l  = face%idomain_l
        ielement_l = face%ielement_l

        idonor_domain_l  = seed%idomain_l
        idonor_element_l = seed%ielement_l


        !
        ! If ielem = ielem_d then the linearization is with respect to the local element. 
        ! So, this is stored in the self%lblks array in the DIAG location, instead of
        ! the self%bc_blks array. In general, the storage location is not important,
        ! but the ILU preconditioner expects the full diagonal contribution to be in 
        ! lblks.
        !
        local_element_linearization = (ielement_l == idonor_element_l)

        if ( local_element_linearization ) then

            call self%store(integral,face,seed,ivar,itime)

        else


            ! Get stored information for the block
            nterms = self%ldata(ielement_l,2)

            ! Find coupled bc densematrix location 
            imat = self%bc_blks(ielement_l,itime)%loc(idonor_domain_l,idonor_element_l)

            ! Store derivatives
            call self%bc_blks(ielement_l,itime)%store(imat,ivar,nterms,integral)


        end if ! check local block.

    end subroutine store_bc
    !******************************************************************************************











    !>  Set all denseblock_t storage to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date   11/14/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(blockmatrix_t),   intent(inout)   :: self

        integer(ik)             :: ielem, itime

        !
        ! For each element
        !
        do ielem = 1,size(self%lblks,1)


            ! Clear local matrices
            do itime = 1,size(self%lblks,2)
                call self%lblks(ielem,itime)%setzero()
            end do  ! itime


            ! Clear chimera matrices
            if (allocated(self%chi_blks)) then
                do itime = 1,size(self%chi_blks,2)
                    call self%chi_blks(ielem,itime)%setzero()
                end do ! itime
            end if


            ! Clear bc matrices
            if (allocated(self%bc_blks)) then
                do itime = 1,size(self%bc_blks,2)
                    call self%bc_blks(ielem,itime)%setzero()
                end do ! itime
            end if


        end do ! ielem


    end subroutine clear
    !******************************************************************************************




    !>
    !!
    !!---------------------------------------------
    subroutine destructor(self)
        type(blockmatrix_t), intent(inout) :: self

    end subroutine
    !**********************************************

end module type_blockmatrix
