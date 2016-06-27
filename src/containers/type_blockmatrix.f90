!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockmatrix
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: DIAG, ZERO, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, &
                                      NFACES, CHIMERA
    use type_mesh,              only: mesh_t
    use type_densematrix,       only: densematrix_t
    use type_face_info,         only: face_info_t
    use type_seed,              only: seed_t
    use type_bcset_coupling,    only: bcset_coupling_t
    use mod_chidg_mpi,          only: IRANK
    use DNAD_D
    implicit none


    !> Container for storing denseblock linearizations that make up a blockmatrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
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
    !-------------------------------------------------------------------------------------------------------------------------------
    type, public :: blockmatrix_t

        type(densematrix_t), allocatable :: lblks(:,:)                      !< Local domain blocks  (nelem, NBLK)
        integer(ik),         allocatable :: ldata(:,:)                      !< Block-local  data    (ielem, 1) -> nvars, (ielem, 2) -> nterms (nvars, nterms)

        !
        ! These may be better located in an array of densematrix_vector containers instead of just an allocatable array.
        !
        type(densematrix_t), allocatable :: chi_blks(:,:)                    !< Chimera inter-domain blocks         (nelem, MaxDonors)
        type(densematrix_t), allocatable :: bc_blks(:,:)                     !< Boundary condition coupling blocks  (nelem, Max coupled elems)

    contains
        ! Initializers
        generic,   public  :: init => initialize_linearization              !< Initialize full linearization matrix
        procedure, private :: initialize_linearization


        ! Setters
        procedure :: store                                                  !< Store linearization data for local blocks
        procedure :: store_chimera                                          !< Store linearization data for chimera blocks
        procedure :: store_bc                                               !< Store linearization data for boundary condition blocks
        procedure :: clear                                                  !< Zero all data storage



        final :: destructor
    end type blockmatrix_t
    !*******************************************************************************************************************************



    private
contains



    !> Subroutine for initializing local linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    mesh_t containing arrays of elements and faces
    !!  @param[in]  mtype   character string indicating the type of matrix to be initialized (ie. Full, Lower-Diagonal, Upper-Diagonal
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine initialize_linearization(self,mesh,bcset_coupling,mtype)
        class(blockmatrix_t),   intent(inout)             :: self
        class(mesh_t),          intent(in)                :: mesh
        type(bcset_coupling_t), intent(in), optional      :: bcset_coupling
        character(*),           intent(in)                :: mtype

        integer(ik), allocatable    :: blocks(:)
        integer(ik)                 :: nelem, nblk, ierr, ielem, iblk, size1d, parent, block_index, neqns, nterms_s
        integer(ik)                 :: nchimera_elements, maxdonors, idonor, iface
        integer(ik)                 :: eparent_l, dparent_l
        integer(ik)                 :: iopen, ChiID, ndonors, max_coupled_elems, ncoupled_elems, icoupled_elem, icoupled_elem_bc, ielem_bc, ibc
        logical                     :: new_elements
        logical                     :: chimera_face
        logical                     :: more_donors
        logical                     :: donor_already_called
        logical                     :: contains_chimera_face
        logical                     :: block_initialized
        logical                     :: init_chimera = .false.
        logical                     :: init_bc      = .false.


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

            case default
                call chidg_signal(FATAL,'blockmatrix%init: unrecognized matrix type')

        end select



        nelem = mesh%nelem      ! Number of elements in the local block
        nblk  = 7               ! Number of potential blocks in the linearization for a given element (1D => 3, 2D => 5, 3D => 7)


        !
        ! Check to make sure the mesh numerics were initialized
        !
        if (.not. mesh%solInitialized) call chidg_signal(FATAL,'blockmatrix_t%initialize_linearization: Incoming mesh_t was not initialized. Make sure to call mesh%init_sol')


        !------------------------------------------------------------------------------
        !
        ! Allocation for 'local blocks'
        !
        !------------------------------------------------------------------------------
        ! If matrix was already allocated, deallocate and then reallocate matrix size
        ! Reallocation would take place if the number of elements were changed
        !
        if (allocated(self%lblks)) then

            !
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            !
            new_elements = (mesh%nelem /= size(self%lblks,1))
            if (new_elements) then
                deallocate(self%lblks, self%ldata)
                allocate(self%lblks(nelem,nblk), self%ldata(nelem,2), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        else

            allocate(self%lblks(nelem,nblk), self%ldata(nelem,2), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if




        !------------------------------------------------------------------------------
        !
        ! Allocation for 'chimera blocks'
        !
        !------------------------------------------------------------------------------
        if (allocated(self%chi_blks)) deallocate(self%chi_blks)
        !
        ! Assemble some Chimera data
        !
        
        !
        ! Get maximum number of donor elements to a given element. May include donors from multiple faces
        !
        maxdonors = 0
        nchimera_elements = 0
        do ielem = 1,mesh%nelem

            ndonors = 0
            contains_chimera_face = .false.
            do iface = 1,NFACES

                !
                ! Check for Chimera face type and add to element donor count
                !
                if (mesh%faces(ielem,iface)%ftype == CHIMERA) then
                    ChiID = mesh%faces(ielem,iface)%ChiID
                    ndonors = ndonors + mesh%chimera%recv%data(ChiID)%ndonors
                    contains_chimera_face = .true.
                end if

            end do


            !
            ! If current face has more donors than the current max, update the max number of donors
            !
            more_donors = ( ndonors > maxdonors )
            if (more_donors) then
                maxdonors = ndonors
            end if

            !
            ! Increment number of chimera elements
            !
            if (contains_chimera_face) then
                nchimera_elements = nchimera_elements + 1
            end if

        end do



        !
        ! Allocate Chimera blocks
        !
        if (init_chimera) then
            if (maxdonors > 0) then

                allocate(self%chi_blks(nelem, maxdonors), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if
        end if





        !------------------------------------------------------------------------------
        !
        ! Allocation for 'boundary condition blocks'
        !
        ! TODO: If the only 'coupled' element to the local face flux is the local element,
        !       then a block is still allocated for it here. However, the linearization of the
        !       local element is stored in self%lblks so the block here is not used. It
        !       would be wasted compute time in the matrix-vector product because it is 
        !       just zeros.
        !
        !------------------------------------------------------------------------------
        if ( init_bc .and. present(bcset_coupling) ) then
            if (allocated(self%bc_blks)) deallocate(self%bc_blks)

            !
            ! Get maximum number of coupled elements across all boundary conditions.
            !
            max_coupled_elems = 0
            do ibc = 1,size(bcset_coupling%bc)

                !
                ! Loop through bc elems and test number of coupled elements against current maximum
                !
                do ielem = 1,size(bcset_coupling%bc(ibc)%elems)
                    ncoupled_elems = bcset_coupling%bc(ibc)%coupled_elems(ielem)%size()
                
                    if ( ncoupled_elems > max_coupled_elems ) then
                        max_coupled_elems = ncoupled_elems
                    end if

                end do ! ielem

            end do ! ibc


            !
            ! Allocate boundary condition blocks
            !
            allocate(self%bc_blks(nelem,max_coupled_elems), stat=ierr)
            if (ierr /= 0) call AllocationError

        end if










        !
        ! Loop through elements and call initialization for 'local', 'chimera', and 'boundary condition' denseblock matrices
        !
        do ielem = 1,mesh%nelem

            
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
!                dparent_l = mesh%idomain_l
!                if (iblk == DIAG) then
!                    eparent_l = mesh%elems(ielem)%ielement_l
!                else
!                    eparent_l = mesh%faces(ielem,iblk)%get_neighbor_element_l()
!                end if

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
                    dparent_g   = mesh%faces(ielem,iblk)%ineighbor_element_g
                    dparent_l   = mesh%faces(ielem,iblk)%ineighbor_element_l
                    parent_proc = mesh%faces(ielem,iblk)%ineighbor_proc
                end if





                !
                ! Call initialization procedure if parent is not 0 (0 meaning there is no parent for that block, probably a boundary)
                !
                if (eparent_l /= 0) then
                    !call self%lblks(ielem,iblk)%init(size1d,dparent_l,eparent_l)
                    call self%lblks(ielem,iblk)%init(size1d,size1d,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc)

                    ! Store data about number of equations and number of terms in solution expansion
                    self%ldata(ielem,1) = mesh%elems(ielem)%neqns
                    self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
                end if

            end do ! init local
            !********************************************


            !--------------------------------------------
            !
            ! Initialization  --  'chimera blocks'
            !
            !--------------------------------------------
            if (init_chimera) then
                do iface = 1,NFACES

                    !
                    ! If facetype is CHIMERA
                    !
                    chimera_face = ( mesh%faces(ielem,iface)%ftype == CHIMERA )
                    if (chimera_face) then

                        !
                        ! Get ChiID and number of donor elements
                        !
                        ChiID = mesh%faces(ielem,iface)%ChiID
                        ndonors = mesh%chimera%recv%data(ChiID)%ndonors

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
                            do iblk = 1,maxdonors
                                donor_already_called = ( dparent_g == self%chi_blks(ielem,iblk)%dparent_g() .and. &
                                                         dparent_l == self%chi_blks(ielem,iblk)%dparent_l() .and. &
                                                         eparent_g == self%chi_blks(ielem,iblk)%eparent_g() .and. &
                                                         eparent_l == self%chi_blks(ielem,iblk)%eparent_l() )
                                if (donor_already_called) exit
                            end do

                        
                            !
                            ! If a block for the donor element hasn't yet been initialized, call initialization procedure
                            !
                            if (.not. donor_already_called) then

                                !
                                ! Find next open block to initialize for the current element
                                !
                                do iblk = 1,maxdonors
                                    if (.not. allocated(self%chi_blks(ielem,iblk)%mat) ) then
                                        iopen = iblk
                                        exit
                                    end if
                                end do

                                !
                                ! Call block initialization
                                !
                                !call self%chi_blks(ielem,iopen)%init(size1d,dparent_l,eparent_l)
                                call self%chi_blks(ielem,iopen)%init(size1d,size1d,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc)

                            end if

                        end do ! idonor

                    end if

                end do ! iface

            end if  ! init_chimera
            !********************************************

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



    end subroutine initialize_linearization
    !*********************************************************************************************************************















    !>  Stores derivative data to the linearization matrix
    !!
    !!      -- Given the integral data computed from the spatial discretization,
    !!         store the derivative values from the AD data types
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  ielem       Element for which the linearization was computed
    !!  @param[in]  iblk        Index of a block for the linearization of the given element
    !!  @param[in]  ivar        Index of the variable
    !!
    !----------------------------------------------------------------------------------------------------------------------------------
    subroutine store(self,integral,ielem,iblk,ivar)
        class(blockmatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        integer(ik),            intent(in)      :: ielem, iblk, ivar

        integer(ik) :: iarray, neqns, nterms, irow_start, irow


        !
        ! Get stored information for the block
        !
        neqns  = self%ldata(ielem,1)
        nterms = self%ldata(ielem,2)


        !
        ! Compute correct row offset for ivar
        !
        irow_start = ( (ivar - 1)  *  nterms)


        !
        ! Loop through integral values, for each value store its derivatives.
        ! The integral values here should be components of the RHS vector. An array of partial derivatives from an AD_D variable
        ! should be stored as a row in the block matrix.
        !
        do iarray = 1,size(integral)

            !
            ! Do a += operation to add derivatives to any that are currently stored
            !
            irow = irow_start + iarray
            self%lblks(ielem,iblk)%mat(irow,:) = self%lblks(ielem,iblk)%mat(irow,:) + integral(iarray)%xp_ad_

        end do


    end subroutine store
    !*******************************************************************************************************************************













    !>  Stores derivative data from Chimera faces to the linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine store_chimera(self,integral,face,seed,ivar)
        class(blockmatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar

        integer(ik) :: idomain_l, ielement_l
        integer(ik) :: idonor_domain_l, idonor_element_l
        integer(ik) :: irow, irow_start, donorblk, i, iblk, iarray
        integer(ik) :: neqns, nterms
        logical     :: block_match    = .false.
        logical     :: no_donor_block = .false.

        idomain_l  = face%idomain_l
        ielement_l = face%ielement_l

        idonor_domain_l  = seed%idomain_l
        idonor_element_l = seed%ielement_l


        !
        ! Get stored information for the block
        !
        neqns  = self%ldata(ielement_l,1)
        nterms = self%ldata(ielement_l,2)

        !
        ! Compute correct row offset for ivar
        !
        irow_start = ( (ivar - 1) * nterms )


        !
        ! Find donor block location 
        !
        donorblk = 0
        do iblk = 1,size(self%chi_blks,2)
            block_match = ( (idonor_domain_l  == self%chi_blks(ielement_l,iblk)%dparent_l()) .and. &
                            (idonor_element_l == self%chi_blks(ielement_l,iblk)%eparent_l()) )

            if ( block_match ) then
                donorblk = iblk
                exit
            end if
        end do

        no_donor_block = (donorblk == 0)
        if (no_donor_block) call chidg_signal(MSG,'blockmatrix%store_chimera: no donor block found to store derivative')



        !
        ! Store derivatives
        !
        do iarray = 1,size(integral)
            ! Do a += operation to add derivatives to any that are currently stored
            irow = irow_start + iarray
            self%chi_blks(ielement_l,donorblk)%mat(irow,:) = self%chi_blks(ielement_l,donorblk)%mat(irow,:) + integral(iarray)%xp_ad_
        end do



    end subroutine store_chimera
    !*******************************************************************************************************************************










    !>  Stores derivative data from boundary condition coupling to the linearization matrix.
    !!
    !!  The linearization of a flux wrt the local element is stored in self%lblks in the DIAG location.
    !!  The linearization of a flux wrt other elements on the boundary is stored in self%bc_blks.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine store_bc(self,integral,face,seed,ivar)
        class(blockmatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_info_t),          intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar

        integer(ik) :: idomain_l, ielement_l
        integer(ik) :: idonor_domain_l, idonor_element_l
        integer(ik) :: irow, irow_start, i, iblk, iarray
        integer(ik) :: neqns, nterms, bcblk
        logical     :: block_match = .false.
        logical     :: no_bc_block = .false.
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

            call self%store(integral,ielement_l,DIAG,ivar)

        else

            !
            ! Get stored information for the block
            !
            neqns  = self%ldata(ielement_l,1)
            nterms = self%ldata(ielement_l,2)

            !
            ! Compute correct row offset for ivar
            !
            irow_start = ( (ivar - 1) * nterms )


            !
            ! Find coupled bc block location 
            !
            bcblk = 0
            do iblk = 1,size(self%bc_blks,2)
                block_match = ( (idonor_domain_l  == self%bc_blks(ielement_l,iblk)%dparent_l()) .and. &
                                (idonor_element_l == self%bc_blks(ielement_l,iblk)%eparent_l()) )

                if ( block_match ) then
                    bcblk = iblk
                    exit
                end if
            end do

            no_bc_block = (bcblk == 0)
            if (no_bc_block) call chidg_signal(FATAL,"blockmatrix%store_bc: no bc block found to store derivatives.")



            !
            ! Store derivatives
            !
            do iarray = 1,size(integral)
                ! Do a += operation to add derivatives to any that are currently stored
                irow = irow_start + iarray
                self%bc_blks(ielement_l,bcblk)%mat(irow,:) = self%bc_blks(ielement_l,bcblk)%mat(irow,:) + integral(iarray)%xp_ad_
            end do


        end if ! check local block.

    end subroutine store_bc
    !*******************************************************************************************************************************

























    !>  Set all denseblock_t storage to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(blockmatrix_t),   intent(inout)   :: self

        integer(ik) :: ielem, iblk  ! do-loop counters

        !
        ! For each element
        !
        do ielem = 1,size(self%lblks,1)



            !
            ! For each local block linearization for the current element
            !
            do iblk = 1,size(self%lblks,2)

                !
                ! Check if the block storage is actually allocated
                !
                if (allocated(self%lblks(ielem,iblk)%mat)) then
                    !
                    ! If so, set to ZERO
                    !
                    self%lblks(ielem,iblk)%mat = ZERO
                end if

            end do ! iblk



            !
            ! For each Chimera block linearization for the current element
            !
            if (allocated(self%chi_blks)) then
                do iblk = 1,size(self%chi_blks,2)

                    !
                    ! Check if the block storage is actually allocated
                    !
                    if (allocated(self%chi_blks(ielem,iblk)%mat)) then

                        !
                        ! If so, set to ZERO
                        !
                        self%chi_blks(ielem,iblk)%mat = ZERO

                    end if

                end do ! iblk
            end if




            !
            ! For each boundary condition block linearization for the current element
            !
            if (allocated(self%bc_blks)) then
                do iblk = 1,size(self%bc_blks,2)

                    !
                    ! Check if the block storage is actually allocated
                    !
                    if (allocated(self%bc_blks(ielem,iblk)%mat)) then

                        !
                        ! If so, set to ZERO
                        !
                        self%bc_blks(ielem,iblk)%mat = ZERO

                    end if

                end do ! iblk
            end if







        end do ! ielem


    end subroutine clear
    !*******************************************************************************************************************************












    !>
    !!
    !!---------------------------------------------
    subroutine destructor(self)
        type(blockmatrix_t), intent(inout) :: self

    end subroutine
    !**********************************************

end module type_blockmatrix
