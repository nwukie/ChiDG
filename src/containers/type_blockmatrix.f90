!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockmatrix
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: DIAG, ZERO, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, &
                                      NFACES, CHIMERA
    use type_mesh,              only: mesh_t
    use type_densematrix,       only: densematrix_t
    use type_face_location,     only: face_location_t
    use type_element_location,  only: element_location_t
    use type_seed,              only: seed_t
    use DNAD_D
    implicit none


    !> Container for storing denseblock linearizations that make up a blockmatrix
    !!
    !!
    !!
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

        type(densematrix_t), allocatable :: lblks(:,:)                      !< Local domain blocks
        integer(ik),         allocatable :: ldata(:,:)                      !< Local block data     (nvars, nterms)

        type(densematrix_t), allocatable :: chiblks(:,:)                    !< Chimera inter-domain blocks (nChiElems, MaxDonors)


    contains
        ! Initializers
        generic,   public  :: init => initialize_linearization              !< Initialize full linearization matrix
        procedure, private :: initialize_linearization


        ! Setters
        procedure :: store                                                  !< Store linearization data for local blocks
        procedure :: store_chimera                                          !< Store linearization data for chimera blocks
        procedure :: clear                                                  !< Zero all data storage



        final :: destructor
    end type blockmatrix_t
    !*******************************************************************************************************************************



    private
contains



    !> Subroutine for initializing local linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    mesh_t containing arrays of elements and faces
    !!  @param[in]  mtype   character string indicating the type of matrix to be initialized (ie. Full, Lower-Diagonal, Upper-Diagonal
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine initialize_linearization(self,mesh,mtype)
        class(blockmatrix_t),   intent(inout)   :: self
        class(mesh_t),          intent(in)      :: mesh
        character(*),           intent(in)      :: mtype

        integer(ik), allocatable    :: blocks(:)
        integer(ik)                 :: nelem, nblk, ierr, ielem, iblk, size1d, parent, block_index, neqns, nterms_s
        integer(ik)                 :: nchimera_elements, maxdonors, idonor, iface, eparent, dparent
        integer(ik)                 :: iopen, ChiID, ndonors
        logical                     :: new_elements
        logical                     :: chimera_face
        logical                     :: more_donors
        logical                     :: donor_already_called
        logical                     :: contains_chimera_face
        logical                     :: init_chimera = .false.


        !
        ! Select matrix blocks to initialize 
        !
        select case (trim(mtype))
            case ('full','Full','FULL')
                blocks = [XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG]
                init_chimera = .true.

            case ('L','l','Lower','lower')
                blocks = [XI_MIN,ETA_MIN,ZETA_MIN]
                init_chimera = .false.

            case ('U','u','Upper','upper')
                blocks = [XI_MAX,ETA_MAX,ZETA_MAX]
                init_chimera = .false.

            case ('LD','ld','LowerDiagonal','lowerdiagonal')
                blocks = [XI_MIN,ETA_MIN,ZETA_MIN,DIAG]
                init_chimera = .false.
                
            case ('UD','ud','UpperDiagonal','upperdiagonal')
                blocks = [XI_MAX,ETA_MAX,ZETA_MAX,DIAG]
                init_chimera = .false.

            case default
                call chidg_signal(FATAL,'blockmatrix%init: unrecognized matrix type')

        end select



        nelem = mesh%nelem      ! Number of elements in the local block
        nblk  = 7               ! Number of potential blocks in the linearization for a given element (1D => 3, 2D => 5, 3D => 7)


        !
        ! Check to make sure the mesh numerics were initialized
        !
        if (.not. mesh%solInitialized) call chidg_signal(FATAL,'blockmatrix_t%initialize_linearization: Incoming mesh_t was not initialized. Make sure to call mesh%init_sol')


        !
        ! Allocate for 'localblocks'
        !
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




        !
        ! Assemble some Chimera data
        !
        if (allocated(self%chiblks)) deallocate(self%chiblks)
        
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
                ! Check for Chimera face type and add donor count
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

                allocate(self%chiblks(nelem, maxdonors), stat=ierr)
                if (ierr /= 0) call AllocationError

            end if
        end if












        !
        ! Loop through elements and call initialization for linearization denseblock matrices
        !
        do ielem = 1,mesh%nelem

            !
            ! Loop through 'blocks' and call initialization. localblocks
            !
            do block_index = 1,size(blocks)
                iblk = blocks(block_index)
                size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s

                !
                ! Parent is the element with respect to which the linearization is computed
                !
                dparent = mesh%idomain
                if (iblk == DIAG) then
                    eparent = mesh%elems(ielem)%ielem
                else
                    eparent = mesh%faces(ielem,iblk)%ineighbor
                end if


                !
                ! Call initialization procedure if parent is not 0 (0 meaning there is no parent for that block, probably a boundary)
                !
                if (eparent /= 0) then
                    call self%lblks(ielem,iblk)%init(size1d,dparent,eparent)

                    ! Store data about number of equations and number of terms in solution expansion
                    self%ldata(ielem,1) = mesh%elems(ielem)%neqns
                    self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
                end if

            end do


            !
            ! Call initialization for Chimera blocks
            !
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
                            neqns    = mesh%chimera%recv%data(ChiID)%donor_neqns%at(idonor)
                            nterms_s = mesh%chimera%recv%data(ChiID)%donor_nterms_s%at(idonor)
                            dparent  = mesh%chimera%recv%data(ChiID)%donor_domain%at(idonor)
                            eparent  = mesh%chimera%recv%data(ChiID)%donor_element%at(idonor)

                            size1d = neqns * nterms_s

                            !
                            ! Check if block initialization was already called for current donor
                            !
                            do iblk = 1,maxdonors
                                donor_already_called = ( dparent == self%chiblks(ielem,iblk)%dparent() .and. &
                                                         eparent == self%chiblks(ielem,iblk)%eparent() )
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
                                    if (.not. allocated(self%chiblks(ielem,iblk)%mat) ) then
                                        iopen = iblk
                                        exit
                                    end if
                                end do

                                !
                                ! Call block initialization
                                !
                                call self%chiblks(ielem,iopen)%init(size1d,dparent,eparent)

                            end if

                        end do ! idonor

                    end if

                end do ! iface

            end if  ! init_chimera



        end do ! ielem






    end subroutine initialize_linearization
    !*********************************************************************************************************************















    !>  Stores derivative data to the linearization matrix
    !!
    !!      -- Given the integral data computed from the spatial discretization,
    !!         store the derivative values from the AD data types
    !!
    !!  @author Nathan A. Wukie
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
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  face        face_location_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!
    !--------------------------------------------------------------------------------------------------------------------------------
    subroutine store_chimera(self,integral,face,seed,ivar)
        class(blockmatrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(face_location_t),      intent(in)      :: face
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar

        integer(ik) :: idom, idom_d, ielem, ielem_d, iblk, iarray
        integer(ik) :: irow, irow_start, donorblk, i
        integer(ik) :: neqns, nterms
        logical     :: block_match    = .false.
        logical     :: no_donor_block = .false.

        idom  = face%idomain
        ielem = face%ielement

        idom_d  = seed%idom
        ielem_d = seed%ielem


        !
        ! Get stored information for the block
        !
        neqns  = self%ldata(ielem,1)
        nterms = self%ldata(ielem,2)

        !
        ! Compute correct row offset for ivar
        !
        irow_start = ( (ivar - 1) * nterms )


        !
        ! Find donor block location 
        !
        donorblk = 0
        do iblk = 1,size(self%chiblks,2)
            block_match = ( (idom_d  == self%chiblks(ielem,iblk)%dparent()) .and. &
                            (ielem_d == self%chiblks(ielem,iblk)%eparent()) )

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
            self%chiblks(ielem,donorblk)%mat(irow,:) = self%chiblks(ielem,donorblk)%mat(irow,:) + integral(iarray)%xp_ad_
        end do



    end subroutine store_chimera
    !*******************************************************************************************************************************












    !>  Set all denseblock_t storage to zero
    !!
    !!  @author Nathan A. Wukie
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
            if (allocated(self%chiblks)) then
                do iblk = 1,size(self%chiblks,2)

                    !
                    ! Check if the block storage is actually allocated
                    !
                    if (allocated(self%chiblks(ielem,iblk)%mat)) then

                        !
                        ! If so, set to ZERO
                        !
                        self%chiblks(ielem,iblk)%mat = ZERO

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
