module type_domain_matrix
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: DIAG, ZERO, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, &
                                      NFACES, CHIMERA, NO_INTERIOR_NEIGHBOR
    use type_mesh,              only: mesh_t
    use type_densematrix,       only: densematrix_t
    use type_densematrix_vector,only: densematrix_vector_t
    use type_element_info,      only: element_info_t
    use type_seed,              only: seed_t
    use type_ivector,           only: ivector_t
    use mod_chidg_mpi,          only: IRANK
    use DNAD_D
    implicit none


    !> Container for storing denseblock linearizations that make up a domain_matrix
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
    type, public :: domain_matrix_t

        !
        ! Primary storage, (nelem,ntime)
        !
        type(densematrix_vector_t), allocatable :: lblks(:,:)       
        type(densematrix_vector_t), allocatable :: chi_blks(:,:)       
        type(densematrix_vector_t), allocatable :: bc_blks(:,:)       
        type(densematrix_vector_t), allocatable :: hb_blks(:,:)


        !
        ! Supporting data
        !
        integer(ik),            allocatable :: ldata(:,:)               ! Block-local data, (nelem,3) nvars, nterms, ntime.
        !integer(ik),            allocatable :: local_transpose(:,:)    ! Block index of the transposed location (nelem,6)
        !type(ivector_t),        allocatable :: local_transpose(:)
        type(ivector_t),        allocatable :: local_lower_blocks(:,:)    ! For each element, which blocks (1-6) are lower blocks
        type(ivector_t),        allocatable :: local_upper_blocks(:,:)    ! For each element, which blocks (1-6) are upper blocks


    contains

        ! Initializers
        generic,   public  :: init => initialize_linearization          ! Initialize full linearization matrix
        procedure, private :: initialize_linearization


        ! Setters
        procedure :: store              ! Store linearization data for local blocks
        procedure :: store_chimera      ! Store linearization data for chimera blocks
        procedure :: store_bc           ! Store linearization data for boundary condition blocks
        procedure :: store_hb           ! Store linearization data for coupling across harmonic balance levels
        procedure :: store_hb_element   ! Make direct contribution to matrix entry instead of from AD variables
        procedure :: clear              ! Zero all data storage

        ! Processors
        procedure :: restrict

        ! Auxiliary info
        procedure :: nelements
        procedure :: ntime



        final :: destructor

    end type domain_matrix_t
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
    subroutine initialize_linearization(self,mesh,idom,mtype)
        class(domain_matrix_t), intent(inout)   :: self
        class(mesh_t),          intent(in)      :: mesh
        integer(ik),            intent(in)      :: idom
        character(*),           intent(in)      :: mtype

        character(:),   allocatable :: user_msg
        integer(ik),    allocatable :: blocks(:)
        integer(ik)                 :: nelem, ierr, ielem, iblk, parent,                &
                                       block_index, nfields, nterms_s, ntime,             &
                                       nchimera_elements, maxdonors, idonor, iface,     &
                                       itime, dparent_g, dparent_l, eparent_g,          &
                                       eparent_l, tparent, parent_proc, eparent_l_trans,&
                                       imat_trans, ChiID, ndonors, max_coupled_elems,   &
                                       ncoupled_elems, icoupled_elem, icoupled_elem_bc, &
                                       ielem_bc, ibc, imat, group_ID, patch_ID,         &
                                       face_ID, elem_ID, idomain_l, ielement_l, itime_couple
        logical                     :: new_elements, chimera_face, more_donors,         &
                                       already_added, contains_chimera_face,            &
                                       block_initialized, lower_block, upper_block,     &
                                       transposed_block, domain_has_face, add_hb_coupling
        logical                     :: init_chimera = .false.
        logical                     :: init_bc      = .false.
        logical                     :: init_hb      = .false.

        type(densematrix_t)         :: temp_blk, temp1
        !
        ! Select matrix blocks to initialize 
        !
        select case (trim(mtype))
            case ('full','Full','FULL')
                blocks       = [XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG]
                init_chimera = .true.
                init_bc      = .true.
                init_hb      = .true.

            case ('L','l','Lower','lower')
                blocks       = [XI_MIN,ETA_MIN,ZETA_MIN]
                init_chimera = .false.
                init_bc      = .false.
                init_hb      = .false.

            case ('U','u','Upper','upper')
                blocks       = [XI_MAX,ETA_MAX,ZETA_MAX]
                init_chimera = .false.
                init_bc      = .false.
                init_hb      = .false.

            case ('LD','ld','LowerDiagonal','lowerdiagonal')
                blocks       = [XI_MIN,ETA_MIN,ZETA_MIN,DIAG]
                init_chimera = .false.
                init_bc      = .false.
                init_hb      = .false.
                
            case ('UD','ud','UpperDiagonal','upperdiagonal')
                blocks       = [XI_MAX,ETA_MAX,ZETA_MAX,DIAG]
                init_chimera = .false.
                init_bc      = .false.
                init_hb      = .false.

            case ('D', 'd', 'Diagonal', 'diagonal')
                blocks       = [DIAG]
                init_chimera = .false.
                init_bc      = .false.
                init_hb      = .false.

            case default
                call chidg_signal(FATAL,'domain_matrix%init: unrecognized matrix type')

        end select



        ! Check to make sure the domain numerics were initialized
        user_msg = "domain_matrix%initialize_linearization: Incoming domain_t was not &
                    initialized. Make sure to call domain%init_sol"
        if (.not. mesh%domain(idom)%solInitialized) call chidg_signal(FATAL,user_msg)



        nelem = mesh%domain(idom)%nelem      ! Number of elements in the local block
        ntime = mesh%domain(idom)%ntime      ! Number of time levels
        !
        ! Allocation for 'local blocks'
        !
        if (allocated(self%lblks)) deallocate(self%lblks)
        allocate(self%lblks(nelem,ntime),           &
                 self%ldata(nelem,3),               &
                 self%local_lower_blocks(nelem,ntime),    &
                 self%local_upper_blocks(nelem,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Allocation for 'chimera', 'bc' coupling
        !
        if (allocated(self%chi_blks)) deallocate(self%chi_blks)
        if (allocated(self%bc_blks))  deallocate(self%bc_blks)
        if (allocated(self%hb_blks))  deallocate(self%hb_blks)

        if (init_chimera) allocate(self%chi_blks(nelem,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        if (init_bc) allocate(self%bc_blks(nelem,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        if (init_hb) allocate(self%hb_blks(nelem,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError





        !
        ! Loop through elements and call initialization for INTERIOR, CHIMERA, and 
        ! denseblock matrices.
        !
        do ielem = 1,mesh%domain(idom)%nelem
            do itime = 1,mesh%domain(idom)%ntime 
                imat = 1
                
                ! Set the element indices that the densematrix_vector is associated with.
                self%lblks(ielem,itime)%idomain_g  = mesh%domain(idom)%elems(ielem)%idomain_g
                self%lblks(ielem,itime)%idomain_l  = mesh%domain(idom)%elems(ielem)%idomain_l
                self%lblks(ielem,itime)%ielement_g = mesh%domain(idom)%elems(ielem)%ielement_g
                self%lblks(ielem,itime)%ielement_l = mesh%domain(idom)%elems(ielem)%ielement_l
                self%lblks(ielem,itime)%mass       = mesh%domain(idom)%elems(ielem)%mass


                !--------------------------------------------
                !
                ! Initialization  --  INTERIOR coupling
                !
                !--------------------------------------------
                do block_index = 1,size(blocks)
                    iblk = blocks(block_index)
                    nterms_s = mesh%domain(idom)%elems(ielem)%nterms_s
                    nfields  = mesh%domain(idom)%elems(ielem)%nfields

                    !
                    ! Parent is the element with respect to which the linearization is computed
                    !
                    dparent_l = mesh%domain(idom)%idomain_l
                    if (iblk == DIAG) then
                        dparent_g   = mesh%domain(idom)%elems(ielem)%idomain_g
                        dparent_l   = mesh%domain(idom)%elems(ielem)%idomain_l
                        eparent_g   = mesh%domain(idom)%elems(ielem)%ielement_g
                        eparent_l   = mesh%domain(idom)%elems(ielem)%ielement_l
                        parent_proc = IRANK
                    else
                        dparent_g   = mesh%domain(idom)%faces(ielem,iblk)%ineighbor_domain_g
                        dparent_l   = mesh%domain(idom)%faces(ielem,iblk)%ineighbor_domain_l
                        eparent_g   = mesh%domain(idom)%faces(ielem,iblk)%ineighbor_element_g
                        eparent_l   = mesh%domain(idom)%faces(ielem,iblk)%ineighbor_element_l
                        parent_proc = mesh%domain(idom)%faces(ielem,iblk)%ineighbor_proc
                    end if



                    !
                    ! Call initialization procedure if parent is not 0 (0 meaning there is no parent, probably a boundary)
                    !
                    if (eparent_l /= NO_INTERIOR_NEIGHBOR) then

                        ! Initialize dense block
                        call temp_blk%init(nterms_s,nfields,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc,itime)
                        call self%lblks(ielem,itime)%push_back(temp_blk)

                        ! Store data about number of equations and number of terms in solution expansion
                        self%ldata(ielem,1) = mesh%domain(idom)%elems(ielem)%nfields
                        self%ldata(ielem,2) = mesh%domain(idom)%elems(ielem)%nterms_s
                        self%ldata(ielem,3) = mesh%domain(idom)%elems(ielem)%ntime


                        !
                        ! If off-diagonal, store block index as 'upper' or 'lower'
                        !
                        ! TODO: Add consideration for ntime and how that affects upper/lower status.
                        !
                        lower_block = ( (parent_proc == IRANK .and. eparent_l < ielem) .or. (parent_proc < IRANK) )
                        upper_block = ( (parent_proc == IRANK .and. eparent_l > ielem) .or. (parent_proc > IRANK) )

                        if ( lower_block ) then
                            call self%local_lower_blocks(ielem,itime)%push_back(imat)
                        else if ( upper_block ) then
                            call self%local_upper_blocks(ielem,itime)%push_back(imat)
                        end if

                        imat = imat + 1
                    end if

                end do ! init local
                !********************************************
                



                !--------------------------------------------
                !
                ! Initialization  --  CHIMERA coupling
                !
                !--------------------------------------------
                if (init_chimera) then
                    ! Set the element indices that the densematrix_vector is associated with.
                    self%chi_blks(ielem,itime)%idomain_g  = mesh%domain(idom)%elems(ielem)%idomain_g
                    self%chi_blks(ielem,itime)%idomain_l  = mesh%domain(idom)%elems(ielem)%idomain_l
                    self%chi_blks(ielem,itime)%ielement_g = mesh%domain(idom)%elems(ielem)%ielement_g
                    self%chi_blks(ielem,itime)%ielement_l = mesh%domain(idom)%elems(ielem)%ielement_l
                    self%chi_blks(ielem,itime)%mass       = mesh%domain(idom)%elems(ielem)%mass

                    do iface = 1,NFACES

                        !
                        ! If facetype is CHIMERA
                        !
                        chimera_face = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )
                        if (chimera_face) then

                            !
                            ! Get ChiID and number of donor elements
                            !
                            ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
                            ndonors = mesh%domain(idom)%chimera%recv(ChiID)%ndonors()

                            !
                            ! Call block initialization for each Chimera donor
                            !
                            do idonor = 1,ndonors
                                dparent_g   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_g
                                dparent_l   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_l
                                eparent_g   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_g
                                eparent_l   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_l
                                parent_proc = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%iproc
                                nfields     = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%nfields
                                nterms_s    = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%nterms_s


                                !
                                ! Check if block initialization was already called for current donor
                                !
                                already_added = .false.
                                do imat = 1,self%chi_blks(ielem,itime)%size()
                                    
                                    ! dummy densematrix to get a specific densematrix inside the chi_blks
                                    ! densematrix_vector temporary variable to access densematrix routine
                                    temp1 = self%chi_blks(ielem,itime)%at(imat)
                                    
                                    already_added = ( dparent_g == temp1%dparent_g() .and. &
                                                      dparent_l == temp1%dparent_l() .and. &
                                                      eparent_g == temp1%eparent_g() .and. &
                                                      eparent_l == temp1%eparent_l() )
                                    if (already_added) exit
                                end do

                            
                                !
                                ! If a block for the donor element hasn't yet been initialized, call initialization procedure
                                !
                                if (.not. already_added) then
                                    call temp_blk%init(nterms_s,nfields,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc,itime)
                                    call self%chi_blks(ielem,itime)%push_back(temp_blk)

                                end if

                            end do ! idonor
                        end if ! if chimera
                    end do ! iface
                end if  ! init_chimera
                !********************************************
            
            end do  ! itime
        end do ! ielem



        !--------------------------------------------
        !
        ! Initialization  --  BOUNDARY coupling
        !
        !--------------------------------------------
        if (init_bc) then
            do itime = 1,mesh%domain(idom)%ntime 

                do group_ID = 1,mesh%nbc_patch_groups()
                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
                        do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                            ! Get indices of the local element to determine if it is on 'idom'
                            idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                            ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                            domain_has_face = (idom == idomain_l)


                            !
                            ! If domain contains current face, add all coupling information.
                            !   ASSUMPTIONS:
                            !       - all coupled elements use the same number of equations
                            !       - all coupled elements run with the same nterms_s
                            !
                            if (domain_has_face) then
                                do elem_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)

                                    ! Set the element indices that the densematrix_vector is associated with.
                                    self%bc_blks(ielement_l,itime)%idomain_g  = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g
                                    self%bc_blks(ielement_l,itime)%idomain_l  = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l
                                    self%bc_blks(ielement_l,itime)%ielement_g = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g
                                    self%bc_blks(ielement_l,itime)%ielement_l = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l
                                    self%bc_blks(ielement_l,itime)%mass       = mesh%domain(idomain_l)%elems(ielement_l)%mass



                                    !
                                    ! Compute size of coupling matrix
                                    !
                                    nfields  = mesh%domain(idomain_l)%elems(ielement_l)%nfields
                                    nterms_s = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s


                                    !
                                    ! Get indices for coupled element
                                    !
                                    dparent_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(elem_ID)
                                    dparent_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(elem_ID)
                                    eparent_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(elem_ID)
                                    eparent_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(elem_ID)
                                    parent_proc = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(elem_ID)


                                    !
                                    ! We need to check that we dont try to add the local element more than once.
                                    ! For example, if an element had two different boundary conditions on different
                                    ! faces, they would both try to add the local element.
                                    !
                                    already_added = .false.
                                    do imat = 1,self%bc_blks(ielement_l,itime)%size()
                                        
                                        ! dummy densematrix to get a specific densematrix inside the bc_blks
                                        ! densematrix_vector temporary variable to access densematrix routine
                                        temp1 = self%bc_blks(ielement_l,itime)%at(imat)
                                        
                                        already_added = ( dparent_g == temp1%dparent_g() .and. &
                                                          dparent_l == temp1%dparent_l() .and. &
                                                          eparent_g == temp1%eparent_g() .and. &
                                                          eparent_l == temp1%eparent_l() )
                                        if (already_added) exit
                                    end do



                                    !
                                    ! Call initialization, store initialized matrix to bc_blks
                                    !
                                    if (.not. already_added) then
                                        call temp_blk%init(nterms_s,nfields,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc,itime)
                                        call self%bc_blks(ielement_l,itime)%push_back(temp_blk)
                                    end if

                                end do !elem_ID, coupling
                            end if !domain_has_face

                        end do !face_ID
                    end do !patch_ID
                end do !group_ID
            end do !itime

        end if !init_bc




        !---------------------------------------------------------------
        !
        ! Initialization  --  HARMONIC BALANCE cross time-level coupling
        !
        !---------------------------------------------------------------
        if (init_hb) then

            ! Default coupling for each element with all its other time levels
            do ielem = 1,mesh%domain(idom)%nelem
                do itime = 1,mesh%domain(idom)%ntime 
                    ! Set the element indices that the densematrix_vector is associated with.
                    self%hb_blks(ielem,itime)%idomain_g  = mesh%domain(idom)%elems(ielem)%idomain_g
                    self%hb_blks(ielem,itime)%idomain_l  = mesh%domain(idom)%elems(ielem)%idomain_l
                    self%hb_blks(ielem,itime)%ielement_g = mesh%domain(idom)%elems(ielem)%ielement_g
                    self%hb_blks(ielem,itime)%ielement_l = mesh%domain(idom)%elems(ielem)%ielement_l
                    self%hb_blks(ielem,itime)%mass       = mesh%domain(idom)%elems(ielem)%mass

                    ! Get problem size for a given time-level
                    nfields  = mesh%domain(idom)%elems(ielem)%nfields
                    nterms_s = mesh%domain(idom)%elems(ielem)%nterms_s

                    ! Get indices for coupled element. Actually coupled with itself, but at different time-level
                    dparent_g   = mesh%domain(idom)%elems(ielem)%idomain_g
                    dparent_l   = mesh%domain(idom)%elems(ielem)%idomain_l
                    eparent_g   = mesh%domain(idom)%elems(ielem)%ielement_g
                    eparent_l   = mesh%domain(idom)%elems(ielem)%ielement_l
                    parent_proc = IRANK

                    ! Call initialization, store initialized matrix to hb_blks
                    do itime_couple = 1,mesh%domain(idom)%ntime
                        if (itime_couple /= itime) then
                            call temp_blk%init(nterms_s,nfields,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc,itime_couple)
                            call self%hb_blks(ielem,itime)%push_back(temp_blk)
                        end if
                    end do !itime_couple

                end do !itime
            end do !ielem


            ! Additional cross time-level coupling between other elements due to boundary condition
            do itime = 1,mesh%domain(idom)%ntime 
                do group_ID = 1,mesh%nbc_patch_groups()
                    do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()

                        ! Check if the patch requires cross-timelevel coupling
                        add_hb_coupling = (mesh%bc_patch_group(group_ID)%patch(patch_ID)%temporal_coupling == 'Global')

                        if (add_hb_coupling) then
                            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                                ! Get indices of the local element to determine if it is on 'idom'
                                idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                                ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                                domain_has_face = (idom == idomain_l)


                                ! If domain contains current face, add all coupling information.
                                !   ASSUMPTIONS:
                                !       - all coupled elements use the same number of equations
                                !       - all coupled elements run with the same nterms_s
                                if (domain_has_face) then
                                    do elem_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%ncoupled_elements(face_ID)
                                        ! Add coupling with all time-levels of spatially-coupled elements 
                                        ! on the boundary.
                                        do itime_couple = 1,mesh%domain(idom)%ntime
                                            if (itime_couple /= itime) then

                                                ! Compute size of coupling matrix
                                                nfields  = mesh%domain(idomain_l)%elems(ielement_l)%nfields
                                                nterms_s = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s

                                                ! Get indices for coupled element
                                                dparent_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g(elem_ID)
                                                dparent_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l(elem_ID)
                                                eparent_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(elem_ID)
                                                eparent_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(elem_ID)
                                                parent_proc = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(elem_ID)

                                                ! We need to check that we dont try to add the local element more than once.
                                                ! For example, if an element had two different boundary conditions on different
                                                ! faces, they would both try to add the local element.
                                                already_added = .false.
                                                do imat = 1,self%hb_blks(ielement_l,itime)%size()
                                                    
                                                    ! dummy densematrix to get a specific densematrix inside the hb_blks
                                                    ! densematrix_vector temporary variable to access densematrix routine
                                                    temp1 = self%hb_blks(ielement_l,itime)%at(imat)
                                                    
                                                    already_added = ( dparent_g    == temp1%dparent_g() .and. &
                                                                      dparent_l    == temp1%dparent_l() .and. &
                                                                      eparent_g    == temp1%eparent_g() .and. &
                                                                      eparent_l    == temp1%eparent_l() .and. &
                                                                      itime_couple == temp1%tparent() )
                                                    if (already_added) exit
                                                end do

                                                ! Call initialization, store initialized matrix to hb_blks
                                                if (.not. already_added) then
                                                    call temp_blk%init(nterms_s,nfields,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc,itime_couple)
                                                    call self%hb_blks(ielement_l,itime)%push_back(temp_blk)
                                                end if

                                            end if
                                        end do !itime_couple
                                    end do !elem_ID, coupling
                                end if !domain_has_face

                            end do !face_ID
                        end if !add_hb_coupling
                    end do !patch_ID
                end do !group_ID
            end do !itime


        end if !init_hb







        ! 
        ! Initialize transpose data
        !
        select case (trim(mtype))
            case ('full','Full','FULL')

                do ielem = 1,mesh%domain(idom)%nelem
                    do itime = 1,mesh%domain(idom)%ntime
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
                                             (transposed_block .eqv. .false.) ) call chidg_signal(FATAL,"domain_matrix%init: no transposed element found")
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
    subroutine store(self,integral,element_info,seed,ivar,itime)
        class(domain_matrix_t), intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        type(element_info_t),   intent(in)      :: element_info
        type(seed_t),           intent(in)      :: seed
        integer(ik),            intent(in)      :: ivar
        integer(ik),            intent(in)      :: itime

        integer(ik) :: nterms, imat, ielement_l

        ielement_l = element_info%ielement_l

        ! Get stored information for the block
        nterms = self%ldata(ielement_l,2)  

        ! Find donor densematrix location 
        imat = self%lblks(ielement_l,itime)%loc(seed%idomain_g,seed%ielement_g,seed%itime)
        if (imat == 0) call chidg_signal_three(FATAL,"domain_matrix%store: no allocation found to store data.",seed%idomain_g,seed%ielement_g,seed%itime)

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
    subroutine store_chimera(self,integral,element_info,seed,ivar,itime)
        class(domain_matrix_t),       intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(element_info_t),          intent(in)      :: element_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar
        integer(ik),                intent(in)      :: itime

        integer(ik) :: nterms, imat, ielement_l

        ielement_l = element_info%ielement_l

        ! Get stored information for the block
        nterms = self%ldata(ielement_l,2)

        ! Find donor densematrix location 
        imat = self%chi_blks(ielement_l,itime)%loc(seed%idomain_g,seed%ielement_g,seed%itime)

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
    !!  @param[in]  integral    Array of modes with embedded partial derivatives 
    !!  @param[in]  face        element_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!  @param[in]  itime       Index of a time level for the linearization of the given element
    !!
    !------------------------------------------------------------------------------------
    subroutine store_bc(self,integral,element_info,seed,ivar,itime)
        class(domain_matrix_t),     intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(element_info_t),       intent(in)      :: element_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar
        integer(ik),                intent(in)      :: itime

        integer(ik) :: imat, nterms
        logical     :: local_coupling = .false.

        ! If cross-timelevel harmonic balance term store separately
        if (itime /= seed%itime) then
            call self%store_hb(integral,element_info,seed,ivar,itime)
        else

            ! If ielem = ielem_d then the linearization is with respect to the local element. 
            ! So, this is stored in the self%lblks array in the DIAG location, instead of
            ! the self%bc_blks array. In general, the storage location is not important,
            ! but the ILU preconditioner expects the full diagonal contribution to be in 
            ! lblks.
            local_coupling = (element_info%idomain_g == seed%idomain_g) .and. (element_info%ielement_g == seed%ielement_g)

            if ( local_coupling ) then
                call self%store(integral,element_info,seed,ivar,itime)
            else

                ! Get stored information for the block
                nterms = self%ldata(element_info%ielement_l,2)

                ! Find coupled bc densematrix location 
                imat = self%bc_blks(element_info%ielement_l,itime)%loc(seed%idomain_g,seed%ielement_g,seed%itime)

                ! Store derivatives
                call self%bc_blks(element_info%ielement_l,itime)%store(imat,ivar,nterms,integral)

            end if ! check local block.

        end if

    end subroutine store_bc
    !********************************************************************************





    !>  Stores data due to harmonic balance cross-timelevel coupling from boundary conditions.
    !!
    !!  This means seed%itime /= itime
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2018
    !!
    !!  @param[in]  integral        Array of modes with embedded partial derivatives 
    !!  @param[in]  element_info    element_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed            seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar            Index of the variable
    !!  @param[in]  itime           Index of a time level for the linearization of the given element
    !!
    !------------------------------------------------------------------------------------
    subroutine store_hb(self,integral,element_info,seed,ivar,itime)
        class(domain_matrix_t),     intent(inout)   :: self
        type(AD_D),                 intent(in)      :: integral(:)
        type(element_info_t),       intent(in)      :: element_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: ivar
        integer(ik),                intent(in)      :: itime

        integer(ik) :: ielement_l, imat, nterms

        ! Get stored information for the block
        ielement_l = element_info%ielement_l
        nterms = self%ldata(ielement_l,2)

        ! Find coupled hb densematrix location 
        imat = self%hb_blks(ielement_l,itime)%loc(seed%idomain_g,seed%ielement_g,seed%itime)

        ! Store derivatives
        call self%hb_blks(ielement_l,itime)%store(imat,ivar,nterms,integral)

    end subroutine store_hb
    !********************************************************************************







    !>  Stores data due to harmonic balance cross-timelevel coupling from boundary conditions.
    !!
    !!  This means seed%itime /= itime
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2018
    !!
    !!  @param[in]  integral    Array of modes with embedded partial derivatives 
    !!  @param[in]  face        face_info_t containing indices for the location of the face being linearized.
    !!  @param[in]  seed        seed_t containing indices of the element against which the linearization was computed.
    !!  @param[in]  ivar        Index of the variable
    !!  @param[in]  itime       Index of a time level for the linearization of the given element
    !!
    !------------------------------------------------------------------------------------
    subroutine store_hb_element(self,contribution,element_info,seed,itime)
        class(domain_matrix_t),     intent(inout)   :: self
        real(rk),                   intent(in)      :: contribution(:,:)
        type(element_info_t),       intent(in)      :: element_info
        type(seed_t),               intent(in)      :: seed
        integer(ik),                intent(in)      :: itime

        integer(ik) :: ielement_l, imat

        ! Get stored information for the block
        ielement_l = element_info%ielement_l

        ! Find coupled hb densematrix location 
        imat = self%hb_blks(ielement_l,itime)%loc(seed%idomain_g,seed%ielement_g,seed%itime)

        ! Store derivatives
        call self%hb_blks(ielement_l,itime)%store_element(imat,contribution)

    end subroutine store_hb_element
    !********************************************************************************









    !>  Set all denseblock_t storage to zero
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date   11/14/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine clear(self)
        class(domain_matrix_t),   intent(inout)   :: self

        integer(ik) :: ielem, itime

        ! For each element
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

            ! Clear hb matrices
            if (allocated(self%hb_blks)) then
                do itime = 1,size(self%hb_blks,2)
                    call self%hb_blks(ielem,itime)%setzero()
                end do ! itime
            end if

        end do ! ielem

    end subroutine clear
    !*******************************************************************************





    !>  Return a domain_matrix instance that is restricted to the specified
    !!  polynomial space.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/23/2017
    !!
    !!
    !------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(domain_matrix_t), intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(domain_matrix_t)   :: restricted
        integer(ik)             :: ierr, ielem, itime


        !
        ! Copy auxiliary data
        !
        restricted%ldata = self%ldata
        restricted%local_lower_blocks = self%local_lower_blocks
        restricted%local_upper_blocks = self%local_upper_blocks



        !
        ! Allocate storage for primary data on restricted object
        !
        if (allocated(self%lblks)) then
            allocate(restricted%lblks(size(self%lblks,1),size(self%lblks,2)), stat=ierr)
            if (ierr /= 0) call AllocationError

            do ielem = 1,size(self%lblks,1)
                do itime = 1,size(self%lblks,2)
                    restricted%lblks(ielem,itime) = self%lblks(ielem,itime)%restrict(nterms_r)
                end do !itime
            end do !ielem
        end if


        if (allocated(self%chi_blks)) then
            allocate(restricted%chi_blks(size(self%chi_blks,1),size(self%chi_blks,2)), stat=ierr)
            if (ierr /= 0) call AllocationError

            do ielem = 1,size(self%chi_blks,1)
                do itime = 1,size(self%chi_blks,2)
                    restricted%chi_blks(ielem,itime) = self%chi_blks(ielem,itime)%restrict(nterms_r)
                end do !itime
            end do !ielem
        end if


        if (allocated(self%bc_blks)) then
            allocate(restricted%bc_blks(size(self%bc_blks,1),size(self%bc_blks,2)), stat=ierr)
            if (ierr /= 0) call AllocationError

            do ielem = 1,size(self%bc_blks,1)
                do itime = 1,size(self%bc_blks,2)
                    restricted%bc_blks(ielem,itime) = self%bc_blks(ielem,itime)%restrict(nterms_r)
                end do !itime
            end do !ielem
        end if


        if (allocated(self%hb_blks)) then
            allocate(restricted%hb_blks(size(self%hb_blks,1),size(self%hb_blks,2)), stat=ierr)
            if (ierr /= 0) call AllocationError

            do ielem = 1,size(self%hb_blks,1)
                do itime = 1,size(self%hb_blks,2)
                    restricted%hb_blks(ielem,itime) = self%hb_blks(ielem,itime)%restrict(nterms_r)
                end do !itime
            end do !ielem
        end if


        !
        ! Update ldata with correct number of terms
        !
        do ielem = 1,size(self%ldata,1)
            restricted%ldata(ielem,2) = nterms_r
        end do !ielem



    end function restrict
    !******************************************************************************




    !>  Return number of elements in the domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/27/2017
    !!
    !--------------------------------------------------------------------------------
    function nelements(self) result(nelements_)
        class(domain_matrix_t), intent(in)  :: self

        integer(ik) :: nelements_

        nelements_ = size(self%lblks,1)

    end function nelements
    !********************************************************************************





    !>  Return number of elements in the domain.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/27/2017
    !!
    !--------------------------------------------------------------------------------
    function ntime(self) result(ntime_)
        class(domain_matrix_t), intent(in)  :: self

        integer(ik) :: ntime_

        ntime_ = size(self%lblks,2)

    end function ntime
    !********************************************************************************





    !>
    !!
    !!---------------------------------------------
    subroutine destructor(self)
        type(domain_matrix_t), intent(inout) :: self

    end subroutine
    !**********************************************





end module type_domain_matrix
