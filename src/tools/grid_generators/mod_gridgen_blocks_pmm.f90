module mod_gridgen_blocks_pmm
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, TWO, THREE, FOUR, SIX, PI
    use mod_string,             only: string_t
    use mod_bc,                 only: create_bc
    use mod_plot3d_utilities,   only: get_block_points_plot3d, &
                                      get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d
    use mod_hdf_utilities,      only: initialize_file_hdf, add_domain_hdf, &
                                      open_file_hdf, close_file_hdf, &
                                      open_domain_hdf, close_domain_hdf, &
                                      set_bc_patch_hdf, add_bc_state_hdf, &
                                      set_contains_grid_hdf, close_hdf, open_hdf, &
                                      create_bc_group_hdf, open_bc_group_hdf, close_bc_group_hdf, &
                                      set_bc_patch_group_hdf,                       &
                                      create_pmm_group_hdf, set_pmm_domain_group_hdf

    use type_point,             only: point_t
    use type_bc_group,          only: bc_group_t
    use type_bc_state,          only: bc_state_t
    use type_bc_state_wrapper,  only: bc_state_wrapper_t
    use hdf5
    implicit none





contains

    !-------------------------------------------------------------------------------------
    !!
    !!
    !!  Create File: Grid + BC's
    !!  -----------------------------
    !!  create_mesh_file__singleblock
    !!  create_mesh_file__multiblock
    !!  create_mesh_file__D2E8M1
    !!
    !!
    !!  Generate grid: point arrays
    !!  --------------------------
    !!  meshgen_1x1x1_linear
    !!  meshgen_1x1x1_unit_linear
    !!  meshgen_2x2x2_linear
    !!  meshgen_2x2x1_linear
    !!  meshgen_3x3x3_linear
    !!  meshgen_3x3x3_unit_linear
    !!  meshgen_3x3x1_linear
    !!  meshgen_4x1x1_linear
    !!  meshgen_4x2x2_linear
    !!  meshgen_3x1x1_linear
    !!  meshgen_2x1x1_linear
    !!  meshgen_40x15x1_linear
    !!  meshgen_15x15x1_linear
    !!  meshgen_15x15x2_linear
    !!  meshgen_15x15x3_linear
    !!
    !!
    !!
    !**************************************************************************************



    !>  Write a ChiDG-formatted grid file consisting of:
    !!
    !!      - one block domain, D1
    !!      - Linear element mapping, M1
    !!      - boundary conditions initialized to Scalar Extrapolate.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file__pmm__singleblock(filename,equation_sets,group_names,bc_groups,nelem_xi,nelem_eta,nelem_zeta,clusterx)
        character(*),               intent(in)              :: filename
        type(string_t),             intent(in), optional    :: equation_sets(:)
        type(string_t),             intent(in), optional    :: group_names(:,:)
        type(bc_group_t),           intent(in), optional    :: bc_groups(:)
        integer(ik),                intent(in)              :: nelem_xi
        integer(ik),                intent(in)              :: nelem_eta
        integer(ik),                intent(in)              :: nelem_zeta
        integer(ik),                intent(in), optional    :: clusterx

        character(:),                   allocatable :: user_msg
        class(bc_state_t),              allocatable :: bc_state
        character(len=10)                           :: patch_names(6)
        integer(HID_T)                              :: file_id, dom_id, patch_id, bcgroup_id
        integer(ik)                                 :: spacedim, mapping, bcface, ierr, igroup, istate
        type(point_t),                  allocatable :: nodes(:)
        integer(ik),                    allocatable :: elements(:,:) 
        integer(ik),                    allocatable :: faces(:,:)
        real(rk),   dimension(:,:,:),   allocatable :: xcoords, ycoords, zcoords


        ! Create/initialize file
        call initialize_file_hdf(filename)
        file_id = open_file_hdf(filename)
        


        ! Generate coordinates for first block
        call meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords,clusterx)



        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,idomain=1)


        !
        ! Add domains
        !
        spacedim = 3

        if ( present(equation_sets) ) then
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian',equation_sets(1)%get(),spacedim)
        else
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian','Scalar Advection',spacedim)
        end if


        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,'01')

        do bcface = 1,size(patch_names)
            
            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom_id,faces,bcface)

        end do !bcface


        !
        ! Create bc_state, 'Scalar Extrapolate'
        !
        call create_bc('Scalar Extrapolate', bc_state)


        !
        ! Add bc_group's
        !
        if (present(bc_groups)) then
            do igroup = 1,size(bc_groups)
                call create_bc_group_hdf(file_id,bc_groups(igroup)%name)

                bcgroup_id = open_bc_group_hdf(file_id,bc_groups(igroup)%name)

                do istate = 1,bc_groups(igroup)%bc_states%size()
                    call add_bc_state_hdf(bcgroup_id, bc_groups(igroup)%bc_states%at(istate))
                end do
                call close_bc_group_hdf(bcgroup_id)
            end do
        else
            call create_bc_group_hdf(file_id,'Default')

            bcgroup_id = open_bc_group_hdf(file_id,'Default')
            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)

        end if



        !
        ! Set boundary condition groups for each patch
        !
        patch_names = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(patch_names)

            call h5gopen_f(dom_id,'BoundaryConditions/'//trim(adjustl(patch_names(bcface))),patch_id,ierr)


            ! Set bc_group
            if (present(group_names)) then
                call set_bc_patch_group_hdf(patch_id,group_names(1,bcface)%get())
            else
                call set_bc_patch_group_hdf(patch_id,'Default')
            end if


            call h5gclose_f(patch_id,ierr)

        end do

        !
        !   Add prescribed mesh motion
        !

        !Add pmm group
        call create_pmm_group_hdf(file_id,'static')

        !Assign pmm to domain
        call set_pmm_domain_group_hdf(dom_id,'static')

        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')

        ! Close file
        call close_domain_hdf(dom_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__pmm__singleblock
    !*************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine create_mesh_file__multiblock(filename,equation_sets,group_names,bc_groups,     &
                                                          nelem_xi,  nelem_eta,  nelem_zeta,  &
                                                          xmax,ymax,zmax,clusterx)
        character(*),       intent(in)              :: filename
        type(string_t),     intent(in), optional    :: equation_sets(:)
        type(string_t),     intent(in), optional    :: group_names(:,:)
        type(bc_group_t),   intent(in), optional    :: bc_groups(:)
        integer(ik),        intent(in)              :: nelem_xi,  nelem_eta,  nelem_zeta
        real(rk),           intent(in), optional    :: xmax, ymax, zmax
        integer(ik),        intent(in), optional    :: clusterx

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(6)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface1_id, bcface2_id, bcgroup_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr, igroup, istate, &
                                                           nxi_max, neta_max, nzeta_max,xi_mid
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, ycoords1, zcoords1, &
                                                           xcoords2, ycoords2, zcoords2
        real(rk)                                        :: xmax_block1,xmin_block2, xmax_current,       &
                                                           xmax_new, ymax_new, zmax_new, ymax_current,  &
                                                           zmax_current


        !
        ! Create/initialize file
        !
        call initialize_file_hdf(filename)
        file_id = open_file_hdf(filename)
        

        !
        ! Generate coordinates
        !
        call meshgen_NxNxN_linear(nelem_xi, nelem_eta, nelem_zeta, xcoords1,ycoords1,zcoords1,clusterx)



        !
        ! Split the block:
        !   - Take the second half of the block, store to block 2
        !   - Take the first half of the block, reset to block 1
        !
        nxi_max   = size(xcoords1,1)
        neta_max  = size(xcoords1,2)
        nzeta_max = size(xcoords1,3)
        xi_mid    = int(nint(real(nxi_max)/2.))

        xcoords2 = xcoords1(xi_mid:nxi_max,1:neta_max,1:nzeta_max)
        xcoords1 = xcoords1(1:xi_mid,      1:neta_max,1:nzeta_max)

        ycoords2 = ycoords1(xi_mid:nxi_max,1:neta_max,1:nzeta_max)
        ycoords1 = ycoords1(1:xi_mid,      1:neta_max,1:nzeta_max)

        zcoords2 = zcoords1(xi_mid:nxi_max,1:neta_max,1:nzeta_max)
        zcoords1 = zcoords1(1:xi_mid,      1:neta_max,1:nzeta_max)


!        !
!        ! Translate block2 to end of block1, in x-direction
!        !
!        xmax_block1 = maxval(xcoords1)
!        xmin_block2 = minval(xcoords2)
!        xcoords2 = xcoords2 + (xmax_block1-xmin_block2)



        !
        ! Scale block coordinates to have xmax,ymax,zmax
        !   - normalize current maximum to 1
        !   - then multiply by desired extent
        !   - default extents are (1,1,1)
        !
        xmax_new = 1._rk
        ymax_new = 1._rk
        zmax_new = 1._rk
        if (present(xmax)) xmax_new = xmax
        if (present(ymax)) ymax_new = ymax
        if (present(zmax)) zmax_new = zmax


        xmax_current = maxval(xcoords2)
        xcoords1 = xcoords1/xmax_current
        xcoords2 = xcoords2/xmax_current

        xcoords1 = xmax_new*xcoords1
        xcoords2 = xmax_new*xcoords2


        ymax_current = maxval(ycoords2)
        ycoords1 = ycoords1/ymax_current
        ycoords2 = ycoords2/ymax_current

        ycoords1 = ymax_new*ycoords1
        ycoords2 = ymax_new*ycoords2

        zmax_current = maxval(zcoords2)
        zcoords1 = zcoords1/zmax_current
        zcoords2 = zcoords2/zmax_current

        zcoords1 = zmax_new*zcoords1
        zcoords2 = zmax_new*zcoords2



        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords1,zcoords1)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords2,zcoords2)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords1,zcoords1,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords2,zcoords2,mapping,idomain=2)


        !
        ! Add domains
        !
        spacedim = 3
        if ( present(equation_sets) ) then
            call add_domain_hdf(file_id,'01',nodes1,elements1,'Cartesian',equation_sets(1)%get(),spacedim)
            call add_domain_hdf(file_id,'02',nodes2,elements2,'Cartesian',equation_sets(2)%get(),spacedim)
        else
            call add_domain_hdf(file_id,'01',nodes1,elements1,'Cartesian','Scalar Advection',spacedim)
            call add_domain_hdf(file_id,'02',nodes2,elements2,'Cartesian','Scalar Advection',spacedim)
        end if






        !
        ! Set boundary conditions patch connectivities
        !
        dom1_id = open_domain_hdf(file_id,'01')
        dom2_id = open_domain_hdf(file_id,'02')

        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords1,zcoords1,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords2,zcoords2,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom1_id,faces1,bcface)
            call set_bc_patch_hdf(dom2_id,faces2,bcface)
        end do !bcface


        !
        ! Create bc_state, 'Scalar Extrapolate'
        !
        call create_bc('Scalar Extrapolate', bc_state)




        !
        ! Add bc_group's
        !
        if (present(bc_groups)) then
            do igroup = 1,size(bc_groups)
                call create_bc_group_hdf(file_id,bc_groups(igroup)%name)

                bcgroup_id = open_bc_group_hdf(file_id,bc_groups(igroup)%name)

                do istate = 1,bc_groups(igroup)%bc_states%size()
                    call add_bc_state_hdf(bcgroup_id, bc_groups(igroup)%bc_states%at(istate))
                end do
                call close_bc_group_hdf(bcgroup_id)
            end do
        else
            call create_bc_group_hdf(file_id,'Default')

            bcgroup_id = open_bc_group_hdf(file_id,'Default')
            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)

        end if




        !
        ! Assign groups to boundary condition patches
        !
        faces = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(faces)
            call h5gopen_f(dom1_id,'BoundaryConditions/'//trim(adjustl(faces(bcface))),bcface1_id,ierr)
            call h5gopen_f(dom2_id,'BoundaryConditions/'//trim(adjustl(faces(bcface))),bcface2_id,ierr)

            if (present(group_names)) then
                call set_bc_patch_group_hdf(bcface1_id,group_names(1,bcface)%get())
                call set_bc_patch_group_hdf(bcface2_id,group_names(2,bcface)%get())
            else
                call set_bc_patch_group_hdf(bcface1_id,'Default')
                call set_bc_patch_group_hdf(bcface2_id,'Default')
            end if

            call h5gclose_f(bcface1_id,ierr)
            call h5gclose_f(bcface2_id,ierr)
        end do



        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')

        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__multiblock
    !*************************************************************************************

























    !>  Write a ChiDG-formatted grid file consisting of:
    !!
    !!      - two block domains
    !!      - overlapping slightly
    !!      - boundary conditions initialized to Scalar Extrapolate. Interior boundaries
    !!        not set so they are detected as Chimera.
    !!
    !!  Incoming Parameter, 'matching' specifies if the elements should overlap with
    !!  a single, or potentially multiple elements
    !!
    !!     Block 1           Block 2 : matching=True        Block 2 : matching=False
    !!  .-----.-----.             .-----.-----.                  .-----.-----.
    !!  |     |     |             |     |     |                  |     |     |
    !!  |     |     |             |     |     |                  |     |     |
    !!  .-----.-----.             .-----.-----.                  |     |     |
    !!  |     |     |             |     |     |                  .-----.-----.
    !!  |     |     |             |     |     |                  |     |     |
    !!  .-----.-----.             .-----.-----.                  .-----.-----.
    !!
    !!  Abutting/Matching combinations:
    !!  -------------------------------
    !!
    !!       abutting = .true.        abutting = .false.
    !!       matching = .true.        matching = .true.
    !!          ----..----               ----.-.----
    !!              ||                       | | 
    !!              ||                       | | 
    !!          ----..----               ----.-.----
    !!              ||                       | |   
    !!              ||                       | | 
    !!          ----..----               ----.-.----
    !!
    !!
    !!       abutting = .true.         abutting = .false.
    !!       matching = .false.        matching = .false.
    !!  
    !!         -----..-----              ----.-.----
    !!              ||                       | | 
    !!              ||                       | | 
    !!         -----.|                   ----|-.   
    !!              |.-----                  .-|----
    !!              ||                       | | 
    !!         -----..-----              ----.-.----
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file__D2E8M1(filename,abutting,matching,equation_sets,group_names,bc_groups)
        character(*),       intent(in)              :: filename
        logical,            intent(in)              :: abutting
        logical,            intent(in)              :: matching
        type(string_t),     intent(in), optional    :: equation_sets(:)
        type(string_t),     intent(in), optional    :: group_names(:,:)
        type(bc_group_t),   intent(in), optional    :: bc_groups(:)

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(5)
        character(:),   allocatable                     :: user_msg
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface_id, bcgroup_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr, igroup, istate
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, xcoords2, ycoords1, ycoords2, zcoords1, zcoords2
        real(rk)                                        :: xmax,ymax, xmax_current, ymax_current, zmax_current


        ! Create/initialize file
        call initialize_file_hdf(filename)
        file_id = open_file_hdf(filename)
        

        ! Generate coordinates for first block
        call meshgen_2x2x2_linear(xcoords1,ycoords1,zcoords1)


        !
        ! Create second block by copying and translating first block.
        !
        xmax = maxval(xcoords1)
        ymax = maxval(ycoords1)


        !
        ! If abutting=.true., Create block2 by copying block1 and translating
        ! it by xmax of block1.
        !
        ! If abutting=.false., only translate by a fraction of xmax so there is overlap
        !
        if (abutting) then
            xcoords2 = xcoords1 + xmax
        else
            xcoords2 = xcoords1 + 0.90_rk*xmax
        end if
        ycoords2 = ycoords1


        !
        ! If matching=false, shift center plane of points so that the overlapping
        ! faces between blocks to not match exactly, and might have multiple Chimera donors
        !
        if (.not. matching) then
            ycoords2(:,2,:) = ycoords2(:,2,:) - 0.2*ymax
        end if


        !
        ! Set zcoords
        !
        zcoords2 = zcoords1


        !
        ! Scale coordinates to xmax, ymax, zmax
        !
        xmax_current = maxval(xcoords2)
        ymax_current = maxval(ycoords2)
        zmax_current = maxval(zcoords2)

        xcoords1 = (1.0_rk)*xcoords1/xmax_current
        ycoords1 = (1.0_rk)*ycoords1/ymax_current
        zcoords1 = (1.0_rk)*zcoords1/zmax_current

        xcoords2 = (1.0_rk)*xcoords2/xmax_current
        ycoords2 = (1.0_rk)*ycoords2/ymax_current
        zcoords2 = (1.0_rk)*zcoords2/zmax_current



        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords1,zcoords1)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords2,zcoords2)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords1,zcoords1,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords2,zcoords2,mapping,idomain=2)


        !
        ! Add domains
        !
        spacedim = 3
        call add_domain_hdf(file_id,'01',nodes1,elements1,'Cartesian','Scalar Advection',spacedim)
        call add_domain_hdf(file_id,'02',nodes2,elements2,'Cartesian','Scalar Advection',spacedim)


        !
        ! Set boundary conditions patch connectivities
        !
        dom1_id = open_domain_hdf(file_id,'01')
        dom2_id = open_domain_hdf(file_id,'02')

        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords1,zcoords1,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords2,zcoords2,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom1_id,faces1,bcface)
            call set_bc_patch_hdf(dom2_id,faces2,bcface)
        end do !bcface





        !
        ! Create bc_state, 'Scalar Extrapolate'
        !
        call create_bc('Scalar Extrapolate', bc_state)



        !
        ! Add bc_group's
        !
        if (present(bc_groups)) then
            user_msg = 'create_mesh_file__D2E8M1: Not quite ready to accept custom bc_group sets.'
            call chidg_signal(FATAL,user_msg)
        else
            call create_bc_group_hdf(file_id,'Default')

            bcgroup_id = open_bc_group_hdf(file_id,'Default')
            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)

        end if






        !
        ! Set boundary condition states for Domain 1: Leave XI_MAX empty
        !
        faces = ['XI_MIN  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(faces)
            call h5gopen_f(dom1_id,'BoundaryConditions/'//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call set_bc_patch_group_hdf(bcface_id,'Default')
            call h5gclose_f(bcface_id,ierr)
        end do

        !
        ! Set boundary condition states for Domain 2: Leave XI_MIN empty
        !
        faces = ['XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(faces)
            call h5gopen_f(dom2_id,'BoundaryConditions/'//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call set_bc_patch_group_hdf(bcface_id,'Default')
            call h5gclose_f(bcface_id,ierr)
        end do


        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')

        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__D2E8M1
    !*************************************************************************************










    !>  Generate a set of points defining a:
    !!      - nelem_xi by nelem_eta by nelem_zeta element, single-block, mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/20/2016
    !!
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !!
    !--------------------------------------------------------------------------------------
    subroutine meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords,clusterx)
        integer(ik)             :: nelem_xi
        integer(ik)             :: nelem_eta
        integer(ik)             :: nelem_zeta
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)
        integer(ik), optional   :: clusterx

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, &
                       npts_xi, npts_eta, npts_zeta
        real(rk)    :: x,y,z


        npts_xi   = nelem_xi   + 1
        npts_eta  = nelem_eta  + 1
        npts_zeta = nelem_zeta + 1


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi

                    if (present(clusterx)) then
                        if ( clusterx == -1 ) then
                            x = ONE - tanh( (PI/TWO)*(ONE - real(ipt_xi-1,rk)/real(npts_xi-1,rk) ) )/tanh(PI/TWO)
                        else if ( clusterx == 1 ) then
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: 'clusterx'=1 not yet implemented.")
                        else
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: Invalid value for 'clusterx'. -1,1.")
                        end if
                    else
                        x = real(ipt_xi-1,rk)/real(npts_xi-1,rk)
                    end if

                    if (ipt_xi == npts_xi) then
                        x = ONE
                    end if

                    y = real(ipt_eta-1,rk)/real(npts_eta-1,rk)
                    z = real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                end do
            end do
        end do


    end subroutine meshgen_NxNxN_linear
    !**************************************************************************************











    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2016
    !!
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, npts_xi, &
                       npts_eta, npts_zeta
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_xi   = 2
        npts_eta  = 2
        npts_zeta = 2

        dx = 1._rk
        dy = 1._rk
        dz = 1._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi


                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z


                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_1x1x1_linear
    !**************************************************************************************














    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_unit_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, npts_xi, &
                       npts_eta, npts_zeta
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear
        npts_xi   = 2
        npts_eta  = 2
        npts_zeta = 2

        dx = 2._rk
        dy = 2._rk
        dz = 2._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = -ONE
        do ipt_zeta = 1,npts_zeta
            y = -ONE
            do ipt_eta = 1,npts_eta
                x = -ONE
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_1x1x1_unit_linear
    !**************************************************************************************












    !> Generate a set of points defining a 2x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 27
        real(rk)                    :: x,y,z

        ! elements (2x2x2) - linear
        !
        !          *-------*-------*
        !         /       /       /|
        !        *-------*-------* |
        !       /       /       /| *
        !      *-------*-------* |/|
        !      |       |       | * |
        !      |       |       |/| *
        !      *-------*-------* |/
        !      |       |       | *
        !      |       |       |/
        !      *-------*-------*
        !
        !
        npts_xi   = 3
        npts_eta  = 3
        npts_zeta = 3



        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    x = ONE*real(ipt_xi  -1,rk)/real(npts_xi  -1,rk)
                    y = ONE*real(ipt_eta -1,rk)/real(npts_eta -1,rk)
                    z = ONE*real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_2x2x2_linear
    !**************************************************************************************













    !> Generate a set of points defining a 2x2x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 18
        real(rk), dimension(npt)    :: x,y,z

        ! elements (2x2x1) - linear
        !
        !          *-------*
        !         /       /|
        !        *-------* |
        !       /       /| * 
        !      *-------* |/|
        !      |       | * |
        !      |       |/| * 
        !      *-------* |/
        !      |       | * 
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 3
        npts_eta  = 3
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_2x2x1_linear
    !***************************************************************************************















    !> Generate a set of points defining a 3x3x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x3_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 64
        !real(rk), dimension(npt)    :: x,y,z
        real(rk)                    :: x,y,z

        ! elements (3x3x3) - linear
        !
        !            *-------*-------*-------*
        !           /       /       /       /|
        !          *-------*-------*-------* |
        !         /       /       /       /| *
        !        *-------*-------*-------* |/|
        !       /       /       /       /| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/|
        !      |       |       |       |/| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/
        !      |       |       |       |/| *
        !      *-------*-------*-------* |/
        !      |       |       |       | *
        !      |       |       |       |/
        !      *-------*-------*-------*
        !
        !
!        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE]
!
!        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE]
!
!        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
!             ZERO, ZERO, ZERO, ZERO, ZERO, &
!             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, & 
!             ONE, ONE, ONE, &
!             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
!             TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, &
!             THREE, THREE, THREE, THREE, THREE, THREE]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 4


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    x = (real(ipt_xi  -1,rk))/real(npts_xi  -1,rk)
                    y = (real(ipt_eta -1,rk))/real(npts_eta -1,rk)
                    z = (real(ipt_zeta-1,rk))/real(npts_zeta-1,rk)

                    !xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    !ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    !zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x3x3_linear
    !***************************************************************************************
















    !> Generate a set of points defining a 3x3x3 unit-element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x3_unit_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 64
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x3) - unit linear
        !
        !            *-------*-------*-------*
        !           /       /       /       /|
        !          *-------*-------*-------* |
        !         /       /       /       /| *
        !        *-------*-------*-------* |/|
        !       /       /       /       /| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/|
        !      |       |       |       |/| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/
        !      |       |       |       |/| *
        !      *-------*-------*-------* |/
        !      |       |       |       | *
        !      |       |       |       |/
        !      *-------*-------*-------*
        !
        !
        x = [ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX]

        y = [ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, &
             SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 4


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_3x3x3_unit_linear
    !***************************************************************************************












    !> Generate a set of points defining a 3x3x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 32
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x1) - linear
        !
        !            *-------*
        !           /       /| 
        !          *-------* |
        !         /       /| *   
        !        *-------* |/|
        !       /       /| * | 
        !      *-------* |/| *
        !      |       | * |/| 
        !      |       |/| * |
        !      *-------* |/| *
        !      |       | * |/  
        !      |       |/| *  
        !      *-------* |/
        !      |       | *
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x3x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_4x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 20
        real(rk), dimension(npt)    :: x,y,z

        ! elements (4x1x1) - linear
        !
        !      *------*------*------*------*
        !      |      |      |      |      | 
        !      |      |      |      |      | 
        !      *------*------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR]
             

        y = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 5
        npts_eta  = 2
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_4x1x1_linear
    !***************************************************************************************










    !> Generate a set of points defining a 4x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_4x2x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        real(rk)    :: x,y,z
        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, &
                       npts_xi, npts_eta, npts_zeta


        ! elements (4x2x2) - linear
        !  *---*---*---*---*
        !  |\   \   \   \   \
        !  * *---*---*---*---*
        !  |\|\   \   \   \   \
        !  * * *---*---*---*---*
        !   \|\|   |   |   |   | 
        !    * *---*---*---*---*
        !     \|   |   |   |   | 
        !      *---*---*---*---*
        !


        npts_xi   = 5
        npts_eta  = 3
        npts_zeta = 3


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    x = TWO*real(ipt_xi  -1,rk)/real(npts_xi  -1,rk)
                    y = ONE*real(ipt_eta -1,rk)/real(npts_eta -1,rk)
                    z = ONE*real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                end do
            end do
        end do



    end subroutine meshgen_4x2x2_linear
    !***************************************************************************************











    !> Generate a set of points defining a 2x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_2x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 12
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*
        !      |      |      | 
        !      |      |      | 
        !      *------*------*
        !



        x = [ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ONE, ONE, ONE]


        npts_xi   = 3
        npts_eta  = 2
        npts_zeta = 2

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_2x1x1_linear
    !**************************************************************************************












    !> Generate a set of points defining a 3x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 16
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*------*
        !      |      |      |      | 
        !      |      |      |      | 
        !      *------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE]
             

        y = [ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE]


        npts_xi   = 4
        npts_eta  = 2
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x1x1_linear
    !**************************************************************************************



















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_40x15x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (40x15x1) - linear

        npts_xi   = 41
        npts_eta  = 16
        npts_zeta = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do


    end subroutine meshgen_40x15x1_linear
    !**************************************************************************************













    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z,dx,dy,dz

        ! elements (15x15x1) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1.0_rk


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do


    end subroutine meshgen_15x15x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 15x15x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi,npts_eta,npts_zeta
        real(rk)    :: x,y,z,dx,dy,dz

        ! elements (15x15x2) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 3

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do




    end subroutine meshgen_15x15x2_linear
    !***************************************************************************************












    !> Generate a set of points defining a 15x15x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x3_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x3) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 4

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_15x15x3_linear
    !**************************************************************************************






end module mod_gridgen_blocks_pmm
