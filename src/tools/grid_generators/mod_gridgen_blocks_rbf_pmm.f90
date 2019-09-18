module mod_gridgen_blocks_rbf_pmm
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, TWO, THREE, FOUR, SIX, PI
    use mod_string,             only: string_t
    use mod_bc,                 only: create_bc
    use mod_plot3d_utilities,   only: get_block_points_plot3d, &
                                      get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d

    use type_point,             only: point_t
    use type_bc_state_group,    only: bc_state_group_t
    use type_bc_state,          only: bc_state_t
    use type_bc_state_wrapper,  only: bc_state_wrapper_t

    use mod_hdf_utilities
    use mod_prescribed_mesh_motion_function,    only: create_prescribed_mesh_motion_function   
    use type_prescribed_mesh_motion_function,   only: prescribed_mesh_motion_function_t
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
    subroutine create_mesh_file__rbf_pmm__singleblock(filename,equation_sets,group_names,bc_state_groups,nelem_xi,nelem_eta,nelem_zeta,clusterx)
        character(*),               intent(in)              :: filename
        type(string_t),             intent(in), optional    :: equation_sets(:)
        type(string_t),             intent(in), optional    :: group_names(:,:)
        type(bc_state_group_t),     intent(in), optional    :: bc_state_groups(:)
        integer(ik),                intent(in)              :: nelem_xi
        integer(ik),                intent(in)              :: nelem_eta
        integer(ik),                intent(in)              :: nelem_zeta
        integer(ik),                intent(in), optional    :: clusterx

        character(:),                   allocatable :: user_msg
        class(bc_state_t),              allocatable :: bc_state
        character(len=10)                           :: patch_names(6)
        integer(HID_T)                              :: file_id, dom_id, patch_id, bcgroup_id, mmgroup_id, rbffcn_id, srcgroup_id, driver_id
        integer(ik)                                 :: mapping, bcface, ierr, igroup, istate
        real(rk),                       allocatable :: nodes(:,:)
        integer(ik),                    allocatable :: elements(:,:) 
        integer(ik),                    allocatable :: faces(:,:)
        real(rk),   dimension(:,:,:),   allocatable :: xcoords, ycoords, zcoords
        character(len=8)                            :: bc_face_strings(6)
        character(:),   allocatable                 :: bc_face_string
        class(prescribed_mesh_motion_function_t), allocatable    :: pmmf

        real(rk)                                    :: radius(3)
    
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
        if ( present(equation_sets) ) then
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian',equation_sets(1)%get())
        else
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian','Scalar Advection')
        end if


        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,'01')

        bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]
        do bcface = 1,size(patch_names)
            
            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            bc_face_string  = trim(bc_face_strings(bcface))
            patch_id = create_patch_hdf(dom_id,bc_face_string)
            call set_patch_hdf(patch_id,faces)
            call close_patch_hdf(patch_id)



        end do !bcface


        !
        ! Create bc_state, 'Scalar Extrapolate'
        !
        call create_bc('Scalar Extrapolate', bc_state)


        !
        ! Add bc_group's
        !
        if (present(bc_state_groups)) then
            do igroup = 1,size(bc_state_groups)
                call create_bc_state_group_hdf(file_id,bc_state_groups(igroup)%name)

                bcgroup_id = open_bc_state_group_hdf(file_id,bc_state_groups(igroup)%name)

                do istate = 1,bc_state_groups(igroup)%nbc_states()
                    call add_bc_state_hdf(bcgroup_id, bc_state_groups(igroup)%bc_state(istate)%state)
                end do
                call close_bc_state_group_hdf(bcgroup_id)
            end do
        else
            call create_bc_state_group_hdf(file_id,'Default')

            bcgroup_id = open_bc_state_group_hdf(file_id,'Default')
            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_state_group_hdf(bcgroup_id)

        end if



        !
        ! Set boundary condition groups for each patch
        !
        patch_names = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(patch_names)


            patch_id = open_patch_hdf(dom_id,trim(patch_names(bcface)))

            ! Set bc_group
            if (present(group_names)) then
                call set_patch_group_hdf(patch_id,group_names(1,bcface)%get())
            else
                call set_patch_group_hdf(patch_id,'Default')
            end if


            call close_patch_hdf(patch_id)

        end do

        !
        !   Add prescribed mesh motion
        !

        !Add pmm group
        call create_mm_group_hdf(file_id,'static','RBF')
        mmgroup_id = open_mm_group_hdf(file_id, 'static')
        ! Add RBF function groups
        rbffcn_id = create_rbf_function_hdf(mmgroup_id)
        call set_rbf_type_hdf(rbffcn_id,'wc6_rbf')
        radius = ONE
        call set_rbf_base_radius_hdf(rbffcn_id, radius)
        ! Add an RBF source
        call create_rbf_src_group_hdf(mmgroup_id, 'mysrc')
        srcgroup_id = open_rbf_src_hdf(mmgroup_id, 'mysrc')
        if (present(group_names)) then
            call set_rbf_src_patchname_hdf(srcgroup_id, group_names(1,bcface)%get())
        else
            call set_rbf_src_patchname_hdf(srcgroup_id, 'Default')
        end if

        ! Add MM driver to the source
        call create_rbf_mm_driver_hdf(srcgroup_id, 'PMM')
        driver_id = open_rbf_mm_driver_hdf(srcgroup_id)

        call create_prescribed_mesh_motion_function(pmmf, 'static')
        call add_pmmf_hdf(driver_id, pmmf)
        call close_rbf_mm_driver_hdf(driver_id)
        call close_rbf_src_hdf(srcgroup_id)
        call close_mm_group_hdf(mmgroup_id)


        !Assign pmm to domain
        call set_mm_domain_group_hdf(dom_id,'static')

        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')

        ! Close file
        call close_domain_hdf(dom_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__rbf_pmm__singleblock
    !*************************************************************************************



    
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
    subroutine create_mesh_file__rbf_pmm__multiscale__singleblock(filename,equation_sets,group_names,bc_state_groups,nelem_xi,nelem_eta,nelem_zeta,clusterx)
        character(*),               intent(in)              :: filename
        type(string_t),             intent(in), optional    :: equation_sets(:)
        type(string_t),             intent(in), optional    :: group_names(:,:)
        type(bc_state_group_t),     intent(in), optional    :: bc_state_groups(:)
        integer(ik),                intent(in)              :: nelem_xi
        integer(ik),                intent(in)              :: nelem_eta
        integer(ik),                intent(in)              :: nelem_zeta
        integer(ik),                intent(in), optional    :: clusterx

        character(:),                   allocatable :: user_msg
        class(bc_state_t),              allocatable :: bc_state
        character(len=10)                           :: patch_names(6)
        integer(HID_T)                              :: file_id, dom_id, patch_id, bcgroup_id, mmgroup_id, rbffcn_id, srcgroup_id, driver_id
        integer(ik)                                 :: mapping, bcface, ierr, igroup, istate
        real(rk),                       allocatable :: nodes(:,:)
        integer(ik),                    allocatable :: elements(:,:) 
        integer(ik),                    allocatable :: faces(:,:)
        real(rk),   dimension(:,:,:),   allocatable :: xcoords, ycoords, zcoords
        character(len=8)                            :: bc_face_strings(6)
        character(:),   allocatable                 :: bc_face_string
        class(prescribed_mesh_motion_function_t), allocatable    :: pmmf

        real(rk)                                    :: radius(3)
    
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
        if ( present(equation_sets) ) then
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian',equation_sets(1)%get())
        else
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian','Scalar Advection')
        end if


        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,'01')

        bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]
        do bcface = 1,size(patch_names)
            
            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            bc_face_string  = trim(bc_face_strings(bcface))
            patch_id = create_patch_hdf(dom_id,bc_face_string)
            call set_patch_hdf(patch_id,faces)
            call close_patch_hdf(patch_id)



        end do !bcface


        !
        ! Create bc_state, 'Scalar Extrapolate'
        !
        call create_bc('Scalar Extrapolate', bc_state)


        !
        ! Add bc_group's
        !
        if (present(bc_state_groups)) then
            do igroup = 1,size(bc_state_groups)
                call create_bc_state_group_hdf(file_id,bc_state_groups(igroup)%name)

                bcgroup_id = open_bc_state_group_hdf(file_id,bc_state_groups(igroup)%name)

                do istate = 1,bc_state_groups(igroup)%nbc_states()
                    call add_bc_state_hdf(bcgroup_id, bc_state_groups(igroup)%bc_state(istate)%state)
                end do
                call close_bc_state_group_hdf(bcgroup_id)
            end do
        else
            call create_bc_state_group_hdf(file_id,'Default')

            bcgroup_id = open_bc_state_group_hdf(file_id,'Default')
            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_state_group_hdf(bcgroup_id)

        end if



        !
        ! Set boundary condition groups for each patch
        !
        patch_names = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(patch_names)


            patch_id = open_patch_hdf(dom_id,trim(patch_names(bcface)))

            ! Set bc_group
            if (present(group_names)) then
                call set_patch_group_hdf(patch_id,group_names(1,bcface)%get())
            else
                call set_patch_group_hdf(patch_id,'Default')
            end if


            call close_patch_hdf(patch_id)

        end do

        !
        !   Add prescribed mesh motion
        !

        !Add pmm group
        call create_mm_group_hdf(file_id,'static','RBF')
        mmgroup_id = open_mm_group_hdf(file_id, 'static')
        ! Add RBF function groups
        rbffcn_id = create_rbf_function_hdf(mmgroup_id)
        call set_rbf_type_hdf(rbffcn_id,'wc6_rbf')
        radius = ONE
        call set_rbf_base_radius_hdf(rbffcn_id, radius)
        call close_rbf_function_hdf(rbffcn_id)
        ! Add explicit RBF function group to enable MS-RBF algorithm
        rbffcn_id = create_rbf_function_exp_hdf(mmgroup_id)
        call set_rbf_type_explicit_hdf(rbffcn_id,'wc6_rbf')
        call set_rbf_base_fraction_hdf(rbffcn_id, 0.05_rk)
        call close_rbf_function_exp_hdf(rbffcn_id)
        ! Add an RBF source
        call create_rbf_src_group_hdf(mmgroup_id, 'mysrc')
        srcgroup_id = open_rbf_src_hdf(mmgroup_id, 'mysrc')
        if (present(group_names)) then
            call set_rbf_src_patchname_hdf(srcgroup_id, group_names(1,bcface)%get())
        else
            call set_rbf_src_patchname_hdf(srcgroup_id, 'Default')
        end if

        ! Add MM driver to the source
        call create_rbf_mm_driver_hdf(srcgroup_id, 'PMM')
        driver_id = open_rbf_mm_driver_hdf(srcgroup_id)

        call create_prescribed_mesh_motion_function(pmmf, 'static')
        call add_pmmf_hdf(driver_id, pmmf)
        call close_rbf_mm_driver_hdf(driver_id)
        call close_rbf_src_hdf(srcgroup_id)
        call close_mm_group_hdf(mmgroup_id)


        !Assign pmm to domain
        call set_mm_domain_group_hdf(dom_id,'static')

        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')

        ! Close file
        call close_domain_hdf(dom_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__rbf_pmm__multiscale__singleblock
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



























    




end module mod_gridgen_blocks_rbf_pmm
