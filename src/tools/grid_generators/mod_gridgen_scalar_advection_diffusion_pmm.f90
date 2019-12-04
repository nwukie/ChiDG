module mod_gridgen_scalar_advection_diffusion_pmm
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: PI, ZERO, ONE, TWO, THREE, FOUR, HALF, &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_string,             only: string_t
    use mod_bc,                 only: create_bc
    use mod_plot3d_utilities,   only: get_block_points_plot3d, get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d
    use mod_hdf_utilities
    use hdf5

    use type_bc_state,          only: bc_state_t
    use type_bc_state_group,    only: bc_state_group_t
    implicit none








contains





    !>  Generate smooth bump grid file + boundary conditions for Euler equations.
    !!
    !!
    !!  x = -1.5              x = 1.5
    !!     |                     |
    !!     v                     v
    !!
    !!     .---------------------.      y = 0.8
    !!     |                     |
    !!     |         ---         |
    !!     .--------/   \--------.      y = 0.0625 e^(-25*x*x)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/19/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    subroutine create_mesh_file__scalar_advection_diffusion_pmm(filename,nelem_xi,nelem_eta,nelem_zeta,equation_sets,group_names,bc_state_groups)
        character(*),           intent(in)              :: filename
        integer(ik),            intent(in)              :: nelem_xi
        integer(ik),            intent(in)              :: nelem_eta
        integer(ik),            intent(in)              :: nelem_zeta
        type(string_t),         intent(in), optional    :: equation_sets(:)
        type(string_t),         intent(in), optional    :: group_names(:,:)
        type(bc_state_group_t), intent(in), optional    :: bc_state_groups(:)

        class(bc_state_t),  allocatable             :: bc_state
        integer(HID_T)                              :: file_id, dom_id, bcface_id, bcgroup_id, patch_id, mmgroup_id
        integer(ik)                                 :: ierr, order, bcface, igroup, istate
        character(8)                                :: patch_names(6)
        real(rk),           allocatable             :: nodes(:,:)
        integer(ik),        allocatable             :: elements(:,:), faces(:,:)
        class(bc_state_t),  allocatable             :: inlet, outlet, wall
        real(rk),           allocatable, dimension(:,:,:)   :: xcoords, ycoords, zcoords

        character(len=8)                            :: bc_face_strings(6)
        character(:),   allocatable                 :: bc_face_string
        integer(ik)                                 :: npoints(3)
        class(prescribed_mesh_motion_function_t), allocatable   :: pmmf


        !
        ! Generate coordinates
        !
        call meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords)

    

        ! Create/initialize file
        call initialize_file_hdf(filename)
        file_id = open_file_hdf(filename)


        !
        ! Get nodes/elements
        !
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,order=1,idomain=1)


        !
        ! Get npoints in each direction for each block
        !
        npoints(1) = size(xcoords,1)
        npoints(2) = size(xcoords,2)
        npoints(3) = size(xcoords,3)


        !
        ! Add domains
        !
        if (present(equation_sets)) then
            call add_domain_hdf(file_id,'01',npoints,nodes,elements,'Cartesian',equation_sets(1)%get())
        else
            call add_domain_hdf(file_id,'01',npoints,nodes,elements,'Cartesian','Scalar Advection Diffusion ALE')
        end if



        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,'01')

        bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]
        do bcface = 1,6

            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping=1,bcface=bcface)

            ! Create/Set patch face indices
            bc_face_string  = trim(bc_face_strings(bcface))
            patch_id = create_patch_hdf(dom_id,bc_face_string)
            call set_patch_hdf(patch_id,faces)
            call close_patch_hdf(patch_id)

        end do !bcface







        !
        ! Add bc_group's
        !
        ! BCs should be defined in the integration test
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

!            ! Create Inlet boundary condition group
!            call create_bc_state_group_hdf(file_id,'Inlet')
!            bcgroup_id = open_bc_group_hdf(file_id,'Inlet')
!
!            call create_bc('Inlet - Total', bc_state)
!            call bc_state%set_fcn_option('Total Pressure',   'val',110000._rk)
!            call bc_state%set_fcn_option('Total Temperature','val',300._rk   )
!
!            call add_bc_state_hdf(bcgroup_id,bc_state)
!            call close_bc_group_hdf(bcgroup_id)
!
!
!            ! Create Outlet boundary condition group
!            call create_bc_state_group_hdf(file_id,'Outlet')
!            bcgroup_id = open_bc_group_hdf(file_id,'Outlet')
!
!            call create_bc('Outlet - Constant Pressure', bc_state)
!            call bc_state%set_fcn_option('Static Pressure','val',100000._rk)
!
!            call add_bc_state_hdf(bcgroup_id,bc_state)
!            call close_bc_group_hdf(bcgroup_id)
!
!
!            ! Create Walls boundary condition group
!            call create_bc_state_group_hdf(file_id,'Walls')
!            bcgroup_id = open_bc_group_hdf(file_id,'Walls')
!
!            call create_bc('Wall', bc_state)
!
!            call add_bc_state_hdf(bcgroup_id,bc_state)
!            call close_bc_group_hdf(bcgroup_id)
!


        end if






        !
        ! Set boundary condition groups for each patch
        !
        patch_names = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(patch_names)

            !call h5gopen_f(dom_id,'BoundaryConditions/'//trim(adjustl(patch_names(bcface))),patch_id,ierr)
            patch_id = open_patch_hdf(dom_id,trim(patch_names(bcface)))


            ! Set bc_group
            if (present(group_names)) then
                call set_patch_group_hdf(patch_id,group_names(1,bcface)%get())
            else
                if (trim(adjustl(patch_names(bcface))) == 'XI_MIN') then
                    call set_patch_group_hdf(patch_id,'Inlet')
                else if (trim(adjustl(patch_names(bcface))) == 'XI_MAX') then
                    call set_patch_group_hdf(patch_id,'Outlet')
                else
                    call set_patch_group_hdf(patch_id,'Walls')
                end if
            end if


            call h5gclose_f(patch_id,ierr)

        end do



        !
        ! Define PMM
        !
        call create_mm_group_hdf(file_id,'sin_uf','PMM')
        mmgroup_id = open_mm_group_hdf(file_id, 'sin_uf')
        call create_prescribed_mesh_motion_function(pmmf, 'sinusoidal_uniform_flow')
        call add_pmmf_hdf(mmgroup_id, pmmf)
        call close_mm_group_hdf(mmgroup_id)

        !Assign pmm to domain
        call set_mm_domain_group_hdf(dom_id,'sin_uf')


        call close_domain_hdf(dom_id)


        ! Set 'Contains Grid', close file
        call set_contains_grid_hdf(file_id,'True')
        call close_file_hdf(file_id)
        call close_hdf()




    end subroutine create_mesh_file__scalar_advection_diffusion_pmm
    !******************************************************************************

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine create_mesh_file__scalar_advection_diffusion_pmm__multiblock(filename,equation_sets,group_names,bc_state_groups,     &
                                                          nelem_xi,  nelem_eta,  nelem_zeta,  &
                                                          xmax,ymax,zmax,clusterx)
        character(*),           intent(in)              :: filename
        type(string_t),         intent(in), optional    :: equation_sets(:)
        type(string_t),         intent(in), optional    :: group_names(:,:)
        type(bc_state_group_t), intent(in), optional    :: bc_state_groups(:)
        integer(ik),            intent(in)              :: nelem_xi,  nelem_eta,  nelem_zeta
        real(rk),               intent(in), optional    :: xmax, ymax, zmax
        integer(ik),            intent(in), optional    :: clusterx

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(6)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface1_id, bcface2_id, &
                                                           bcgroup_id, patch1_id, patch2_id, mmgroup_id
        integer(ik)                                     :: order, bcface, ierr, igroup, istate, &
                                                           nxi_max, neta_max, nzeta_max,xi_mid
        real(rk),       allocatable                     :: nodes1(:,:), nodes2(:,:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, ycoords1, zcoords1, &
                                                           xcoords2, ycoords2, zcoords2
        real(rk)                                        :: xmax_block1,xmin_block2, xmax_current,       &
                                                           xmax_new, ymax_new, zmax_new, ymax_current,  &
                                                           zmax_current
        character(len=8)                                :: bc_face_strings(6)
        character(:),   allocatable                     :: bc_face_string
        integer(ik)                                     :: npoints_1(3), npoints_2(3)
        class(prescribed_mesh_motion_function_t), allocatable   :: pmmf

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
        order = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords1,zcoords1)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords2,zcoords2)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords1,zcoords1,order,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords2,zcoords2,order,idomain=2)


        !
        ! Get npoints in each direction for each block
        !
        npoints_1(1) = size(xcoords1,1)
        npoints_1(2) = size(xcoords1,2)
        npoints_1(3) = size(xcoords1,3)
        npoints_2(1) = size(xcoords2,1)
        npoints_2(2) = size(xcoords2,2)
        npoints_2(3) = size(xcoords2,3)


        !
        ! Add domains
        !
        if ( present(equation_sets) ) then
            call add_domain_hdf(file_id,'01',npoints_1,nodes1,elements1,'Cartesian',equation_sets(1)%get())
            call add_domain_hdf(file_id,'02',npoints_2,nodes2,elements2,'Cartesian',equation_sets(2)%get())
        else
            call add_domain_hdf(file_id,'01',npoints_1,nodes1,elements1,'Cartesian','Scalar Advection Diffusion ALE')
            call add_domain_hdf(file_id,'02',npoints_2,nodes2,elements2,'Cartesian','Scalar Advection Diffusion ALE')
        end if






        !
        ! Set boundary conditions patch connectivities
        !
        dom1_id = open_domain_hdf(file_id,'01')
        dom2_id = open_domain_hdf(file_id,'02')

        bc_face_strings = ["XI_MIN  ","XI_MAX  ","ETA_MIN ","ETA_MAX ","ZETA_MIN","ZETA_MAX"]
        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords1,zcoords1,order,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords2,zcoords2,order,bcface)


            ! Set bc patch face indices
            bc_face_string  = trim(bc_face_strings(bcface))
            patch1_id = create_patch_hdf(dom1_id,bc_face_string)
            patch2_id = create_patch_hdf(dom2_id,bc_face_string)

            call set_patch_hdf(patch1_id,faces1)
            call set_patch_hdf(patch2_id,faces2)

            call close_patch_hdf(patch1_id)
            call close_patch_hdf(patch2_id)

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
        ! Assign groups to boundary condition patches
        !
        faces = ['XI_MIN  ','XI_MAX  ', 'ETA_MIN ', 'ETA_MAX ', 'ZETA_MIN', 'ZETA_MAX']
        do bcface = 1,size(faces)
            patch1_id = open_patch_hdf(dom1_id,trim(adjustl(faces(bcface))))
            patch2_id = open_patch_hdf(dom2_id,trim(adjustl(faces(bcface))))

            if (present(group_names)) then
                call set_patch_group_hdf(patch1_id,group_names(1,bcface)%get())
                call set_patch_group_hdf(patch2_id,group_names(2,bcface)%get())
            else
                call set_patch_group_hdf(patch1_id,'Default')
                call set_patch_group_hdf(patch2_id,'Default')
            end if

            call close_patch_hdf(patch1_id)
            call close_patch_hdf(patch2_id)
        end do



        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,'True')
!
        ! Define PMM
        !
        call create_mm_group_hdf(file_id,'sin_uf','PMM')
        mmgroup_id = open_mm_group_hdf(file_id, 'sin_uf')
        call create_prescribed_mesh_motion_function(pmmf, 'sinusoidal_uniform_flow')
        call add_pmmf_hdf(mmgroup_id, pmmf)
        call close_mm_group_hdf(mmgroup_id)

        !Assign pmm to domain
        call set_mm_domain_group_hdf(dom1_id,'sin_uf')
        call set_mm_domain_group_hdf(dom2_id,'sin_uf')



        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__scalar_advection_diffusion_pmm__multiblock
    !*************************************************************************************






    subroutine meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords,clusterx,x_max_in,x_min_in)
        integer(ik),    intent(in)                  :: nelem_xi
        integer(ik),    intent(in)                  :: nelem_eta
        integer(ik),    intent(in)                  :: nelem_zeta
        real(rk),       intent(inout),  allocatable :: xcoords(:,:,:)
        real(rk),       intent(inout),  allocatable :: ycoords(:,:,:)
        real(rk),       intent(inout),  allocatable :: zcoords(:,:,:)
        integer(ik),    intent(in),     optional    :: clusterx
        real(rk),       intent(in),     optional    :: x_max_in
        real(rk),       intent(in),     optional    :: x_min_in

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, &
                       npts_xi, npts_eta, npts_zeta
        real(rk)    :: x,y,z, x_max, x_min


        npts_xi   = nelem_xi   + 1
        npts_eta  = nelem_eta  + 1
        npts_zeta = nelem_zeta + 1


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Set max/min X-Coordinate
        !
        if (present(x_max_in)) then
            x_max = x_max_in
        else
            x_max = ONE
        end if

        if (present(x_min_in)) then
            x_min = x_min_in
        else
            x_min = ZERO
        end if




        !
        ! Generate points
        !
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi

                    if (present(clusterx)) then
                        if ( clusterx == -1 ) then
                            !x = ONE - tanh( (PI/TWO)*(ONE - real(ipt_xi-1,rk)/real(npts_xi-1,rk) ) )/tanh(PI/TWO)
                            x = x_max - (x_max-x_min)*tanh( (PI/TWO)*(ONE - real(ipt_xi-1,rk)/real(npts_xi-1,rk) ) )/tanh(PI/TWO)
                        else if ( clusterx == 1 ) then
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: 'clusterx'=1 not yet implemented.")
                        else
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: Invalid value for 'clusterx'. -1,1.")
                        end if
                    else
                        !x = real(ipt_xi-1,rk)/real(npts_xi-1,rk)
                        x = x_min + (x_max - x_min)*real(ipt_xi-1,rk)/real(npts_xi-1,rk)
                    end if

                    if (ipt_xi == npts_xi) then
                        !x = ONE
                        x = x_max
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




end module mod_gridgen_scalar_advection_diffusion_pmm
