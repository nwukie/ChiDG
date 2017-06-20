module mod_gridgen_smoothbump
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: PI, ZERO, ONE, TWO, THREE, HALF, &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_string,             only: string_t
    use mod_bc,                 only: create_bc
    use mod_plot3d_utilities,   only: get_block_points_plot3d, get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d
    use mod_hdf_utilities,      only: add_domain_hdf, initialize_file_hdf, open_domain_hdf,         &
                                      close_domain_hdf, add_bc_state_hdf, close_file_hdf,           &
                                      close_hdf, set_patch_hdf, set_contains_grid_hdf,              &
                                      open_bc_group_hdf, close_bc_group_hdf, create_bc_state_group_hdf,   &
                                      set_patch_group_hdf, open_file_hdf, create_patch_hdf,         &
                                      open_patch_hdf, close_patch_hdf
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
    subroutine create_mesh_file__smoothbump(filename,nelem_xi,nelem_eta,nelem_zeta,equation_sets,group_names,bc_state_groups)
        character(*),           intent(in)              :: filename
        integer(ik),            intent(in)              :: nelem_xi
        integer(ik),            intent(in)              :: nelem_eta
        integer(ik),            intent(in)              :: nelem_zeta
        type(string_t),         intent(in), optional    :: equation_sets(:)
        type(string_t),         intent(in), optional    :: group_names(:,:)
        type(bc_state_group_t), intent(in), optional    :: bc_state_groups(:)

        class(bc_state_t),  allocatable             :: bc_state
        integer(HID_T)                              :: file_id, dom_id, bcface_id, bcgroup_id, patch_id
        integer(ik)                                 :: ierr, mapping, bcface, igroup, istate
        character(8)                                :: patch_names(6)
        real(rk),           allocatable             :: nodes(:,:)
        integer(ik),        allocatable             :: elements(:,:), faces(:,:)
        class(bc_state_t),  allocatable             :: inlet, outlet, wall
        real(rk),           allocatable, dimension(:,:,:)   :: xcoords, ycoords, zcoords



        !
        ! Generate coordinates
        !
        call meshgen_smoothbump_quartic(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords)

    

        ! Create/initialize file
        call initialize_file_hdf(filename)
        file_id = open_file_hdf(filename)


        !
        ! Get nodes/elements
        !
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping=4,idomain=1)



        !
        ! Add domains
        !
        if (present(equation_sets)) then
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian',equation_sets(1)%get())
        else
            call add_domain_hdf(file_id,'01',nodes,elements,'Cartesian','Euler')
        end if



        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,'01')

        do bcface = 1,6

            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping=4,bcface=bcface)

            ! Create/Set patch face indices
            patch_id = create_patch_hdf(dom_id,bcface)
            call set_patch_hdf(patch_id,faces)
            call close_patch_hdf(patch_id)

        end do !bcface







        !
        ! Add bc_group's
        !
        if (present(bc_state_groups)) then
            do igroup = 1,size(bc_state_groups)
                call create_bc_state_group_hdf(file_id,bc_state_groups(igroup)%name)

                bcgroup_id = open_bc_group_hdf(file_id,bc_state_groups(igroup)%name)

                do istate = 1,bc_state_groups(igroup)%nbc_states()
                    call add_bc_state_hdf(bcgroup_id, bc_state_groups(igroup)%bc_state(istate)%state)
                end do
                call close_bc_group_hdf(bcgroup_id)
            end do
        else

            ! Create Inlet boundary condition group
            call create_bc_state_group_hdf(file_id,'Inlet')
            bcgroup_id = open_bc_group_hdf(file_id,'Inlet')

            call create_bc('Inlet - Total', bc_state)
            call bc_state%set_fcn_option('Total Pressure',   'val',110000._rk)
            call bc_state%set_fcn_option('Total Temperature','val',300._rk   )

            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)


            ! Create Outlet boundary condition group
            call create_bc_state_group_hdf(file_id,'Outlet')
            bcgroup_id = open_bc_group_hdf(file_id,'Outlet')

            call create_bc('Outlet - Constant Pressure', bc_state)
            call bc_state%set_fcn_option('Static Pressure','val',100000._rk)

            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)


            ! Create Walls boundary condition group
            call create_bc_state_group_hdf(file_id,'Walls')
            bcgroup_id = open_bc_group_hdf(file_id,'Walls')

            call create_bc('Wall', bc_state)

            call add_bc_state_hdf(bcgroup_id,bc_state)
            call close_bc_group_hdf(bcgroup_id)



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




        call close_domain_hdf(dom_id)


        ! Set 'Contains Grid', close file
        call set_contains_grid_hdf(file_id,'True')
        call close_file_hdf(file_id)
        call close_hdf()




    end subroutine create_mesh_file__smoothbump
    !******************************************************************************



















    !>  Generate a single-block grid representing a smooth bump in a channel.
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
    !!
    !---------------------------------------------------------------------------------
    subroutine meshgen_smoothbump_quartic(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords)
        integer(ik)             :: nelem_xi
        integer(ik)             :: nelem_eta
        integer(ik)             :: nelem_zeta
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik)             :: npt_x, npt_y, npt_z, &
                                   ierr, i, j, k, n, &
                                   ipt_x, ipt_y, ipt_z
        real(rk)                :: x,y,z,alpha


        !
        ! Compute number of points in each direction, quartic elements
        !
        npt_x = nelem_xi*4   + 1
        npt_y = nelem_eta*4  + 1
        npt_z = nelem_zeta*4 + 1 


        !
        ! Allocate coordinate array
        !
        allocate(xcoords(npt_x,npt_y,npt_z), &
                 ycoords(npt_x,npt_y,npt_z), &
                 zcoords(npt_x,npt_y,npt_z), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Generate coordinates
        !
        do ipt_z = 1,npt_z
            do ipt_y = 1,npt_y
                do ipt_x = 1,npt_x

                    x = -1.5_rk + real(ipt_x-1,kind=rk)*(THREE / real(npt_x-1,kind=rk))
                    alpha = 0.0625_rk * exp(-25._rk * (x**TWO))
                    y = alpha + real(ipt_y-1,kind=rk)*(0.8_rk - alpha)/real(npt_y-1,kind=rk)
                    z = ZERO + real(ipt_z-1,kind=rk)*(ONE / real(npt_z-1,kind=rk))

                    xcoords(ipt_x,ipt_y,ipt_z) = x
                    ycoords(ipt_x,ipt_y,ipt_z) = y
                    zcoords(ipt_x,ipt_y,ipt_z) = z

                end do
            end do
        end do





    end subroutine meshgen_smoothbump_quartic
    !**********************************************************************************






end module mod_gridgen_smoothbump
