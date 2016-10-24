module mod_gridgen_smoothbump
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: PI, ZERO, ONE, TWO, THREE, HALF, &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use mod_string,             only: string_t
    use mod_bc,                 only: create_bc
    use type_bc_state_wrapper,  only: bc_state_wrapper_t
    use mod_plot3d_utilities,   only: get_block_points_plot3d, get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d
    use mod_hdf_utilities,      only: add_domain_hdf, initialize_file_hdf, open_domain_hdf, &
                                      close_domain_hdf, add_bc_state_hdf, close_file_hdf, &
                                      close_hdf, set_bc_patch_hdf, set_contains_grid_hdf
    use hdf5

    use type_point,             only: point_t
    use type_bc_state,          only: bc_state_t
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
    subroutine create_mesh_file__smoothbump(filename,nelem_xi,nelem_eta,nelem_zeta,equation_sets,bc_states)
        character(*),               intent(in)              :: filename
        integer(ik),                intent(in)              :: nelem_xi
        integer(ik),                intent(in)              :: nelem_eta
        integer(ik),                intent(in)              :: nelem_zeta
        type(string_t),             intent(in), optional    :: equation_sets(:)
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states(:,:)

        type(bc_state_wrapper_t)                    :: bc_states_set(6)
        integer(HID_T)                              :: file_id, dom_id, bcface_id
        integer(ik)                                 :: ierr, spacedim, mapping, bcface
        character(8)                                :: face_strings(6)
        type(point_t),      allocatable             :: nodes(:)
        integer(ik),        allocatable             :: elements(:,:), faces(:,:)
        class(bc_state_t),  allocatable             :: inlet, outlet, wall
        real(rk),           allocatable, dimension(:,:,:)   :: xcoords, ycoords, zcoords



        !
        ! Generate coordinates
        !
        call meshgen_smoothbump_quartic(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords)

    

        ! Create/initialize file
        file_id = initialize_file_hdf(filename)


        !
        ! Get nodes/elements
        !
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping=4,idomain=1)



        !
        ! Add domains
        !
        spacedim = 3
        if (present(equation_sets)) then
            call add_domain_hdf(file_id,"01",nodes,elements,equation_sets(1)%get(),spacedim)
        else
            call add_domain_hdf(file_id,"01",nodes,elements,"Euler",spacedim)
        end if



        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,"01")

        do bcface = 1,6

            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping=4,bcface=bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom_id,faces,bcface)

        end do !bcface




        if (present(bc_states)) then

            bc_states_set(XI_MIN)   = bc_states(1,XI_MIN)
            bc_states_set(XI_MAX)   = bc_states(1,XI_MAX)
            bc_states_set(ETA_MIN)  = bc_states(1,ETA_MIN)
            bc_states_set(ETA_MAX)  = bc_states(1,ETA_MAX)
            bc_states_set(ZETA_MIN) = bc_states(1,ZETA_MIN)
            bc_states_set(ZETA_MAX) = bc_states(1,ZETA_MAX)

        else

            ! Create boundary conditions
            call create_bc("Total Inlet",     bc_states_set(XI_MIN)%state)
            call create_bc("Pressure Outlet", bc_states_set(XI_MAX)%state)
            call create_bc("Wall", bc_states_set(ETA_MIN)%state )
            call create_bc("Wall", bc_states_set(ETA_MAX)%state )
            call create_bc("Wall", bc_states_set(ZETA_MIN)%state)
            call create_bc("Wall", bc_states_set(ZETA_MAX)%state)


            ! Set Inlet bc parameters
            call bc_states_set(XI_MIN)%state%set_fcn_option("TotalPressure","val",110000._rk)
            call bc_states_set(XI_MIN)%state%set_fcn_option("TotalTemperature","val",300._rk)

            ! Set Outlet bc parameter
            call bc_states_set(XI_MAX)%state%set_fcn_option("Static Pressure","val",100000._rk)

        end if



        !
        ! Set all boundary conditions to walls, inlet, outlet...
        !
        !
        face_strings = ["XI_MIN  ","XI_MAX  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(face_strings)

            call h5gopen_f(dom_id,"BoundaryConditions/"//trim(adjustl(face_strings(bcface))),bcface_id,ierr)

            call add_bc_state_hdf(bcface_id,bc_states_set(bcface)%state)

            call h5gclose_f(bcface_id,ierr)

        end do


        call close_domain_hdf(dom_id)


        ! Set 'Contains Grid', close file
        call set_contains_grid_hdf(file_id,"True")
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
        if (ierr /= 0) stop "Allocation Error"



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
