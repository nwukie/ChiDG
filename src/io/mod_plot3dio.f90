module mod_plot3dio
#include <messenger.h>
    use mod_kinds,                  only: rk,ik,rdouble
    use mod_constants,              only: ZERO, NFACES, TWO_DIM, THREE_DIM, NO_PROC, &
                                          CARTESIAN, CYLINDRICAL
    use type_chidg_data,            only: chidg_data_t
    use type_mr4vector,             only: mr4vector_t
    use type_nvector,               only: nvector_t
    implicit none


contains


  
   
    !>  Write plot3d files for mesh sensitivities, these are the files written out:
    !!  
    !!      - grid file (.x) containing the nodes' coordinates
    !!      - function file (.q) containing the sensitivities
    !!      - name file (.nam) containg the name of the variables
    !!
    !!  The name file is meant to be used for visualization purposes for example in TecPlot.
    !!  If read in together with the grid and function file, it overwrites the variables
    !!  names and Tecplot will show the variable names like they are listed in the name file.
    !!  
    !!  @author Matteo Ugolotti
    !!  @date   8/26/2018
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine write_x_sensitivities_plot3d(data,ifunc,filename)
        type(chidg_data_t),      intent(inout)  :: data
        integer(ik),             intent(in)     :: ifunc
        character(*),            intent(in)     :: filename

        character(:),   allocatable     :: filename_grd, filename_fun, filename_nam
        integer(ik)                     :: nblocks, iblock, ipt_x, ipt_y, ipt_z, inode, &
                                           i, j, k, idir, ivar, ierr, n, nnodes_grid,   &
                                           nnodes_stored
        integer(ik),    allocatable     :: nvars(:,:), npoints(:,:)
        logical                         :: node_count_match
        real(rk)                        :: r, theta, z, dJdr, dJdtheta, dJdz
        real(rk),       allocatable     :: blk_coords(:,:,:,:), blk_sensitivities(:,:,:,:)
        type(mr4vector_t)               :: coords, sensitivities
        type(nvector_t),allocatable     :: node_sensitivities(:)


        ! File names
        filename_grd = trim(filename) // ".x"
        filename_fun = trim(filename) // ".q"
        filename_nam = trim(filename) // ".nam"

         
        ! Get overall number of blocks
        npoints = data%mesh%npoints
        nblocks = size(npoints,1)
        

        ! Define matrix for .q file variables
        allocate(nvars(nblocks,4),stat=ierr)
        if (ierr/=0) call AllocationError
        nvars(:,1:3) = npoints
        nvars(:,4)  = 3 ! number of variables

         
        ! Get sensitivities for the ifunc
        node_sensitivities = data%sdata%adjointx%node_sensitivities(ifunc,:)


        ! Retrieve data and set it in plot3d friendly-format
        do iblock = 1,nblocks


            ! For sanity check, check that the number of points correspond to the original number of points in the grid
            nnodes_stored = node_sensitivities(iblock)%size()
            nnodes_grid   = npoints(iblock,1)*npoints(iblock,2)*npoints(iblock,3)
            node_count_match = ( nnodes_stored == nnodes_grid )
            if (.not. node_count_match) then
                call chidg_signal_two(FATAL,"write_x_sensitivities: mismatch of node count in block ",iblock,".")
            end if
           

            ! Re/allocate matrices based on block dimensions
            if (allocated(blk_coords))        deallocate(blk_coords)
            if (allocated(blk_sensitivities)) deallocate(blk_sensitivities)
            allocate(blk_coords(npoints(iblock,1),npoints(iblock,2),npoints(iblock,3),3), stat=ierr)
            if (ierr/=0) call AllocationError
            allocate(blk_sensitivities(npoints(iblock,1),npoints(iblock,2),npoints(iblock,3),3), stat=ierr)
            if (ierr/=0) call AllocationError

            
            ! Reset node count
            inode = 1


            do ipt_z = 1,npoints(iblock,3)
                do ipt_y = 1,npoints(iblock,2)
                    do ipt_x = 1,npoints(iblock,1)
                        
                        ! Simply get sensitivities and coords since they are already in cartesian coordinates
                        if ( node_sensitivities(iblock)%data(inode)%coordinate_system == CARTESIAN ) then
                            blk_coords(ipt_x,ipt_y,ipt_z,1)        = node_sensitivities(iblock)%data(inode)%coords(1)
                            blk_coords(ipt_x,ipt_y,ipt_z,2)        = node_sensitivities(iblock)%data(inode)%coords(2)
                            blk_coords(ipt_x,ipt_y,ipt_z,3)        = node_sensitivities(iblock)%data(inode)%coords(3)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,1) = node_sensitivities(iblock)%data(inode)%sensitivities(1)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,2) = node_sensitivities(iblock)%data(inode)%sensitivities(2)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,3) = node_sensitivities(iblock)%data(inode)%sensitivities(3)
                        end if

                        ! Tranform cylindrical coordinates back to cartesian and rotate sensitivity vector
                        if ( node_sensitivities(iblock)%data(inode)%coordinate_system == CYLINDRICAL ) then
                            r     = node_sensitivities(iblock)%data(inode)%coords(1)
                            theta = node_sensitivities(iblock)%data(inode)%coords(2)
                            z     = node_sensitivities(iblock)%data(inode)%coords(3)
                            blk_coords(ipt_x,ipt_y,ipt_z,1) = r*cos(theta) 
                            blk_coords(ipt_x,ipt_y,ipt_z,2) = r*sin(theta)
                            blk_coords(ipt_x,ipt_y,ipt_z,3) = z
                            dJdr     = node_sensitivities(iblock)%data(inode)%sensitivities(1)
                            dJdtheta = node_sensitivities(iblock)%data(inode)%sensitivities(2)
                            dJdz     = node_sensitivities(iblock)%data(inode)%sensitivities(3)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,1) = dJdr*cos(theta) - dJdtheta*sin(theta)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,2) = dJdr*sin(theta) + dJdtheta*cos(theta)
                            blk_sensitivities(ipt_x,ipt_y,ipt_z,3) = dJdz
                        end if
                        
                        inode = inode + 1

                    end do !ipt_x

                end do !ipt_y

            end do !ipt_z
            
            ! Push-back matrix for iblock
            call coords%push_back(blk_coords)
            call sensitivities%push_back(blk_sensitivities)


        end do !iblock


        ! Write grid file .x
        open(unit=7, file=trim(filename_grd), form='unformatted')
        write(7) nblocks
        write(7) (( npoints(iblock,idir), idir=1,3), iblock=1,nblocks)
        do iblock = 1,nblocks
            write(7) (((( coords%data(iblock)%mat(i,j,k,n), i=1,npoints(iblock,1)), j=1,npoints(iblock,2)), k=1,npoints(iblock,3)), n=1,3)
        end do
        close(7)


        ! Write function file .q
        open(unit=8, file=trim(filename_fun), form='unformatted')
        write(8) nblocks
        write(8) (( nvars(iblock,ivar), ivar=1,4), iblock=1,nblocks)
        do iblock = 1,nblocks
            write(8) (((( sensitivities%data(iblock)%mat(i,j,k,n), i=1,npoints(iblock,1)), j=1,npoints(iblock,2)), k=1,npoints(iblock,3)), n=1,3)
        end do
        close(8)


        ! Write name file .nam
        open(unit=9, file=trim(filename_nam) )
        write(9,'(A)') "dJ/dx" 
        write(9,'(A)') "dJ/dy" 
        write(9,'(A)') "dJ/dz" 
        close(9)

    end subroutine write_x_sensitivities_plot3d
    !****************************************************************************************
  
  


end module mod_plot3dio
