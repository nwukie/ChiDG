module mod_hdfio
    use mod_kinds,      only: rk,ik
    use type_point,     only: point_t
    use type_domain,    only: domain_t
    use hdf5
    use h5lt
    implicit none



contains


    subroutine read_grid_hdf5(filename, domains)
        character(*),                 intent(in)    :: filename
        type(domain_t), allocatable, intent(inout) :: domains(:)

        integer(HID_T)   :: fid, gid, sid, did_x, did_y, did_z      !> Identifiers
        integer(HSIZE_T) :: dims(3), maxdims(3)                     !> Dataspace dimensions

        type(c_ptr)                             :: pts
        type(point_t), allocatable              :: points(:,:,:)
        real(rk), dimension(:,:,:), allocatable :: xpts, ypts, zpts
        character(10)                           :: gname
        integer                                 :: nmembers, type, ierr, ndomains, igrp,    &
                                                   npts, izeta, ieta, ixi, idom, nterms_1d, &
                                                   nterms_c, mapping
        integer, dimension(1)                   :: buf
        logical                                 :: file_exists

        !>  Check file exists
        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            print*, "Found input file: ", filename
        else
            print*, "Error: read_grid_hdf5 - file not found: ", filename
            stop
        end if


        !>  Initialize Fortran interface.
        call h5open_f(ierr)
        if (ierr /= 0) stop "Error: read_grid_hdf5 - h5open_f"



        !>  Open input file using default properties.
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) stop "Error: read_grid_hdf5 - h5fopen_f"



        !>  Get number of domains from attribute 'ndomains' in file root
        call h5ltget_attribute_int_f(fid, "/", 'ndomains', buf, ierr)
        ndomains = buf(1)
        if (ierr /= 0) stop "Error: read_grid_hdf5 - h5ltget_attribute_int_f"



        !>  Allocate number of domains
        if (ndomains == 0) then
            stop "Error: read_hdf5 - no Domains found in file."
        else
            allocate(domains(ndomains), stat=ierr)
            if (ierr /= 0) stop "Error - read_grid_hdf5: domain allocation"
        end if


        !>  Get number of groups in the file root
        call h5gn_members_f(fid, "/", nmembers, ierr)

        !>  Loop through groups and read domains
        idom = 1
        do igrp = 0,nmembers-1
            call h5gget_obj_info_idx_f(fid,"/", igrp, gname, type, ierr)

            if (gname(1:2) == 'D_') then

                !> Open the Domain/Grid group
                call h5gopen_f(fid, trim(gname)//"/Grid", gid, H5P_DEFAULT_F)


                !>  Get number of terms
                call h5ltget_attribute_int_f(fid, "/", 'mapping', buf, ierr)
                mapping = buf(1)
                if (ierr /= 0) stop "Error: read_grid_hdf5 - h5ltget_attribute_int_f"
                nterms_1d = (mapping + 1)
                nterms_c = nterms_1d*nterms_1d*nterms_1d


                !>  Open the Coordinate datasets
                call h5dopen_f(gid, "CoordinateX", did_x, ierr)
                call h5dopen_f(gid, "CoordinateY", did_y, ierr, H5P_DEFAULT_F)
                call h5dopen_f(gid, "CoordinateZ", did_z, ierr, H5P_DEFAULT_F)


                !>  Get the dataspace id and dimensions
                call h5dget_space_f(did_x, sid, ierr)
                call h5sget_simple_extent_dims_f(sid, dims, maxdims, ierr)

                !>  Read x-points
                allocate(xpts(dims(1),dims(2),dims(3)))
                call h5dread_f(did_x, H5T_NATIVE_DOUBLE, xpts, dims, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"

                allocate(ypts, source=xpts)
                call h5dread_f(did_y, H5T_NATIVE_DOUBLE, ypts, dims, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"

                allocate(zpts, source=xpts)
                call h5dread_f(did_z, H5T_NATIVE_DOUBLE, zpts, dims, ierr)
                if (ierr /= 0) stop "Error: read_grid_hdf5 -- h5dread_f"


                !>  Accumulate pts into a single points_t matrix to initialize domain
                npts = dims(1)*dims(2)*dims(3)
                allocate(points(dims(1),dims(2),dims(3)))
                do izeta = 1,dims(3)
                    do ieta = 1,dims(2)
                        do ixi = 1,dims(1)
                            call points(ixi,ieta,izeta)%set(xpts(ixi,ieta,izeta),ypts(ixi,ieta,izeta),zpts(ixi,ieta,izeta))
                        end do
                    end do
                end do


                !> Call domain geometry initialization
                call domains(idom)%init_geom(nterms_c,points)



                !> Close the Coordinate datasets
                call h5dclose_f(did_x,ierr)
                call h5dclose_f(did_y,ierr)
                call h5dclose_f(did_z,ierr)


                !> Close the Domain/Grid group
                call h5gclose_f(gid,ierr)

                !> Deallocate points for the current domain
                deallocate(zpts,ypts,xpts,points)
                idom = idom + 1
            end if
        end do

        !>  Close file
        call h5fclose_f(fid, ierr)

        !>  Close Fortran interface
        call h5close_f(ierr)


    end subroutine







end module mod_hdfio
