module mod_vtk_file_unstr

! Module containing subroutines and functions for writing hdf5 data to a xml vtk format


! For unstructured grids

    use mod_kinds,          only: rk,ik,rdouble
    use mod_constants,      only: OUTPUT_RES

    use mod_vtk_calc_func,  only: get_cons_var, get_piece_nums, get_piece_coord, get_piece_data, get_piece_connectivity_data

    use type_chidg_data,    only: chidg_data_t

    implicit none



contains



    !>  Open vtk file
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/26/2016
    !!
    !!
    !!
    !!
    !
    !-------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine open_vtk_unstr(file_name,byte_order)

        character(len = 100),intent(in)         ::  file_name
        integer(ik)         ,intent(in)         ::  byte_order  ! 0 - LittleEndian, 1 - BigEndian


        character(len = 100)                    ::  grid_type_name, byte_order_name
        integer(ik),parameter                   ::  funit = 10
        logical                                 ::  exist


        grid_type_name = 'UnstructuredGrid'


        select case(byte_order)
            case(0)
                byte_order_name = 'LittleEndian'
            case(1)
                byte_order_name = 'BigEndian'
        end select


        !
        ! Open new vtk file and write the part common to all xml vtk files
        !
        inquire(file=trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', action = 'write')
        else
            open(funit, file = trim(file_name), status = 'new', action = 'write')
        end if
            write(funit, '(A)') '<?xml version="1.0"?>'
            write(funit, '(5A)') '<VTKFile type="',trim(grid_type_name),'" version="0.1" byte_order="',trim(byte_order_name),'">'
            write(funit, '(3A)') '  <',trim(grid_type_name),'>'
            close(funit)


    end subroutine open_vtk_unstr
    !***************************************************************************************************************************************************



    !>  Close vtk file
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/27/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine close_vtk_unstr(file_name)

        character(len = 100),intent(in)         ::  file_name


        character(len = 100)                    ::  grid_type_name
        integer(ik)                             ::  funit = 10
        logical                                 ::  exist


        grid_type_name = 'UnstructuredGrid'


        !
        ! Open existing vtk file (show an error if file doesn't exist) and write out the closing statements common to all
        ! xml vtk files
        !
        inquire(file = trim(file_name), exist=exist)
        if(exist) then
            open(funit, file=trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit,'(3A)') '  </',trim(grid_type_name),'>'
            write(funit,'(A)') '</VTKFile>'
            close(funit)


    end subroutine close_vtk_unstr
    !***************************************************************************************************************************************************



    !>  Initialize a unstructrued grid piece (domain in the given geometry) in the vtk file
    !!  A block may contain more than one pieces
    !!  Currently the converter is set up such that each block contains only 1 piece
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/27/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine init_piece_unstr(file_name,num_pts,num_cells)
        
        character(len = 100),intent(in)         ::  file_name
        integer(ik),         intent(in)         ::  num_pts, num_cells


        integer(ik)                             ::  funit = 10
        logical                                 ::  exist  


        !
        ! Open existing vtk file (show error if file doesn't exist) and initialize a new unstructured grid piece 
        ! num_pts = Number of points in the given piece (/in a block) 
        ! num_cells = Number of cells in the given piece (/in a block)
        ! num_pts and num_cells determined in mod_vtk_calc_func.f90
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit,'(A,I10,A,I10,A)' ) '    <Piece NumberOfPoints="',num_pts,'" NumberOfCells="',num_cells,'">'
            close(funit)


    end subroutine init_piece_unstr
    !***************************************************************************************************************************************************



    !>  End an initialized unstructured grid piece in the vtk file 
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/27/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine end_piece_unstr(file_name)

        character(len = 100),intent(in)         ::  file_name


        integer(ik)                             ::  funit = 10
        logical                                 ::  exist


        !
        ! Open existing vtk file (show error if file doesn't exist) and close current piece
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit, '(A)') '    </Piece>'
            close(funit)


    end subroutine end_piece_unstr 
    !***************************************************************************************************************************************************



    !>  Write point coordinates of an unstructured grid piece 
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/27/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine write_piece_coord(file_name,num_pts,X_coord,Y_coord,Z_coord)
        
        character(len = 100), intent(in)        ::  file_name
        integer(ik),          intent(in)        ::  num_pts
        real(rk),dimension(:),intent(in)        ::  X_coord,Y_coord,Z_coord


        character(len = 100)                    :: data_type
        integer(ik)                             :: funit = 10, ipt
        logical                                 :: exist


        data_type = 'Float32'


        !
        ! Open existing vtk file (show error if file doesn't exist)
        ! Initialize the points field (<Points>) for the current piece and write out x,y,z coordinate values
        ! Size of the X,Y,Z arrays (along with the arrays) is calculated in mod_vtk_calc_func
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status= 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit, '(A)') '      <Points>'
            write(funit, '(3A)') '        <DataArray type="',trim(data_type),'" NumberOfComponents="3" format="ascii">'
            do ipt = 1, num_pts

                write(funit, '(A,3F20.16)') '        ',X_coord(ipt),Y_coord(ipt),Z_coord(ipt)

            end do
            write(funit, '(A)') '        </DataArray>'
            write(funit, '(A)') '      </Points>'
            close(funit)


    end subroutine write_piece_coord
    !***************************************************************************************************************************************************



    !>  Write point data for conservative variables for an unstructured grid piece
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/27/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine write_piece_data(file_name,noeq,num_pts,cons_var,cons_var_val)

        character(Len = 100),             intent(in)        ::  file_name
        integer(ik),                      intent(in)        ::  num_pts,noeq
        character(len = 100),dimension(:),intent(in)        ::  cons_var        ! Array containing conservative  variable names
        real(rk),dimension(:,:),          intent(in)        ::  cons_var_val    ! Array containing values of conservative variables


        character(len = 100)                                :: data_type, varstring
        integer(ik)                                         :: funit = 10, ieq, ival, ivar
        logical                                             :: exist


        !
        ! Get the variable string as required to open the point data field <PointData...>
        !
        varstring = trim(cons_var(1))
        do ivar = 1,size(cons_var) - 1

            varstring = trim(varstring)//" "//trim(cons_var(ivar + 1))

        end do


        !
        ! TODO: Classify variables as scalrs or vectore according to data read from code execution
        !


        data_type = 'Float32'


        !
        ! Open existing vtk file (show error if file doesn't exist)
        ! Initialize point data field (<PointData...>) for the current piece, currently only for scalar quantities
        ! Write out conservative variable values for current place
        ! Name of conservative variables (along with their values) obtained from mod_vtk_calc_func
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit, '(3A)') '      <PointData Scalars="',trim(varstring),'">'
            do ieq = 1,noeq

                write(funit, '(5A)') '        <DataArray type="',trim(data_type),'" Name="',trim(cons_var(ieq)),'" format="ascii">'
                do ival = 1,num_pts

                    write(funit,*) '        ',cons_var_val(ieq,ival)

                end do
                write(funit, '(A)') '        </DataArray>'

            end do
            write(funit, '(A)') '      </PointData>'
            close(funit)


    end subroutine write_piece_data
    !***************************************************************************************************************************************************



    !>  Initialize a cell field (<Cells>) associated with a given piece 
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/28/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine init_cell(file_name)

        character(len = 100),intent(in)         :: file_name


        integer(ik)                             :: funit = 10
        logical                                 :: exist

        !
        ! Open existing vtk file (show error if file doesn't exist)
        ! Initialize the cells field (<Cells>) associated with a given piece
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit,'(A)' ) '      <Cells>'
            close(funit)


    end subroutine init_cell
    !***************************************************************************************************************************************************



    !>  End an initialized  cells field (</Cells>) associated with a given piece 
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/28/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine end_cell(file_name)

        character(len = 100),intent(in)         :: file_name


        integer(ik)                             :: funit = 10
        logical                                 :: exist


        !
        ! Open an existing vtk file (show error if file doesn't exist)
        ! End existing cells field
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit,'(A)') '      </Cells>'
            close(funit)


    end subroutine end_cell
    !***************************************************************************************************************************************************



    !>  Write piece connectivity, offsets and types fields 
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/28/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine write_piece_connectivity_data(file_name,num_cells,connectivity,offsets,types)

        character(len = 100),                  intent(in)           ::  file_name
        integer(ik),                           intent(in)           ::  num_cells
        integer(ik),dimension(:,:),allocatable,intent(in)           ::  connectivity
        integer(ik),dimension(:),allocatable  ,intent(in)           ::  offsets,types

 
        character(len = 100)                                        ::  type_1, type_2
        integer(ik)                                                 ::  icell_row,icell_col
        integer(ik)                                                 ::  funit = 10
        logical                                                     ::  exist


        type_1 = 'Int32'; type_2 = 'UInt8'

        !
        ! Open existing vtk file (show error if file doesn't exist)
        ! Write connectivity, offsets and types data for a given piece to the vtk file
        !
        inquire(file = trim(file_name), exist=exist)
        if (exist) then
            open(funit, file = trim(file_name), status = 'old', position = 'append', action = 'write')
        else
            print *, 'Error: vtk file does not exist'
        end if
            write(funit,'(3A)' ) '        <DataArray type="',trim(type_1),'" Name="connectivity" format="ascii">'
            do icell_row = 1,num_cells
 
                write(funit,*) '        ',connectivity(icell_row,:)

            end do
            write(funit,'(A)') '        </DataArray>'
            write(funit,'(3A)' ) '        <DataArray type="',trim(type_1),'" Name="offsets" format="ascii">'
            do icell_row = 1,num_cells

                write(funit,*) '        ',offsets(icell_row)

            end do
            write(funit,'(A)') '        </DataArray>'
            write(funit,'(3A)' ) '        <DataArray type="',trim(type_2),'" Name="types" format="ascii">'
            do icell_row = 1,num_cells

                write(funit,*) '        ',types(icell_row)

            end do
            write(funit,'(A)') '        </DataArray>'
            close(funit)


    end subroutine write_piece_connectivity_data
    !***************************************************************************************************************************************************



    !>  Write final multi block .pvd file for the entire simulation
    !!  The .pvd file calls upon separate blockwise .vtu files
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/28/2016
    !!
    !!
    !!
    !!
    !
    !----------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine write_pvd_final(data,pvd_file,file_arr,timestep)

        type(chidg_data_t),               intent(inout)         ::  data
        character(len = *),               intent(in)            ::  pvd_file
        character(len = 100),dimension(:),intent(in)            ::  file_arr
        integer(ik),                      intent(in)            ::  timestep ! No. of timesteps in the solution


        character(len = 100)                                    ::  type_name, compressor
        integer(ik)                                             ::  funit = 10, ifile, itime, d
        logical                                                 ::  exist


        !
        ! .pvd file type, collection of blockwise .vtu files
        !
        type_name = 'Collection'


        !
        ! Compressor to be used for any compressed data
        !
        compressor = 'vtkZLibDataCompressor'


        !
        ! Open new .pvd file or overwrite an exisiting file
        ! Write all blockwise .vtu filenames in the .pvd format
        ! Close .pvd file after writing output
        !
        inquire(file = trim(pvd_file), exist=exist)
        if (exist) then
            open(funit, file = trim(pvd_file), status = 'old', action = 'write')
        else
            open(funit, file = trim(pvd_file), status = 'new', action = 'write')
        end if
            write(funit,'(A)') '<?xml version="1.0"?>'
            write(funit,'(5A)') '<VTKFile type="',trim(type_name),'" version="0.1" byte_order="LittleEndian" Compressor="',trim(compressor),'">'
            write(funit,'(3A)') '   <',trim(type_name),'>'
            !
            ! Write individual .vtu file names to the .pvd file
            !
            do itime = 1,timestep

                d = (itime - 1)*data%ndomains()

                do ifile = 1,data%ndomains()

                    write(funit,'(A,I5,A,I5,3A)') '        <DataSet timestep="',itime - 1,'" part="',ifile - 1,'" file="',trim(file_arr(d + ifile)),'"/>'
                    
                end do
                
            end do
            write(funit,'(3A)') '   </',trim(type_name),'>'
            write(funit,'(A)') '</VTKFile>'
            close(funit)
            
            
    end subroutine write_pvd_final
    !***************************************************************************************************************************************************



















     
end module mod_vtk_file_unstr
