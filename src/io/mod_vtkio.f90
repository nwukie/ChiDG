module mod_vtkio

! Module containing subroutines for writing vtk output for a given ChiDG run

    use mod_kinds,              only: rk,ik,rdouble
    use mod_constants,          only: ONE,HALF,TWO,OUTPUT_RES

    use mod_vtk_file_unstr,     only: open_vtk_unstr,close_vtk_unstr,init_piece_unstr,end_piece_unstr,write_piece_coord,write_piece_data,init_cell,end_cell,write_piece_connectivity_data,write_pvd_final
    use mod_vtk_calc_func,      only: get_cons_var,get_piece_nums,get_piece_coord,get_piece_data,get_piece_connectivity_data

    use type_element,           only: element_t
    use type_blockvector,       only: blockvector_t
    use type_solverdata,        only: solverdata_t
    use type_chidg_data,        only: chidg_data_t

    implicit none



contains



    !>  Subroutine for writing separate .vtu files for each domain in the given geometry
    !!  Also writes final Paraview file (.pvd) as a multiblock collection
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date 10/28/2016
    !!
    !!
    !
    !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
    subroutine write_vtk_file(data,pvd_filename)

        type(chidg_data_t),intent(inout)                        ::  data
        character(len = *),intent(in)                           ::  pvd_filename


        integer(ik),parameter                                   ::  bo_type = 0_ik
        integer(ik)                                             ::  idom,ielem,nelem,noeq,s,num_pts,num_cells,ntime
        character(len = 100),dimension(:),allocatable           ::  cons_var
        real(rk),dimension(:),allocatable                       ::  X,Y,Z
        real(rk),dimension(:,:),allocatable                     ::  cons_var_val
        integer(ik),dimension(:,:),allocatable                  ::  connectivity,connectivity_A
        integer(ik),dimension(:),allocatable                    ::  offsets,types
        character(len = 100),dimension(:),allocatable           ::  file_arr
        integer(ik)                                             ::  itime,d ! Counters for outer time loop and file name array


        ntime = 1   ! No. of time steps in the solution file (1 for steady cases)


        !
        ! Allocate array for storing individual .vtu file names for each block over all time steps
        ! Each block corresponds to a ChiDG domain
        !
        allocate(file_arr(data%ndomains()*ntime))


        !
        ! Loop through all time steps of a time-dependent solution to get data
        !
        do itime = 1,ntime


            d = (itime - 1)*data%ndomains()

            !
            ! Loop through all domains to get data
            !
            do idom = 1,data%ndomains()

                !
                ! Get the file names for the individual vtk files for individual domains
                !
                write(file_arr(d + idom), fmt = '(a,i0,a,i0,a)') 'chidg_results_',itime - 1,'_',idom - 1,'.vtu'

                !
                ! Get number of elements in the current block
                !
                nelem = data%mesh(idom)%nelem

                !
                ! Get number of equations in the equation set
                !
<<<<<<< HEAD
                noeq = data%eqnset(1)%prop%nequations()
=======
                noeq = data%eqnset(idom)%prop%nequations()
>>>>>>> dev_vtk

                !
                ! Open individual .vtu files for individual domains
                !
                call open_vtk_unstr(file_arr(d + idom),bo_type)
 
                !
                ! Get conservative variable names
                !
                call get_cons_var(data,idom,noeq,cons_var)
                
                !
                ! Get number of points and cells created after sub-sampling in the current block
                !
                call get_piece_nums(data,nelem,num_pts,num_cells)
                
                !
                ! Initialize an unstructured grid piece (block) in the current .vtu file
                !
                call init_piece_unstr(file_arr(d + idom),num_pts,num_cells)
                
                !
                ! Get x,y and z coordinates for the current piece
                !
                call get_piece_coord(data,idom,nelem,num_pts,X,Y,Z)

                !
                ! Write x,y and z coordinates to the current .vtu file
                !
                call write_piece_coord(file_arr(d + idom),num_pts,X,Y,Z)

                !
                ! Get conservative variable values for the current piece
                !
                call get_piece_data(data,idom,nelem,num_pts,noeq,cons_var_val)

                !
                ! Write conservative variable values to the current .vtu file
                !
                call write_piece_data(file_arr(d + idom),noeq,num_pts,cons_var,cons_var_val)

                !
                ! Initialize cells field (<Cells>) in the current .vtu file
                !
                call init_cell(file_arr(d + idom))

                !
                ! Get connectivity, offsets and element types for the current piece
                !
                call get_piece_connectivity_data(data,idom,nelem,num_pts,num_cells,connectivity,offsets,types)

                !
                ! Write connectivity, offsets and element types to the current  .vtu file
                !
                call write_piece_connectivity_data(file_arr(d + idom),num_cells,connectivity,offsets,types)

                !
                ! End initialized cells field in the current .vtu file
                !
                call end_cell(file_arr(d + idom))

                !
                ! End current initialized piece
                !
                call end_piece_unstr(file_arr(d + idom))
                 
                !
                ! Close current .vtu file
                !
                call close_vtk_unstr(file_arr(d + idom))


            end do ! idom


        end do ! itime

        !
        ! Write the final Paraview data file (.pvd)
        ! This file is a multi block collection file for entire geometry
        !
        call write_pvd_final(data,pvd_filename,file_arr,ntime)


    end subroutine write_vtk_file
   !**************************************************************************************************************************************************************************************************************
 



















 end module mod_vtkio
