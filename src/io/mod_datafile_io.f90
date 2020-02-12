module mod_datafile_io
#include <messenger.h>

    use mod_chidg_mpi,          only: IRANK, GLOBAL_MASTER
    use type_chidg_data,        only: chidg_data_t
    use mod_string
    use mod_file_utilities
    implicit none

contains


    !>  Write computed functionals to files
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/20/2017
    !!
    !!  Restructured file
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/8/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine write_functionals_file(data)
        type(chidg_data_t), intent(in)  :: data

        character(len=50)       :: filename, func_name, ref_geom, aux_geom
        logical                 :: stat, func_stored
        integer(ik)             :: nfunc, nstep, unit_file, ifunc
        integer(ik)             :: ierr, igeom, istep, step_
        real(rk)                :: time_, value_(1)
        type(string_t)          :: noblank_name
        character(len=4)        :: file_format = '.txt'
        !character(len=20)       :: FMT = '(ES20.10E3,I6,F12.9)'
        character(len=12)       :: FMT = '(*,I6,F12.9)'

        ! Only Global Master proc writes the file
        if (IRANK == GLOBAL_MASTER) then

            ! Get number of functionals and ntime
            nfunc = data%sdata%functional%nfunc()
            nstep = data%sdata%functional%nstep()
        
            ! Write functionals
            do ifunc = 1,nfunc
            
                associate ( functional => data%functional_group%fcl_entities(ifunc)%func )

                    ! Retrive functional name
                    func_name = functional%get_name()
                    ! Retrieve functional reference geometries
                    ref_geom = functional%get_ref_geom()
                    ! Retrieve auxiliary reference geometries
                    aux_geom = functional%get_aux_geom()

                    ! Define the file name for the functional (replace blank spaces with '_')
                    noblank_name = functional%name%replace(' ','_')
                    filename = noblank_name%get() // file_format
                    
                    ! Display action
                    call write_line("   writing functional to    :", filename, ltrim=.false., io_proc=GLOBAL_MASTER)

                    ! Check if exists, if so delete it and recreate it
                    inquire(file=filename, exist=stat)
                    call write_file_header(filename,stat,func_name,ref_geom,aux_geom)

                    ! Open file 
                    open(newunit=unit_file, file=filename, status='old', position='append', action='write', iostat=ierr)
                    if ( ierr /= 0 ) call chidg_signal(FATAL,"write_functionals_file: error opening target file.")

                    ! Wirte data
                    do istep = 1,nstep
                        value_ = data%sdata%functional%get_func(ifunc,istep)
                        step_  = data%sdata%functional%get_step(istep)
                        time_  = data%sdata%functional%get_time(istep)
                        write(unit_file, *) value_(1),step_,time_
                    end do

                    ! Close file
                    close(unit_file)

                end associate
            end do
        end if

    end subroutine write_functionals_file
    !***************************************************************************************    





    !>  Write computed gradient of functionals wrt BC parameters to files
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine write_BC_sensitivities_file(data,ifunc,file_prefix)
        type(chidg_data_t), intent(in)  :: data
        integer(ik),        intent(in)  :: ifunc
        character(*),       intent(in)  :: file_prefix

        character(len=50)       :: filename, func_name, ref_geom, aux_geom
        logical                 :: stat, func_stored
        integer(ik)             :: nfunc, nstep, unit_file
        integer(ik)             :: ierr, igeom, istep, step_
        real(rk)                :: time_, value_(1)
        type(string_t)          :: noblank_name
        character(len=4)        :: file_format = '.txt'
        character(len=18)       :: FMT = '(ES20.14)'

        ! Only Global Master proc writes the file
        if (IRANK == GLOBAL_MASTER) then

            ! Get number of functionals and ntime
            !nfunc = data%sdata%functional%nfunc()
            !nstep = data%sdata%functional%nstep()
            filename = file_prefix // file_format
        
            ! Write functional BC sensitivity
            call write_line("   writing functional gradient to: ", trim(filename), ltrim=.false., io_proc=GLOBAL_MASTER)

            ! Delete old file if exists
            inquire(file=filename, exist=stat)
            if (stat) call delete_file(filename)

            ! Create and open a new one
            open(newunit=unit_file, file=trim(filename), status='new', action='write', iostat=ierr)
            if ( ierr /= 0 ) call chidg_signal(FATAL,"write_BC_sensitivities_file: error opening target file.")

            ! Write data
            write(unit_file, *) data%sdata%adjointbc%Ja_unsteady(ifunc)

            ! Close file
            close(unit_file)

        end if

    end subroutine write_BC_sensitivities_file
    !***************************************************************************************    


end module mod_datafile_io
