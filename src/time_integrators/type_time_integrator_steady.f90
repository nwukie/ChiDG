module type_time_integrator_steady
#include<messenger.h>
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_time_integrator,   only: time_integrator_t
    use hdf5
    use h5lt
    use mod_hdf_utilities,      only: open_file_hdf, close_file_hdf, &
                                      get_ntimes_hdf, set_time_integrator_hdf,  &
                                      get_time_integrator_hdf
    use mpi_f08
    implicit none


    !>  Abstraction for 'Steady' time integrators.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/7/2017
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    type, abstract, extends(time_integrator_t), public  :: time_integrator_steady_t

    contains

        ! Define state initialization for a 'Steady' time integrator.
        procedure   :: initialize_state 
        procedure   :: write_time_options
        procedure   :: read_time_options
        procedure   :: process_data_for_output

    end type time_integrator_steady_t
    !*******************************************************************************



contains




    !>  Initialize the working state vector from input data.
    !!
    !!  Interpret data read from file for use in working time integrator.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2017
    !!
    !!
    !-------------------------------------------------------------------------------
    subroutine initialize_state(self,data)
        class(time_integrator_steady_t),    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data

        associate( q => data%sdata%q, q_in => data%sdata%q_in)

            q = q_in

        end associate

    end subroutine initialize_state
    !*******************************************************************************



    !>  Write time_integrator options to hdf file
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine write_time_options(self,data,filename)
        class(time_integrator_steady_t),    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename

        integer(kind = 8)   :: fid
        integer(ik)         :: ierr, iwrite


        do iwrite = 0,NRANK - 1
            if (iwrite == IRANK) then
                
                !
                ! Assuming file exists, open hdf file
                !
                fid = open_file_hdf(filename)


!                !
!                ! Write time_integrator name to hdf file
!                !
!                call h5ltset_attribute_string_f(fid,"/","time_integrator",trim(data%time_manager%get_name()), ierr)
!                if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltset_attribute_string_f")
                call set_time_integrator_hdf(fid,trim(data%time_manager%get_name()))

                
                call close_file_hdf(fid)

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do


    end subroutine
    !*******************************************************************************



    !>  Read time integrator options from hdf file
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine read_time_options(self,data,filename)
        class(time_integrator_steady_t),    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename

        integer(HID_T)              :: fid
        character(:),   allocatable :: temp_string
        integer(ik)                 :: ierr, ntime

        
        !
        ! Open file
        !
        fid = open_file_hdf(filename)


        !
        ! Read time_integrator name
        !
        !call h5ltget_attribute_string_f(fid,"/","time_integrator",temp_string,ierr)
        !if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltget_attribute_string_f")
        temp_string = get_time_integrator_hdf(fid)


        !
        ! Get ntime
        !
        ntime = get_ntimes_hdf(fid)

        call close_file_hdf(fid)
        

        !
        ! Set time_options in time_manager
        !
        data%time_manager%time_scheme = trim(temp_string)
        data%time_manager%ntime       = ntime


    end subroutine read_time_options
    !*******************************************************************************



    !>  Modify data for post processing
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine process_data_for_output(self,data)
        class(time_integrator_steady_t),    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data


        !
        ! Set q_out: for steady runs no time_specific modifications are necessary
        !
        call data%sdata%q_out%init(data%mesh,data%time_manager%ntime)
        call data%sdata%q_out%set_ntime(data%time_manager%ntime)
        call data%sdata%q_out%clear()

        data%sdata%q_out = data%sdata%q_in


    end subroutine process_data_for_output
    !*******************************************************************************




end module type_time_integrator_steady
