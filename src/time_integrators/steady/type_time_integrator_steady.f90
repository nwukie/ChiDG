module type_time_integrator_steady
#include<messenger.h>
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_time_integrator,   only: time_integrator_t
    use hdf5
    use h5lt
    use mod_hdf_utilities
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

!        data%sdata%q = data%sdata%q_in

        associate( q => data%sdata%q, q_in => data%sdata%q_in)
            q = q_in
            call q%assemble()
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

        integer(HID_T)  :: fid
        integer(ik)     :: ierr, iwrite


        do iwrite = 0,NRANK-1
            if (iwrite == IRANK) then
                
                fid = open_file_hdf(filename)
                call set_time_integrator_hdf(fid, trim(data%time_manager%get_name()))
                call set_times_hdf(          fid, [ZERO]                            )
                call set_time_step_hdf(      fid, ZERO                              )
                call set_nsteps_hdf(         fid, data%time_manager%nsteps          )
                call set_nwrite_hdf(         fid, data%time_manager%nwrite          )
                call close_file_hdf(fid)

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do


    end subroutine write_time_options
    !*******************************************************************************



    !>  Read time integrator options from hdf file
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine read_time_options(self,data,filename,read_type)
        class(time_integrator_steady_t),    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename
        character(*),                       intent(in)      :: read_type

        integer(HID_T)  :: fid

        
        !
        ! Open file
        !
        fid = open_file_hdf(filename)


        !
        ! Get ntime
        !
        select case(trim(read_type))
            case('run')

            case('process')
                data%time_manager%time_scheme = trim(get_time_integrator_hdf(fid))
                data%time_manager%times       = get_times_hdf(fid)
                data%time_manager%nsteps      = get_nsteps_hdf(fid)
                data%time_manager%nwrite      = get_nwrite_hdf(fid)
                data%time_manager%ntime       = size(data%time_manager%times)
                data%time_manager%t           = data%time_manager%times(1)

            case default
                call chidg_signal(FATAL,"time_integrator_steady%read_time_options: Invalid read_type. 'run' or 'process'.")
        end select

        call close_file_hdf(fid)
        

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


        associate( q_out => data%sdata%q_out, q_in => data%sdata%q_in)

        q_out = q_in

        end associate


    end subroutine process_data_for_output
    !*******************************************************************************




end module type_time_integrator_steady
