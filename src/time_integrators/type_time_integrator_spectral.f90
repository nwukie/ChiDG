module type_time_integrator_spectral
#include<messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_time_integrator,   only: time_integrator_t
    use hdf5
    use h5lt
    use mod_hdf_utilities,      only: open_file_hdf, close_file_hdf, get_ntimes_hdf, &
                                      set_time_integrator_hdf, get_time_integrator_hdf
    use mod_HB_post,            only: get_post_processing_data
    use mpi_f08
    implicit none


    !>  Abstraction for time integrators.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!  @date   2/7/2017
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !---------------------------------------------------------------------------------------
    type, abstract, extends(time_integrator_t), public  :: time_integrator_spectral_t

    contains

        ! Define deferred state initialization for a 'Spectral' time integrator.
        procedure   :: initialize_state
        procedure   :: write_time_options
        procedure   :: read_time_options
        procedure   :: process_data_for_output

    end type time_integrator_spectral_t
    !*****************************************************************************************



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
        class(time_integrator_spectral_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data

        associate( q => data%sdata%q, q_in => data%sdata%q_in)

            q = q_in

        end associate

    end subroutine initialize_state
    !*******************************************************************************



    !>  Write time integrator options to hdf file
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine write_time_options(self,data,filename)
        class(time_integrator_spectral_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename

        integer(HID_T)          :: fid
        integer(kind = 8)       :: nfreq, ntime, SIZE_ONE = 1
        real(rk),   allocatable :: freq(:), time_lev(:)
        integer                 :: ierr, iwrite


        do iwrite = 0,NRANK - 1
            if (iwrite == IRANK) then
                
                !
                ! Assuming that the file exists, open the hdf file
                !
                fid = open_file_hdf(filename)


                !
                ! Set number of frequencies and time levels
                ! Also set frequency and time level data
                !
                nfreq = int(data%time_manager%freq_data%size(),8)
                ntime = int(data%time_manager%time_lev%size(),8)

                if (allocated(freq) .and. allocated(time_lev)) deallocate(freq,time_lev)
                allocate(freq(nfreq), time_lev(ntime), stat=ierr)
                if (ierr /= 0) call AllocationError

                freq     = data%time_manager%freq_data%data()
                time_lev = data%time_manager%time_lev%data()


                !
                ! Write time_integrator name to hdf file
                !
                !call h5ltset_attribute_string_f(fid,"/","time_integrator",trim(data%time_manager%get_name()), ierr)
                !if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltset_attribute_string_f")
                call set_time_integrator_hdf(fid,trim(data%time_manager%get_name()))

                call h5ltset_attribute_int_f(fid,"/","nsteps",[data%time_manager%nsteps],SIZE_ONE,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltset_attribute_int_f")

                call h5ltset_attribute_int_f(fid,"/","nwrite",[data%time_manager%nwrite],SIZE_ONE,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltset_attribute_int_f")
                    
            
                !
                ! Write frequencies and time levels to hdf file
                !
                call h5ltset_attribute_double_f(fid,"/","frequencies",freq,nfreq,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error_h5ltset_attribute_double_f")

                call h5ltset_attribute_double_f(fid,"/","time_levels",time_lev,ntime,ierr)
                if (ierr /= 0) call chidg_signal(FATAL,"Error_h5ltset_attribute_double_f")

                
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
    subroutine read_time_options(self,data,filename)
        class(time_integrator_spectral_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename

        integer(HID_T)                  :: fid
        integer(kind = 8)               :: nfreq, ntime
        character(:),   allocatable     :: temp_string
        real(rk),       allocatable     :: freq(:), time_lev(:)
        integer(ik)                     :: ierr, ifreq, itime 
        integer(ik),    dimension(1)    :: nsteps, nwrite

        
        !
        ! Open hdf file
        !  
        fid = open_file_hdf(filename)


        !
        ! Set no. of time levels and frequencies
        !
        ntime = get_ntimes_hdf(fid)
        nfreq = (ntime - 1)/2


        !
        ! Allocate frequency and time level data arrays
        !
        if (allocated(freq) .and. allocated(time_lev)) deallocate(freq,time_lev)
        allocate(freq(nfreq), time_lev(ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        
        !
        ! Read time integrator name
        !
        !call h5ltget_attribute_string_f(fid,"/","time_integrator",temp_string,ierr)
        !if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltget_attribute_string_f")
        temp_string = get_time_integrator_hdf(fid)

        call h5ltget_attribute_int_f(fid,"/","nsteps",nsteps,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltget_attribute_int_f")

        call h5ltget_attribute_int_f(fid,"/","nwrite",nwrite,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"Error h5ltget_attribute_int_f")


        !
        ! Read frequencies and time levels
        !
        call h5ltget_attribute_double_f(fid,"/","frequencies",freq,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"Error_h5ltget_attribute_double_f")

        call h5ltget_attribute_double_f(fid,"/","time_levels",time_lev,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"Error_h5ltget_attribute_double_f")

        call close_file_hdf(fid)


        !
        ! Set time_options in time_manager
        !
        data%time_manager%time_scheme = trim(temp_string)
        data%time_manager%ntime       = ntime
        data%time_manager%nsteps      = nsteps(1)
        data%time_manager%nwrite      = nwrite(1)

        do ifreq = 1,int(nfreq,ik)

            call data%time_manager%freq_data%push_back(freq(ifreq))

        end do

        do itime = 1,int(ntime,ik)

            call data%time_manager%time_lev%push_back(time_lev(itime))

        end do


    end subroutine read_time_options
    !*******************************************************************************



    !>  Modify data for post processing
    !!  
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine process_data_for_output(self,data)
        class(time_integrator_spectral_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data


        !
        ! Set q_out: all subroutines defined in mod_HB_post
        !
        call get_post_processing_data(data)


    end subroutine process_data_for_output
    !*******************************************************************************




end module type_time_integrator_spectral
