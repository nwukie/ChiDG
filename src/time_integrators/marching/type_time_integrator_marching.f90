module type_time_integrator_marching
#include<messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_time_integrator,   only: time_integrator_t
    use hdf5
    use h5lt
    use mod_hdf_utilities
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
    type, abstract, extends(time_integrator_t), public  :: time_integrator_marching_t

    contains

        ! Define deferred state initialization for a 'Marching' time integrator.
        procedure   :: initialize_state
        procedure   :: write_time_options
        procedure   :: read_time_options
        procedure   :: process_data_for_output

    end type time_integrator_marching_t
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
        class(time_integrator_marching_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data

        integer(ik)                 :: idom,ielem,ivar,ierr,eqn_ID
        real(rk),   allocatable     :: temp(:)

        associate( q          => data%sdata%q,             &
                   q_in       => data%sdata%q_in,          &
                   ntime_q    => data%sdata%q%get_ntime(), &
                   ntime_q_in => data%sdata%q_in%get_ntime())

            if (ntime_q == ntime_q_in) then
                q = q_in
            else if (ntime_q .ne. ntime_q_in .and. ntime_q == 1) then

                    do idom = 1,data%mesh%ndomains()

                        eqn_ID = data%mesh%domain(idom)%eqn_ID

                        if (allocated(temp)) deallocate(temp)
                        allocate(temp(data%mesh%domain(idom)%nterms_s), stat=ierr)
                        if (ierr /= 0) call AllocationError

                        do ielem = 1,data%mesh%domain(idom)%nelem
                            do ivar = 1,data%eqnset(idom)%prop%nprimary_fields()

                                temp = q_in%dom(idom)%vecs(ielem)%getvar(ivar,ntime_q_in)
                                call q%dom(idom)%vecs(ielem)%setvar(ivar,1,temp)

                            end do
                        end do

                    end do
            else
                call chidg_signal(FATAL, 'Initialization array incompatible with solution array')
            end if

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
        class(time_integrator_marching_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename

        integer(HID_T)                  :: fid
        integer(ik)                     :: ierr, iwrite
        

        do iwrite = 0,NRANK -1
            if (iwrite == IRANK) then

                !
                ! Write dt, no. of time steps and nwrite to hdf file
                !
                fid = open_file_hdf(filename)
                call set_time_integrator_hdf(fid, trim(data%time_manager%get_name()))
                call set_time_step_hdf(      fid, data%time_manager%dt              )
                call set_times_hdf(          fid, [data%time_manager%t]             )
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
        class(time_integrator_marching_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        character(*),                       intent(in)      :: filename
        character(*),                       intent(in)      :: read_type

        integer(HID_T)  :: fid


        !
        ! Open hdf file
        !
        fid = open_file_hdf(filename)


        !
        ! Read dt, no. of time steps and nwrite
        !
        select case(trim(read_type))
            case('run')
                ! For running a time-marching case, 
                data%time_manager%times       = get_times_hdf(fid)
                data%time_manager%t           = data%time_manager%times(1)
                data%time_manager%ntime       = size(data%time_manager%times)

            case('process')
                data%time_manager%time_scheme = trim(get_time_integrator_hdf(fid))
                data%time_manager%dt          = get_time_step_hdf(fid)
                data%time_manager%nsteps      = get_nsteps_hdf(fid)
                data%time_manager%nwrite      = get_nwrite_hdf(fid)
                data%time_manager%times       = get_times_hdf(fid)
                data%time_manager%t           = data%time_manager%times(1)
                data%time_manager%ntime       = size(data%time_manager%times)
            case default
                call chidg_signal(FATAL,"time_integrator_marching%read_time_options: Invalid read_type. 'run' or 'process'.")
        end select


        !
        ! Close file
        !
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
        class(time_integrator_marching_t),  intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data


        !
        ! Set q_out
        ! TODO: This part will change after implementation of time marching integrators
        !
        call data%sdata%q_out%init(data%mesh,data%time_manager%ntime)
        call data%sdata%q_out%set_ntime(data%time_manager%ntime)
        call data%sdata%q_out%clear()

        data%sdata%q_out = data%sdata%q_in


    end subroutine process_data_for_output
    !*******************************************************************************




end module type_time_integrator_marching
