module type_time_manager
#include <messenger.h>
    
    use mod_kinds,       only: rk,ik
    use mod_constants,   only: PI,ZERO,ONE,TWO
    use type_rvector,    only: rvector_t
    use mod_HB_matrices, only: calc_pseudo_spectral_operator
    use mod_io

    implicit none


    !>  Time manager data type
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/22/2016
    !!
    !!
    !------------------------------------------------------------------------------------------

    type, public        :: time_manager_t
        
        !Time scheme
        character(len=100)      :: time_scheme  ! 'steady'

        ! Unsteady time parameter
        real(rk)                :: t            ! Current time
        real(rk),   allocatable :: times(:)     ! List of times: times(ntime)
        real(rk),   allocatable :: freqs(:)     ! List of frequencies: freqs(ntime)
        real(rk)                :: dt           ! Time interval

        integer(ik)             :: itime  = 1   ! Current time index
        integer(ik)             :: ntime  = 1   ! Number of time levels in HB (=1 for steady)
        integer(ik)             :: nsteps = 1   ! Number of time steps in time_marching solution
        integer(ik)             :: nwrite = 0
        
        ! Harmonic Balance pseudo-spectral matrices
        real(rk), allocatable   :: D(:,:)
        real(rk), allocatable   :: E(:,:)

    contains

        procedure   :: init         ! Initialization procedure to store all the time information needed
        procedure   :: set_name     ! Procedure to set the name of the time_integrator used
        procedure   :: get_name     ! Procedure to get the neme of the time_integrator used
    

    end type time_manager_t
    !------------------------------------------------------------------------------------------


contains



    !> Time Manager initialization
    !!
    !!  read data from namelist (mod_io.f90) and save time informations needed
    !!  
    !!  
    !!  @author Matteo Ugolotti
    !!  @date   12/25/2016
    !!
    !!  
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(time_manager_t),  intent(inout)   :: self

        character(:),   allocatable     :: user_msg, dev_msg
        integer(ik)                     :: i, nfreq, ierr
        
        
        !
        ! Check if the time scheme typed in belongs to the options available
        ! if so, initialize the time-Manager with appropriate attributes
        !
        !
        ! to allow different number of frequencies we can use a vector of inputs in mod_io.f90
        ! see here: http://docs.cray.com/books/S-3693-51/html-S-3693-51/i5lylchri.html 
        !
        select case (trim(time_integrator))
            
            !
            ! TIME: Steady
            !
            case ('steady', 'Steady') 
                
                call self%set_name(time_integrator)

                self%dt     = 0
                self%ntime  = 1
                self%nsteps = 1
                self%nwrite = 0     ! don't write intermediate file
                self%times  = [ZERO]


            !
            ! TIME: Time-Marching
            !
            case ('Forward_Euler', 'Forward Euler', 'forward euler', 'forward_euler',    &
                  'Second Order Runge-Kutta', 'Explicit Midpoint', 'Second Order RK',    &
                  'Modified Euler', 'Second Order Heun Method',                          &
                  'Ralston Method', 'Second Order Ralston Method',                       &
                  'Third Order Runge-Kutta', 'Third Order Kutta', 'Third Order RK',      &
                  'Runge-Kutta Method', 'Fourth Runge-Kutta Method',                     &
                  'Fourth Order RK Method', 'RK4',                                       &
                  'Three-Eighth Rule', 'Fourth Order Kutta',                             &
                  'Backward_Euler', 'Backward Euler', 'backward euler', 'backward_euler',&
                  'DIRK', 'DIRK_coupled_oscillator')

                call self%set_name(time_integrator)

                !
                ! add dt, ntimes, nsteps and nwrite to the time_manager
                !
                
                self%dt     = dt
                self%nsteps = time_steps
                self%ntime  = 1
                self%nwrite = nwrite
                self%times  = [ZERO]


            
            !
            ! TIME: Time-Spectral
            !
            case ('Harmonic Balance', 'Harmonic_Balance', 'harmonic balance',   &
                  'harmonic_balance', 'HB')
                
                call self%set_name(time_integrator)
                self%nsteps     = 1
                self%nwrite     = 0     ! don't write intermediate file
                !
                ! Verify that at least one frequency has been passed in
                !
                user_msg = "time_integrator%init: The time scheme set needs frequencies &
                            passed in along with it. Please define at least one frequency &
                            different than 0 in chidg.nml"
                if ( (maxval(frequencies) == ZERO) .and. (minval(frequencies) == ZERO) ) call chidg_signal_one(FATAL,user_msg,trim(time_integrator))


                !
                ! Determine number of frequencies and time levels
                !   : determine number of frequencies, by number of nonzero entries in IO variable 'frequencies'
                !   : number of time levels = 2*nfreq + 1
                !
                nfreq = 0
                do i = 1,size(frequencies)
                    if ( frequencies(i) /= ZERO ) then
                        nfreq = nfreq + 1
                    end if
                end do
                self%ntime = 2*nfreq + 1

                
                !
                ! Allocate times(:),freqs(:) storage
                !
                if (allocated(self%times)) deallocate(self%times)
                if (allocated(self%freqs)) deallocate(self%freqs)
                allocate(self%times(self%ntime), self%freqs(nfreq), stat=ierr)
                if (ierr /= 0) call AllocationError

                
                !
                ! Store input frequencies
                !
                do i = 1,size(frequencies)
                    if ( frequencies(i) /= ZERO ) then
                        self%freqs(i) = frequencies(i)
                    end if
                end do

                !
                ! Compute, store time levels
                !
                do i = 1,self%ntime
                    !self%times(i) = ((TWO*PI)/minval(self%freqs)) * (real(i)/real(self%ntime))
                    self%times(i) = ((TWO*PI)/minval(abs(self%freqs))) * (real(i)/real(self%ntime))
                end do


                !
                ! Compute the pseudo spectral operator when the HB time integrator is specified
                ! in the namelist file
                !
                self%D = calc_pseudo_spectral_operator(self%freqs,self%times)


            case default
                user_msg = "We can't seem to find a time integrator that matches the input &
                            string. Maybe check that the time integrator string in the input &
                            file or drive script is valid"

                dev_msg  = "Check that the time integrator is registered properly in type_time_manager."

                call chidg_signal_two(OOPS, user_msg,trim(time_integrator), dev_msg=dev_msg)

        end select


    end subroutine init
    !------------------------------------------------------------------------------------------




    
    !> Set Name of the time integrator based on the input string from namelist
    !!
    !!  @author Matteo Ugolotti
    !!  @date   01/18/2017
    !!
    !!  @param[in]     >name of the time integrator
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_name(self,name_string)
        class(time_manager_t),                 intent(inout)       :: self
        character(100),                        intent(in)          :: name_string
        
        self%time_scheme = trim(name_string)

    end subroutine set_name
    !-----------------------------------------------------------------------------------------


    

    
    !> Get Name of the time integrator based on the input string from namelist
    !!
    !!  @author Matteo Ugolotti
    !!  @date   01/18/2017
    !!
    !!  @param[out]     >name of the time integrator
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_name(self) result(res)
        class(time_manager_t),  intent(inout)   :: self
        
        character(:),   allocatable :: res

        res = self%time_scheme

    end function get_name
    !-----------------------------------------------------------------------------------------


    


end module type_time_manager
