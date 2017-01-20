module type_time_manager
#include <messenger.h>
    
    use mod_kinds,      only: rk,ik
    use mod_constants, only: PI,ZERO,ONE,TWO
    use type_rvector,   only: rvector_t
    use mod_io

    implicit none


    !>  Time manager data type
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/22/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------

    type, public        :: time_manager_t
        
        !Time scheme
        character(len=100)      :: time_scheme != 'steady'

        ! Unsteady time parameter
        real(rk)                :: dt          != 0.001_rk
        integer(ik)             :: time_steps  != 100
        integer(ik)             :: nwrite      != 10
        
        ! HB time parameter
        type(rvector_t)         :: freq_data   !> we have a limit of 100 freq's based on the array size in the namelist file
        type(rvector_t)         :: time_lev
    
    contains

        procedure   :: init             !< Initialization procedure to store all the time information needed
        procedure   :: push_freq        !< Procedure to push a new frequency unto freq_data checking for duplicate frequencies first
        procedure   :: loc_freq         !< Procedure to check if a frequency is already stored in freq_data
        procedure   :: set_name         !< Procedure to set the name of the time_integrator used
        procedure   :: get_name         !< Procedure to get the neme of the time_integrator used
    

    end type time_manager_t
    !--------------------------------------------------------------------------------------------------------


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
    !--------------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(time_manager_t),  intent(inout)   :: self

        character(len=:),   allocatable  :: user_msg, dev_msg
        integer(ik)                      :: i
        integer(ik)                      :: n_times
        
        !
        ! Check if the time scheme typed in belongs to the options available
        ! if so, initialize the time-Manager with appropriate attributes
        !
        !
        ! to allow different number of frequencies we can use a vector of inputs in mod_io.f90
        ! see here: http://docs.cray.com/books/S-3693-51/html-S-3693-51/i5lylchri.html 
        !


        select case (trim(time_integrator))
            
            case ('steady', 'Steady') 
                
                call self%set_name(time_integrator)



            case ('Forward_Euler', 'Forward Euler', 'forward euler', 'forward_euler', &
                  'Second Order Runge-Kutta', 'Explicit Midpoint', 'Second Order RK', &
                  'Modified Euler', 'Second Order Heun Method', &
                  'Ralston Method', 'Second Order Ralston Method', &
                  'Third Order Runge-Kutta', 'Third Order Kutta', 'Third Order RK', &
                  'Runge-Kutta Method', 'Fourth Runge-Kutta Method', 'Fourth Order RK Method', 'RK4', &
                  'Three-Eighth Rule', 'Fourth Order Kutta', &
                  'Backward_Euler', 'Backward Euler', 'backward euler', 'backward_euler')

                call self%set_name(time_integrator)

                !
                ! add dt, time_steps and nwrite to the time_manager
                !
                
                self%dt         = dt
                self%time_steps = time_steps
                self%nwrite     = nwrite


            
            case ('Harmonic Balance', 'Harmonic_Balance', 'harmonic balance', 'harmonic_balance', 'HB')
                
                call self%set_name(time_integrator)
                !
                ! Verify that at least one frequency has been passed in
                !
                user_msg = "time_integrator%init: The time scheme set needs frequencies passed in along with it. &
                            Please define at least one frequency different than 0 in chidg.nml"
                if ( (maxval(frequencies) == ZERO) .and. (minval(frequencies) == ZERO) ) call chidg_signal_one(FATAL,user_msg,trim(time_integrator))

                
                !
                ! add number of frequencies and frequencies to time_manager
                !

                do i = 1, size(frequencies)
                    
                    if ( frequencies(i) /= ZERO ) then
                        call self%push_freq(frequencies(i))
                    end if

                end do


                !
                ! Calculate and store time levels
                !
                
                ! n = 2*(number of frequencies) + 1
                n_times = TWO*self%freq_data%size() + ONE

                do i = 1, n_times

                   call  self%time_lev%push_back( (TWO*PI)/minval(self%freq_data%data()) * (i/n_times) )

                end do

            case default
                user_msg = "We can't seem to find a time integrator that matches the input string. &
                            Maybe check that the time integrator string in the input file or drive &
                            script is valid"

                dev_msg  = "Check that the time integrator is registered properly in type_time_manager."

                call chidg_signal_two(OOPS, user_msg,trim(time_integrator), dev_msg=dev_msg)

        end select


    end subroutine init
    !-----------------------------------------------------------------------------------------------------             





    !> Push a new frequency to self%freq
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/25/2016
    !!
    !!  @param[in] frequency
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine push_freq(self,frequency)
        class(time_manager_t),  intent(inout)       :: self
        real(rk),               intent(in)          :: frequency
        
        logical     :: test

        character(len=:),   allocatable    :: user_msg

        test = self%loc_freq(frequency)

        !
        ! Verify if the frequency is already stored, if not it is added
        !
        if (test) then 

            user_msg = "time_manager%push_freq: A HB frequency has been typed in twice, only one will be considered"
        
            call chidg_signal(WARN, user_msg)

        else

            call self%freq_data%push_back(frequency)

        end if

    end subroutine push_freq
    !------------------------------------------------------------------------------------------------------









    !> locate a frequency in self%freq_data, return TRUE if the frequency is already stored
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/25/2016
    !!
    !!  @param[in] frequency    >frequency we are looking for
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function loc_freq(self,frequency) result(val)
        class(time_manager_t),  intent(in)       :: self
        real(rk),               intent(in)       :: frequency
        
        logical     :: val
        integer(ik) :: i

        val = .false.

        do i = 1, self%freq_data%size()
            
            if ( self%freq_data%data_(i) == frequency ) then
                val = .true.
                exit
            end if
            
        end do  

    end function loc_freq
    !------------------------------------------------------------------------------------------------------





    
    !> Set Name of the time integrator based on the input string from namelist
    !!
    !!  @author Matteo Ugolotti
    !!  @date   01/18/2017
    !!
    !!  @param[in]     >name of the time integrator
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine set_name(self,name_string)
        class(time_manager_t),                 intent(inout)       :: self
        character(100),                        intent(in)          :: name_string
        
        self%time_scheme = trim(name_string)

    end subroutine set_name
    !------------------------------------------------------------------------------------------------------


    

    
    !> Get Name of the time integrator based on the input string from namelist
    !!
    !!  @author Matteo Ugolotti
    !!  @date   01/18/2017
    !!
    !!  @param[out]     >name of the time integrator
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function get_name(self) result(res)
        class(time_manager_t),  intent(inout)       :: self
        
        character(len=:),   allocatable        :: res

        res = self%time_scheme

    end function get_name
    !------------------------------------------------------------------------------------------------------


    


end module type_time_manager
