module type_chidg_functional
#include<messenger.h>
    
    use mod_kinds,      only: ik,rk
    use mod_constants,  only: ZERO
    use type_rvector,   only: rvector_t
    use type_ivector,   only: ivector_t
    use mod_io,         only: verbosity
    use type_svector,   only: svector_t
    use mod_string
    implicit none


    !>  Storage for real functionals (computed with or w/o adjoint)
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    type,   public  ::chidg_functional_t
    
        ! Functional storage
        type(rvector_t),    allocatable :: func(:) ! Real values of the functionals at each nwrite
        type(ivector_t)                 :: step
        type(rvector_t)                 :: time

        ! Initialization completed
        logical :: functional_storage_initialized = .false.

    contains

        procedure   :: init
        procedure   :: check_functional_stored
        procedure   :: nfunc
        procedure   :: nstep
        procedure   :: get_func
        procedure   :: get_step
        procedure   :: get_time
        procedure   :: release
        procedure   :: report
         
    end type chidg_functional_t
    !******************************************************************************




contains





    
    !>  Storage for real functionals (computed with or w/o adjoint)
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    subroutine init(self,nfunc)
        class(chidg_functional_t),      intent(inout)   :: self
        integer(ik),                    intent(in)      :: nfunc

        integer(ik)     :: ierr

        if (allocated(self%func)) deallocate (self%func)
        allocate (self%func(nfunc), stat = ierr)
        if (ierr/=0) call AllocationError
        
        self%functional_storage_initialized = .true.

    end subroutine init
    !******************************************************************************







    !>  Check if any functional has been stored. Return .true. if at least one
    !!  functional is stored. This is used in write_solution subroutine
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    function check_functional_stored(self) result(stat)
        class(chidg_functional_t),      intent(in)   :: self

        logical     :: stat, functional_stored

        stat = .false.

        if (self%functional_storage_initialized) then
            functional_stored = (self%func(1)%size() /= 0)
            if (functional_stored) then
                stat = .true.
            end if
        end if
        
    end function check_functional_stored
    !******************************************************************************

    
    
    
    
    !>  Return the number of functionals stored 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2017
    !!
    !------------------------------------------------------------------------------
    function nfunc(self) result(func_number)
        class(chidg_functional_t),      intent(in)   :: self

        integer(ik) :: func_number

        func_number = size(self%func)

    end function nfunc
    !******************************************************************************




    
    !>  Return the number of steps (times) at which the functionala have computed and
    !!  stored (for unsteady/HB)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2017
    !!
    !------------------------------------------------------------------------------
    function nstep(self) result(func_times)
        class(chidg_functional_t),      intent(in)   :: self

        integer(ik) :: func_times

        func_times = self%step%size()

    end function nstep
    !******************************************************************************





    !>  Return the functional data stored
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/20/2017
    !!
    !------------------------------------------------------------------------------
    function get_func(self,ifunc,istep) result(data_)
        class(chidg_functional_t),      intent(in)   :: self
        integer(ik),                    intent(in)   :: ifunc
        integer(ik),    optional,       intent(in)   :: istep

        real(rk)    ,allocatable  :: data_(:)
        integer(ik)               :: index_
        
        if (present(istep)) then 
            ! Check for index out of range
            if (istep > self%nstep()) call chidg_signal(FATAL,'chidg_functional_t%get_func: index out of range.')
            
            if (.not. allocated(data_)) allocate(data_(1))
            data_(1) = self%func(ifunc)%at(istep)
        else
            data_ = self%func(ifunc)%data()
        end if

    end function get_func
    !******************************************************************************





    !>  Return the step stored
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/20/2017
    !!
    !------------------------------------------------------------------------------
    function get_step(self,istep) result(step_)
        class(chidg_functional_t),      intent(in)   :: self
        integer(ik),                    intent(in)   :: istep

        integer(ik)               :: index_, step_
            
        if (istep > self%nstep()) then
            call chidg_signal(FATAL,'chidg_functional_t%get_step: index out of range.')
        else
            step_ = self%step%at(istep)
        end if

    end function get_step
    !******************************************************************************





    !>  Return the time stored
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/20/2017
    !!
    !------------------------------------------------------------------------------
    function get_time(self,istep) result(time_)
        class(chidg_functional_t),      intent(in)   :: self
        integer(ik),                    intent(in)   :: istep

        integer(ik)     :: index_
        real(rk)        :: time_

        if (istep > self%nstep()) then
            call chidg_signal(FATAL,'chidg_functional_t%get_time: index out of range.')
        else
            time_ = self%time%at(istep)
        end if
        
    end function get_time
    !******************************************************************************




    
    
    !>  Write report for functional, if any functional is registered
    !!
    !!  @author Matteo Ugolotti
    !!  @date   1/30/2019
    !!
    !!  TODO: improve it
    !!
    !------------------------------------------------------------------------------------
    subroutine report(self,funcs_name)
        class(chidg_functional_t),  intent(inout)   :: self
        type(svector_t),            intent(in)      :: funcs_name


        integer(ik)          :: step, ifunc
        real(rk)             :: time, func 
        logical              :: print_report
        
        
        ! If no functional is stored do not print the report
        print_report = self%check_functional_stored()

        if (print_report) then

            call write_line(' ',silence=(verbosity<2))
            call write_line('---------------------------------      Functional Report     ----------------------------------',silence=(verbosity<2))

            do ifunc = 1,size(self%func)
                call write_line(' ',silence=(verbosity<2))
                step = self%step%at(self%step%size())
                time = self%time%at(self%time%size())
                func = self%func(ifunc)%at(self%step%size())
                call write_line( string_to_upper(funcs_name%data_(ifunc)%get()), func , columns=.True., column_width=30,silence=(verbosity<2))
            end do !ifunc

            call write_line(' ',silence=(verbosity<2))
            call write_line('-----------------------------------------------------------------------------------------------',silence=(verbosity<2))

        end if ! print_report

    end subroutine report
    !************************************************************************************



    !>  Release memory 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/13/2017
    !!
    !------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_functional_t),      intent(inout)   :: self

        if (allocated(self%func)) deallocate (self%func)

    end subroutine release
    !******************************************************************************





end module type_chidg_functional
