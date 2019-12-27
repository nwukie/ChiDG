module type_storage_flags
#include <messenger.h>

    implicit none

    !>  This data type contains the flags for specialized initialization of the 
    !!  storage required for the specific solver
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/3/2018
    !!
    !-------------------------------------------------------------------------
    type, public ::    storage_flags_t
        

        ! Primal solver storage
        !------------------------
        logical     :: q 
        logical     :: q_in 
        logical     :: q_out
        logical     :: dq
        logical     :: rhs
        logical     :: lhs
        logical     :: dt
        logical     :: function_status

        ! Adjoint solver storage
        !------------------------
        logical     :: v 
        logical     :: q_time 
        logical     :: Jq
        logical     :: lhs_trans
        logical     :: Rd
        logical     :: Rd_trans
        logical     :: vRd
        logical     :: solver_iter
        logical     :: solver_time
        logical     :: solver_err

        ! AdjointX solver storage
        !------------------------
        logical     :: vRx 
        logical     :: Rx 
        logical     :: Jx 
        logical     :: Jx_unsteady
        logical     :: Rx_trans 
        logical     :: wAx 

        ! AdjointBC solver storage
        !------------------------
        logical     :: vRa 
        logical     :: wAa 
        logical     :: Ra 
        logical     :: Ja_unsteady
        
        ! Check functional input 
        !------------------------
        logical     :: func_check               ! Check if user input at least a functional
                                                ! if not, adjoint or adjointx not doable
        
        contains

        procedure   :: set
        procedure   :: clear

    end type storage_flags_t
    !*************************************************************************


contains






    !>  Set all flags to false 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/3/2018
    !!
    !-------------------------------------------------------------------------
    subroutine clear(self)
        class(storage_flags_t),  intent(inout)   :: self


        self%q                  = .false. 
        self%q_in               = .false.
        self%q_out              = .false.
        self%dq                 = .false.
        self%rhs                = .false.
        self%lhs                = .false.
        self%dt                 = .false.
        self%function_status    = .false.

        self%v                  = .false.
        self%q_time             = .false.
        self%Jq                 = .false.
        self%lhs_trans          = .false.
        self%vRd                = .false.
        self%Rd                 = .false.
        self%Rd_trans           = .false.
        self%solver_iter        = .false.
        self%solver_time        = .false.
        self%solver_err         = .false.

        self%vRx                = .false.
        self%Rx                 = .false.
        self%Jx                 = .false.
        self%Jx_unsteady        = .false.
        self%Rx_trans           = .false.
        self%wAx                = .false.

        self%vRa                = .false.
        self%wAa                = .false.
        self%Ra                 = .false.
        self%Ja_unsteady        = .false.

        
        self%func_check         = .false.



    end subroutine clear
    !*************************************************************************







    !>  Set flags based on solver type
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/3/2018
    !!
    !-------------------------------------------------------------------------
    subroutine set(self,storage_type,nauxiliary_fields)
        class(storage_flags_t),  intent(inout)  :: self
        character(*),            intent(in)     :: storage_type
        integer(ik),             intent(in)     :: nauxiliary_fields

        call self%clear()

        select case(trim(storage_type))
            case('primal storage')

                self%q                  = .true. 
                self%q_in               = .true. 
                self%q_out              = .true. 
                self%dq                 = .true. 
                self%rhs                = .true. 
                self%lhs                = .true. 
                self%dt                 = .true. 
                self%function_status    = .true. 


            case('adjoint storage')

                self%q                  = .true. 
                self%q_in               = .true. 
                self%rhs                = .true. 
                self%lhs                = .true. 
                self%function_status    = .true. 

                self%v                  = .true.
                self%q_time             = .true.
                self%Jq                 = .true.
                self%lhs_trans          = .true.
                self%solver_iter        = .true.
                self%solver_time        = .true.
                self%solver_err         = .true.

                self%func_check         = .true.

                !if (compute_auxiliary) then
                if (nauxiliary_fields > 0) then
                    self%Rd             = .true.
                    self%Rd_trans       = .true.
                    self%vRd            = .true.
                end if

            case('adjointx storage')

                self%q                  = .true. 
                self%q_in               = .true. 
                self%rhs                = .true. 
                self%lhs                = .true. 
                self%function_status    = .true. 

                self%v                  = .true.
                self%q_time             = .true.

                self%vRx                = .true.
                self%Rx                 = .true.
                self%Jx                 = .true.
                self%Jx_unsteady        = .true.
                self%Rx_trans           = .true.

                self%func_check         = .true.
                
                !if (compute_auxiliary) then
                if (nauxiliary_fields > 0) then
                    self%wAx             = .true.
                end if

            case('adjointbc storage')

                self%q                  = .true. 
                self%q_in               = .true. 
                self%function_status    = .true. 

                self%v                  = .true.
                self%q_time             = .true.

                self%vRa                = .true.
                self%wAa                = .true.
                self%Ra                 = .true.
                self%Ja_unsteady        = .true.

                self%func_check         = .true.

            case default

                call chidg_signal_one(FATAL,'storage_flags_t%set: storage type input argument not recognized.', storage_type)

        end select


    end subroutine set
    !*************************************************************************



end module type_storage_flags

