module type_steady
    use messenger,              only: write_line
    use mod_kinds,              only: rk,ik
    use type_time_integrator,   only: time_integrator_t
    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_vector

    use mod_spatial,    only: update_space

    use mod_entropy,    only: compute_entropy_error
    implicit none
    private



    !>  Solve to steady-state.
    !!
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !----------------------------------------------------------------------------------------
    type, extends(time_integrator_t), public :: steady_t



    contains

        procedure   :: iterate

        final :: destructor

    end type steady_t
    !****************************************************************************************










contains





    !>  Steady equation set. No time.
    !!
    !!  Solving R(Q) = 0
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(steady_t),                        intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner


        !
        ! Simply solve the nonlinear system. No iteration in time.
        !
        call nonlinear_solver%solve(data,linear_solver,preconditioner)



    end subroutine iterate
    !******************************************************************************************







    
    subroutine destructor(self)
        type(steady_t),      intent(in) :: self

    end subroutine




end module type_steady
