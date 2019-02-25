module mod_nonlinear_solver
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_dict,              only: dict_t

    ! Import solverdata types
    use type_newton,            only: newton_t
    use type_petsc_nonlinear,   only: petsc_nonlinear_t
    use type_jfnk,              only: jfnk_t
    use type_quasi_newton_sst,  only: quasi_newton_sst_t
    use type_quasi_newton_rs,   only: quasi_newton_rs_t
    implicit none


    ! Instantiate solver types for sourcing
    type(newton_t)              :: NEWTON
    type(petsc_nonlinear_t)     :: PETSC_NONLINEAR
    type(jfnk_t)                :: JFNK 
    type(quasi_newton_sst_t)    :: QUASI_NEWTON_SST
    type(quasi_newton_rs_t)     :: QUASI_NEWTON_RS

    logical :: initialized = .false.



contains




    !>  Create a concrete nonlinear solver.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine create_nonlinear_solver(string,instance,options)
        character(*),                           intent(in)      :: string
        class(nonlinear_solver_t), allocatable, intent(inout)   :: instance
        type(dict_t), optional,                 intent(inout)   :: options

        character(:),   allocatable :: user_msg, dev_msg

        select case (trim(string))

            case ('newton','Newton','NEWTON')
                allocate(instance, source=NEWTON)

            case ('petsc','PETSC')
                allocate(instance, source=PETSC_NONLINEAR)

            case ('jfnk','Jfnk','JFNK')
                allocate(instance, source=JFNK)

            case ('quasi_newton_sst','Quasi_Newton_SST','quasi-newton-sst','Quasi-Newton-SST')
                allocate(instance, source=QUASI_NEWTON_SST)

            case ('quasi_newton_rs','Quasi_Newton_RS','quasi-newton-rs','Quasi-Newton-RS')
                allocate(instance, source=QUASI_NEWTON_RS)

            case default
                user_msg = "We can't seem to find a nonlinear solver that matches the input &
                            string in chidg.nml. Maybe check that the nonlinear solver string &
                            in the input file or driver script is valid."
                dev_msg  = "Check that the nonlinear solver is registered properly in &
                            create_nonlinear_solver."
                call chidg_signal_two(OOPS, user_msg, trim(string), dev_msg=dev_msg)

        end select


        ! Call options initialization if present
        if (present(options)) call instance%set(options)


        ! Make sure the solver was allocated
        user_msg = "create_nonlinear_solver: solver was not allocated. Check that the &
                    desired solver was registered and instantiated in the &
                    mod_nonlinear_solver module"
        if (.not. allocated(instance)) call chidg_signal(FATAL,user_msg)


    end subroutine create_nonlinear_solver
    !*****************************************************************************************








end module mod_nonlinear_solver
