module type_time_integrator_marching
    use mod_kinds,              only: rk,ik
    use type_chidg_data,        only: chidg_data_t
    use type_time_integrator,   only: time_integrator_t
    implicit none


    !>  Abstraction for time integrators.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!  @date   2/7/2017
    !!
    !---------------------------------------------------------------------------------------
    type, abstract, extends(time_integrator_t), public  :: time_integrator_marching_t

    contains

        ! Define deferred state initialization for a 'Marching' time integrator.
        procedure   :: initialize_state

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

        associate( q => data%sdata%q, q_in => data%sdata%q_in)

            q = q_in

        end associate

    end subroutine initialize_state
    !*******************************************************************************


end module type_time_integrator_marching
