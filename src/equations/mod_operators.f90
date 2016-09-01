module mod_operators
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_operator,  only: operator_t
    use type_ovector,   only: ovector_t

    ! Linear Advection Equations
    use LA_volume_advective_flux,               only: LA_volume_advective_flux_t
    use LA_boundary_average_advective_flux,     only: LA_boundary_average_advective_flux_t
    use LA_LaxFriedrichs_flux,                  only: LA_LaxFriedrichs_flux_t

    ! Dual Linear Advection Equations
    use DLA_volume_advective_flux,              only: DLA_volume_advective_flux_t
    use DLA_boundary_average_advective_flux,    only: DLA_boundary_average_advective_flux_t
    use DLA_LaxFriedrichs_flux,                 only: DLA_LaxFriedrichs_flux_t

    ! Euler Equations
    use euler_volume_operator,                  only: euler_volume_operator_t
    use euler_boundary_average_operator,        only: euler_boundary_average_operator_t
    use euler_roe_operator,                     only: euler_roe_operator_t
    use euler_laxfriedrichs_operator,           only: euler_laxfriedrichs_operator_t
    implicit none



    type(ovector_t) :: registered_operators
    logical         :: operators_initialized = .false.




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine register_operators()
        integer(ik) :: iop

        ! Linear Advection Equations
        type(LA_volume_advective_flux_t)            :: LA_volume_operator
        type(LA_boundary_average_advective_flux_t)  :: LA_average_operator
        type(LA_LaxFriedrichs_flux_t)               :: LA_laxfriedrichs_operator
        
        ! Dual Linear Advection Equations
        type(DLA_volume_advective_flux_t)           :: DLA_volume_operator
        type(DLA_boundary_average_advective_flux_t) :: DLA_average_operator
        type(DLA_LaxFriedrichs_flux_t)              :: DLA_laxfriedrichs_operator

        ! Euler Equations
        type(euler_volume_operator_t)               :: euler_volume_operator
        type(euler_boundary_average_operator_t)     :: euler_average_operator
        type(euler_roe_operator_t)                  :: euler_roe_operator
        type(euler_laxfriedrichs_operator_t)        :: euler_laxfriedrichs_operator







        if (.not. operators_initialized) then

            ! Register Linear Advection
            call registered_operators%push_back(LA_volume_operator)
            call registered_operators%push_back(LA_average_operator)
            call registered_operators%push_back(LA_laxfriedrichs_operator)


            ! Register Dual Linear Advection
            call registered_operators%push_back(DLA_volume_operator)
            call registered_operators%push_back(DLA_average_operator)
            call registered_operators%push_back(DLA_laxfriedrichs_operator)


            ! Register Euler
            call registered_operators%push_back(euler_volume_operator)
            call registered_operators%push_back(euler_average_operator)
            call registered_operators%push_back(euler_roe_operator)
            call registered_operators%push_back(euler_laxfriedrichs_operator)


            ! Initialize all operators
            do iop = 1,registered_operators%size()
                call registered_operators%data(iop)%op%init()
            end do
            

            operators_initialized = .true.

        end if


    end subroutine register_operators
    !*****************************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function build_operator(string) result(op)
        character(len=*),   intent(in)  :: string

        integer(ik)                     :: oindex, ierr
        class(operator_t),  allocatable :: op

        !
        ! Find equation set in 'available_equations' vector
        !
        oindex = registered_operators%index_by_name(string)


        !
        ! Check equationset was found in 'available_equations'
        !
        if (oindex == 0) call chidg_signal_one(FATAL,"build_operator: We couldn't find the operator string in the list of registered operators.", trim(string))


        !
        ! Get equation set builder
        !
        allocate(op, source=registered_operators%at(oindex), stat=ierr)
        if (ierr /= 0) call AllocationError


        if (.not. allocated(op)) call chidg_signal(FATAL,"build_operator: For some reason, the operator didn't get allocated.")

    end function build_operator
    !*****************************************************************************************************













end module mod_operators
