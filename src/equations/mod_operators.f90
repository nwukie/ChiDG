module mod_operators
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_operator,  only: operator_t
    use type_ovector,   only: ovector_t

    use euler_volume_operator,              only: euler_volume_operator_t
    use euler_boundary_average_operator,    only: euler_boundary_average_operator_t
    use euler_roe_operator,                 only: euler_roe_operator_t
    use euler_laxfriedrichs_operator,       only: euler_laxfriedrichs_operator_t
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

        type(euler_volume_operator_t)           :: euler_volume_operator
        type(euler_boundary_average_operator_t) :: euler_average_operator
        type(euler_roe_operator_t)              :: euler_roe_operator
        type(euler_laxfriedrichs_operator_t)    :: euler_laxfriedrichs_operator







        if (.not. operators_initialized) then

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
