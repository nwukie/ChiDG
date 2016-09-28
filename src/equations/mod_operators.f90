module mod_operators
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_operator,  only: operator_t
    use type_ovector,   only: ovector_t

    ! Linear Advection Operators
    use SA_volume_advective_operator,               only: SA_volume_advective_operator_t
    use SA_boundary_average_advective_operator,     only: SA_boundary_average_advective_operator_t
    use SA_LaxFriedrichs_operator,                  only: SA_LaxFriedrichs_operator_t
    use SA_bc_operator,                             only: SA_bc_operator_t

    ! Dual Linear Advection Operators
    use DLA_volume_advective_flux,                  only: DLA_volume_advective_flux_t
    use DLA_boundary_average_advective_flux,        only: DLA_boundary_average_advective_flux_t
    use DLA_LaxFriedrichs_flux,                     only: DLA_LaxFriedrichs_flux_t

    ! Scalar Diffusion Operators
    use SD_volume_operator,                         only: SD_volume_operator_t
    use SD_boundary_operator,                       only: SD_boundary_operator_t
    use SD_volume_source,                           only: SD_volume_source_t
    use SD_bc_operator,                             only: SD_bc_operator_t

    ! Fluid Inviscid Operators
    use euler_volume_operator,                      only: euler_volume_operator_t
    use euler_boundary_average_operator,            only: euler_boundary_average_operator_t
    use euler_roe_operator,                         only: euler_roe_operator_t
    use euler_laxfriedrichs_operator,               only: euler_laxfriedrichs_operator_t
    use euler_bc_operator,                          only: euler_bc_operator_t

    ! Fluid Viscous Operators
    use fluid_viscous_volume_operator,              only: fluid_viscous_volume_operator_t
    use fluid_viscous_boundary_average_operator,    only: fluid_viscous_boundary_average_operator_t
    use fluid_viscous_bc_operator,                  only: fluid_viscous_bc_operator_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    type, public :: operator_factory_t

        type(ovector_t) :: operators

    contains

        procedure   :: register
        procedure   :: produce

    end type operator_factory_t
    !**************************************************************************************************



    type(operator_factory_t)    :: operator_factory
    logical                     :: operators_initialized = .false.


contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine register(self,operator_instance)
        class(operator_factory_t),  intent(inout)   :: self
        class(operator_t),          intent(inout)   :: operator_instance

        call self%operators%push_back(operator_instance)

    end subroutine register
    !**************************************************************************************************



    !>  Build an operator based on an incoming string specification. Initialize the operator,
    !!  and return it to the caller.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function produce(self,string) result(op)
        class(operator_factory_t),  intent(inout)   :: self
        character(len=*),           intent(in)      :: string

        integer(ik)                     :: oindex, ierr
        class(operator_t),  allocatable :: op

        !
        ! Find equation set in 'available_equations' vector
        !
        oindex = self%operators%index_by_name(string)


        !
        ! Check equationset was found in 'available_equations'
        !
        if (oindex == 0) call chidg_signal_one(FATAL,"build_operator: We couldn't find the operator string in the list of registered operators.", trim(string))


        !
        ! Get equation set builder
        !
        allocate(op, source=self%operators%at(oindex), stat=ierr)
        if (ierr /= 0) call AllocationError


        if (.not. allocated(op)) call chidg_signal(FATAL,"build_operator: For some reason, the operator didn't get allocated.")

    end function produce
    !*****************************************************************************************************












    !>  Register new operators in the module vector, registered_operators.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine register_operators()
        integer(ik) :: iop

        ! Linear Advection Operators
        type(SA_volume_advective_operator_t)            :: SA_volume_operator
        type(SA_boundary_average_advective_operator_t)  :: SA_average_operator
        type(SA_LaxFriedrichs_operator_t)               :: SA_laxfriedrichs_operator
        type(SA_bc_operator_t)                          :: SA_bc_operator
        
        ! Linear Diffusion Operators
        type(SD_volume_operator_t)                      :: SD_volume_operator
        type(SD_boundary_operator_t)                    :: SD_boundary_operator
        type(SD_bc_operator_t)                          :: SD_bc_operator
        type(SD_volume_source_t)                        :: SD_volume_source

        ! Dual Linear Advection Operators
        type(DLA_volume_advective_flux_t)               :: DLA_volume_operator
        type(DLA_boundary_average_advective_flux_t)     :: DLA_average_operator
        type(DLA_LaxFriedrichs_flux_t)                  :: DLA_laxfriedrichs_operator

        ! Fluid Inviscid Operators
        type(euler_volume_operator_t)                   :: euler_volume_operator
        type(euler_boundary_average_operator_t)         :: euler_average_operator
        type(euler_roe_operator_t)                      :: euler_roe_operator
        type(euler_laxfriedrichs_operator_t)            :: euler_laxfriedrichs_operator
        type(euler_bc_operator_t)                       :: euler_bc_operator

        ! Fluid Viscous Operators
        type(fluid_viscous_volume_operator_t)           :: fluid_viscous_volume_operator
        type(fluid_viscous_boundary_average_operator_t) :: fluid_viscous_boundary_average_operator
        type(fluid_viscous_bc_operator_t)               :: fluid_viscous_bc_operator






        if (.not. operators_initialized) then

            ! Register Linear Advection
            call operator_factory%register(SA_volume_operator)
            call operator_factory%register(SA_average_operator)
            call operator_factory%register(SA_laxfriedrichs_operator)
            call operator_factory%register(SA_bc_operator)

            ! Register Linear Diffusion
            call operator_factory%register(SD_volume_operator)
            call operator_factory%register(SD_boundary_operator)
            call operator_factory%register(SD_volume_source)
            call operator_factory%register(SD_bc_operator)


            ! Register Dual Linear Advection
            call operator_factory%register(DLA_volume_operator)
            call operator_factory%register(DLA_average_operator)
            call operator_factory%register(DLA_laxfriedrichs_operator)


            ! Register Fluid Inviscid
            call operator_factory%register(euler_volume_operator)
            call operator_factory%register(euler_average_operator)
            call operator_factory%register(euler_roe_operator)
            call operator_factory%register(euler_laxfriedrichs_operator)
            call operator_factory%register(euler_bc_operator)

            ! Register Fluid Viscous
            call operator_factory%register(fluid_viscous_volume_operator)
            call operator_factory%register(fluid_viscous_boundary_average_operator)
            call operator_factory%register(fluid_viscous_bc_operator)


            ! Initialize all operators
            do iop = 1,operator_factory%operators%size()
                call operator_factory%operators%data(iop)%op%init()
            end do
            

            operators_initialized = .true.

        end if


    end subroutine register_operators
    !*****************************************************************************************************






















end module mod_operators
