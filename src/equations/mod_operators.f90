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
    use SD_bc_operator,                             only: SD_bc_operator_t

    ! Mesh Motion Diffusion Operators
    use MMD_volume_operator,                         only: MMD_volume_operator_t
    use MMD_boundary_operator,                       only: MMD_boundary_operator_t
    use MMD_bc_operator,                             only: MMD_bc_operator_t

    ! Mesh Motion Linear Elasticity Operators
    use MMLE_volume_operator,                         only: MMLE_volume_operator_t
    use MMLE_boundary_operator,                       only: MMLE_boundary_operator_t
    use MMLE_bc_operator,                             only: MMLE_bc_operator_t


    ! Fluid Inviscid Operators
    use euler_volume_operator,                      only: euler_volume_operator_t
    use euler_volume_cylindrical_source,            only: euler_volume_cylindrical_source_t
    use euler_boundary_average_operator,            only: euler_boundary_average_operator_t
    use euler_roe_operator,                         only: euler_roe_operator_t
    use euler_laxfriedrichs_operator,               only: euler_laxfriedrichs_operator_t
    use euler_bc_operator,                          only: euler_bc_operator_t

    ! Fluid Viscous Operators
    use fluid_viscous_volume_operator,              only: fluid_viscous_volume_operator_t
    use fluid_viscous_boundary_average_operator,    only: fluid_viscous_boundary_average_operator_t
    use fluid_viscous_bc_operator,                  only: fluid_viscous_bc_operator_t
    use fluid_viscous_volume_cylindrical_source,    only: fluid_viscous_volume_cylindrical_source_t

    ! Fluid Turbulence Operators
    use spalart_allmaras_source,                    only: spalart_allmaras_source_operator_t
    use spalart_allmaras_laxfriedrichs,             only: spalart_allmaras_laxfriedrichs_operator_t
    use spalart_allmaras_volume_advection,          only: spalart_allmaras_volume_advection_operator_t
    use spalart_allmaras_bc_advection,              only: spalart_allmaras_bc_advection_operator_t
    use spalart_allmaras_boundary_diffusion,        only: spalart_allmaras_boundary_diffusion_operator_t
    use spalart_allmaras_volume_diffusion,          only: spalart_allmaras_volume_diffusion_operator_t
    use spalart_allmaras_bc_diffusion,              only: spalart_allmaras_bc_diffusion_operator_t


    ! RANS Low-Cache operators
    use RANS_volume_advection,                      only: RANS_volume_advection_t
    use RANS_boundary_advection,                    only: RANS_boundary_advection_t
    use RANS_bc_advection,                          only: RANS_bc_advection_t
    use RANS_volume_diffusion,                      only: RANS_volume_diffusion_t
    use RANS_boundary_diffusion,                    only: RANS_boundary_diffusion_t
    use RANS_bc_diffusion,                          only: RANS_bc_diffusion_t
    use RANS_source,                                only: RANS_source_t

    ! Artificial Viscosity Operators
    use artificial_viscosity_boundary_average_operator, only: artificial_viscosity_boundary_average_operator_t
    use artificial_viscosity_volume_operator,           only: artificial_viscosity_volume_operator_t
    use artificial_viscosity_bc_operator,               only: artificial_viscosity_bc_operator_t
    use artificial_viscosity_source,                    only: artificial_viscosity_source_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public :: operator_factory_t

        type(ovector_t) :: operators

    contains

        procedure   :: register
        procedure   :: produce

    end type operator_factory_t
    !***************************************************************************************



    type(operator_factory_t)    :: operator_factory
    logical                     :: operators_initialized = .false.


contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/19/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine register(self,operator_instance)
        class(operator_factory_t),  intent(inout)   :: self
        class(operator_t),          intent(inout)   :: operator_instance

        ! Initialize the new operator
        call operator_instance%init()

        ! Add to the list of registered operators
        call self%operators%push_back(operator_instance)

    end subroutine register
    !***************************************************************************************



    !>  Build an operator based on an incoming string specification. Initialize the operator,
    !!  and return it to the caller.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function produce(self,string) result(op)
        class(operator_factory_t),  intent(inout)   :: self
        character(*),               intent(in)      :: string

        integer(ik)                     :: oindex, ierr
        character(:),       allocatable :: user_msg
        class(operator_t),  allocatable :: op

        !
        ! Find equation set in 'available_equations' vector
        !
        oindex = self%operators%index_by_name(string)


        !
        ! Check equationset was found in 'available_equations'
        !
        user_msg = "operator_factory%produce: We couldn't find the operator string in &
                    the list of registered operators. Make sure the operator was registered &
                    in the operator factory."
        if (oindex == 0) call chidg_signal_one(FATAL,user_msg,trim(string))


        !
        ! Get equation set builder
        !
        allocate(op, source=self%operators%at(oindex), stat=ierr)
        if (ierr /= 0) call AllocationError


        user_msg = "operator_factory%produce: For some reason, the operator didn't get allocated"
        if (.not. allocated(op)) call chidg_signal(FATAL,user_msg)

    end function produce
    !***************************************************************************************












    !>  Register new operators in the module vector, registered_operators.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
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
        
        ! Mesh Motion Diffusion Operators
        type(MMD_volume_operator_t)                      :: MMD_volume_operator
        type(MMD_boundary_operator_t)                    :: MMD_boundary_operator
        type(MMD_bc_operator_t)                          :: MMD_bc_operator

        ! Mesh Motion Linear Elasticity Operators
        type(MMLE_volume_operator_t)                      :: MMLE_volume_operator
        type(MMLE_boundary_operator_t)                    :: MMLE_boundary_operator
        type(MMLE_bc_operator_t)                          :: MMLE_bc_operator



        ! Dual Linear Advection Operators
        type(DLA_volume_advective_flux_t)               :: DLA_volume_operator
        type(DLA_boundary_average_advective_flux_t)     :: DLA_average_operator
        type(DLA_LaxFriedrichs_flux_t)                  :: DLA_laxfriedrichs_operator

        ! Fluid Inviscid Operators
        type(euler_volume_operator_t)                   :: euler_volume_operator
        type(euler_volume_cylindrical_source_t)         :: euler_volume_cylindrical_source
        type(euler_boundary_average_operator_t)         :: euler_average_operator
        type(euler_roe_operator_t)                      :: euler_roe_operator
        type(euler_laxfriedrichs_operator_t)            :: euler_laxfriedrichs_operator
        type(euler_bc_operator_t)                       :: euler_bc_operator

        ! Fluid Viscous Operators
        type(fluid_viscous_volume_operator_t)           :: fluid_viscous_volume_operator
        type(fluid_viscous_boundary_average_operator_t) :: fluid_viscous_boundary_average_operator
        type(fluid_viscous_bc_operator_t)               :: fluid_viscous_bc_operator
        type(fluid_viscous_volume_cylindrical_source_t) :: fluid_viscous_volume_cylindrical_source


        ! Fluid Turbulence Operators
        type(spalart_allmaras_source_operator_t)                :: spalart_allmaras_source_operator
        type(spalart_allmaras_laxfriedrichs_operator_t)         :: spalart_allmaras_laxfriedrichs_operator
        type(spalart_allmaras_volume_advection_operator_t)      :: spalart_allmaras_volume_advection_operator
        type(spalart_allmaras_bc_advection_operator_t)          :: spalart_allmaras_bc_advection_operator
        type(spalart_allmaras_boundary_diffusion_operator_t)    :: spalart_allmaras_boundary_diffusion_operator
        type(spalart_allmaras_volume_diffusion_operator_t)      :: spalart_allmaras_volume_diffusion_operator
        type(spalart_allmaras_bc_diffusion_operator_t)          :: spalart_allmaras_bc_diffusion_operator

        type(RANS_volume_advection_t)                           :: rans_volume_advection
        type(RANS_boundary_advection_t)                         :: rans_boundary_advection
        type(RANS_bc_advection_t)                               :: rans_bc_advection
        type(RANS_volume_diffusion_t)                           :: rans_volume_diffusion
        type(RANS_boundary_diffusion_t)                         :: rans_boundary_diffusion
        type(RANS_bc_diffusion_t)                               :: rans_bc_diffusion
        type(RANS_source_t)                                     :: rans_source


        ! Artificial Viscosity Operators
        type(artificial_viscosity_boundary_average_operator_t)  :: artificial_viscosity_boundary_average_operator
        type(artificial_viscosity_volume_operator_t)            :: artificial_viscosity_volume_operator
        type(artificial_viscosity_bc_operator_t)                :: artificial_viscosity_bc_operator
        type(artificial_viscosity_source_t)                     :: artificial_viscosity_source


        if (.not. operators_initialized) then

            ! Register Linear Advection
            call operator_factory%register(SA_volume_operator)
            call operator_factory%register(SA_average_operator)
            call operator_factory%register(SA_laxfriedrichs_operator)
            call operator_factory%register(SA_bc_operator)

            ! Register Linear Diffusion
            call operator_factory%register(SD_volume_operator)
            call operator_factory%register(SD_boundary_operator)
            call operator_factory%register(SD_bc_operator)

            ! Register Mesh Motion Diffusion
            call operator_factory%register(MMD_volume_operator)
            call operator_factory%register(MMD_boundary_operator)
            call operator_factory%register(MMD_bc_operator)

            ! Register Mesh Motion  Linear Elasticity
            call operator_factory%register(MMLE_volume_operator)
            call operator_factory%register(MMLE_boundary_operator)
            call operator_factory%register(MMLE_bc_operator)


            ! Register Dual Linear Advection
            call operator_factory%register(DLA_volume_operator)
            call operator_factory%register(DLA_average_operator)
            call operator_factory%register(DLA_laxfriedrichs_operator)


            ! Register Fluid Inviscid
            call operator_factory%register(euler_volume_operator)
            call operator_factory%register(euler_volume_cylindrical_source)
            call operator_factory%register(euler_average_operator)
            call operator_factory%register(euler_roe_operator)
            call operator_factory%register(euler_laxfriedrichs_operator)
            call operator_factory%register(euler_bc_operator)

            ! Register Fluid Viscous
            call operator_factory%register(fluid_viscous_boundary_average_operator)
            call operator_factory%register(fluid_viscous_bc_operator)
            call operator_factory%register(fluid_viscous_volume_operator)
            call operator_factory%register(fluid_viscous_volume_cylindrical_source)

            ! Register Fluid Turbulence
            call operator_factory%register(spalart_allmaras_source_operator)
            call operator_factory%register(spalart_allmaras_laxfriedrichs_operator)
            call operator_factory%register(spalart_allmaras_volume_advection_operator)
            call operator_factory%register(spalart_allmaras_bc_advection_operator)
            call operator_factory%register(spalart_allmaras_boundary_diffusion_operator)
            call operator_factory%register(spalart_allmaras_volume_diffusion_operator)
            call operator_factory%register(spalart_allmaras_bc_diffusion_operator)


            ! Register RANS Low-Cache operators
            call operator_factory%register(rans_volume_advection)
            call operator_factory%register(rans_volume_diffusion)
            call operator_factory%register(rans_boundary_advection)
            call operator_factory%register(rans_boundary_diffusion)
            call operator_factory%register(rans_bc_advection)
            call operator_factory%register(rans_bc_diffusion)
            call operator_factory%register(rans_source)



            ! Register Artificial Viscosity
            call operator_factory%register(artificial_viscosity_boundary_average_operator)
            call operator_factory%register(artificial_viscosity_volume_operator)
            call operator_factory%register(artificial_viscosity_bc_operator)
            call operator_factory%register(artificial_viscosity_source)


            operators_initialized = .true.

        end if


    end subroutine register_operators
    !**************************************************************************************










end module mod_operators
