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
    use DLA_bc_operator,                            only: DLA_bc_operator_t

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
    use euler_boundary_average_operator,            only: euler_boundary_average_operator_t
    use euler_roe_operator,                         only: euler_roe_operator_t
    use euler_bc_operator,                          only: euler_bc_operator_t
    use euler_volume_cylindrical_source,            only: euler_volume_cylindrical_source_t
    !use euler_laxfriedrichs_operator,               only: euler_laxfriedrichs_operator_t

    ! Fluid Viscous Operators
    use fluid_viscous_volume_operator,              only: fluid_viscous_volume_operator_t
    use fluid_viscous_boundary_average_operator,    only: fluid_viscous_boundary_average_operator_t
    use fluid_viscous_bc_operator,                  only: fluid_viscous_bc_operator_t
    use fluid_viscous_volume_cylindrical_source,    only: fluid_viscous_volume_cylindrical_source_t

    ! Fluid Laplacian AV Operators
    use fluid_laplacian_av_volume_operator,              only: fluid_laplacian_av_volume_operator_t
    use fluid_laplacian_av_boundary_average_operator,    only: fluid_laplacian_av_boundary_average_operator_t
    use fluid_laplacian_av_bc_operator,                  only: fluid_laplacian_av_bc_operator_t

    ! Fluid Laplacian Anisotropic AV Operators
    use fluid_laplacian_anisotropic_av_volume_operator,              only: fluid_laplacian_anisotropic_av_volume_operator_t
    use fluid_laplacian_anisotropic_av_boundary_average_operator,    only: fluid_laplacian_anisotropic_av_boundary_average_operator_t
    use fluid_laplacian_anisotropic_av_bc_operator,                  only: fluid_laplacian_anisotropic_av_bc_operator_t





    ! Fluid Turbulence Operators
    ! Spalart-Allmaras
    use spalart_allmaras_source,                            only: spalart_allmaras_source_operator_t
    use spalart_allmaras_advection_boundary_average,        only: spalart_allmaras_advection_boundary_average_operator_t
    use spalart_allmaras_laxfriedrichs,                     only: spalart_allmaras_laxfriedrichs_operator_t
    use spalart_allmaras_volume_advection,                  only: spalart_allmaras_volume_advection_operator_t
    use spalart_allmaras_bc_advection,                      only: spalart_allmaras_bc_advection_operator_t
    use spalart_allmaras_boundary_diffusion,                only: spalart_allmaras_boundary_diffusion_operator_t
    use spalart_allmaras_volume_diffusion,                  only: spalart_allmaras_volume_diffusion_operator_t
    use spalart_allmaras_bc_diffusion,                      only: spalart_allmaras_bc_diffusion_operator_t

    !! SST
    use sst_source,                                         only: sst_source_operator_t
    use sst_roe_operator,                                   only: sst_roe_operator_t
    use sst_advection_boundary_average,                     only: sst_advection_boundary_average_operator_t
    use sst_laxfriedrichs,                                  only: sst_laxfriedrichs_operator_t
    use sst_volume_advection,                               only: sst_volume_advection_operator_t
    use sst_bc_advection,                                   only: sst_bc_advection_operator_t
    use sst_boundary_diffusion,                             only: sst_boundary_diffusion_operator_t
    use sst_volume_diffusion,                               only: sst_volume_diffusion_operator_t
    use sst_bc_diffusion,                                   only: sst_bc_diffusion_operator_t
    use sst_artificial_viscosity_operator,                  only: sst_artificial_viscosity_operator_t
    use sst_artificial_viscosity_bc_operator,               only: sst_artificial_viscosity_bc_operator_t
    use sst_artificial_viscosity_boundary_average_operator, only: sst_artificial_viscosity_boundary_average_operator_t

    !! RSTM SSG-LRR-w operators
    use rstm_ssglrrw_source,                        only: rstm_ssglrrw_source_operator_t
    use rstm_ssglrrw_advection_boundary_average,    only: rstm_ssglrrw_advection_boundary_average_operator_t
    use rstm_ssglrrw_laxfriedrichs,                 only: rstm_ssglrrw_laxfriedrichs_operator_t
    use rstm_ssglrrw_volume_advection,              only: rstm_ssglrrw_volume_advection_operator_t
    use rstm_ssglrrw_bc_advection,                  only: rstm_ssglrrw_bc_advection_operator_t
    use rstm_ssglrrw_boundary_diffusion,            only: rstm_ssglrrw_boundary_diffusion_operator_t
    use rstm_ssglrrw_volume_diffusion,              only: rstm_ssglrrw_volume_diffusion_operator_t
    use rstm_ssglrrw_bc_diffusion,                  only: rstm_ssglrrw_bc_diffusion_operator_t
    use rstm_ssglrrw_artificial_viscosity_operator,                          only: rstm_ssglrrw_artificial_viscosity_operator_t
    use rstm_ssglrrw_artificial_viscosity_bc_operator,                       only: rstm_ssglrrw_artificial_viscosity_bc_operator_t
    use rstm_ssglrrw_artificial_viscosity_boundary_average_operator,         only: rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t


    ! RANS Efficient operators
    use rans_bc_advection,                          only: rans_bc_advection_t
    use rans_bc_diffusion,                          only: rans_bc_diffusion_t
    use rans_boundary_average_advection,            only: rans_boundary_average_advection_t
    use rans_boundary_average_diffusion,            only: rans_boundary_average_diffusion_t
    use rans_volume_advection,                      only: rans_volume_advection_t
    use rans_volume_diffusion,                      only: rans_volume_diffusion_t
    use rans_volume_source,                         only: rans_volume_source_t
    use rans_upwind_operator,                       only: rans_upwind_operator_t

    ! Geometric Conservation Operators
    use GCL_volume_advective_operator,               only: GCL_volume_advective_operator_t
    use GCL_boundary_average_advective_operator,     only: GCL_boundary_average_advective_operator_t
    use GCL_bc_operator,                             only: GCL_bc_operator_t


    ! Primitive Variable Time-Linearized Euler
    use PRIMLINEULER_LaxFriedrichs,             only: PRIMLINEULER_LaxFriedrichs_t
    use PRIMLINEULER_boundary_average,          only: PRIMLINEULER_boundary_average_t
    use PRIMLINEULER_bc,                        only: PRIMLINEULER_bc_t
    use PRIMLINEULER_volume_advection,          only: PRIMLINEULER_volume_advection_t
    use PRIMLINEULER_temporal_source,           only: PRIMLINEULER_temporal_source_t
    use PRIMLINEULER_axial_source,              only: PRIMLINEULER_axial_source_t
    use PRIMLINEULER_circumferential_source,    only: PRIMLINEULER_circumferential_source_t
    use PRIMLINEULER_equation_source,           only: PRIMLINEULER_equation_source_t
    use PRIMLINEULER_divergence_source,         only: PRIMLINEULER_divergence_source_t


    ! Hyperbolized Poisson
    use HP_LaxFriedrichs,                       only: HP_LaxFriedrichs_t
    use HP_boundary_average,                    only: HP_boundary_average_t
    use HP_volume,                              only: HP_volume_t
    use HP_bc,                                  only: HP_bc_t
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

        ! Find equation set in 'available_equations' vector
        oindex = self%operators%index_by_name(string)

        ! Check equationset was found in 'available_equations'
        user_msg = "operator_factory%produce: We couldn't find the operator string in &
                    the list of registered operators. Make sure the operator was registered &
                    in the operator factory."
        if (oindex == 0) call chidg_signal_one(FATAL,user_msg,trim(string))

        ! Allocate operator
        allocate(op, source=self%operators%at(oindex), stat=ierr)
        if (ierr /= 0) call AllocationError

        user_msg = "operator_factory%produce: For some reason, the operator didn't get allocated"
        if (.not. allocated(op)) call chidg_signal(FATAL,user_msg)

    end function produce
    !***************************************************************************************









    !>  Register new operators in the module vector, registered_operators.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine register_operators()
        integer(ik) :: iop

        ! Linear Advection Operators
        type(SA_volume_advective_operator_t)                :: SA_volume_operator
        type(SA_boundary_average_advective_operator_t)      :: SA_average_operator
        type(SA_LaxFriedrichs_operator_t)                   :: SA_laxfriedrichs_operator
        type(SA_bc_operator_t)                              :: SA_bc_operator

        ! Linear Diffusion Operators
        type(SD_volume_operator_t)                          :: SD_volume_operator
        type(SD_boundary_operator_t)                        :: SD_boundary_operator
        type(SD_bc_operator_t)                              :: SD_bc_operator

        ! Mesh Motion Diffusion Operators
        type(MMD_volume_operator_t)                         :: MMD_volume_operator
        type(MMD_boundary_operator_t)                       :: MMD_boundary_operator
        type(MMD_bc_operator_t)                             :: MMD_bc_operator

        ! Mesh Motion Linear Elasticity Operators
        type(MMLE_volume_operator_t)                        :: MMLE_volume_operator
        type(MMLE_boundary_operator_t)                      :: MMLE_boundary_operator
        type(MMLE_bc_operator_t)                            :: MMLE_bc_operator

        ! Dual Linear Advection Operators
        type(DLA_volume_advective_flux_t)                   :: DLA_volume_operator
        type(DLA_boundary_average_advective_flux_t)         :: DLA_average_operator
        type(DLA_LaxFriedrichs_flux_t)                      :: DLA_laxfriedrichs_operator
        type(DLA_bc_operator_t)                             :: DLA_bc_operator

        ! Fluid Inviscid Operators
        type(euler_volume_operator_t)                       :: euler_volume_operator
        type(euler_boundary_average_operator_t)             :: euler_average_operator
        type(euler_roe_operator_t)                          :: euler_roe_operator
        type(euler_bc_operator_t)                           :: euler_bc_operator
        type(euler_volume_cylindrical_source_t)             :: euler_volume_cylindrical_source
        !type(euler_laxfriedrichs_operator_t)                :: euler_laxfriedrichs_operator

        ! Fluid Viscous Operators
        type(fluid_viscous_volume_operator_t)               :: fluid_viscous_volume_operator
        type(fluid_viscous_boundary_average_operator_t)     :: fluid_viscous_boundary_average_operator
        type(fluid_viscous_bc_operator_t)                   :: fluid_viscous_bc_operator
        type(fluid_viscous_volume_cylindrical_source_t)     :: fluid_viscous_volume_cylindrical_source

        ! Fluid Laplacian AV Operators
        type(fluid_laplacian_av_volume_operator_t)               :: fluid_laplacian_av_volume_operator
        type(fluid_laplacian_av_boundary_average_operator_t)     :: fluid_laplacian_av_boundary_average_operator
        type(fluid_laplacian_av_bc_operator_t)                   :: fluid_laplacian_av_bc_operator
        
        ! Fluid Laplacian Anisotropic AV Operators
        type(fluid_laplacian_anisotropic_av_volume_operator_t)               :: fluid_laplacian_anisotropic_av_volume_operator
        type(fluid_laplacian_anisotropic_av_boundary_average_operator_t)     :: fluid_laplacian_anisotropic_av_boundary_average_operator
        type(fluid_laplacian_anisotropic_av_bc_operator_t)                   :: fluid_laplacian_anisotropic_av_bc_operator




        ! Fluid Turbulence Operators
        !! Spalart-Allmaras
        type(spalart_allmaras_source_operator_t)                        :: spalart_allmaras_source_operator
        type(spalart_allmaras_advection_boundary_average_operator_t)    :: spalart_allmaras_advection_boundary_average_operator
        type(spalart_allmaras_laxfriedrichs_operator_t)                 :: spalart_allmaras_laxfriedrichs_operator
        type(spalart_allmaras_volume_advection_operator_t)              :: spalart_allmaras_volume_advection_operator
        type(spalart_allmaras_bc_advection_operator_t)                  :: spalart_allmaras_bc_advection_operator
        type(spalart_allmaras_boundary_diffusion_operator_t)            :: spalart_allmaras_boundary_diffusion_operator
        type(spalart_allmaras_volume_diffusion_operator_t)              :: spalart_allmaras_volume_diffusion_operator
        type(spalart_allmaras_bc_diffusion_operator_t)                  :: spalart_allmaras_bc_diffusion_operator


        !! SST
        type(sst_source_operator_t)                        :: sst_source_operator
        type(sst_advection_boundary_average_operator_t)    :: sst_advection_boundary_average_operator
        type(sst_roe_operator_t)                           :: sst_roe_operator
        type(sst_laxfriedrichs_operator_t)                 :: sst_laxfriedrichs_operator
        type(sst_volume_advection_operator_t)              :: sst_volume_advection_operator
        type(sst_bc_advection_operator_t)                  :: sst_bc_advection_operator
        type(sst_boundary_diffusion_operator_t)            :: sst_boundary_diffusion_operator
        type(sst_volume_diffusion_operator_t)              :: sst_volume_diffusion_operator
        type(sst_bc_diffusion_operator_t)                  :: sst_bc_diffusion_operator
        type(sst_artificial_viscosity_operator_t)          :: sst_artificial_viscosity_operator
        type(sst_artificial_viscosity_bc_operator_t)          :: sst_artificial_viscosity_bc_operator
        type(sst_artificial_viscosity_boundary_average_operator_t)          :: sst_artificial_viscosity_boundary_average_operator


        !! Reynolds-Stress 
        type(rstm_ssglrrw_source_operator_t)                        :: rstm_ssglrrw_source_operator
        type(rstm_ssglrrw_advection_boundary_average_operator_t)    :: rstm_ssglrrw_advection_boundary_average_operator
        type(rstm_ssglrrw_laxfriedrichs_operator_t)                 :: rstm_ssglrrw_laxfriedrichs_operator
        type(rstm_ssglrrw_volume_advection_operator_t)              :: rstm_ssglrrw_volume_advection_operator
        type(rstm_ssglrrw_bc_advection_operator_t)                  :: rstm_ssglrrw_bc_advection_operator
        type(rstm_ssglrrw_boundary_diffusion_operator_t)            :: rstm_ssglrrw_boundary_diffusion_operator
        type(rstm_ssglrrw_volume_diffusion_operator_t)              :: rstm_ssglrrw_volume_diffusion_operator
        type(rstm_ssglrrw_bc_diffusion_operator_t)                  :: rstm_ssglrrw_bc_diffusion_operator
        type(rstm_ssglrrw_artificial_viscosity_operator_t)          :: rstm_ssglrrw_artificial_viscosity_operator
        type(rstm_ssglrrw_artificial_viscosity_bc_operator_t)          :: rstm_ssglrrw_artificial_viscosity_bc_operator
        type(rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t)          :: rstm_ssglrrw_artificial_viscosity_boundary_average_operator



        ! RANS Efficient
        type(rans_bc_advection_t)                       :: RANS_BC_ADVECTION
        type(rans_bc_diffusion_t)                       :: RANS_BC_DIFFUSION
        type(rans_boundary_average_advection_t)         :: RANS_BOUNDARY_AVERAGE_ADVECTION
        type(rans_boundary_average_diffusion_t)         :: RANS_BOUNDARY_AVERAGE_DIFFUSION
        type(rans_volume_advection_t)                   :: RANS_VOLUME_ADVECTION
        type(rans_volume_diffusion_t)                   :: RANS_VOLUME_DIFFUSION
        type(rans_volume_source_t)                      :: RANS_VOLUME_SOURCE
        type(rans_upwind_operator_t)                    :: RANS_UPWIND_OPERATOR


        ! Geometric Conservation Operators
        type(GCL_volume_advective_operator_t)                   :: GCL_volume_operator
        type(GCL_boundary_average_advective_operator_t)         :: GCL_average_operator
        type(GCL_bc_operator_t)                                 :: GCL_bc_operator


        ! Linearized Euler Eigen
        type(PRIMLINEULER_boundary_average_t)           :: PRIMLINEULER_boundary_average
        type(PRIMLINEULER_volume_advection_t)           :: PRIMLINEULER_volume_advection
        type(PRIMLINEULER_bc_t)                         :: PRIMLINEULER_bc
        type(PRIMLINEULER_LaxFriedrichs_t)              :: PRIMLINEULER_LaxFriedrichs
        type(PRIMLINEULER_temporal_source_t)            :: PRIMLINEULER_temporal_source
        type(PRIMLINEULER_axial_source_t)               :: PRIMLINEULER_axial_source
        type(PRIMLINEULER_circumferential_source_t)     :: PRIMLINEULER_circumferential_source
        type(PRIMLINEULER_equation_source_t)            :: PRIMLINEULER_equation_source
        type(PRIMLINEULER_divergence_source_t)          :: PRIMLINEULER_divergence_source

        ! Hyperbolized Poisson
        type(HP_LaxFriedrichs_t)                        :: HP_LaxFriedrichs_operator
        type(HP_volume_t)                               :: HP_volume_operator
        type(HP_bc_t)                                   :: HP_bc_operator
        type(HP_boundary_average_t)                     :: HP_boundary_average_operator


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
            call operator_factory%register(DLA_bc_operator)


            ! Register Fluid Inviscid
            call operator_factory%register(euler_volume_operator)
            call operator_factory%register(euler_average_operator)
            call operator_factory%register(euler_roe_operator)
            call operator_factory%register(euler_bc_operator)
            call operator_factory%register(euler_volume_cylindrical_source)
            !call operator_factory%register(euler_laxfriedrichs_operator)
            

            ! Register Fluid Viscous
            call operator_factory%register(fluid_viscous_boundary_average_operator)
            call operator_factory%register(fluid_viscous_bc_operator)
            call operator_factory%register(fluid_viscous_volume_operator)
            call operator_factory%register(fluid_viscous_volume_cylindrical_source)

            ! Register Fluid Laplacian AV
            call operator_factory%register(fluid_laplacian_av_boundary_average_operator)
            call operator_factory%register(fluid_laplacian_av_bc_operator)
            call operator_factory%register(fluid_laplacian_av_volume_operator)

            ! Register Fluid Laplacian Anisotropic AV
            call operator_factory%register(fluid_laplacian_anisotropic_av_boundary_average_operator)
            call operator_factory%register(fluid_laplacian_anisotropic_av_bc_operator)
            call operator_factory%register(fluid_laplacian_anisotropic_av_volume_operator)


            ! Register Fluid Turbulence
            ! Spalart Allmaras
            call operator_factory%register(spalart_allmaras_source_operator)
            call operator_factory%register(spalart_allmaras_advection_boundary_average_operator)
            call operator_factory%register(spalart_allmaras_laxfriedrichs_operator)
            call operator_factory%register(spalart_allmaras_volume_advection_operator)
            call operator_factory%register(spalart_allmaras_bc_advection_operator)
            call operator_factory%register(spalart_allmaras_boundary_diffusion_operator)
            call operator_factory%register(spalart_allmaras_volume_diffusion_operator)
            call operator_factory%register(spalart_allmaras_bc_diffusion_operator)

            !! SST
            call operator_factory%register(sst_source_operator)
            call operator_factory%register(sst_roe_operator)
            call operator_factory%register(sst_advection_boundary_average_operator)
            call operator_factory%register(sst_laxfriedrichs_operator)
            call operator_factory%register(sst_volume_advection_operator)
            call operator_factory%register(sst_bc_advection_operator)
            call operator_factory%register(sst_boundary_diffusion_operator)
            call operator_factory%register(sst_volume_diffusion_operator)
            call operator_factory%register(sst_bc_diffusion_operator)
            call operator_factory%register(sst_artificial_viscosity_operator)
            call operator_factory%register(sst_artificial_viscosity_bc_operator)
            call operator_factory%register(sst_artificial_viscosity_boundary_average_operator)


            !! Reynolds-Stress
            call operator_factory%register(rstm_ssglrrw_source_operator)
            call operator_factory%register(rstm_ssglrrw_advection_boundary_average_operator)
            call operator_factory%register(rstm_ssglrrw_laxfriedrichs_operator)
            call operator_factory%register(rstm_ssglrrw_volume_advection_operator)
            call operator_factory%register(rstm_ssglrrw_bc_advection_operator)
            call operator_factory%register(rstm_ssglrrw_boundary_diffusion_operator)
            call operator_factory%register(rstm_ssglrrw_volume_diffusion_operator)
            call operator_factory%register(rstm_ssglrrw_bc_diffusion_operator)
            call operator_factory%register(rstm_ssglrrw_artificial_viscosity_operator)
            call operator_factory%register(rstm_ssglrrw_artificial_viscosity_bc_operator)
            call operator_factory%register(rstm_ssglrrw_artificial_viscosity_boundary_average_operator)


            ! Register RANS Efficient operators
            call operator_factory%register(RANS_BC_ADVECTION)
            call operator_factory%register(RANS_BC_DIFFUSION)
            call operator_factory%register(RANS_BOUNDARY_AVERAGE_ADVECTION)
            call operator_factory%register(RANS_BOUNDARY_AVERAGE_DIFFUSION)
            call operator_factory%register(RANS_VOLUME_ADVECTION)
            call operator_factory%register(RANS_VOLUME_DIFFUSION)
            call operator_factory%register(RANS_VOLUME_SOURCE)
            call operator_factory%register(RANS_UPWIND_OPERATOR)

            ! Register Geometric Conservation
            call operator_factory%register(GCL_volume_operator)
            call operator_factory%register(GCL_average_operator)
            call operator_factory%register(GCL_bc_operator)

            ! Linearized Euler Eigen
            call operator_factory%register(PRIMLINEULER_boundary_average)
            call operator_factory%register(PRIMLINEULER_volume_advection)
            call operator_factory%register(PRIMLINEULER_bc)
            call operator_factory%register(PRIMLINEULER_LaxFriedrichs)
            call operator_factory%register(PRIMLINEULER_temporal_source)
            call operator_factory%register(PRIMLINEULER_axial_source)
            call operator_factory%register(PRIMLINEULER_circumferential_source)
            call operator_factory%register(PRIMLINEULER_equation_source)
            call operator_factory%register(PRIMLINEULER_divergence_source)

            ! Hyperbolized Poisson
            call operator_factory%register(HP_LaxFriedrichs_operator)
            call operator_factory%register(HP_volume_operator)
            call operator_factory%register(HP_boundary_average_operator)
            call operator_factory%register(HP_bc_operator)



            operators_initialized = .true.

        end if


    end subroutine register_operators
    !**************************************************************************************










end module mod_operators
