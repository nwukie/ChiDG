module mod_models
#include <messenger.h>
    use mod_kinds,              only: ik
    use mod_string,             only: string_to_upper
    use type_model_wrapper,     only: model_wrapper_t
    use type_model,             only: model_t

    use model_ideal_gas,                                only: ideal_gas_t
    use model_ideal_gas_sst,                            only: ideal_gas_sst_t
    use model_ideal_gas_rstm,                           only: ideal_gas_rstm_t

    use model_reynolds_shear_stress,                    only: reynolds_shear_stress_t
    use model_shear_stress,                             only: shear_stress_t
    use model_vorticity,                                only: vorticity_t
    use model_temperature_gradient,                     only: temperature_gradient_t
    use model_pressure_gradient,                        only: pressure_gradient_t
    use model_velocity_gradients,                       only: velocity_gradients_t
    use model_sutherlands_law,                          only: sutherlands_law_t
    use model_constant_viscosity,                       only: constant_viscosity_t
    use model_stokes_hypothesis,                        only: stokes_hypothesis_t
    use model_reynolds_analogy,                         only: reynolds_analogy_t
    use model_zero_turbulent_model_fields,              only: zero_turbulent_model_fields_t
    use model_zero_reynolds_stress,                     only: zero_reynolds_stress_t
    use model_critical_sound_speed,                     only: critical_sound_speed_t

    use model_spalart_allmaras_turbulent_model_fields,  only: spalart_allmaras_turbulent_model_fields_t
    use model_fluid_wave_speed,                         only: fluid_wave_speed_t
    use model_wall_distance,                            only: wall_distance_m


    use model_strain_rate,                              only: strain_rate_t
    use model_rotation_rate,                            only: rotation_rate_t
    use model_velocity_gradient,                        only: velocity_gradient_t
    use model_velocity_div_curl,                        only: velocity_div_curl_t

    use model_modified_ducros_sensor,                   only: modified_ducros_sensor_t
    use model_artificial_viscosity,                     only: artificial_viscosity_t
    use model_zero_artificial_viscosity,                only: zero_artificial_viscosity_t
    use model_unsmoothed_artificial_viscosity,          only: unsmoothed_artificial_viscosity_t
    use model_rbf_smoothed_artificial_viscosity,        only: rbf_smoothed_artificial_viscosity_t
    use model_vertex_smoothed_artificial_viscosity,     only: vertex_smoothed_artificial_viscosity_t
    use model_pde_smoothed_artificial_viscosity,        only: pde_smoothed_artificial_viscosity_t

    use model_mnp_shock_sensor,                         only: mnp_shock_sensor_t
    use model_mnp_artificial_viscosity,                 only: mnp_artificial_viscosity_t
    use model_vertex_smoothed_mnp_artificial_viscosity, only: vertex_smoothed_mnp_artificial_viscosity_t

    use model_mnph_shock_sensor,                        only: mnph_shock_sensor_t
    use model_mnph_artificial_viscosity,                only: mnph_artificial_viscosity_t
    
    use model_mnpha_artificial_viscosity,               only: mnpha_artificial_viscosity_t

    use model_sst_turbulence_kinetic_energy,            only: sst_turbulence_kinetic_energy_t
    use model_sst_turbulence_quantities,                only: sst_turbulence_quantities_t
    use model_sst_blended_coefficients,                 only: sst_blended_coefficients_t
    use model_sst_source_terms,                         only: sst_source_terms_t

    use model_rstm_ssglrrw_blended_coefficients,        only: rstm_ssglrrw_blended_coefficients_t
    use model_rstm_ssglrrw_generalized_diffusion,       only: rstm_ssglrrw_generalized_diffusion_t
    use model_rstm_ssglrrw_simple_diffusion,            only: rstm_ssglrrw_simple_diffusion_t
    use model_rstm_ssglrrw_isotropic_dissipation,       only: rstm_ssglrrw_isotropic_dissipation_t
    use model_rstm_ssglrrw_pressure_strain_correlation, only: rstm_ssglrrw_pressure_strain_correlation_t
    use model_rstm_ssglrrw_production,                  only: rstm_ssglrrw_production_t
!    use model_rstm_ssglrrw_reynolds_stress,             only: rstm_ssglrrw_reynolds_stress_t
    use model_rstm_ssglrrw_turbulence_quantities,       only: rstm_ssglrrw_turbulence_quantities_t
!    use model_rstm_ssglrrw_turbulence_kinetic_energy,       only: rstm_ssglrrw_turbulence_kinetic_energy_t
    implicit none




    !>  A model_t factory for initializing, holding, and producing model's on-demand.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    type, private :: model_factory_t

        type(model_wrapper_t),  allocatable :: models(:)

    contains

        procedure   :: register
        procedure   :: produce
        procedure   :: nmodels

        procedure, private  :: index_by_name

    end type model_factory_t
    !**************************************************************************************



    type(model_factory_t)   :: model_factory
    logical                 :: models_initialized = .false.



contains




    !>  Register a new model in the model factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine register(self,model)
        class(model_factory_t), intent(inout)   :: self
        class(model_t),         intent(inout)   :: model

        integer(ik)                         :: ierr, imodel
        type(model_wrapper_t),  allocatable :: temp(:)

        ! Initialize the incoming model
        call model%init()

        ! Extend storage
        if (allocated(self%models)) then

            allocate(temp(size(self%models) + 1), stat=ierr)
            if (ierr /= 0) call AllocationError

            do imodel = 1,size(self%models)
                allocate(temp(imodel)%model, source=self%models(imodel)%model, stat=ierr)
                if (ierr /= 0) call AllocationError
            end do

        else
            allocate(temp(1), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if

        ! Store new model
        allocate(temp(size(temp))%model, source=model, stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Move allocation
        call move_alloc(from=temp, to=self%models)

    end subroutine register
    !**************************************************************************************





    !>  Given a string, find a model with a name matching the string
    !!  and return a model allocated with that concrete type.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !--------------------------------------------------------------------------------------
    function produce(self,string) result(model)
        class(model_factory_t), intent(in)  :: self
        character(*),           intent(in)  :: string

        character(:),       allocatable :: user_msg
        class(model_t),     allocatable :: model
        integer(ik)                     :: imodel, ierr

        ! Find location of model
        imodel = self%index_by_name(string)

        ! Check model was found
        user_msg = "model_factory%produce: We couldn't find the model string in &
                    the list of registered models. Make sure the model was registered &
                    in the model factory."
        if (imodel == 0) call chidg_signal_one(FATAL,user_msg,trim(string))

        ! Allocate model to be returned
        allocate(model, source=self%models(imodel)%model, stat=ierr)
        if (ierr /= 0) call AllocationError

        user_msg = "model_factory%produce: For some reason, the model didn't get allocated"
        if (.not. allocated(model)) call chidg_signal(FATAL,user_msg)

    end function produce
    !**************************************************************************************





    !>  Return number of models registered in the model factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !-------------------------------------------------------------------------------------
    function nmodels(self) result(nmodels_)
        class(model_factory_t), intent(in)  :: self

        integer(ik) :: nmodels_

        if (allocated(self%models)) then
            nmodels_ = size(self%models)
        else
            nmodels_ = 0
        end if

    end function nmodels
    !*************************************************************************************





    !>  Given a string indicating a model, return the index of the model in the factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !------------------------------------------------------------------------------------
    function index_by_name(self,string) result(ind)
        class(model_factory_t), intent(in)  :: self
        character(*),           intent(in)  :: string

        integer(ik)                 :: ind, imodel
        character(:),   allocatable :: model_name
        logical                     :: found

        ind = 0
        do imodel = 1,self%nmodels()

            model_name = self%models(imodel)%model%get_name()
            found      = string_to_upper(trim(string)) == string_to_upper(trim(model_name))

            if (found) then
                ind = imodel
                exit
            end if

        end do

    end function index_by_name
    !************************************************************************************





    !>  Register model_t's in the model_factory.
    !!
    !!  This should be getting called at start-up by the ChiDG
    !!  framework.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine register_models()
        integer(ik) :: imodel

        type(ideal_gas_t)                               :: IDEAL_GAS
        type(ideal_gas_sst_t)                           :: IDEAL_GAS_SST
        type(ideal_gas_rstm_t)                          :: IDEAL_GAS_RSTM
        type(reynolds_shear_stress_t)                   :: REYNOLDS_SHEAR_STRESS

        type(shear_stress_t)                            :: SHEAR_STRESS
        type(vorticity_t)                               :: VORTICITY
        type(temperature_gradient_t)                    :: TEMPERATURE_GRADIENT
        type(pressure_gradient_t)                       :: PRESSURE_GRADIENT
        type(velocity_gradients_t)                      :: VELOCITY_GRADIENTS
        type(sutherlands_law_t)                         :: SUTHERLANDS_LAW
        type(constant_viscosity_t)                      :: CONSTANT_VISCOSITY
        type(stokes_hypothesis_t)                       :: STOKES_HYPOTHESIS
        type(reynolds_analogy_t)                        :: REYNOLDS_ANALOGY
        type(zero_turbulent_model_fields_t)             :: ZERO_TURBULENT_MODEL_FIELDS
        type(zero_reynolds_stress_t)                    :: ZERO_REYNOLDS_STRESS

        type(spalart_allmaras_turbulent_model_fields_t) :: SPALART_ALLMARAS_TURBULENT_MODEL_FIELDS
        type(fluid_wave_speed_t)                        :: FLUID_WAVE_SPEED
        type(critical_sound_speed_t)                    :: CRITICAL_SOUND_SPEED

        type(strain_rate_t)                             :: STRAIN_RATE
        type(rotation_rate_t)                           :: ROTATION_RATE

        type(velocity_div_curl_t)                       :: VEL_DIV_CURL
        type(velocity_gradient_t)                       :: VEL_GRAD

        type(modified_ducros_sensor_t)                  :: MOD_DUCROS
        type(artificial_viscosity_t)                    :: ART_VISC
        type(zero_artificial_viscosity_t)               :: ZERO_ART_VISC
        type(unsmoothed_artificial_viscosity_t)         :: UNS_ART_VISC
        type(rbf_smoothed_artificial_viscosity_t)       :: RBF_ART_VISC
        type(vertex_smoothed_artificial_viscosity_t)    :: VERTEX_ART_VISC
        type(pde_smoothed_artificial_viscosity_t)       :: PDE_ART_VISC

        type(mnp_shock_sensor_t)                        :: MNP_SHOCK_SENSOR
        type(mnp_artificial_viscosity_t)                :: MNP_ARTIFICIAL_VISCOSITY
        type(vertex_smoothed_mnp_artificial_viscosity_t)    :: VERTEX_SMOOTHED_MNP_ARTIFICIAL_VISCOSITY

        type(mnph_shock_sensor_t)                        :: MNPH_SHOCK_SENSOR
        type(mnph_artificial_viscosity_t)                :: MNPH_ARTIFICIAL_VISCOSITY
        type(mnpha_artificial_viscosity_t)               :: MNPHA_ARTIFICIAL_VISCOSITY

        type(sst_turbulence_kinetic_energy_t)            :: SST_TKE
        type(sst_turbulence_quantities_t)                :: SST_TQ
        type(sst_blended_coefficients_t)                 :: SST_COEFF
        type(sst_source_terms_t)                        :: SST_SRC


        type(rstm_ssglrrw_blended_coefficients_t)       :: RSTM_SSGLRRW_BLENDED_COEFFICIENTS
        type(rstm_ssglrrw_generalized_diffusion_t)      :: RSTM_SSGLRRW_GENERALIZED_DIFFUSION
        type(rstm_ssglrrw_simple_diffusion_t)           :: RSTM_SSGLRRW_SIMPLE_DIFFUSION
        type(rstm_ssglrrw_isotropic_dissipation_t)      :: RSTM_SSGLRRW_ISOTROPIC_DISSIPATION
        type(rstm_ssglrrw_pressure_strain_correlation_t)  :: RSTM_SSGLRRW_PRESSURE_STRAIN_CORRELATION
        type(rstm_ssglrrw_production_t)                 :: RSTM_SSGLRRW_PRODUCTION
        !type(rstm_ssglrrw_reynolds_stress_t)            :: RSTM_SSGLRRW_REYNOLDS_STRESS
        type(rstm_ssglrrw_turbulence_quantities_t)      :: RSTM_SSGLRRW_TURBULENCE_QUANTITIES
!        type(rstm_ssglrrw_turbulence_kinetic_energy_t)      :: RSTM_SSGLRRW_TURBULENCE_KE


        type(wall_distance_m)                           :: WALL_DISTANCE_NORMALIZATION


        if (.not. models_initialized) then

            call model_factory%register(IDEAL_GAS)
            call model_factory%register(IDEAL_GAS_SST)
            call model_factory%register(IDEAL_GAS_RSTM)
            call model_factory%register(REYNOLDS_SHEAR_STRESS)

            call model_factory%register(SHEAR_STRESS)
            call model_factory%register(VORTICITY)
            call model_factory%register(TEMPERATURE_GRADIENT)
            call model_factory%register(PRESSURE_GRADIENT)
            call model_factory%register(VELOCITY_GRADIENTS)
            call model_factory%register(SUTHERLANDS_LAW)
            call model_factory%register(CONSTANT_VISCOSITY)
            call model_factory%register(STOKES_HYPOTHESIS)
            call model_factory%register(REYNOLDS_ANALOGY)
            call model_factory%register(ZERO_TURBULENT_MODEL_FIELDS)
            call model_factory%register(ZERO_REYNOLDS_STRESS)
            call model_factory%register(SPALART_ALLMARAS_TURBULENT_MODEL_FIELDS)
            call model_factory%register(WALL_DISTANCE_NORMALIZATION)
            call model_factory%register(FLUID_WAVE_SPEED)
            call model_factory%register(CRITICAL_SOUND_SPEED)

            call model_factory%register(STRAIN_RATE)
            call model_factory%register(ROTATION_RATE)
            call model_factory%register(VEL_DIV_CURL)
            call model_factory%register(VEL_GRAD)


            call model_factory%register(MOD_DUCROS)
            call model_factory%register(ART_VISC)
            call model_factory%register(ZERO_ART_VISC)
            call model_factory%register(UNS_ART_VISC)
            call model_factory%register(RBF_ART_VISC)
            call model_factory%register(PDE_ART_VISC)
            call model_factory%register(VERTEX_ART_VISC)

            call model_factory%register(MNP_SHOCK_SENSOR)
            call model_factory%register(MNP_ARTIFICIAL_VISCOSITY)
            call model_factory%register(VERTEX_SMOOTHED_MNP_ARTIFICIAL_VISCOSITY)

            call model_factory%register(MNPH_SHOCK_SENSOR)
            call model_factory%register(MNPH_ARTIFICIAL_VISCOSITY)
            call model_factory%register(MNPHA_ARTIFICIAL_VISCOSITY)

            call model_factory%register(SST_TKE)
            call model_factory%register(SST_TQ)
            call model_factory%register(SST_COEFF)
            call model_factory%register(SST_SRC)

            call model_factory%register(RSTM_SSGLRRW_BLENDED_COEFFICIENTS)
            call model_factory%register(RSTM_SSGLRRW_GENERALIZED_DIFFUSION)
            call model_factory%register(RSTM_SSGLRRW_SIMPLE_DIFFUSION)
            call model_factory%register(RSTM_SSGLRRW_ISOTROPIC_DISSIPATION)
            call model_factory%register(RSTM_SSGLRRW_PRESSURE_STRAIN_CORRELATION)
            call model_factory%register(RSTM_SSGLRRW_PRODUCTION)
            !call model_factory%register(RSTM_SSGLRRW_REYNOLDS_STRESS)
            call model_factory%register(RSTM_SSGLRRW_TURBULENCE_QUANTITIES)
!            call model_factory%register(RSTM_SSGLRRW_TURBULENCE_KE)

            models_initialized = .true.

        end if


    end subroutine register_models
    !*************************************************************************************


end module mod_models
