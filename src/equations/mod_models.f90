module mod_models
#include <messenger.h>
    use mod_kinds,              only: ik
    use mod_string,             only: string_to_upper
    use type_model_wrapper,     only: model_wrapper_t
    use type_model,             only: model_t

    use type_ideal_gas,                                 only: ideal_gas_t
    use type_sutherlands_law,                           only: sutherlands_law_t
    use type_stokes_hypothesis,                         only: stokes_hypothesis_t
    use type_reynolds_analogy,                          only: reynolds_analogy_t
    use type_zero_turbulent_model_fields,               only: zero_turbulent_model_fields_t
    use type_spalart_allmaras_turbulent_model_fields,   only: spalart_allmaras_turbulent_model_fields_t
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

        
        !
        ! Initialize the incoming model
        !
        call model%init()


        !
        ! Extend storage
        !
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



        !
        ! Store new model
        !
        allocate(temp(size(temp))%model, source=model, stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Move allocation
        !
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

        !
        ! Find location of model
        !
        imodel = self%index_by_name(string)


        !
        ! Check model was found
        !
        user_msg = "model_factory%produce: We couldn't find the model string in &
                    the list of registered models. Make sure the model was registered &
                    in the model factory."
        if (imodel == 0) call chidg_signal_one(FATAL,user_msg,trim(string))


        !
        ! Allocate model to be returned
        !
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
        type(sutherlands_law_t)                         :: SUTHERLANDS_LAW
        type(stokes_hypothesis_t)                       :: STOKES_HYPOTHESIS
        type(reynolds_analogy_t)                        :: REYNOLDS_ANALOGY
        type(zero_turbulent_model_fields_t)             :: ZERO_TURBULENT_MODEL_FIELDS
        type(spalart_allmaras_turbulent_model_fields_t) :: SPALART_ALLMARAS_TURBULENT_MODEL_FIELDS


        if (.not. models_initialized) then

            call model_factory%register(IDEAL_GAS)
            call model_factory%register(SUTHERLANDS_LAW)
            call model_factory%register(STOKES_HYPOTHESIS)
            call model_factory%register(REYNOLDS_ANALOGY)
            call model_factory%register(ZERO_TURBULENT_MODEL_FIELDS)
            call model_factory%register(SPALART_ALLMARAS_TURBULENT_MODEL_FIELDS)

            models_initialized = .true.

        end if


    end subroutine register_models
    !*************************************************************************************


end module mod_models
