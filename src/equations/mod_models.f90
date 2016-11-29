module mod_models
#include <messenger.h>
    implicit none




    !>  A model_t factory for initializing, holding, and producing
    !!  model_t's on-demand.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !---------------------------------------------------------------
    type, public :: model_factory_t

        type(model_wrapper_t),  allocatable :: models(:)

    contains

        procedure   :: register
        procedure   :: produce

    end type model_factory_t
    !***************************************************************



    type(model_factory_t)   :: model_factory
    logical                 :: models_initialized = .false.



contains




    !>
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------
    subroutine register(self,model)
        class(model_factory_t), intent(inout)   :: self
        class(model_t),         intent(in)      :: model







    end subroutine register
    !****************************************************************





    !>  Given a string, find a model with a name matching the string
    !!  and return a model allocated with that concrete type.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !----------------------------------------------------------------
    function produce(self,string) result(model)
        class(model_factory_t), intent(in)  :: self
        character(*),           intent(in)  :: string

        class(model_t), allocatable :: model




    end function produce
    !****************************************************************











    !>  Register model_t's in the model_factory.
    !!
    !!  This should be getting called at start-up by the ChiDG
    !!  framework.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !----------------------------------------------------------------
    subroutine register_models()
        integer(ik) :: imodel


        if (.not. models_initialized) then



            ! Initialize all models
            do imodel = 1,model_factory%nmodels()
                call model_factory%models(imodel)%model%init()
            end do

            models_initialized = .true.


        end if





    end subroutine register_models
    !****************************************************************


end module mod_models
