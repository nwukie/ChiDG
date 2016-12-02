module mock_model
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    implicit none



    !>  A mock model implementation for testing purposes.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !--------------------------------------------------------------
    type, extends(model_t), public :: mock_model_t

    contains

        procedure   :: init
        procedure   :: compute

    end type mock_model_t
    !**************************************************************




contains



    !>
    !!
    !!
    !-------------------------------------------------------------
    subroutine init(self)
        class(mock_model_t),    intent(inout)   :: self

        ! Set model
        call self%set_name('Mock Model')

        ! Add model parameters
        call self%add_model_field('Parameter One')
        call self%add_model_field('Parameter Two')

    end subroutine init
    !*************************************************************









    !>
    !!
    !!
    !!
    !-------------------------------------------------------------
    subroutine compute(self,worker)
        class(mock_model_t),    intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker





    end subroutine compute
    !**************************************************************


end module mock_model
