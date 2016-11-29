module type_model
#include <messenger.h>
    use mod_kinds,          only: ik
    use type_chidg_worker,  only: chidg_worker_t
    use mod_string,         only: string_t
    implicit none




    !>  A class for implementing model parameters.
    !!
    !!
    !!  Examples: 
    !!  ---------
    !!  Equation of state:
    !!      name       = "Perfect Gas Law"
    !!      parameters = "Pressure" , "Temperature"
    !!
    !!  Viscosity model:
    !!      name       = "Sutherland's Law"
    !!      parameters = "Viscosity"
    !!
    !!
    !!  Usage:
    !!  ---------
    !!  Every new model is required to implement the 'init' and 'compute'
    !!  procedures. 'init' is executed by the framework at start-up 
    !!  and informs the model of what parameters it is contributing to.
    !!  'compute' is responsible for computing a contribution to the 
    !!  parameter.
    !!
    !!
    !!  Note:
    !!  ---------
    !!  Multiple models can contribute to the same parameter. Examples
    !!  of this would be viscosity models - laminar + turbulent eddy 
    !!  viscosity. Also the equation of state can have source terms
    !!  that could be accounted for by including extra models.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !------------------------------------------------------------------
    type, public, abstract :: model_t

        character(:),   allocatable :: name
        type(string_t), allocatable :: parameters(:)

    contains

        procedure(self_interface),      deferred :: init
        procedure(compute_interface),   deferred :: compute

        procedure   :: set_name
        procedure   :: get_name

        procedure   :: add_parameter
        procedure   :: get_parameter
        procedure   :: nparameters

    end type model_t
    !******************************************************************




    abstract interface
        subroutine self_interface(self)
            import model_t
            class(model_t), intent(inout)   :: self
        end subroutine
    end interface




    abstract interface
        subroutine compute_interface(self,worker)
            import model_t
            import chidg_worker_t
            class(model_t),         intent(in)      :: self
            type(chidg_worker_t),   intent(inout)   :: worker
        end subroutine
    end interface





contains



    !>  Set the 'name' component of the model_t.
    !!
    !!  This should be called from 'init' in a model implementation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !-------------------------------------------------------------------
    subroutine set_name(self,string)
        class(model_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        self%name = string

    end subroutine set_name
    !*******************************************************************




    !>  Return the name of the model_t.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !------------------------------------------------------------------
    function get_name(self) result(string)
        class(model_t), intent(in)  :: self

        character(:),   allocatable :: string

        string = self%name

    end function get_name
    !******************************************************************





    !>  Add a parameter that the model is providing a definition for.
    !!
    !!  This should be called from 'init' in a model implementation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !------------------------------------------------------------------
    subroutine add_parameter(self,string)
        class(model_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        integer(ik)                 :: ierr, ieq
        type(string_t), allocatable :: temp(:)


        !
        ! Extend fields if necessary
        !
        if (allocated(self%parameters)) then
            
            allocate(temp(size(self%parameters) + 1), stat=ierr)
            do ieq = 1,size(self%parameters)
                temp(ieq) = self%parameters(ieq)
            end do
            
        else
            allocate(temp(1), stat=ierr)
        end if
        if (ierr /= 0) call AllocationError


        !
        ! Set new variable
        !
        call temp(size(temp))%set(string)


        !
        ! Copy temp back to self
        !
        self%parameters = temp


    end subroutine add_parameter
    !******************************************************************







    !>  Return the parameter from an index location in self%parameters.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !------------------------------------------------------------------
    function get_parameter(self,iparam) result(string)
        class(model_t), intent(in)  :: self
        integer(ik),    intent(in)  :: iparam

        character(:),   allocatable :: string

        string = self%parameters(iparam)%get()

    end function get_parameter
    !******************************************************************









    !>  Return number of parameters the model is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !------------------------------------------------------------------
    function nparameters(self) result(nparams)
        class(model_t), intent(in)  :: self

        integer(ik) :: nparams

        nparams = size(self%parameters)

    end function nparameters
    !******************************************************************





end module type_model
