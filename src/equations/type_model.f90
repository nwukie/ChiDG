module type_model
#include <messenger.h>
    use mod_kinds,          only: ik
    use type_chidg_worker,  only: chidg_worker_t
    use mod_string,         only: string_t
    implicit none




    !>  A class for implementing models.
    !!
    !!
    !!  Examples: 
    !!  ---------
    !!  Equation of state:
    !!      name         = "Perfect Gas Law"
    !!      model fields = "Pressure" , "Temperature"
    !!
    !!  Viscosity model:
    !!      name         = "Sutherland's Law"
    !!      model fields = "Viscosity"
    !!
    !!
    !!  Usage:
    !!  ---------
    !!  Every new model is required to implement the 'init' and 'compute'
    !!  procedures. 'init' is executed by the framework at start-up 
    !!  and informs the model of what fields it is contributing to.
    !!  'compute' is responsible for computing a contribution to the 
    !!  field.
    !!
    !!
    !!  Note:
    !!  ---------
    !!  Multiple models can contribute to the same field. Examples
    !!  of this would be viscosity models - laminar + turbulent eddy 
    !!  viscosity. Also the equation of state can have source terms
    !!  that could be accounted for by including extra models.
    !!
    !!  Model dependencies:
    !!  -------------------
    !!  Models can have different dependency strings, indicating to the infrastructure
    !!  what the model needs differentiated with respect to when they are being computed.
    !!  The model initialization routine shall set the dependency string for the model
    !!  being implemented:
    !!      Options:
    !!          call self%set_dependency('f(Q-)')
    !!          call self%set_dependency('f(Q-,Q+)')
    !!          call self%set_dependency('f(Grad(Q))')
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !------------------------------------------------------------------------------------
    type, public, abstract :: model_t

        character(:),   allocatable :: name
        character(:),   allocatable :: dependency
        type(string_t), allocatable :: model_fields(:)

    contains

        procedure(self_interface),      deferred :: init
        procedure(compute_interface),   deferred :: compute

        procedure   :: set_name
        procedure   :: get_name

        procedure   :: set_dependency
        procedure   :: get_dependency

        procedure   :: add_model_field
        procedure   :: get_model_field
        procedure   :: nmodel_fields

    end type model_t
    !*************************************************************************************




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
    !-------------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(model_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        self%name = string

    end subroutine set_name
    !*************************************************************************************




    !>  Return the name of the model_t.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !-------------------------------------------------------------------------------------
    function get_name(self) result(string)
        class(model_t), intent(in)  :: self

        character(:),   allocatable :: string

        if (allocated(self%name)) then
            string = self%name
        else
            string = ' '
        end if

    end function get_name
    !*************************************************************************************












    !>  Set the 'dependency' component of the model_t.
    !!
    !!  This should be called from 'init' in a model implementation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/21/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine set_dependency(self,string)
        class(model_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        character(:),   allocatable :: user_msg

        if ( (trim(string) == 'f(Q-)')     .or.    &
             (trim(string) == 'f(Q-,Q+)')  .or.    &
             (trim(string) == 'f(Grad(Q))') ) then
            self%dependency = string
        else
            user_msg = "model%set_dependency: Invalid model dependency string. Valid &
                        options are 'f(Q-)', 'f(Q-,Q+)', and 'f(Grad(Q))'."
            call chidg_signal_one(FATAL,user_msg,self%get_name())
        end if

    end subroutine set_dependency
    !*************************************************************************************




    !>  Return the 'dependency' of the model_t.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/21/2016
    !!
    !-------------------------------------------------------------------------------------
    function get_dependency(self) result(string)
        class(model_t), intent(in)  :: self

        character(:),   allocatable :: string, user_msg

        if (allocated(self%dependency)) then
            string = self%dependency
        else
            user_msg = "model%set_dependency: A model dependency string was never set &
                        for the current model. Make sure self%set_dependency is getting &
                        called in the model initialization routine, model%init()"
            call chidg_signal_one(FATAL,user_msg,self%get_name())
        end if


    end function get_dependency
    !*************************************************************************************


























    !>  Add a model field that the model is providing a definition for.
    !!
    !!  This should be called from 'init' in a model implementation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine add_model_field(self,string)
        class(model_t), intent(inout)   :: self
        character(*),   intent(in)      :: string

        integer(ik)                 :: ierr, ieq
        type(string_t), allocatable :: temp(:)


        !
        ! Extend fields if necessary
        !
        if (allocated(self%model_fields)) then
            
            allocate(temp(size(self%model_fields) + 1), stat=ierr)
            do ieq = 1,size(self%model_fields)
                temp(ieq) = self%model_fields(ieq)
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
        self%model_fields = temp


    end subroutine add_model_field
    !*************************************************************************************







    !>  Return the field from an index location in self%model_fields.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !-------------------------------------------------------------------------------------
    function get_model_field(self,ifield) result(string)
        class(model_t), intent(in)  :: self
        integer(ik),    intent(in)  :: ifield

        character(:),   allocatable :: string, user_msg

        ! Check bounds
        user_msg = "get_model_field: field index is out of bounds."
        if (ifield > self%nmodel_fields()) call chidg_signal(FATAL,user_msg)


        string = self%model_fields(ifield)%get()

    end function get_model_field
    !*************************************************************************************








    !>  Return number of model fields the model is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/29/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    function nmodel_fields(self) result(nfields)
        class(model_t), intent(in)  :: self

        integer(ik) :: nfields

        if (allocated(self%model_fields)) then
            nfields = size(self%model_fields)
        else
            nfields = 0
        end if


    end function nmodel_fields
    !*************************************************************************************





end module type_model
