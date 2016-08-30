module type_operator
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: BOUNDARY_ADVECTIVE_FLUX, VOLUME_ADVECTIVE_FLUX, &
                                  BOUNDARY_DIFFUSIVE_FLUX, VOLUME_DIFFUSIVE_FLUX
    use mod_string,         only: string_t, string_to_upper
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  Valid Operator Types are:
    !!      BOUNDARY_ADVECTIVE_FLUX
    !!      VOLUME_ADVECTIVE_FLUX
    !!
    !----------------------------------------------------------------------------------
    type, abstract, public :: operator_t

        integer(ik)                         :: operator_type
        character(len=:),       allocatable :: name
        type(string_t),         allocatable :: eqns(:)

    contains

        procedure(self_interface),      deferred :: init
        procedure(compute_interface),   deferred :: compute

        procedure   :: set_name
        procedure   :: get_name

        procedure   :: set_operator_type
        procedure   :: get_operator_type
        procedure   :: set_equation

    end type operator_t
    !**********************************************************************************




    abstract interface
        subroutine compute_interface(self,worker,prop)
            import operator_t
            import chidg_worker_t
            import properties_t

            class(operator_t),      intent(in)      :: self
            type(chidg_worker_t),   intent(inout)   :: worker
            class(properties_t),    intent(inout)   :: prop
        end subroutine
    end interface




    abstract interface
        subroutine self_interface(self)
            import operator_t

            class(operator_t),      intent(inout)   :: self
        end subroutine
    end interface






contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine set_name(self, string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        self%name = trim(string)

    end subroutine set_name
    !****************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_name(self) result(string)
        class(operator_t),  intent(inout)   :: self

        character(len=:), allocatable :: string

        string = trim(self%name)

    end function get_name
    !****************************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine set_operator_type(self, string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        integer(ik)                     :: ind
        character(len=:),   allocatable :: upper
        logical                         :: volume, boundary, advective, diffusive
        


        ! Convert to upper case
        upper = string_to_upper(string)


        ! Get operator string key words
        ind = index(upper,"VOLUME")
        volume    = (ind /= 0)

        ind = index(upper,"BOUNDARY")
        boundary  = (ind /= 0)

        ind = index(upper,"ADVECTIVE")
        advective = (ind /= 0)

        ind = index(upper,"DIFFUSIVE")
        diffusive = (ind /= 0)

        
        ! Set operator type
        if (volume .and. advective) then
            self%operator_type = VOLUME_ADVECTIVE_FLUX
        else if (volume .and. diffusive) then
            self%operator_type = VOLUME_DIFFUSIVE_FLUX
        else if (boundary .and. advective) then
            self%operator_type = BOUNDARY_ADVECTIVE_FLUX
        else if (boundary .and. diffusive) then
            self%operator_type = BOUNDARY_DIFFUSIVE_FLUX
        else
            call chidg_signal(FATAL,"operator%set_operator_type: Did not recognize a valid operator type")
        end if


    end subroutine set_operator_type
    !**************************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    function get_operator_type(self) result(otype)
        class(operator_t),  intent(in)  :: self

        integer(ik) :: otype

        otype = self%operator_type

    end function get_operator_type
    !***************************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine set_equation(self,string)
        class(operator_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        integer(ik)     :: ierr, ieq
        type(string_t), allocatable  :: temp(:)



        !
        ! Extend eqns if necessary
        !
        if (allocated(self%eqns)) then
            
            allocate(temp(size(self%eqns) + 1), stat=ierr)

            do ieq = 1,size(self%eqns)
                temp(ieq) = self%eqns(ieq)
            end do
            
        else

            allocate(temp(1), stat=ierr)
        
        end if
        if (ierr /= 0) call AllocationError


        !
        ! Set new variable
        !
        temp(size(temp))%str = string


        !
        ! Copy temp back to self
        !
        self%eqns = temp


    end subroutine set_equation
    !***************************************************************************************************

















end module type_operator
