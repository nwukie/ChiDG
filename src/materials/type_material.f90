module type_material
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public, abstract :: material_t

        character(len=:),   allocatable :: name

!        ! List of available equation names for use by the operators
!        type(string_t),     allocatable :: equation_names

    contains


        procedure   :: set_name
        procedure   :: get_name

!        procedure(self_interface), deferred :: init
!        procedure   :: contribute_equation
!        procedure   :: access_equation

    end type material_t
    !***************************************************************************************



    abstract interface
        subroutine self_interface(self)
            import material_t

            class(material_t), intent(inout)    :: self
        end subroutine
    end interface




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_name(self,string)
        class(material_t),  intent(inout)   :: self
        character(len=*),   intent(in)      :: string

        self%name = string

    end subroutine set_name
    !****************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_name(self) result(string)
        class(material_t),  intent(inout)   :: self

        character(len=:), allocatable :: string

        string = self%name

    end function get_name
    !****************************************************************************************








!    !>
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   8/30/2016
!    !!
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine contribute_equation(self,string)
!        class(material_t),  intent(inout)   :: self
!        character(len=*),   intent(in)      :: string
!
!        integer(ik) :: neq, ieq
!        type(string_t), allocatable :: temp(:)
!        
!
!        neq = self%nequations()
!
!
!        ! Allocate temp storage to be bigger by one
!        allocate(temp(neq + 1), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!
!        ! Copy all old entries if necessary
!        if (allocated(self%equation_names)) then
!            do ieq = size(self%equation_names)
!                temp(ieq) = self%equation_names(ieq)
!            end do
!        end if
!
!
!
!        ! Set new variable string
!        temp(size(temp)) = string
!
!
!        ! Copy back to material
!        self%equation_names = temp
!
!    end subroutine contribute_equation
!    !******************************************************************************************
!
!
!
!
!
!
!
!
!    !>
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   8/30/2016
!    !!
!    !!
!    !!
!    !-------------------------------------------------------------------------------------------
!    function nequations(self) result(neq)
!        class(material_t),  intent(in)  :: self
!
!        integer(ik) :: neq
!
!        if (allocated(self%equation_names)) then
!            neq = size(self%equation_names)
!        else
!            neq = 0
!        end if
!
!    end function nequations
!    !*******************************************************************************************



end module type_material
