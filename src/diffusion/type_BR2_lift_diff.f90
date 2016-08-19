module type_BR2_lift_diff
#include <messenger.h>
    use mod_kinds,              only: ik
    use type_mesh,              only: mesh_t
    use type_element_info,      only: element_info_t
    use type_chidgVector,       only: chidgVector_t
    use type_BR2_lift,          only: BR2_lift_t
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public :: BR2_lift_diff_t

        type(BR2_lift_t),  allocatable :: eqn(:)

    contains

        procedure   :: init
        procedure   :: update

    end type BR2_lift_diff_t
    !****************************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self,mesh,elem_info,iface,idiff)
        class(BR2_lift_diff_t), intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: idiff

        integer(ik) :: neqns, ieqn, ierr
        logical     :: allocate_equations, reallocate_equations


        !
        ! (Re)Allocate storage for each equation
        !
        neqns = mesh(elem_info%idomain_l)%neqns


        allocate_equations = (.not. allocated(self%eqn))

        if (allocate_equations) then

            allocate(self%eqn(neqns), stat=ierr)
            if (ierr /= 0) call AllocationError

        
        else 

            reallocate_equations = (neqns /= size(self%eqn))

            if (reallocate_equations) then
                deallocate(self%eqn)
                allocate(self%eqn(neqns), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        end if


        !
        ! Call initialization for each equation
        !
        do ieqn = 1,size(self%eqn)
            call self%eqn(ieqn)%init(mesh,elem_info,iface,idiff)
        end do




    end subroutine init
    !*****************************************************************************************















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine update(self,mesh,elem_info,iface,idiff,q,BR2_TYPE)
        class(BR2_lift_diff_t),     intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(element_info_t),       intent(in)      :: elem_info
        integer(ik),                intent(in)      :: iface
        integer(ik),                intent(in)      :: idiff
        type(chidgVector_t),        intent(in)      :: q
        integer(ik),                intent(in)      :: BR2_TYPE

        integer(ik) :: ieqn


        !
        ! Compute the lifting operators for each equations
        !
        do ieqn = 1,size(self%eqn)
            call self%eqn(ieqn)%update(mesh,elem_info,iface,idiff,ieqn,q,BR2_TYPE)
        end do

    end subroutine update
    !********************************************************************************************
























end module type_BR2_lift_diff
