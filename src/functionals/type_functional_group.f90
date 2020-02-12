module type_functional_group
#include <messenger.h>
    use mod_kinds,                  only: ik, rk
    use type_functional_wrapper,    only: functional_wrapper_t
    use type_evaluator,             only: evaluator_t
    use type_svector,               only: svector_t
    implicit none

    !>  Group of functional that needs to be optimized
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/14/2017
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: functional_group_t

        class(functional_wrapper_t), allocatable :: fcl_entities(:)

        logical :: compute_functionals = .false.

    contains

        procedure   :: init
        procedure   :: n_functionals
        procedure   :: functionals_name
        procedure   :: add_functional
        procedure   :: release

    end type functional_group_t
    !******************************************************************************************



contains




    !> Initialize functional_group by allocating it with the input number of functionals
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/14/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self,nfunc)
        class(functional_group_t),  intent(inout)   :: self
        integer(ik),                intent(in)      :: nfunc
        
        integer(ik)     :: ierr

        ! Allocate fcl_entities
        if (allocated(self%fcl_entities)) deallocate(self%fcl_entities)
        allocate(self%fcl_entities(nfunc), stat=ierr)
        if (ierr/=0) call chidg_signal(FATAL,"type_functional_group%init: error allocating the functional group")

        ! Confirm functionals calculation
        self%compute_functionals = .true.

    end subroutine init
    !******************************************************************************************






    !>  Add functional to group 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/4/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine add_functional(self,ifunc)
        class(functional_group_t),          intent(inout)   :: self
        class(evaluator_t), allocatable,    intent(in)      :: ifunc
        
        type(functional_wrapper_t), allocatable     :: temp_group(:)

        integer(ik)     :: nfunc, new_size, ierr

        ! Get number of functionals already stored
        if (allocated(self%fcl_entities)) then
            nfunc = self%n_functionals()
        else
            nfunc = 0
        end if

        ! Allocate fcl_entities
        select case (nfunc)
            ! No functional stored yet
            case (0)
            
                allocate(self%fcl_entities(1), stat=ierr)
                if (ierr/=0) call chidg_signal(FATAL,"type_functional_group%init: error allocating the functional group")
                allocate(self%fcl_entities(1)%func, source = ifunc, stat=ierr)
                if (ierr/=0) call AllocationError
            
            ! Other functionals already stored
            case default

                new_size = nfunc + 1

                allocate (temp_group(new_size),stat=ierr)
                if (ierr/=0) call AllocationError

                !Copy functional from current to temporaru vector
                if (allocated(self%fcl_entities)) then
                    temp_group(lbound(self%fcl_entities,1):ubound(self%fcl_entities,1)) = self%fcl_entities
                end if

                ! Move allocation data back to self%fcl_entities and dellocate temp_group
                call move_alloc(FROM=temp_group,TO=self%fcl_entities)
                
                allocate(self%fcl_entities(new_size)%func, source = ifunc, stat=ierr)
                if (ierr/=0) call AllocationError

        end select

        ! Confirm functionals calculation
        self%compute_functionals = .true.

    end subroutine add_functional
    !******************************************************************************************





    !>  Overall number of functional entered
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/30/2017
    !!
    !-----------------------------------------------------------------------------------------
    function n_functionals(self) result(n_entities)
        class(functional_group_t),  intent(inout)   :: self

        integer(ik)                 :: n_entities
        character(:),   allocatable :: user_msg

        if (allocated(self%fcl_entities)) then
            n_entities = size(self%fcl_entities)
        else
            n_entities = 0
            !user_msg = "functional_group%n_functionals: no functional was requested in 'chidg edit'"
            !call chidg_signal(FATAL,user_msg)
        end if

    end function n_functionals
    !******************************************************************************************





    !>  Return functionals' names
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/30/2018
    !!
    !-----------------------------------------------------------------------------------------
    function functionals_name(self) result(names)
        class(functional_group_t),  intent(inout)   :: self

        type(svector_t)     :: names
        integer(ik)         :: nfuncs, ifunc 
        
        nfuncs = self%n_functionals()
        
        do ifunc = 1,nfuncs
            call names%push_back(self%fcl_entities(ifunc)%func%name)
        end do 

    end function functionals_name
    !******************************************************************************************





    !>  Release allocated memory
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/13/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine release(self)
        class(functional_group_t),  intent(inout)   :: self

        if (allocated (self%fcl_entities)) deallocate (self%fcl_entities)

    end subroutine release
    !******************************************************************************************


end module type_functional_group
