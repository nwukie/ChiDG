module mod_matrixsolver
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use atype_matrixsolver, only: matrixsolver_t


    ! IMPORT MATRIX SOLVERS
    use type_directsolver,  only: directsolver_t
    
    



    type(directsolver_t)    :: DIRECT






contains





    !>  Factory method for creating matrixsolver objects
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mstring     Character string used to select the appropriate matrixsolver for allocation
    !!  @param[inout]   msolver     matrixsolver_t that will be allocated to a concrete type.
    !------------------------------------------------------
    subroutine create_matrixsolver(mstring,msolver)
        character(len=*),                    intent(in)      :: mstring
        class(matrixsolver_t), allocatable,  intent(inout)   :: msolver

        integer(ik) :: ierr


        select case (trim(mstring))
            case ('direct','Direct')
                allocate(msolver, source=DIRECT, stat=ierr)
                if (ierr /= 0) call AllocationError

            case default
                call signal(FATAL,"create_matrixsolver: matrix solver string did not match any valid type")

        end select


    end subroutine create_matrixsolver




end module mod_matrixsolver
