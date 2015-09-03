module mod_matrixsolver
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use atype_matrixsolver, only: matrixsolver_t
    use type_dict,          only: dict_t


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
    subroutine create_matrixsolver(mstring,msolver,options)
        character(len=*),                    intent(in)      :: mstring
        class(matrixsolver_t), allocatable,  intent(inout)   :: msolver
        type(dict_t), optional,              intent(inout)   :: options

        integer(ik) :: ierr


        select case (trim(mstring))
            case ('direct','Direct')
                allocate(msolver, source=DIRECT, stat=ierr)
                if (ierr /= 0) call AllocationError

            case default
                call signal(FATAL,"create_matrixsolver: matrix solver string did not match any valid type")

        end select




        !
        ! Call options initialization if present
        !
        if (present(options)) then
            call msolver%set(options)
        end if

        


        !
        ! Make sure the solver was allocated
        !
        if (.not. allocated(msolver)) call signal(FATAL,"create_matrixsolver: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_matrixsolver module")












    end subroutine create_matrixsolver




end module mod_matrixsolver
