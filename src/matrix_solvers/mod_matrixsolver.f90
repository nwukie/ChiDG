module mod_matrixsolver
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use atype_matrixsolver, only: matrixsolver_t
    use type_dict,          only: dict_t


    ! IMPORT MATRIX SOLVERS
!    use type_directsolver,  only: directsolver_t
!    use type_blockjacobi,   only: blockjacobi_t
!    use type_gaussseidel,   only: gaussseidel_t
!    use type_sor,           only: sor_t
!    use type_gmres,         only: gmres_t
    use type_fgmres,        only: fgmres_t
    

    



!    type(directsolver_t)    :: DIRECT
!    type(blockjacobi_t)     :: BLOCKJACOBI
!    type(gaussseidel_t)     :: GAUSSSEIDEL
!    type(sor_t)             :: SOR
!    type(gmres_t)           :: GMRES
    type(fgmres_t)          :: FGMRES






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
!            case ('direct','Direct')
!                allocate(msolver, source=DIRECT, stat=ierr)
!
!            case ('blockjacobi','BlockJacobi')
!                allocate(msolver, source=BLOCKJACOBI, stat=ierr)
!
!            case ('gaussseidel','GaussSeidel')
!                allocate(msolver, source=GAUSSSEIDEL, stat=ierr)
!
!            case ('sor','SOR')
!                allocate(msolver, source=SOR, stat=ierr)
!
!            case ('gmres','GMRES')
!                allocate(msolver, source=GMRES, stat=ierr)
!
            case ('fgmres','FGMRES')
                allocate(msolver, source=FGMRES, stat=ierr)


            case default
                call chidg_signal(FATAL,"create_matrixsolver: matrix solver string did not match any valid type")

        end select
        if (ierr /= 0) call AllocationError




        !
        ! Call options initialization if present
        !
        if (present(options)) then
            call msolver%set(options)
        end if

        


        !
        ! Make sure the solver was allocated
        !
        if (.not. allocated(msolver)) call chidg_signal(FATAL,"create_matrixsolver: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_matrixsolver module")












    end subroutine create_matrixsolver




end module mod_matrixsolver
