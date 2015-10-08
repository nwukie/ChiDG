module mod_preconditioner
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_preconditioner,    only: preconditioner_t



    ! Import preconditioner types
    use precon_identity,            only: precon_identity_t
    use precon_jacobi,              only: precon_jacobi_t
    use precon_ILU0,                only: precon_ILU0_t
    implicit none



    ! Instantiate preconditioner types for sourcing
    type(precon_identity_t)             :: IDENTITY
    type(precon_jacobi_t)               :: BLOCKJACOBI
    type(precon_ILU0_t)                 :: ILU0



contains


    subroutine create_preconditioner(precon_string,instance)
        character(*),                           intent(in)      :: precon_string
        class(preconditioner_t), allocatable,   intent(inout)   :: instance

        select case (trim(precon_string))
            case ('identity','Identity','IDENTITY')
                allocate(instance, source=IDENTITY)


            case ('jacobi','Jacobi','blockjacobi','BlockJacobi')
                allocate(instance, source=BLOCKJACOBI)

            case('ilu0','ILU0')
                allocate(instance, source=ILU0)

            case default
                call signal(FATAL,'create_preconditioner -- preconditioner string not recognized')

        end select




        !
        ! Make sure the preconditioner was allocated
        !
        if (.not. allocated(instance)) call signal(FATAL,"create_preconditioner: preconditioner was not allocated. Check that the desired solver was registered and instantiated in the mod_preconditioner module")




    end subroutine create_preconditioner







end module mod_preconditioner
