module mod_preconditioner
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_preconditioner,    only: preconditioner_t



    ! Import preconditioner types
    use precon_identity,            only: precon_identity_t
    use precon_jacobi,              only: precon_jacobi_t
    use precon_ILU0,                only: precon_ILU0_t
!    use precon_RASILU0,             only: precon_RASILU0_t
    implicit none



    ! Instantiate preconditioner types for sourcing
    type(precon_identity_t)             :: IDENTITY
    type(precon_jacobi_t)               :: BLOCKJACOBI
    type(precon_ILU0_t)                 :: ILU0
!    type(precon_RASILU0_t)              :: RASILU0



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/22/2016
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine create_preconditioner(precon_string,instance)
        character(*),                           intent(in)      :: precon_string
        class(preconditioner_t), allocatable,   intent(inout)   :: instance

        character(:),   allocatable :: user_msg

        select case (trim(precon_string))
            case ('identity','Identity','IDENTITY')
                allocate(instance, source=IDENTITY)


            case ('jacobi','Jacobi','blockjacobi','BlockJacobi')
                allocate(instance, source=BLOCKJACOBI)

            case('ilu0','ILU0')
                allocate(instance, source=ILU0)

!            case('rasilu0','RASILU0', 'ras-ilu0', 'RAS-ILU0')
!                allocate(instance, source=RASILU0)

            case default
                user_msg = "create_preconditioner: It seems like we can't find the preconditioner &
                            that was specified in the input string. You might check the spelling &
                            and also make sure the preconditioner is registered in the file &
                            mod_preconditioner.f90."
                call chidg_signal_one(FATAL,user_msg,trim(precon_string))

        end select


        !
        ! Make sure the preconditioner was allocated
        !
        user_msg = "create_preconditioner: preconditioner was not allocated. &
                    Check that the desired solver was registered and instantiated &
                    in the mod_preconditioner module."
        if (.not. allocated(instance)) call chidg_signal_one(FATAL,user_msg,trim(precon_string))

    end subroutine create_preconditioner
    !***********************************************************************************







end module mod_preconditioner
