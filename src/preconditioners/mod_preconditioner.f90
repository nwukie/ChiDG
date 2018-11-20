module mod_preconditioner
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_preconditioner,    only: preconditioner_t



    ! Import preconditioner types
    use precon_identity,        only: precon_identity_t
    use precon_jacobi,          only: precon_jacobi_t
    use precon_ILU0,            only: precon_ILU0_t
    use precon_ILU0_overset,    only: precon_ILU0_overset_t
    use precon_HB,              only: precon_HB_t
    use precon_RASILU0,         only: precon_RASILU0_t
    use precon_line,            only: precon_line_t
    use precon_schur_element,   only: precon_schur_element_t
    implicit none



    ! Instantiate preconditioner types for sourcing
    type(precon_identity_t)         :: P_IDENTITY
    type(precon_jacobi_t)           :: P_BLOCKJACOBI
    type(precon_ILU0_t)             :: P_ILU0
    type(precon_ILU0_overset_t)     :: P_ILU0_OVERSET
    type(precon_HB_t)               :: P_HB
    type(precon_RASILU0_t)          :: P_RASILU0
    type(precon_line_t)             :: P_LINE
    type(precon_schur_element_t)    :: P_SCHUR_ELEMENT



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
                allocate(instance, source=P_IDENTITY)


            case ('jacobi','Jacobi','blockjacobi','BlockJacobi')
                allocate(instance, source=P_BLOCKJACOBI)

            case('ilu0','ILU0')
                allocate(instance, source=P_ILU0)

            case('ilu0_overset','ILU0_OVERSET')
                allocate(instance, source=P_ILU0_OVERSET)

            case('HB','hb')
                allocate(instance, source=P_HB)

            case('rasilu0','RASILU0', 'ras-ilu0', 'RAS-ILU0', 'RAS+ILU0', 'ras+ilu0', 'RAS ILU0', 'ras ilu0')
                allocate(instance, source=P_RASILU0)

            case('line','Line', 'LINE')
                allocate(instance, source=P_LINE)

            case('Schur Element')
                allocate(instance, source=P_SCHUR_ELEMENT)

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
