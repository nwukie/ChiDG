module mod_solverdata
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use atype_solverdata,   only: solverdata_t



    ! Import solverdata types
    use solverdata_base,        only: base_d
    implicit none


    ! Instantiate solver data structures for sourcing
    type(base_d) :: BASE_DATA




contains



    !>  Factory method for returning specialized solver data structures
    !!      - In general, the base data structure should provide most functionality (q, rhs, lin)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      dataString      String containing the name of the solverdata type to be sourced
    !!  @param[inout]   sdata           Allocatable solverdata_t to be source-allocated
    !-----------------------------------------------------------------------------
    subroutine create_solverdata(dataString,sdata)
        character(*),                       intent(in)      :: dataString
        class(solverdata_t), allocatable,   intent(inout)   :: sdata



        select case (trim(dataString))
            case ('base')
                allocate(sdata, source=BASE_DATA)

            case default
                call signal(FATAL,'create_solverdata -- equation string not recognized')
        end select




    end subroutine






end module mod_solverdata
