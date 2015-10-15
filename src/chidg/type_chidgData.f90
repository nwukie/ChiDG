module type_chidgData
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_domain,            only: domain_t
    use type_solverdata,        only: solverdata_t
    implicit none


    !> solver abstract type definition
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    type, public  :: chidgData_t

        !> Base solver data
        type(domain_t),         allocatable    :: domains(:)         !< Domain storage
        type(solverdata_t)                     :: sdata              !< Storage for matrix/vector data

        logical                         :: solverInitialized = .false.




        !! NOTE: if one wanted to add specialized data, instead of deriving from chidgData, you could add a
        !!       chidgExtension class that could be specialized further which could contain non-standard data 
        !!  class(chidgExtension_t)

    contains



    end type chidgData_t
    !-------------------------------------------------------------------------------------------------------




contains







end module type_chidgData
