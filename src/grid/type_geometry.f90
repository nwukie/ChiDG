module type_geometry
#include <messenger.h>
    implicit none





    !>  A geometry data type for implementing definitions for Hexahedrals,
    !!  Tets, etc.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/15/2017
    !!
    !!
    !------------------------------------------------------------------------
    type, abstract, public :: geometry_t


    contains

        procedure, deferred :: map
        procedure, deferred :: 

    end type geometry_t
    !************************************************************************








contains









end module type_geometry
