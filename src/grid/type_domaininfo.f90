module type_domaininfo


    !> Container for domain information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    type, public :: domaininfo_t

        character(len=:), allocatable   :: name

    end type domaininfo_t
    !************************************************************************









end module type_domaininfo
