module type_file_properties
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    implicit none




    !>
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
    type, public :: file_properties_t

        logical     :: contains_grid        = .false.
        logical     :: contains_solution    = .false.

        integer(ik)                      :: ndomains = 0     !< Number of domains in the file
        character(len=1024), allocatable :: domain_names(:)
        integer(ik),         allocatable :: order_c(:)       !< Coordinate order ( 1st, 2nd, 3rd, etc. )
        integer(ik),         allocatable :: order_s(:)       !< Solution order ( 1st, 2nd, 3rd, etc. )
        integer(ik),         allocatable :: nterms_c(:)      !< Number of terms in the coordinate expansion
        integer(ik),         allocatable :: nterms_s(:)      !< Number of terms in the solution expansion
 
        character(len=1024), allocatable :: eqnset(:)        !< Equation set specified in the file.


    contains

        procedure   :: set_ndomains

    end type file_properties_t
    !##############################################################################################






contains




    !> Initialize the file_properties type by allocating the other components to the right sizes
    !! according to the number of domains detected.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine set_ndomains(self,ndomains)
        class(file_properties_t),   intent(inout)   :: self
        integer(ik),                intent(in)      :: ndomains

        integer :: ierr

        !
        ! Set number of domains in the file
        !
        self%ndomains = ndomains


        !
        ! Allocate storage for other data
        !
        if ( allocated(self%order_c)  ) deallocate( self%order_c  )
        if ( allocated(self%order_s)  ) deallocate( self%order_s  )
        if ( allocated(self%nterms_c) ) deallocate( self%nterms_c )
        if ( allocated(self%nterms_s) ) deallocate( self%nterms_s )
        if ( allocated(self%eqnset)   ) deallocate( self%eqnset   )


        allocate( self%order_c(ndomains),   &
                  self%order_s(ndomains),   &
                  self%nterms_c(ndomains),  &
                  self%nterms_s(ndomains),  &
                  self%eqnset(ndomains), stat=ierr )

        if (ierr /= 0) call AllocationError


    end subroutine
    !################################################################################################





end module type_file_properties
