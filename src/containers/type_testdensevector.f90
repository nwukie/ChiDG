module type_testdensevector
    implicit none


    !
    ! NOTE: uncomment destructor to cause ICE. Cannot find small reproducer
    !



    type, public :: testdensevector_t

        real,  dimension(:), allocatable :: vec

    contains

!        final :: destructor
    end type testdensevector_t



    public operator (*)
    interface operator (*)
        module procedure mult_real_dv   ! real * testdensevector,   ELEMENTAL
    end interface



contains


    elemental function mult_real_dv(left,right) result(res)
        real,                       intent(in)  :: left
        type(testdensevector_t),    intent(in)  :: right
        type(testdensevector_t)                 :: res


        res%vec = left * right%vec
        
    end function





    !subroutine destructor(self)
    !    type(testdensevector_t),    intent(inout)   :: self
!
!    end subroutine

end module type_testdensevector
