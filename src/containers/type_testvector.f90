module type_testvector
    use type_testdensevector
    implicit none


    type, public :: testvector_t

        type(testdensevector_t), allocatable :: lvecs(:)

    end type testvector_t



    !-----------------      OPERATOR INTERFACE       ------------------------
    public operator (*)
    interface operator (*)
        module procedure mult_real_bv   ! real * testvector
    end interface



contains


    !-----------        OPERATOR IMPLEMENTATIONS        ---------------------
    function mult_real_bv(left,right) result(res)
        real,                   intent(in)  :: left
        type(testvector_t),     intent(in)  :: right
        type(testvector_t)      :: res
        
        res%lvecs = left * right%lvecs

    end function





end module type_testvector
