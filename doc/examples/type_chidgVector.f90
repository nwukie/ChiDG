module type_chidgVector
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: ZERO, TWO
    use type_mesh,                  only: mesh_t
    use type_blockvector
    implicit none





    !> Container stores a blockvector_t for each domain_t
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    !> [chidgvector_t]
    type, public :: chidgVector_t

        type(blockvector_t), allocatable    :: dom(:)

    end type chidgVector_t
    !> [chidgvector_t]
    !-------------------------------------------------------------------------------------










    !------------------------       OPERATORS       --------------------------------------

    public operator (*)
    interface operator(*)
        module procedure mult_real_chidgVector          ! real * chidgVector
        module procedure mult_chidgVector_real          ! chidgVector * real
    end interface


    public operator (/)
    interface operator (/)
        module procedure div_real_chidgVector           ! real / chidgVector
        module procedure div_chidgVector_real           ! chidgVector / real
    end interface


    public operator (-)
    interface operator (-)
        module procedure sub_chidgVector_chidgVector    ! chidgVector - chidgVector
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_chidgVector_chidgVector    ! chidgVector + chidgVector
    end interface
















contains





    !> Subroutine for initializing chidgVector_t
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine initialize(self,mesh)
        class(chidgVector_t),   intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)

        integer(ik) :: ierr, ndomains, idom


        !
        ! Allocate blockvector_t for each mesh
        !
        ndomains = size(mesh)
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization procedure for each blockvector_t
        !
        do idom = 1,ndomains
            call self%dom(idom)%init(mesh(idom))
        end do



    end subroutine initialize











    !>
    !!
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidgVector_t),   intent(inout)   :: self

        integer :: idom


        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do


    end subroutine clear







    !> Compute the L2-Norm of the vector
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function norm(self) result(res)
        class(chidgVector_t),   intent(in)   :: self

        real(rk)    :: res
        integer(ik) :: idom, ielem


        res = ZERO

        !
        ! Loop through domain vectors and compute contribution to vecotr L2-Norm
        !
        do idom = 1,size(self%dom)
            do ielem = 1,size(self%dom(idom)%lvecs)
            
                res = res + sum( self%dom(idom)%lvecs(ielem)%vec ** TWO )

            end do ! ielem
        end do ! idom


        !
        ! Take the square root of the result
        !
        res = sqrt(res)

    end function norm






    !> Dump contents of the vector
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine dump(self)
        class(chidgVector_t),   intent(in)   :: self

        integer(ik) :: idom



        !
        ! Loop through domain vectors and compute contribution to vecotr L2-Norm
        !
        do idom = 1,size(self%dom)
            
            call self%dom(idom)%dump()

        end do ! idom



    end subroutine dump






    !---------------------------------------------------------------------------------------------
    !
    !
    !                              Operator Implementations
    !
    !---------------------------------------------------------------------------------------------



    function mult_real_chidgVector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left * right%dom(idom)
        end do

    end function




    function mult_chidgVector_real(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) * right
        end do

    end function





    function div_real_chidgVector(left,right) result(res)
        real(rk),               intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(right%dom)
            res%dom(idom) = left / right%dom(idom)
        end do

    end function




    function div_chidgVector_real(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        real(rk),               intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(left%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) / right
        end do

    end function





    function add_chidgVector_chidgVector(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))

        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) + right%dom(idom)
        end do

    end function


    function sub_chidgVector_chidgVector(left,right) result(res)
        type(chidgVector_t),    intent(in)  :: left
        type(chidgVector_t),    intent(in)  :: right

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ndom

        ndom = size(right%dom)

        allocate(res%dom(ndom))


        do idom = 1,size(left%dom)
            res%dom(idom) = left%dom(idom) - right%dom(idom)
        end do


    end function





































































end module type_chidgVector
