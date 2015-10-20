module type_chidgVector
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use type_blockvector,           only: blockvector_t
    use type_mesh,                  only: mesh_t
    implicit none





    !> Container stores a blockvector_t for each domain_t
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, public :: chidgVector_t

        type(blockvector_t), allocatable    :: dom(:)

    contains
        !> Initializers
        generic, public     :: init => initialize
        procedure, private  :: initialize

        !> Modifiers
        procedure, public   :: clear

        !> Interogators
        !procedure, public   :: norm
        !procedure, public   :: nentries
        !procedure, public   :: ndomains

        !final               :: destructor

    end type chidgVector_t
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
















end module type_chidgVector
