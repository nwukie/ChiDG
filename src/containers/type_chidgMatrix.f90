module type_chidgMatrix
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_blockmatrix,       only: blockmatrix_t
    use type_mesh,              only: mesh_t
    use DNAD_D
    implicit none




    !> ChiDG matrix type. Contains an array of blockmatrix_t types, each corresponding to a domain.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    type, public :: chidgMatrix_t

        type(blockmatrix_t), allocatable    :: dom(:)          !< Array of block-matrices. One for each domain



    contains
        !> Initializers
        generic,    public  :: init => initialize
        procedure,  private :: initialize           !< ChiDGMatrix initialization

        !> Setters
        procedure   :: store                        !< Store linearization data
        procedure   :: clear                        !< Zero matrix-values


        !procedure   :: build                       !< Build full-matrix representation of the block-sparse matrix

        final       :: destructor

    end type chidgMatrix_t
    !-----------------------------------------------------------------------------------------------------------



    private
contains




    !> Subroutine for initializing chidgMatrix_t
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  domains     Array of domain_t instances
    !!  
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine initialize(self,mesh,mtype)
        class(chidgMatrix_t),   intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh(:)
        character(*),           intent(in)      :: mtype

        integer(ik) :: ierr, ndomains, imesh, nmesh


        !
        ! Allocate blockmatrix_t for each domain
        !
        ndomains = size(mesh)
        allocate(self%dom(ndomains), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization procedure for each blockmatrix_t
        !
        do imesh = 1,nmesh
            call self%dom(imesh)%init(mesh(imesh),mtype)
        end do



    end subroutine
    !-----------------------------------------------------------------------------------------------------------








    !> Procedure for storing linearization information
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  idom        Domain index for storing the linearization
    !!  @param[in]  ielem       Element index for which the linearization was computed
    !!  @param[in]  iblk        Index of the block for the linearization of the given elemen
    !!  @param[in]  ivar        Index of the variable, for which the linearization was computed
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine store(self, integral, idom, ielem, iblk, ivar)
        class(chidgMatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        integer(ik),            intent(in)      :: idom, ielem, iblk, ivar

        !
        ! Store linearization in associated domain blockmatrix_t
        !
        call self%dom(idom)%store(integral,ielem,iblk,ivar)

    end subroutine
    !-----------------------------------------------------------------------------------------------------------











    !> Set all ChiDGMatrix matrix-values to zero
    !!
    !!  @author Nathan A. Wukie
    !! 
    !! 
    !! 
    !----------------------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidgMatrix_t),   intent(inout)   :: self

        integer(ik) :: idom
    

        !
        ! Call blockmatrix_t%clear() on all matrices
        !
        do idom = 1,size(self%dom)
           call self%dom(idom)%clear() 
        end do
    
    
    end subroutine 
    !-----------------------------------------------------------------------------------------------------------











    !> ChiDGMatrix destructor.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-----------------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidgMatrix_t),    intent(inout)   :: self

    end subroutine
    !-----------------------------------------------------------------------------------------------------------



end module type_chidgMatrix
