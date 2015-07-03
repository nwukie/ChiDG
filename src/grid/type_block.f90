module type_block
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: HALF, TWO, ONE, DIAG, NFACES, &
                                  XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use mod_io,             only: output_res

    use type_mesh,          only: mesh_t
    use type_point,         only: point_t
    implicit none

    private

    !====================================================
    type, public :: block_t
        type(mesh_t)                    :: mesh

    contains
        procedure       :: init
        procedure       :: write

        final           :: destructor

    end type block_t
    !=====================================================
contains
    
    subroutine init(self,eqnset,nterms_sol,nterms_mesh,points)
        class(block_t),         intent(inout)       :: self
        class(equationset_t),   intent(in), target  :: eqnset
        type(point_t),          intent(in)          :: points(:,:,:)
        integer(ik),            intent(in)          :: nterms_sol
        integer(ik),            intent(in)          :: nterms_mesh

        ! Set domain equation set
        self%eqnset => eqnset

        ! Initialize mesh
        call self%mesh%init(eqnset%neqns,nterms_sol,nterms_mesh,points)

    end subroutine



    !===============================================
    !
    !   Destructor
    !
    !===============================================
    subroutine destructor(self)
        type(block_t), intent(in) :: self
    end subroutine

end module type_block
