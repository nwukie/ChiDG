module operator_chidg_mv
    use mod_kinds,          only: rk, ik
    use type_chidgMatrix,   only: chidgMatrix_t
    use type_chidgVector,   only: chidgVector_t
    implicit none



    public operator(*)
    interface operator(*)
        module procedure MULT_chidgMatrix_chidgVector
    end interface

contains


    !> This function implements the important matrix-vector multiplication 
    !! operation : A*x : for multi-domain configurations, which use the chidg'Container' 
    !! type containers.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    function MULT_chidgMatrix_chidgVector(A,x) results(res)
        type(chidgMatrix_t),    intent(in)  :: A
        type(chidgVector_t),    intent(in)  :: x

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ielem, iblk, iparent


        !
        ! Allocate result and clear
        !
        res = x
        call res%clear




        !
        ! Compute A*x for local blocks
        !
        do idom = 1,size(A%dom)

            !
            ! Routine for local blocks (lblks)
            !
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do iblk = 1,size(A%dom(idom)%lblks,2)
                    
                    if (allocated(A%dom(idom)%lblks(ielem,iblk)%mat)) then
                        iparent = A%dom(idom)%lblks(ielem,iblk)%parent()

                        associate ( resvec => res%dom(idom)%lvecs(ielem)%vec, &
                                    xvec => x%dom(idom)%lvecs(iparent)%vec,   &
                                    Amat => A%dom(idom)%lblks(ielem,iblk)%mat)

                        !res%dom(idom)%lvecs(ielem)%vec = res%dom(idom)%lvecs(ielem)%vec + matmul(A%dom(idom)%lblks(ielem,iblk)%mat,x%dom(idom)%lvecs(iparent)%vec)
                        resvec = resvec + matmul(Amat,xvec)

                        end associate
                    end if

                end do
            end do

            !
            ! TODO: implement routine for off-diagonal, chimera blocks
            !

        end do






    end function


end module operator_chidg_mv
