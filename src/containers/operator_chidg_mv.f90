module operator_chidg_mv
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_chidgMatrix,   only: chidgMatrix_t
    use type_chidgVector
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
    function MULT_chidgMatrix_chidgVector(A,x) result(res)
        type(chidgMatrix_t),    intent(in)  :: A
        type(chidgVector_t),    intent(in)  :: x

        type(chidgVector_t) :: res
        integer(ik)         :: idom, ielem, iblk
        integer(ik)         :: eparent, dparent
        logical             :: nonconforming = .false.


        !
        ! Allocate result and clear
        !
        res = x
        call res%clear




        !
        ! Compute A*x for global matrix-vector product
        !
        do idom = 1,size(A%dom)

            !
            ! Routine for local blocks (lblks)
            !
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do iblk = 1,size(A%dom(idom)%lblks,2)
                    
                    if (allocated(A%dom(idom)%lblks(ielem,iblk)%mat)) then
                        eparent = A%dom(idom)%lblks(ielem,iblk)%eparent()

                        associate ( resvec => res%dom(idom)%lvecs(ielem)%vec, &
                                    xvec => x%dom(idom)%lvecs(eparent)%vec,   &
                                    Amat => A%dom(idom)%lblks(ielem,iblk)%mat)

                            resvec = resvec + matmul(Amat,xvec)

                        end associate
                    end if

                end do
            end do



            !
            ! Routine for off-diagonal, chimera blocks
            !
            if (allocated(A%dom(idom)%chiblks)) then
                do ielem = 1,size(A%dom(idom)%chiblks,1)
                    do iblk = 1,size(A%dom(idom)%chiblks,2)


                        if (allocated(A%dom(idom)%chiblks(ielem,iblk)%mat)) then
                            dparent = A%dom(idom)%chiblks(ielem,iblk)%dparent()
                            eparent = A%dom(idom)%chiblks(ielem,iblk)%eparent()

                            associate ( resvec => res%dom(idom)%lvecs(ielem)%vec, &
                                        xvec => x%dom(dparent)%lvecs(eparent)%vec, &
                                        Amat => A%dom(idom)%chiblks(ielem,iblk)%mat) 

                                !
                                ! Test matrix vector sizes
                                !
                                nonconforming = ( size(Amat,2) /= size(xvec) )
                                if (nonconforming) call chidg_signal(FATAL,"operator_chidg_mv: nonconforming Chimera m-v operation")

                                resvec = resvec + matmul(Amat,xvec)

                            end associate
                        end if

                    end do ! iblk
                end do ! ielem
            end if

        end do ! idom






    end function


end module operator_chidg_mv
