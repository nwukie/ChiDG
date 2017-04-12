module operator_domain_dot
    use mod_kinds,          only: rk, rdouble, ik
    use mod_constants,      only: ZERO
    use type_domain_vector, only: domain_vector_t

    implicit none


    ! Import BLAS dot product
    external DDOT
    real(rdouble) :: DDOT

    


contains





    !> Compute vector-vector dot product from two domain_vector_t types.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    function domain_dot(a,b) result(res)
        type(domain_vector_t),    intent(in)  :: a
        type(domain_vector_t),    intent(in)  :: b

        real(rk)                :: res
        integer(ik)             :: ielem
        integer(4)              :: vecsize

        res = ZERO

        !
        ! Compute vector dot-product
        !
        do ielem = 1,size(a%vecs)

            res = res + dot_product(a%vecs(ielem)%vec, b%vecs(ielem)%vec)

            !vecsize = size(a%vecs(ielem)%vec)
            !res = res + ddot(vecsize, a%vecs(ielem)%vec, 1, b%vecs(ielem)%vec, 1)

        end do  ! ielem
        

    end function domain_dot
    !*****************************************************************************










end module operator_domain_dot

