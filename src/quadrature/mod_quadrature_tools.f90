module mod_quadrature_tools
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    implicit none




contains







    !>  Compute number of quadrature nodes in one dimension, from a set of volume nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------
    function compute_nnodes1d_volume(spacedim,nnodes) result(nnodes1d)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nnodes

        integer(ik) :: nnodes1d

        nnodes1d = 0

        if ( spacedim == 3 ) then

            do while (nnodes1d*nnodes1d*nnodes1d /= nnodes)
                nnodes1d = nnodes1d + 1
            end do
            if (nnodes1d*nnodes1d*nnodes1d > nnodes) call chidg_signal(FATAL, "Error in 1d volume quadrature node count")

        else if ( spacedim == 2 ) then

            do while (nnodes1d*nnodes1d /= nnodes)
                nnodes1d = nnodes1d + 1
            end do
            if (nnodes1d*nnodes1d > nnodes) call chidg_signal(FATAL, "Error in 1d volume quadrature node count")

        else
            call chidg_signal(FATAL,'compute_nnodes1d_volume: Invalid value for spatial dimension - spacedim.')

        end if


    end function compute_nnodes1d_volume
    !*********************************************************************************************************











    !>  Compute number of quadrature nodes in one dimension, from a set of volume nodes.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------
    function compute_nnodes1d_face(spacedim,nnodes) result(nnodes1d)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nnodes

        integer(ik) :: nnodes1d

        nnodes1d = 0

        if ( spacedim == 3 ) then

            do while (nnodes1d*nnodes1d /= nnodes)
                nnodes1d = nnodes1d + 1
            end do
            if (nnodes1d*nnodes1d > nnodes) call chidg_signal(FATAL, "Error in 1d face quadrature node count")

        else if ( spacedim == 2 ) then
            nnodes1d = nnodes
            if (nnodes1d > nnodes) call chidg_signal(FATAL, "Error in 1d face quadrature node count")

        else
            call chidg_signal(FATAL,'compute_nnodes1d_face: Invalid value for spatial dimension - spacedim.')

        end if


    end function compute_nnodes1d_face
    !*********************************************************************************************************















end module mod_quadrature_tools
