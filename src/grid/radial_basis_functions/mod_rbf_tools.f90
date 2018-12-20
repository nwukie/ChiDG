module mod_rbf_tools
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO
    implicit none

contains


    function compute_radius_norm(eval_node, support_node, support_radius) result(radius_norm)
        real(rk),               intent(in)      :: eval_node(3)
        real(rk),               intent(in)      :: support_node(3)
        real(rk),               intent(in)      :: support_radius(3)
        real(rk) :: radius_norm
 
        radius_norm = sqrt(&
            ((eval_node(1)-support_node(1))/support_radius(1))**TWO+&
            ((eval_node(2)-support_node(2))/support_radius(2))**TWO+&
            ((eval_node(3)-support_node(3))/support_radius(3))**TWO)


    end function compute_radius_norm

    function compute_radius_norm_grad(eval_node, support_node, support_radius) result(radius_norm_grad)
        real(rk),               intent(in)      :: eval_node(3)
        real(rk),               intent(in)      :: support_node(3)
        real(rk),               intent(in)      :: support_radius(3)
        real(rk) :: radius_norm_grad(3)
        real(rk) :: radius_norm
 
        radius_norm = compute_radius_norm(eval_node, support_node, support_radius)

        if (radius_norm > 1.0e-16_rk) then
            radius_norm_grad = ONE/sqrt(radius_norm)*(eval_node-support_node)/support_radius**TWO
        else

            ! The RBF should have derivative zero here.
            radius_norm_grad = ZERO 
        end if


    end function compute_radius_norm_grad

    function node_dist(node1, node2) result(ndist)
        real(rk), intent(in) :: node1(3), node2(3)
        real(rk) :: ndist

        integer(ik) :: idir

        ndist = ZERO
        do idir = 1,3
            ndist = ndist+(node1(idir)-node2(idir))**TWO
        end do
        ndist = sqrt(ndist)

    end function node_dist



end module mod_rbf_tools
