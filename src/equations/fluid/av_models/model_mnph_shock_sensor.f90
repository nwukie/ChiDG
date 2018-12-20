module model_mnph_shock_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, PI
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use mod_interpolate,           only: interpolate_from_vertices
    use DNAD_D
    use ieee_arithmetic 
    implicit none


    


    !> Int. J. Numer. Meth. Fluids 2016; 82:398–416
    !! Dilation-based shock capturing for high-order methods
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: mnph_shock_sensor_t


    contains

        procedure   :: init
        procedure   :: compute

    end type mnph_shock_sensor_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(mnph_shock_sensor_t), intent(inout)   :: self

        call self%set_name('MNPH Shock Sensor')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('MNPH Shock Sensor')
    end subroutine init
    !***************************************************************************************






    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(mnph_shock_sensor_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: order
        type(AD_D), allocatable,    dimension(:) ::         &
            div_vel, D, D_smooth, c_star, temp1, temp2

        real(rk),   allocatable,    dimension(:) :: h_ref
        real(rk) :: alpha, beta

        integer(ik)             :: p, ii, nvertex, inode, ivertex, idom, ielem, idom_g, inode_g
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:), nodes(:,:), h_field(:)
        real(rk)                :: eval_node(3), center(3), radius(3), vert_vals_hmin(8)


        idom = worker%element_info%idomain_l
        ielem = worker%element_info%ielement_l

        idom_g = worker%element_info%idomain_g


        if (worker%interpolation_source == 'element') then
            nodes = worker%mesh%domain(idom)%elems(ielem)%basis_s%nodes_elem_

        else

            nodes = worker%mesh%domain(idom)%elems(ielem)%basis_s%nodes_face_(:,:, worker%iface)

        end if

                
        do ivertex = 1, 8

            inode_g = sum(worker%solverdata%nnodes_per_domain(1:idom_g-1)) + worker%mesh%domain(idom)%elems(ielem)%vertex_indices(ivertex)


            vert_vals_hmin(ivertex) = worker%solverdata%avg_mesh_size_vertex(inode_g)

        end do

        allocate(h_field(size(nodes(:,1))))
        do inode = 1, size(nodes(:,1))
            h_field(inode)         = interpolate_from_vertices(vert_vals_hmin, nodes(inode,:))

        end do
 
        div_vel = worker%get_field('Velocity Divergence', 'value')
        c_star = worker%get_field('Critical Sound Speed', 'value')
        order = worker%solution_order('interior')

        if (order == 0 ) order = 1

        D = -(1.5_rk*h_field/real(order, rk))*div_vel/(c_star + 1.0e-15_rk)
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'D is infinite'
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'c_star', c_star(:)%x_ad_

        ! The reference gives the value alpha = 1.0e4, but this has resulted
        ! in values of Infinity for the exponential below.
        ! Use a smaller value of alpha, or a different soft-max variant?
        alpha = 1.0e3_rk !1.0e4_rk
        beta = 0.01_rk

        D_smooth = D
        !D_smooth = log(ONE +  exp(alpha*(D-beta)))/alpha
        D_smooth = max(beta, D)
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) D_smooth = max(beta, D)
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D_smooth is infinite'
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D', D(:)%x_ad_
        !print *, 'D_smooth', D_smooth(:)%x_ad_

        call worker%store_model_field('MNPH Shock Sensor', 'value', D_smooth)

    end subroutine compute
    !***************************************************************************************




end module model_mnph_shock_sensor
