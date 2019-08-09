module model_mnph_shock_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO,HALF, ONE, TWO, THREE, PI
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

        logical :: use_lift = .true.
        logical :: elem_avg = .false.
        logical :: use_pressure_jump_indicator = .false.


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
        integer                                     :: unit, msg
        logical                                     :: file_exists, use_lift, elem_avg, use_pressure_jump_indicator

        namelist /shock_sensor_options/ use_lift, elem_avg, use_pressure_jump_indicator

        call self%set_name('MNPH Shock Sensor')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('MNPH Shock Sensor')

        inquire(file='artificial_viscosity.nml', exist=file_exists)
         if (file_exists) then
             open(newunit=unit,form='formatted',file='artificial_viscosity.nml')
             read(unit,nml=shock_sensor_options,iostat=msg)
             if (msg == 0) self%use_lift                        = use_lift
             if (msg == 0) self%elem_avg                        = elem_avg
             if (msg == 0) self%use_pressure_jump_indicator     = use_pressure_jump_indicator
             close(unit)
         end if

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
            div_vel, D, D_smooth, c_star, temp1, temp2, temp_sensor, pjump

        real(rk),   allocatable,    dimension(:) :: h_ref
        real(rk) :: alpha, beta

        integer(ik)             :: p, ii, nvertex, inode, ivertex, idom, ielem, idom_g, inode_g
        type(AD_D)              :: pjump1
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:), nodes(:,:), h_field(:, :), h_scalar(:)
        real(rk)                :: eval_node(3), center(3), radius(3), vert_vals_hmin(8)

        real(rk), allocatable, dimension(:) :: weights, jinv


        idom = worker%element_info%idomain_l
        ielem = worker%element_info%ielement_l

        idom_g = worker%element_info%idomain_g


        !if (worker%interpolation_source == 'element') then
        !    nodes = worker%mesh%domain(idom)%elems(ielem)%basis_s%nodes_elem_

        !else

        !    nodes = worker%mesh%domain(idom)%elems(ielem)%basis_s%nodes_face_(:,:, worker%iface)

        !end if

        !        
        !do ivertex = 1, 8

        !    inode_g = sum(worker%solverdata%nnodes_per_domain(1:idom_g-1)) + worker%mesh%domain(idom)%elems(ielem)%vertex_indices(ivertex)


        !    vert_vals_hmin(ivertex) = worker%solverdata%avg_mesh_size_vertex(inode_g)

        !end do

        !allocate(h_field(size(nodes(:,1))))
        !do inode = 1, size(nodes(:,1))
        !    h_field(inode)         = interpolate_from_vertices(vert_vals_hmin, nodes(inode,:))

        !end do
 
        h_field = worker%h_smooth()
        h_scalar = (h_field(:,1) + h_field(:,2) + h_field(:,3))/THREE
        !h_scalar = (h_field(:,1)*h_field(:,2)*h_field(:,3))**(ONE/THREE)
        !h_scalar = (h_field(:,1) + h_field(:,2))/TWO
        !print *, 'h_scalar: ', h_scalar
        !h_scalar = ONE
        if (self%use_lift) then
            div_vel = worker%get_field('Velocity Divergence', 'value')
        else
            div_vel = worker%get_field('Velocity Divergence No Lift', 'value')
        end if
        c_star = worker%get_field('Critical Sound Speed', 'value')
        order = worker%solution_order('interior')

        if (order == 0 ) order = 1

        D = -(1.5_rk*h_scalar/real(order, rk))*div_vel/(c_star + 1.0e-11_rk)
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'D is infinite'
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'c_star', c_star(:)%x_ad_

        ! The reference gives the value alpha = 1.0e4, but this has resulted
        ! in values of Infinity for the exponential below.
        ! Use a smaller value of alpha, or a different soft-max variant?
        alpha = 1.0e3_rk !1.0e4_rk
        beta = 0.01_rk

        D_smooth = D
        !D_smooth = log(ONE +  exp(alpha*(D-beta)))/alpha
        D_smooth = D*sin_ramp(D-beta, ZERO, 0.1_rk) !+ beta
        !if (self%elem_avg) then
        !    if (worker%interpolation_source == 'element') then
        !        weights = worker%quadrature_weights('element')
        !        jinv    = worker%inverse_jacobian('element')

        !        temp_sensor = D_smooth 

        !        temp_sensor = sum(weights*jinv*D_smooth)/sum(weights*jinv)

        !        D_smooth = temp_sensor

        !    else
        !        weights = worker%quadrature_weights('face')
        !        jinv    = worker%inverse_jacobian('face')

        !        temp_sensor = D_smooth 

        !        temp_sensor = sum(weights*jinv*D_smooth)/sum(weights*jinv)

        !        D_smooth = temp_sensor

        !    end if

        !end if

        !if (self%use_pressure_jump_indicator) then

        !    if (worker%interpolation_source == 'element') then
        !        !pjump = worker%get_pressure_jump_indicator()
        !        pjump = worker%get_pressure_jump_shock_sensor()

        !        D_smooth = D_smooth+pjump

        !        !D_smooth = D*sin_ramp(D, ZERO, 0.1_rk) !+ beta
        !        !D_smooth = pjump
        !    else
        !        pjump = worker%get_pressure_jump_shock_sensor()
        !        pjump1 = pjump(1)
        !        do inode = 1, size(D_smooth)
        !            D_smooth(inode) = D_smooth(inode)+pjump1
        !        end do
        !    end if
        !end if

        !D_smooth = D*sin_ramp(D-beta, ZERO, 0.1_rk) !+ beta
        !D_smooth = max(beta, D)
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) D_smooth = max(beta, D)
        if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D_smooth is infinite'
        if (any( (ieee_is_nan(D_smooth(:)%x_ad_)))) print *, 'D_smooth is nan'
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D', D(:)%x_ad_
        !print *, 'D_smooth', D_smooth(:)%x_ad_

        call worker%store_model_field('MNPH Shock Sensor', 'value', D_smooth)

    end subroutine compute
    !***************************************************************************************




end module model_mnph_shock_sensor
