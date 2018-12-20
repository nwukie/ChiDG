module model_vertex_smoothed_mnp_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use mod_interpolate,    only: interpolate_from_vertices
    use DNAD_D
    use ieee_arithmetic
    implicit none


    


    !>  Artificial viscosity model based on
    !!  A physics-based shock capturing method for unsteady laminar and turbulent flows
    !!      Fernandez et al, 2018, AIAA SciTech Forum
    !!
    !!  Model Fields:
    !!      - Artifical Bulk Viscosity
    !!      - Artifical Shear Viscosity
    !!      - Artifical Thermal Conductivity 
    !!
    !!  @author Eric M. Wolf
    !!  @date   07/11/2018
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: vertex_smoothed_mnp_artificial_viscosity_t

    contains

        procedure   :: init
        procedure   :: compute

    end type vertex_smoothed_mnp_artificial_viscosity_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(vertex_smoothed_mnp_artificial_viscosity_t), intent(inout)   :: self

        call self%set_name('Vertex Smoothed MNP Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Smoothed Artificial Viscosity')

    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(vertex_smoothed_mnp_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_viscosity

        integer(ik)             :: p, ii, nvertex, inode, ivertex, idom, ielem, idom_g, inode_g
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:), nodes(:,:)
        real(rk)                :: eval_node(3), center(3), radius(3), vert_vals_bv(8)
        

        
        density = worker%get_field('Density','value')

        art_viscosity = ZERO*density

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


            vert_vals_bv(ivertex) = worker%solverdata%vertex_artificial_bulk_viscosity(inode_g)

        end do
        !print *, 'vert vals', vert_vals_bv

        if (any(ieee_is_nan(vert_vals_bv))) print *, 'av vert vals is nan'
        if (any(.not. (ieee_is_finite(vert_vals_bv)))) print *, 'av vert vals is infinite'

        do inode = 1, size(nodes(:,1))
            art_viscosity(inode)%x_ad_         = interpolate_from_vertices(vert_vals_bv, nodes(inode,:))

        end do
               
        if (any(ieee_is_nan(art_viscosity(:)%x_ad_))) print *, 'smoothed av is nan'
        !if (any(ieee_is_nan(art_viscosity(:)%x_ad_))) print *, 'vert vals', vert_vals_bv
        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Smoothed Artificial Viscosity', 'value', art_viscosity)


    end subroutine compute
    !***************************************************************************************




end module model_vertex_smoothed_mnp_artificial_viscosity
