module model_mnpha_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use mod_interpolate,           only: interpolate_from_vertices
    use ieee_arithmetic
    implicit none


    


    !> Int. J. Numer. Meth. Fluids 2016; 82:398–416
    !! Dilation-based shock capturing for high-order methods
    !! Presmoothed h
    !!  Model Fields:
    !!      - Smoothed Artifical Viscosity
    !!
    !!  @author Eric M. Wolf
    !!  @date   07/11/2018
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: mnpha_artificial_viscosity_t

        real(rk) :: av_constant = 1.5_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type mnpha_artificial_viscosity_t
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
        class(mnpha_artificial_viscosity_t), intent(inout)   :: self
        
        real(rk)            :: av_constant
        integer             :: unit, msg
        logical             :: file_exists

        namelist /mnpha_artificial_viscosity/   av_constant



        call self%set_name('MNPHA Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Smoothed Anisotropic Artificial Viscosity - 1')
        call self%add_model_field('Smoothed Anisotropic Artificial Viscosity - 2')
        call self%add_model_field('Smoothed Anisotropic Artificial Viscosity - 3')

        !!
        !! Check if input from 'models.nml' is available.
        !!   1: if available, read and set self%mu
        !!   2: if not available, do nothing and mu retains default value
        !!
        !inquire(file='models.nml', exist=file_exists)
        !if (file_exists) then
        !    open(newunit=unit,form='formatted',file='models.nml')
        !    read(unit,nml=mnpha_artificial_viscosity,iostat=msg)
        !    if (msg == 0) self%av_constant = av_constant
        !    close(unit)
        !end if


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
        class(mnpha_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, sensor, av 

        real(rk), dimension(:)  :: h(3)
        real(rk)                :: hmin
        real(rk)                :: Pr_star  = 0.9_rk     

        integer(ik)             :: p, ii, nvertex, inode, ivertex, idom, ielem, idom_g, inode_g, idir
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:), nodes(:,:), h_field(:,:)
        real(rk)                :: eval_node(3), center(3), radius(3), vert_vals_h(8,3)

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


            vert_vals_h(ivertex,:) = worker%solverdata%mesh_size_vertex(inode_g,:)

        end do
        allocate(h_field(size(nodes(:,1)),3))
        do inode = 1, size(nodes(:,1))
            do idir = 1, 3
                h_field(inode,idir)         = interpolate_from_vertices(vert_vals_h(:,idir), nodes(inode,:))
            end do

        end do
 

        p = worker%solution_order('interior')
        if (p == 0) p = 1
        h_field = h_field/real(p, rk)
        
        
        sensor = worker%get_field('MNPHA Shock Sensor', 'value')
        density = worker%get_field('Density','value')

        vel1 = worker%get_field('Momentum-1','value')/density
        vel2 = worker%get_field('Momentum-2','value')/density
        vel3 = worker%get_field('Momentum-3','value')/density

        vel1 = vel1/density
        vel2 = vel2/density
        vel3 = vel3/density

        c = worker%get_field('Pressure', 'value')

        wave_speed = c
        c = (gam*wave_speed/density)
        wave_speed = sqrt(vel1**TWO+vel2**TWO+vel3**TWO+c)

        
        !
        ! Contribute laminar viscosity
        !
        av = (1.5_rk*h_field(:,1))*wave_speed*sensor
        call worker%store_model_field('Smoothed Anisotropic Artificial Viscosity - 1', 'value', av)


        av = (1.5_rk*h_field(:,2))*wave_speed*sensor
        call worker%store_model_field('Smoothed Anisotropic Artificial Viscosity - 2', 'value', av)


        av = (1.5_rk*h_field(:,3))*wave_speed*sensor
        call worker%store_model_field('Smoothed Anisotropic Artificial Viscosity - 3', 'value', av)
    end subroutine compute
    !***************************************************************************************




end module model_mnpha_artificial_viscosity
