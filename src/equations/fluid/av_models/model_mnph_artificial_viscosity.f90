module model_mnph_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, THREE, TWO, ONE, ZERO
    use mod_io,             only: h_field_dimension
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
    type, extends(model_t)  :: mnph_artificial_viscosity_t

        real(rk) :: av_constant = 1.5_rk
        logical :: elem_avg = .false.

    contains

        procedure   :: init
        procedure   :: compute

    end type mnph_artificial_viscosity_t
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
        class(mnph_artificial_viscosity_t), intent(inout)   :: self

        integer                                     :: unit, msg
        logical                                     :: file_exists, use_lift, elem_avg

        namelist /av_options/ elem_avg




        call self%set_name('MNPH Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Artificial Viscosity')
        call self%add_model_field('Artificial Viscosity - 1')
        call self%add_model_field('Artificial Viscosity - 2')
        call self%add_model_field('Artificial Viscosity - 3')

        !!
        !! Check if input from 'models.nml' is available.
        !!   1: if available, read and set self%mu
        !!   2: if not available, do nothing and mu retains default value
        !!
        !inquire(file='models.nml', exist=file_exists)
        !if (file_exists) then
        !    open(newunit=unit,form='formatted',file='models.nml')
        !    read(unit,nml=mnph_artificial_viscosity_unsmoothed_ani,iostat=msg)
        !    if (msg == 0) self%av_constant = av_constant
        !    close(unit)
        !end if

        inquire(file='artificial_viscosity.nml', exist=file_exists)
         if (file_exists) then
             open(newunit=unit,form='formatted',file='artificial_viscosity.nml')
             read(unit,nml=av_options,iostat=msg)
             if (msg == 0) self%elem_avg  = elem_avg
             close(unit)
         end if


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
        class(mnph_artificial_viscosity_t), intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, sensor, av, av1, av2, av3, avtemp, temp_av

        real(rk), dimension(:)  :: h(3)
        real(rk)                :: hmin
        real(rk)                :: Pr_star  = 0.9_rk     

        integer(ik)             :: p, ii, nvertex, inode, ivertex, idom, ielem, idom_g, inode_g
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:), nodes(:,:), h_field(:,:), h_scalar(:)
        real(rk)                :: eval_node(3), center(3), radius(3), vert_vals_hmin(8)

        real(rk), allocatable, dimension(:) :: weights, jinv
        idom = worker%element_info%idomain_l
        ielem = worker%element_info%ielement_l

        idom_g = worker%element_info%idomain_g


        h_field = worker%h_smooth()

        select case (trim(h_field_dimension))
            case('2D','2d')
                h_scalar = (h_field(:,1) + h_field(:,2))/TWO
            case('3D','3d')
                h_scalar = (h_field(:,1) + h_field(:,2) + h_field(:,3))/THREE
            case default
                call chidg_signal(FATAL,'mnph_artificial_viscosity: invalid input for h_field_dimension (2D,3D).')

        end select
        !h_scalar = (h_field(:,1)*h_field(:,2)*h_field(:,3))**(ONE/THREE)
        !h_scalar = ONE

        p = worker%solution_order('interior')
        if (p == 0) p = 1
        h_scalar = h_scalar/real(p, rk)
        h_field = h_field/real(p,rk)
        
        
        sensor = worker%get_field('MNPH Shock Sensor', 'value')
        density = worker%get_field('Density','value')

        vel1 = worker%get_field('Momentum-1','value')
        vel2 = worker%get_field('Momentum-2','value')
        vel3 = worker%get_field('Momentum-3','value')

        vel1 = vel1/(density)
        vel2 = vel2/(density)
        vel3 = vel3/(density)

        c = worker%get_field('Pressure', 'value')

        wave_speed = c
        c = (gam*wave_speed/(density))
        c = c*sin_ramp(c, ZERO, ONE)
        wave_speed = sqrt(vel1**TWO+vel2**TWO+vel3**TWO+c)

        av = 1.0_rk*(1.5_rk*h_scalar)*wave_speed*sensor
        
        !avtemp = av
        !av = sin_ramp2(avtemp,0.01_rk*1.5_rk*h_scalar*wave_speed, 1.5_rk*h_scalar*wave_speed)
        !if (self%elem_avg) then
        !    if (worker%interpolation_source == 'element') then
        !        weights = worker%quadrature_weights('element')
        !        jinv    = worker%inverse_jacobian('element')

        !        temp_av = av

        !        temp_av = sum(weights*jinv*av)/sum(weights*jinv)

        !        av = temp_av

        !    else
        !        weights = worker%quadrature_weights('face')
        !        jinv    = worker%inverse_jacobian('face')

        !        temp_av = av

        !        temp_av = sum(weights*jinv*av)/sum(weights*jinv)

        !        av = temp_av

        !    end if

        !end if


        if (any(ieee_is_nan(av(:)%x_ad_))) print *, 'unsmoothed av is nan'
        if (any(ieee_is_nan(av(:)%x_ad_))) print *, worker%interpolation_source
        

        av1 = 1.5_rk*(h_field(:,1))*wave_speed*sensor
        av2 = 1.5_rk*(h_field(:,2))*wave_speed*sensor
        av3 = 1.5_rk*(h_field(:,3))*wave_speed*sensor

        !av1 = 1.0_rk*(1.5_rk*h_field(:,1))*wave_speed*sensor
        !av2 = 1.0_rk*(1.5_rk*h_field(:,2))*wave_speed*sensor
        !av3 = 1.0_rk*(1.5_rk*h_field(:,3))*wave_speed*sensor
        !print *, 'unsmoothed av: ', av(1)%x_ad_

        !! Average to improve robustness
        !if (worker%interpolation_source == 'element') then
        !    weights = worker%quadrature_weights('element')
        !    jinv    = worker%inverse_jacobian('element')

        !    av2 = av
        !    av = sum(weights*jinv*av2)/sum(weights*jinv)

        !else 
        !    weights = worker%quadrature_weights('face')
        !    jinv    = worker%inverse_jacobian('face')

        !    av2 = av
        !    av = sum(weights*jinv*av2)/sum(weights*jinv)


        !end if
       
        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Artificial Viscosity', 'value', av)
        call worker%store_model_field('Artificial Viscosity - 1', 'value', av1)
        call worker%store_model_field('Artificial Viscosity - 2', 'value', av2)
        call worker%store_model_field('Artificial Viscosity - 3', 'value', av3)


    end subroutine compute
    !***************************************************************************************




end module model_mnph_artificial_viscosity
