module model_mnp_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none


    


    !> Int. J. Numer. Meth. Fluids 2016; 82:398–416
    !! Dilation-based shock capturing for high-order methods
    !!
    !!  Model Fields:
    !!      - Artifical Viscosity
    !!
    !!  @author Eric M. Wolf
    !!  @date   07/11/2018
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: mnp_artificial_viscosity_t

        real(rk) :: av_constant = 1.5_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type mnp_artificial_viscosity_t
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
        class(mnp_artificial_viscosity_t), intent(inout)   :: self
        
        real(rk)            :: av_constant
        integer             :: unit, msg
        logical             :: file_exists

        namelist /mnp_artificial_viscosity/   av_constant



        call self%set_name('MNP Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Artificial Viscosity')

        !!
        !! Check if input from 'models.nml' is available.
        !!   1: if available, read and set self%mu
        !!   2: if not available, do nothing and mu retains default value
        !!
        !inquire(file='models.nml', exist=file_exists)
        !if (file_exists) then
        !    open(newunit=unit,form='formatted',file='models.nml')
        !    read(unit,nml=mnp_artificial_viscosity,iostat=msg)
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
        class(mnp_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, sensor, av 

        real(rk), dimension(:)  :: h(3)
        real(rk)                :: hmin
        real(rk)                :: Pr_star  = 0.9_rk     
        integer(ik)             :: p, ii
        

        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        if (p == 0) p = 1
        hmin = minval(h)/real(p, rk)
        
        
        sensor = worker%get_field('MNP Shock Sensor', 'value')
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

        av = (1.5_rk*hmin)*wave_speed*sensor
        if (any(ieee_is_nan(av(:)%x_ad_))) print *, 'av is nan'
        if (any(.not. (ieee_is_finite(av(:)%x_ad_)))) print *, 'av is infinity'
        if (any(.not. (ieee_is_finite(av(:)%x_ad_)))) print *, 'density', density(:)%x_ad_
        if (any(.not. (ieee_is_finite(av(:)%x_ad_)))) print *, 'sensor', sensor(:)%x_ad_
        
        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Artificial Viscosity', 'value', av)


    end subroutine compute
    !***************************************************************************************




end module model_mnp_artificial_viscosity
