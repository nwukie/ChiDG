module type_constant_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Constant Viscosity for Laminar Viscosity.
    !!
    !!  Model Fields:
    !!  ------------------------------
    !!      : Laminar Viscosity
    !!
    !!  User-Input:
    !!  ------------------------------
    !!  The user can set a constant viscosity to be used at run-time by
    !!  placing a 'viscosity' namelist with a variable 'mu' in the file 
    !!  'models.nml' in the working directory.
    !!
    !!  Example 'models.nml' contents:
    !!  ------------------------------
    !!
    !!  &viscosity
    !!      mu = 1.e-5
    !!  /
    !!
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/26/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: constant_viscosity_t

        ! Default viscosity value. Can be reset by using 'models.nml'
        real(rk)    :: mu = 1.6e-5_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type constant_viscosity_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/26/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(constant_viscosity_t), intent(inout)   :: self

        real(rk)            :: mu
        integer             :: unit, msg
        logical             :: file_exists

        namelist /viscosity/    mu



        !
        ! Initialize model object
        !
        call self%set_name('Constant Viscosity')
        call self%set_dependency('f(Q-)')
        call self%add_model_field('Laminar Viscosity')


        !
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%mu
        !   2: if not available, do nothing and mu retains default value
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=viscosity,iostat=msg)
            if (msg == 0) self%mu = mu
            close(unit)
        end if


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing a viscosity contribution from Sutherland's Law.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(constant_viscosity_t),    intent(in)      :: self
        type(chidg_worker_t),           intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            viscosity, T

        !
        ! Interpolate solution to quadrature nodes
        !
        T = worker%get_field('Temperature','value')
    

        !
        ! Constant Viscosity for Laminar Viscosity
        !   - initialize derivatives first...
        !
        viscosity = T
        viscosity = self%mu


        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Laminar Viscosity', 'value', viscosity)


    end subroutine compute
    !***************************************************************************************




end module type_constant_viscosity
