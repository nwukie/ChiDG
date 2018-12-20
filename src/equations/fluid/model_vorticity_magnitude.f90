module model_vorticity_magnitude
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use mod_fluid,          only: omega
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing vorticity_magnitude
    !!
    !!  Model Fields:
    !!      : vorticity_magnitude-1  
    !!      : vorticity_magnitude-2
    !!      : vorticity_magnitude-3
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/23/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: vorticity_magnitude_t


    contains

        procedure   :: init
        procedure   :: compute

    end type vorticity_magnitude_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/23/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(vorticity_magnitude_t), intent(inout)   :: self

        call self%set_name('Vorticity Magnitude')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Vorticity Magnitude')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(vorticity_magnitude_t),   intent(in)      :: self
        type(chidg_worker_t),           intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) :: &
            vort1, vort2, vort3, z_sq, z

        real(rk) :: epsz 


        !
        ! Interpolate solution to quadrature nodes
        !
        vort1    = worker%get_field('Vorticity-1', 'value')
        vort2    = worker%get_field('Vorticity-2', 'value')
        vort3    = worker%get_field('Vorticity-3', 'value')

        z_sq = vort1
        z       = vort1

        z_sq = vort1**TWO + vort2**TWO + vort3**TWO

        epsz = 1.0e-6_rk

        ! NEED SCALAR ARGUMENT
!        if (z_sq < epsz) then
!            z = HALF*(epsz + z_sq/epsz)
!
!        else
!            z = sqrt(z_sq)
!
!        end if

        where (z_sq < epsz) 
            z = HALF*(epsz + z_sq/epsz)
        elsewhere
            z = sqrt(z_sq)
        end where
        


        call worker%store_model_field('Vorticity Magnitude', 'value', z)



    end subroutine compute
    !***************************************************************************************




end module model_vorticity_magnitude
