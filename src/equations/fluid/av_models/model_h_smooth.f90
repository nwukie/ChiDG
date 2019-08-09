module model_h_smooth
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


    


    !> A model to enable IO-based visualization of the smoothed h field used in AV.
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    03/01/2019 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: h_smooth_t


    contains

        procedure   :: init
        procedure   :: compute

    end type h_smooth_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    03/01/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(h_smooth_t), intent(inout)   :: self
        
       
        call self%set_name('h Smooth')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Smoothed h Field - 1')
        call self%add_model_field('Smoothed h Field - 2')
        call self%add_model_field('Smoothed h Field - 3')


    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    03/01/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(h_smooth_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, h_model

        real(rk), allocatable   :: h_field(:, :)


 
        h_field = worker%h_smooth()

        
        h_model = worker%get_field('Density','value')

        h_model = ZERO
        h_model%x_ad_ = h_field(:,1)
        call worker%store_model_field('Smoothed h Field - 1', 'value', h_model)

        h_model%x_ad_ = h_field(:,2)
        call worker%store_model_field('Smoothed h Field - 2', 'value', h_model)

        h_model%x_ad_ = h_field(:,3)
        call worker%store_model_field('Smoothed h Field - 3', 'value', h_model)



    end subroutine compute
    !***************************************************************************************




end module model_h_smooth
