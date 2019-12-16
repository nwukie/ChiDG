module mod_functional_operators
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_chidg_worker,      only: chidg_worker_t
    use mod_constants,          only: HALF, ONE, ZERO
    use DNAD_D
    implicit none


    !>  Operations available:
    !!
    !!      1) integrate_surface_mass_weighted              mass weighted integral over a face
    !!      2) integrate_surface_area_weighted              area weighted integral over a face
    !!      3) integrate_volume                             integral over a volume
    !!      4) integrate_surface                            integral over a face
    !!      5) average                                      average of the values of a given 
    !!                                                      primary or compound variable at quadrature
    !!                                                      points of a face or element
    !!
    !-----------------------------------------------------------------------------------------


contains


    !>  Function that computes the mass weighted integral over a face of given input
    !!  quantity.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   06/26/2017
    !!
    !!  param[in]               worker      
    !!  param[in], optional,    integrand   quantity to integrate
    !!
    !!
    !-------------------------------------------------------------------------------
    function  integrate_surface_mass_weighted(worker,integrand) result(integral)
        type(chidg_worker_t),   intent(inout)           :: worker
        type(AD_D), allocatable,intent(in),  optional   :: integrand(:)

        type(AD_D), allocatable, dimension(:) :: density, mom_1, mom_2, mom_3, &
                                                 u_1, u_2, u_3, vn
        real(rk),   allocatable, dimension(:) :: unorm_1, unorm_2, unorm_3,    &
                                                 norm_1, norm_2, norm_3, r,    &
                                                 areas
        real(rk),   allocatable, dimension(:)   :: weights

        type(AD_D) :: integral
        
        ! Get primary fields
        density = worker%get_field('Density'   ,'value','face interior')
        mom_1   = worker%get_field('Momentum-1','value','face interior')
        mom_2   = worker%get_field('Momentum-2','value','face interior')
        mom_3   = worker%get_field('Momentum-3','value','face interior')
        
        ! Consider cylindrical coordinates
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom_2 = mom_2/r
        end if

        ! Compute velocity components
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density

        ! Compute unit normals: negate sign so that if the flow
        ! is entering the element the mass flux will be positive
        ! and viceversa.
        unorm_1 = -worker%unit_normal(1)
        unorm_2 = -worker%unit_normal(2)
        unorm_3 = -worker%unit_normal(3)

        ! Compute unit normals and normals
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        ! Compute normal velocity
        vn = u_1*unorm_1 + u_2*unorm_2 + u_3*unorm_3

        ! Retrieve weights
        weights = worker%quadrature_weights('face')

        ! Compute differential areas
        areas = sqrt( norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)

        ! Compute mass-weighted integral 
        if (present(integrand)) then
            integral = sum( integrand * density * vn * weights * areas )
        else
            integral = sum( density * vn * weights * areas )
        end if

    end function  integrate_surface_mass_weighted
    !*******************************************************************************






!    !>  Function that computes the area weighted integral over a face of given input
!    !!  quantity.
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   06/26/2017
!    !!
!    !!  param[in]               worker      
!    !!  param[in], optional,    integrand   quantity to integrate
!    !!
!    !-------------------------------------------------------------------------------
!    function  integrate_surface_area_weighted(worker,integrand) result(integral)
!        type(chidg_worker_t),   intent(inout)           :: worker
!        type(AD_D), allocatable,intent(in),  optional   :: integrand(:)
!
!        type(AD_D)                              :: integral
!
!        real(rk),   allocatable, dimension(:)   :: weights
!        type(AD_D), allocatable, dimension(:)   :: norm_1, norm_2, norm_3, r,    &
!                                                   areas
!
!        !
!        ! Compute normals
!        !
!        norm_1  = worker%normal(1)
!        norm_2  = worker%normal(2)
!        norm_3  = worker%normal(3)
!
!
!        !
!        ! Retrieve weights
!        !
!        weights = worker%quadrature_weights('face')
!
!        !
!        ! Compute differential areas
!        !
!        areas = sqrt( norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)
!
!
!        !
!        ! Compute area-weighted integral 
!        !
!        if (present(integrand)) then
!            integral = sum( integrand * weights * areas )
!        else
!            !integral = sum( weights * areas )
!            integral = sum( abs(weights * areas) )
!        end if
!
!
!    end function  integrate_surface_area_weighted
!    !*******************************************************************************






    !>  Function that computes the integral over a face of given input quantity.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   06/26/2017
    !!
    !!  param[in]               worker      
    !!  param[in], optional,    integrand   quantity to integrate
    !!
    !-------------------------------------------------------------------------------
    function integrate_surface(worker,integrand) result(integral)
        type(chidg_worker_t),   intent(inout)           :: worker
        type(AD_D), allocatable,intent(in),  optional   :: integrand(:)

        type(AD_D)                              :: integral

        real(rk),   allocatable, dimension(:)   :: weights
        !type(AD_D), allocatable, dimension(:)   :: norm_1, norm_2, norm_3, dareas
        real(rk),   allocatable, dimension(:)   :: norm_1, norm_2, norm_3, dareas
        
        ! Compute normals
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        ! Retrieve weights
        weights = worker%quadrature_weights('face')

        ! Compute differential areas
        dareas = sqrt( norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)
        
        ! Compute area-weighted integral 
        if (present(integrand)) then
            integral = sum( integrand * weights * dareas )
        else
            integral = sum( abs(weights * dareas) )
        end if

    end function  integrate_surface
    !*******************************************************************************





    !>  Function that computes the mass integral over an element of given input
    !!  quantity.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   06/26/2017
    !!
    !!  param[in]               worker      
    !!  param[in], optional,    integrand   quantity to integrate
    !!
    !-------------------------------------------------------------------------------
    function integrate_volume(worker,integrand) result(integral)
        type(chidg_worker_t),       intent(inout)       :: worker
        type(AD_D), allocatable,    intent(in)          :: integrand(:)

        real(rk),   allocatable, dimension(:)   :: weights
        !type(AD_D), allocatable, dimension(:)   :: jinv
        real(rk),   allocatable, dimension(:)   :: jinv
        integer(ik)                             :: idomain_l, ielement_l
        type(AD_D)                              :: integral

        ! Retrieve weights and jinv
        jinv    = worker%inverse_jacobian('element')
        weights = worker%quadrature_weights('element')

        ! Compute volume integral 
        integral = sum( integrand * weights * jinv)
        

    end function integrate_volume
    !*******************************************************************************






    !>  Compute the average of the value at the quadrature points for a given face 
    !!  or element 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   01/22/2018
    !!
    !!  param[in]               worker      
    !!  param[in], optional,    qg_values   values of the function at quadrature nodes
    !!
    !-------------------------------------------------------------------------------
    function average(worker,gq_values) result(ave)
        type(chidg_worker_t),       intent(inout)       :: worker
        type(AD_D), allocatable,    intent(in)          :: gq_values(:)

        type(AD_D)      :: ave
        type(AD_D)      :: temp_sum
        integer(ik)     :: igq
        real(rk)        :: counter

        ! Compute average 
        temp_sum = gq_values(1)
        temp_sum = ZERO
        counter  = ZERO
        do igq = 1, size(gq_values)
            temp_sum = temp_sum + gq_values(igq)
            counter = counter + ONE
        end do

        ave = temp_sum / counter


    end function average
    !*******************************************************************************





    !>  Compute the mass-weighted average of the value at the quadrature points for a given face 
    !!  or element 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/24/2019
    !!
    !!  param[in]               worker      
    !!  param[in], optional,    qg_values   values of the function at quadrature nodes
    !!
    !-------------------------------------------------------------------------------
    function  mass_weighted_average(worker,gq_values) result(ave)
        type(chidg_worker_t),       intent(inout)               :: worker
        type(AD_D), allocatable,    intent(in),     optional    :: gq_values(:)

        type(AD_D)      :: ave
        type(AD_D)      :: temp_sum
        integer(ik)     :: igq
        real(rk)        :: mass_sum
        
        type(AD_D), allocatable, dimension(:) :: density, mom_1, mom_2, mom_3, &
                                                 u_1, u_2, u_3, vn, gq_vals   
        real(rk),   allocatable, dimension(:) :: unorm_1, unorm_2, unorm_3,    &
                                                 norm_1, norm_2, norm_3, r,    &
                                                 areas
        real(rk),   allocatable, dimension(:) :: weights, m_vals
        

        ! Get primary fields
        density = worker%get_field('Density'   ,'value','face interior')
        mom_1   = worker%get_field('Momentum-1','value','face interior')
        mom_2   = worker%get_field('Momentum-2','value','face interior')
        mom_3   = worker%get_field('Momentum-3','value','face interior')
        

        ! Consider cylindrical coordinates
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom_2 = mom_2/r
        end if


        ! Compute velocity components
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density


        ! Compute unit normals: negate sign so that if the flow
        ! is entering the element the mass flux will be positive
        ! and viceversa.
        unorm_1 = -worker%unit_normal(1)
        unorm_2 = -worker%unit_normal(2)
        unorm_3 = -worker%unit_normal(3)

         
        ! Compute unit normals and normals
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

         
        ! Compute normal velocity
        vn = u_1*norm_1 + u_2*norm_2 + u_3*norm_3

         
        ! Compute mass-weighted values at quadrature nodes 
        if (present(gq_values)) then
            gq_vals = gq_values * density%x_ad_ * vn%x_ad_ 
        else
            gq_vals = density * vn * density%x_ad_ * vn%x_ad_ 
        end if
       

        ! Compute mass flux as GQ points that will be use as weights 
        m_vals  = density%x_ad_ * vn%x_ad_ 
        
        
        ! Compute average 
        temp_sum = gq_vals(1)
        temp_sum = ZERO
        mass_sum = ZERO
        do igq = 1, size(gq_vals)
            temp_sum = temp_sum + gq_vals(igq)
            mass_sum = mass_sum + m_vals(igq)
        end do

        ave = temp_sum / mass_sum


    end function mass_weighted_average
    !*******************************************************************************



end module mod_functional_operators
