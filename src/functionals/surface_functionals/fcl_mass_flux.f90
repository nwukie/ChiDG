module fcl_mass_flux
#include <messenger.h>

    use mod_kinds,                  only: ik, rk
    use type_evaluator,             only: evaluator_t
    use type_chidg_worker,          only: chidg_worker_t
    use mod_fluid,                  only: Rgas,cp, gam
    use mod_constants,              only: HALF, ONE, TWO, PI
    use type_functional_cache,      only: functional_cache_t
    use type_chidg_vector,          only: chidg_vector_t
    use mod_functional_operators

    use DNAD_D


    implicit none

    !>  This functional computes the mass flux through a surface
    !!      - reference geom
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: mass_flux_t

        ! Quantities needed for the computation
        ! Used to check completeness of information provided
        integer(ik)     :: min_ref_geom = 1
        integer(ik)     :: min_aux_geom = 0
        integer(ik)     :: max_ref_geom = 100
        integer(ik)     :: max_aux_geom = 0

    contains

        procedure, public   :: init
        procedure, public   :: check
        procedure, public   :: compute_functional
        procedure, public   :: store_value
        procedure, public   :: store_deriv

    end type mass_flux_t
    !******************************************************************************************
    


contains



    !>  Initialize the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(mass_flux_t), intent(inout)   :: self

        call self%set_name("Mass Flux")
        call self%set_eval_type("Functional")
        call self%set_int_type("FACE INTEGRAL")

        call self%add_integral("mass flux")

    end subroutine init
    !******************************************************************************************





    !>  This procedure checks that all the information have been provided by the user to fully 
    !!  compute the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine check(self)
        class(mass_flux_t), intent(inout)   :: self

        integer(ik)                     :: aux_geoms, ref_geoms
        character(len=:), allocatable   :: usr_msg_r, usr_msg_a
        logical                         :: ref_geoms_exceed, aux_geoms_exceed


        ! Check that the functional has all the information needed
        ref_geoms = self%n_ref_geom()
        aux_geoms = self%n_aux_geom()

        ref_geoms_exceed = (ref_geoms < self%min_ref_geom .or. ref_geoms > self%max_ref_geom)
        aux_geoms_exceed = (aux_geoms < self%min_aux_geom .or. aux_geoms > self%max_aux_geom)


        usr_msg_r = "fcl_mass_average_flowangle12: wrong number of reference geometries. The minimum number of reference geometries is 1."
        usr_msg_a = "fcl_mass_average_flowangle12: wrong number of auxiliary geometries. The minimum number of auxiliary geometries is 1."

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check 
    !******************************************************************************************






    !>  Compute mass flux 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !!  param[inout]       worker
    !!  param[inout]       cache        ! Storage for integrals computed in each element/face
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(mass_flux_t),             intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(functional_cache_t),       intent(inout)   :: cache

        type(AD_D)                              :: mass_flux
        type(AD_D), allocatable, dimension(:)   :: density, mom_1, mom_2, mom_3, &
                                                   u_1, u_2, u_3, vn, r, mass,   &
                                                   unorm_1, unorm_2, unorm_3
        real(rk),   allocatable, dimension(:)   :: weights
        
        
        ! Get primary fields
        density = worker%get_field('Density'   ,'value','face interior')
        mom_1   = worker%get_field('Momentum-1','value','face interior')
        mom_2   = worker%get_field('Momentum-2','value','face interior')
        mom_3   = worker%get_field('Momentum-3','value','face interior')
        
        ! Consider cylindrical coordinates
        r = worker%coordinate('1','face interior')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom_2 = mom_2/r
!        else if (worker%coordinate_system() == 'Cartesian') then
!
!        else
!            call chidg_signal(FATAL,'bad coordinate system')
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

        ! Compute normal velocity
        vn = u_1*unorm_1 + u_2*unorm_2 + u_3*unorm_3

        ! Compute mass flux at GQ
        mass = density*vn

        ! Compute mass-weighted integral 
        mass_flux = integrate_surface(worker,mass)
       
        ! Store in cache 
        call cache%set_entity_value(worker%mesh,mass_flux,'mass flux','reference',worker%function_info)

    end subroutine compute_functional
    !******************************************************************************************


    
    
    
    !>  Store the real value of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_value(self,cache) result(res) 
        class(mass_flux_t),             intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        real(rk) :: res

        res = cache%ref_cache%get_real('mass flux')

    end function store_value
    !*********************************************************************************************



    
    
    !>  Store the derivatives of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/23/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_deriv(self,cache) result(res) 
        class(mass_flux_t),             intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        type(chidg_vector_t) :: res
        
        res = cache%ref_cache%get_deriv('mass flux')

    end function store_deriv
    !*********************************************************************************************



end module fcl_mass_flux
