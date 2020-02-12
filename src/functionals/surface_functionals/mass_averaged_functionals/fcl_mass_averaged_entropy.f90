module fcl_mass_averaged_entropy
#include <messenger.h>
    use mod_kinds,                  only: ik, rk
    use type_evaluator,             only: evaluator_t
    use type_chidg_worker,          only: chidg_worker_t
    use mod_fluid,                  only: Rgas,cp, gam
    use mod_constants,              only: HALF, ONE
    use type_functional_cache,      only: functional_cache_t
    use type_chidg_vector,          only: chidg_vector_t
    use mod_functional_operators
    use DNAD_D
    implicit none


    !>  This functional computes the mass-averaged entropy rise between two section:
    !!      - auxiliary geom
    !!      - reference geom
    !!
    !!  Auxiliary geometry is most likely an inlet surface (or in proximity of the inlet)
    !!  Reference geometry is most likely a surface close to the exit of the domain
    !!
    !!
    !!  ds_iexit = [ ( cp * log(T0_iexit/T0_in) ) - ( Rgas * log(P0_iexit/P0_in) ) ]
    !!
    !!
    !!       sum_i(ref_geom) [ ds_iexit ]*[ da_iexit * rho_iexit * vn_iexit ]
    !!  dS = -----------------------------------------------------------------
    !!             sum_i(ref_geom) [ da_iexit * rho_iexit * vn_iexit ]
    !!
    !!  
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: mass_averaged_entropy_t
   
        ! Quantities needed for the computation
        ! Used to check completeness of information provided
        integer(ik)     :: min_ref_geom = 1
        integer(ik)     :: min_aux_geom = 1
        integer(ik)     :: max_ref_geom = 100
        integer(ik)     :: max_aux_geom = 100

    contains

        procedure, public   :: init
        procedure, public   :: check
        procedure, public   :: compute_functional
        procedure, public   :: compute_auxiliary
        procedure, public   :: finalize_functional
        procedure, public   :: finalize_auxiliary
        procedure, public   :: store_value
        procedure, public   :: store_deriv

    end type mass_averaged_entropy_t
    !******************************************************************************************
    


contains



    !>  Initialize the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(mass_averaged_entropy_t), intent(inout)   :: self

        ! Start defining the evaluator information
        call self%set_name("MA Entropy")
        call self%set_eval_type("Functional")
        call self%set_int_type("FACE INTEGRAL")

        call self%add_integral("temperature flux")
        call self%add_integral("pressure flux")
        call self%add_integral("mass flux")
        call self%add_integral("MA total pressure")
        call self%add_integral("MA total temperature")
        call self%add_integral("MA entropy")

    end subroutine init
    !******************************************************************************************






    !>  This procedure checks that all the information have been provided by the user to fully 
    !!  compute the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/03/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine check(self)
        class(mass_averaged_entropy_t), intent(inout)   :: self

        integer(ik)                     :: aux_geoms, ref_geoms
        character(len=:), allocatable   :: usr_msg_r, usr_msg_a
        logical                         :: ref_geoms_exceed, aux_geoms_exceed


        ! Check that the functional has all the information needed
        ref_geoms = self%n_ref_geom()
        aux_geoms = self%n_aux_geom()

        ref_geoms_exceed = (ref_geoms < self%min_ref_geom .or. ref_geoms > self%max_ref_geom)
        aux_geoms_exceed = (aux_geoms < self%min_aux_geom .or. aux_geoms > self%max_aux_geom)

        usr_msg_r = "fcl_mass_average_entropy: wrong number of reference geometries. The minimum number of reference geometries is 1."
        usr_msg_a = "fcl_mass_average_entropy: wrong number of auxiliary geometries. The minimum number of auxiliary geometries is 1."

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check 
    !******************************************************************************************






    !>  This subroutine is meant to compute the auxiliary integral.
    !!  Indeed, we need to compute the mass-averaged total temperature and pressure on a 
    !!  reference geometry (let's say an inlet boundary) and use it later on for computing 
    !!  the entropy production 
    !!  
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       worker
    !!  param[inout]       cache        ! Storage for functional integrals
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_auxiliary(self,worker,cache)
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(functional_cache_t),       intent(inout)   :: cache
        
        type(AD_D), allocatable, dimension(:) :: pressure, temperature,         &
                                                 density, mom_1, mom_2, mom_3,  &
                                                 u_1, u_2, u_3, tot_pressure,   &
                                                 tot_temperature, umag, c, mach,&
                                                 norm_1, norm_2, norm_3,        &
                                                 unorm_1, unorm_2, unorm_3, r   
        type(AD_D)                            :: temperature_flux, pressure_flux,&
                                                 mass_flux
        
        
        ! Get primary fields
        density = worker%get_field('Density'   ,'value','face interior')
        mom_1   = worker%get_field('Momentum-1','value','face interior')
        mom_2   = worker%get_field('Momentum-2','value','face interior')
        mom_3   = worker%get_field('Momentum-3','value','face interior')

        ! Account for cylindrical coordinates
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom_2 = mom_2/r
        end if

        ! Get necessary models
        pressure    = worker%get_field('Pressure'  ,'value','face interior')
        temperature = worker%get_field('Temperature','value','face interior')
        
        ! Compute velocities
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density
        umag = sqrt(u_1*u_1 + u_2*u_2 + u_3*u_3)
        
        ! Compute Mach number
        c = sqrt( gam * Rgas * temperature)
        mach = umag/c
        
        ! Compute total pressure and total temperature
        tot_pressure    = pressure * (ONE + (gam - ONE)*HALF*mach*mach)**(gam/(gam - ONE))
        tot_temperature = temperature * (ONE + (gam - ONE)*HALF*mach*mach)

        ! Send the total pressure and temperature to the averaging subroutine
        temperature_flux = integrate_surface_mass_weighted(worker,tot_temperature)
        pressure_flux    = integrate_surface_mass_weighted(worker,tot_pressure)
        
        ! Call the averaging subroutine for computing pure mass flux through the face
        mass_flux = integrate_surface_mass_weighted(worker)
     
        ! Store in cache 
        call cache%set_entity_value(worker%mesh,temperature_flux,'temperature flux','auxiliary',worker%function_info) 
        call cache%set_entity_value(worker%mesh,pressure_flux,   'pressure flux',   'auxiliary',worker%function_info) 
        call cache%set_entity_value(worker%mesh,mass_flux,       'mass flux',       'auxiliary',worker%function_info) 
    
    end subroutine compute_auxiliary
    !******************************************************************************************

    
    
    
    
    
    
    !>  This subroutine compute the final total temperature and total pressure mass-average
    !!  on the auxiliary boundary.
    !! 
    !!  The quantities stored are:
    !!      - P0aux = P0_flux/mass_flux
    !!      - T0aux = T0_flux/mass_flux
    !!
    !!  NOTE: here the integral contribution from each process has been shared and summed.
    !!        Therefore, we are dealing here with integrals over the entire geometry.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       cache     cache contains the overall iintegrals value at this point,
    !!                               since the parallel communication already happened 
    !!                               
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_auxiliary(self,worker,cache)
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(in)      :: worker
        type(functional_cache_t),       intent(inout)   :: cache

        type(AD_D)        :: temperature_flux, pressure_flux, mass_flux
        
        temperature_flux = cache%get_global_value(worker%mesh,'temperature flux','auxiliary')
        pressure_flux    = cache%get_global_value(worker%mesh,'pressure flux','auxiliary')
        mass_flux        = cache%get_global_value(worker%mesh,'mass flux','auxiliary')

        call cache%set_global_value(worker%mesh,pressure_flux/mass_flux,   'MA total pressure',   'auxiliary')
        call cache%set_global_value(worker%mesh,temperature_flux/mass_flux,'MA total temperature','auxiliary')

    end subroutine finalize_auxiliary
    !*********************************************************************************************






    !>  Compute mass weighted total pressure and temperature 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       worker
    !!  param[inout]       cache        ! Storage for integrals computed in each element/face
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(functional_cache_t),       intent(inout)   :: cache

        integer(ik)             :: idomain_g, ielement_g, iface, ierr
        
        type(AD_D), allocatable, dimension(:) :: pressure, temperature,     &
                                                 deltaS, density,           &
                                                 mom_1, mom_2, mom_3,       &
                                                 u_1, u_2, u_3, norm_vel,   &
                                                 tot_pressure,              &
                                                 tot_temperature, umag, c,  &
                                                 mach, r,                   &
                                                 norm_1, norm_2, norm_3,    &
                                                 unorm_1, unorm_2, unorm_3
        type(AD_D)                            :: mass_flux,                 &
                                                 temperature_flux,          &
                                                 pressure_flux

        
        ! Get primary fields
        density     = worker%get_field('Density'   ,'value','face interior')
        mom_1       = worker%get_field('Momentum-1','value','face interior')
        mom_2       = worker%get_field('Momentum-2','value','face interior')
        mom_3       = worker%get_field('Momentum-3','value','face interior')

        ! Account for cylindrical coordinates
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom_2 = mom_2/r
        end if

        ! Get necessary models
        pressure    = worker%get_field('Pressure'  ,'value','face interior')
        temperature = worker%get_field('Temperature','value','face interior')
        
        ! Compute velocities
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density
        umag = sqrt(u_1*u_1 + u_2*u_2 + u_3*u_3)

        ! Compute Mach number
        c = sqrt( gam * Rgas * temperature)
        mach = umag/c

        ! Compute total pressure and total temperature
        tot_pressure    = pressure * (ONE + (gam - ONE)*HALF*mach*mach)**(gam/(gam - ONE))
        tot_temperature = temperature * (ONE + (gam - ONE)*HALF*mach*mach)
        
        ! Compute mass-weighted face integral
        pressure_flux    = integrate_surface_mass_weighted(worker,tot_pressure)
        temperature_flux = integrate_surface_mass_weighted(worker,tot_temperature)

        ! Compute mass flux
        mass_flux =  integrate_surface_mass_weighted(worker)
        
        ! Store in cache 
        call cache%set_entity_value(worker%mesh,pressure_flux,   'pressure flux',   'reference',worker%function_info)
        call cache%set_entity_value(worker%mesh,temperature_flux,'temperature flux','reference',worker%function_info)
        call cache%set_entity_value(worker%mesh,mass_flux,       'mass flux',       'reference',worker%function_info)

    end subroutine compute_functional
    !******************************************************************************************

    
    
    
    
    
    
    
    !>  Finalize subroutine carries out the final computation before storing the functional
    !!  It comes after the communcation among processore and input argument cache has 
    !!  the over all integrals (in this case 'entropy flux' and 'mass flux') over the reference
    !!  geometry.
    !!
    !!
    !!  NOTE: here the integral contribution from each process has been shared and summed.
    !!        Therefore, we are dealing here with integrals over the entire geometry.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       cache     cache contains the overall iintegrals value at this point,
    !!                               since the parallel communication already happened 
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_functional(self,worker,cache) 
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(in)      :: worker
        type(functional_cache_t),       intent(inout)   :: cache

        type(AD_D)        :: ref_pressure_flux, ref_mass_flux, deltaS, aux_total_P, aux_total_T, &
                             ref_total_P, ref_total_T, ref_temperature_flux
        
        ! Compute mass-averaged total pressure and temperature on reference geometries
        ref_pressure_flux    = cache%get_global_value(worker%mesh,'pressure flux',   'reference')
        ref_temperature_flux = cache%get_global_value(worker%mesh,'temperature flux','reference')
        ref_mass_flux        = cache%get_global_value(worker%mesh,'mass flux',       'reference')

        ref_total_P = ref_pressure_flux/ref_mass_flux
        ref_total_T = ref_temperature_flux/ref_mass_flux

        ! Retrieve mass-averaged total pressure and temperature on auxiliary geoemtries
        aux_total_P = cache%get_global_value(worker%mesh,'MA total pressure',   'auxiliary')
        aux_total_T = cache%get_global_value(worker%mesh,'MA total temperature','auxiliary')

        ! Compute entropy production
        deltaS = ( cp * log(ref_total_T/aux_total_T) ) - ( Rgas * log(ref_total_P/aux_total_P) )

        ! Store entropy 
        call cache%set_global_value(worker%mesh,deltaS,'MA entropy','reference')

    end subroutine finalize_functional
    !*********************************************************************************************


    
    
    
    
    !>  Store the real value of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_value(self,cache) result(res) 
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        real(rk)          :: res

        res = cache%ref_cache%get_real('MA entropy')

    end function store_value
    !*********************************************************************************************



    
    
    !>  Store the derivatives of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_deriv(self,cache) result(res) 
        class(mass_averaged_entropy_t), intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        type(chidg_vector_t)       :: res
        
        res = cache%ref_cache%get_deriv('MA entropy')

    end function store_deriv
    !*********************************************************************************************



end module fcl_mass_averaged_entropy
