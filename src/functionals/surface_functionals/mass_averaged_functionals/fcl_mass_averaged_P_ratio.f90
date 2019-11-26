module fcl_mass_averaged_P_ratio
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


    !>  This functional computes the mass-averaged pressure ratio between two section:
    !!      - auxiliary geom
    !!      - reference geom (derivatives)
    !!
    !!  Auxiliary geometry is most likely an stage inlet surface (or in proximity of the inlet)
    !!  Reference geometry is most likely a surface close to the exit of the domain
    !!
    !!
    !!        P_exit 
    !!  PR = ---------
    !!        P_in 
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
    !!  TODO: to be tested
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: mass_averaged_P_ratio_t
   
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

    end type mass_averaged_P_ratio_t
    !******************************************************************************************
    


contains



    !>  Initialize the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/20/2018
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(mass_averaged_P_ratio_t), intent(inout)   :: self

        call self%set_name("MA Static Pressure Ratio")
        call self%set_eval_type("Functional")
        call self%set_int_type("FACE INTEGRAL")

    end subroutine init
    !******************************************************************************************






    !>  This procedure checks that all the information have been provided by the user to fully 
    !!  compute the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/20/2018
    !!
    !------------------------------------------------------------------------------------------
    subroutine check(self)
        class(mass_averaged_P_ratio_t), intent(inout)   :: self

        integer(ik)                     :: aux_geoms, ref_geoms
        character(len=:), allocatable   :: usr_msg_r, usr_msg_a
        logical                         :: ref_geoms_exceed, aux_geoms_exceed


        ! Check that the functional has all the information needed
        ref_geoms = self%n_ref_geom()
        aux_geoms = self%n_aux_geom()

        ref_geoms_exceed = (ref_geoms < self%min_ref_geom .or. ref_geoms > self%max_ref_geom)
        aux_geoms_exceed = (aux_geoms < self%min_aux_geom .or. aux_geoms > self%max_aux_geom)

        usr_msg_r = "fcl_mass_average_P_ratio: wrong number of reference geometries. The minimum number of reference geometries is 1."
        usr_msg_a = "fcl_mass_average_P_ratio: wrong number of auxiliary geometries. The minimum number of auxiliary geometries is 1."

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check 
    !******************************************************************************************






    !>  This subroutine is meant to compute the auxiliary integral.
    !!  Indeed, we need to compute the mass-averaged pressure on a 
    !!  auxiliary geometry (let's say an inlet boundary) and use it later on for computing 
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
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(functional_cache_t),               intent(inout)   :: cache
        
        type(AD_D), allocatable, dimension(:) :: pressure
        type(AD_D)                            :: pressure_flux, mass_flux
        
        ! Get necessary models
        pressure    = worker%get_field('Pressure'  ,'value','face interior')
        
        ! Send the pressure to the averaging subroutine
        pressure_flux = integrate_surface_mass_weighted(worker,pressure)
        
        ! Call the averaging subroutine for computing pure mass flux through the face
        mass_flux    = integrate_surface_mass_weighted(worker)
     
        ! Store in cache 
        call cache%set_value(pressure_flux,   'pressure flux','auxiliary',worker%function_info) 
        call cache%set_value(mass_flux,       'mass flux',    'auxiliary',worker%function_info) 
    
    end subroutine compute_auxiliary
    !******************************************************************************************

    
    

    
    
    !>  This subroutine compute the final pressure mass-average
    !!  on the auxiliary boundary.
    !! 
    !!  The quantities stored are:
    !!      - Paux = P_flux/mass_flux
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
    subroutine finalize_auxiliary(self,cache)
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(functional_cache_t),               intent(inout)   :: cache

        type(AD_D)        :: pressure_flux, mass_flux
        
        pressure_flux    = cache%get_value('pressure flux','auxiliary')
        mass_flux        = cache%get_value('mass flux',    'auxiliary')

        call cache%set_value(pressure_flux/mass_flux,'MA pressure','auxiliary')

    end subroutine finalize_auxiliary
    !*********************************************************************************************






    !>  Compute pressure flux over each element face.
    !!  The sum of the contributions of the pressure ratio flux of all the faces belonging to the 
    !!  patch will be divided by the mass flux to achieve the mass-averaged pressure ratio
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
    !!  param[inout]       cache     cache contains the overall iintegrals value at this point,
    !!                               since the parallel communication already happened 
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        type(functional_cache_t),               intent(inout)   :: cache

        integer(ik)             :: idomain_g, ielement_g, iface, ierr
        
        type(AD_D), allocatable, dimension(:) :: pressure
        type(AD_D)                            :: pressure_flux, mass_flux


        ! Get necessary models
        pressure    = worker%get_field('Pressure'  ,'value','face interior')
        
        ! Compute mass-weighted face integral
        pressure_flux = integrate_surface_mass_weighted(worker,pressure)

        ! Compute mass flux
        mass_flux =  integrate_surface_mass_weighted(worker)
        
        ! Store in cache 
        call cache%set_value(pressure_flux,'pressure flux','reference',worker%function_info) 
        call cache%set_value(mass_flux,    'mass flux',    'reference',worker%function_info) 

    end subroutine compute_functional
    !******************************************************************************************

    
    
    
    
    
    
    !>  Finalize subroutine carries out the final computation before storing the functional
    !!  It comes after the communcation among processore and input argument cache has 
    !!  the over all integrals (in this case 'pressure flux' and 'mass flux') over the reference
    !!  geometry.
    !!
    !!  At this point the final division 'pressure flux'/'mass flux' can be accomplished
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
    !!  param[inout]       worker
    !!  param[inout]       cache     cache contains the overall iintegrals value at this point,
    !!                               since the parallel communication already happened 
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_functional(self,cache) 
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(functional_cache_t),               intent(inout)   :: cache

        type(AD_D)        :: pressure_flux, mass_flux, MA_ref_pressure, MA_aux_pressure
        
        ! Finilize functiona on reference geometry
        pressure_flux    = cache%get_value('pressure flux','reference')
        mass_flux        = cache%get_value('mass flux',    'reference')
        MA_ref_pressure  = pressure_flux/mass_flux

        ! Compute the final PRESSURE RATIO using auxiliary and reference pressures
        MA_aux_pressure  = cache%get_value('MA pressure',  'auxiliary')
        call cache%set_value(MA_ref_pressure/MA_aux_pressure,'MA pressure ratio','reference')

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
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(functional_cache_t),               intent(inout)   :: cache

        real(rk)          :: res

        res = cache%ref_cache%get_real('MA pressure ratio')

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
        class(mass_averaged_P_ratio_t),  intent(inout)         :: self
        type(functional_cache_t),               intent(inout)   :: cache

        type(chidg_vector_t)       :: res
        
        res = cache%ref_cache%get_deriv('MA pressure ratio')

    end function store_deriv
    !*********************************************************************************************





end module fcl_mass_averaged_P_ratio
