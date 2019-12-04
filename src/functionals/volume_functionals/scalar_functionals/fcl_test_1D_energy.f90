module fcl_test_1D_energy
#include <messenger.h>

    use mod_kinds,                  only: ik, rk
    use type_evaluator,             only: evaluator_t
    use type_chidg_worker,          only: chidg_worker_t
    use mod_constants,              only: HALF, ONE, TWO
    use type_functional_cache,      only: functional_cache_t
    use type_chidg_vector,          only: chidg_vector_t
    use mod_functional_operators
    
    use DNAD_D


    implicit none

    !>  Functional for 1D tests.
    !!  Dummy energy: 1/2*u**2
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/17/2017
    !!
    !!  Restructured functional computation
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: test_1D_energy_t
   
        ! Quantities needed for the computation
        ! Used to check completeness of information provided
        integer(ik)     :: min_ref_geom = 1
        integer(ik)     :: min_aux_geom = 0
        integer(ik)     :: max_ref_geom = 100
        integer(ik)     :: max_aux_geom = 0

    contains

        procedure, public    :: init
        procedure, public    :: check
        procedure, public    :: compute_functional
        procedure, public    :: store_value
        procedure, public    :: store_deriv

    end type test_1D_energy_t
    !******************************************************************************************
    


contains







    !>  Initialize the functional
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   05/11/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(test_1D_energy_t), intent(inout)   :: self


        !
        ! Start defining the evaluator information
        !
        call self%set_name("Test 1D Energy")
        call self%set_eval_type("Functional")
        call self%set_int_type("VOLUME INTEGRAL")
        

    end subroutine init
    !******************************************************************************************








    !>  This procedure checks that all the information have been provided by the user to fully 
    !!  compute the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/17/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine check(self)
        class(test_1D_energy_t), intent(inout)   :: self

        integer(ik)                     :: aux_geoms, ref_geoms
        character(len=:), allocatable   :: usr_msg_r, usr_msg_a
        logical                         :: ref_geoms_exceed, aux_geoms_exceed


        !
        ! Check that the functional has all the information needed
        ! 
        ref_geoms = self%n_ref_geom()
        aux_geoms = self%n_aux_geom()

        ref_geoms_exceed = (ref_geoms < self%min_ref_geom .or. ref_geoms > self%max_ref_geom)
        aux_geoms_exceed = (aux_geoms < self%min_aux_geom .or. aux_geoms > self%max_aux_geom)


        usr_msg_r = "fcl_test_1D_energy: wrong number of reference geometries. The minimum number of reference geometries required is 1"
        usr_msg_a = "fcl_test_1D_energy: This functional does not require any auxiliary geometry. Please, use `chidg edit` to remove the input auxiliary geometry"

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check
    !******************************************************************************************






    !>  Compute Kinetic Energy Integral of a single element.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/17/2017
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
        class(test_1D_energy_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(functional_cache_t),   intent(inout)   :: cache

        integer(ik)             :: idomain_g, ielement_g, iface, ierr
        
        type(AD_D), allocatable, dimension(:) :: u, Ek
        type(AD_D)                            :: energy_integral

        
        ! 
        ! Get primary fields
        !
        u = worker%get_field('u','value','element')
        
        !
        ! Compute Kinetic Energy
        !
        Ek = HALF * u * u

        !
        ! Compute integral over the volume of the functional
        !
        energy_integral = integrate_volume(worker,Ek)

        
        !
        ! Store in cache
        !
        call cache%set_value(worker%mesh,energy_integral,'kinetic energy','reference',worker%function_info) 

    end subroutine compute_functional
    !******************************************************************************************





    
    
    
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
        class(test_1D_energy_t),    intent(inout)   :: self
        type(functional_cache_t),   intent(inout)   :: cache

        real(rk)          :: res

        res = cache%ref_cache%get_real('kinetic energy')

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
        class(test_1D_energy_t),    intent(inout)   :: self
        type(functional_cache_t),   intent(inout)   :: cache

        type(chidg_vector_t)       :: res
        
        
        res = cache%ref_cache%get_deriv('kinetic energy')

    end function store_deriv
    !*********************************************************************************************











end module fcl_test_1D_energy
