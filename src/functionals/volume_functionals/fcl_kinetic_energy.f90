module fcl_kinetic_energy
#include <messenger.h>

    use mod_kinds,                  only: ik, rk
    use type_evaluator,             only: evaluator_t
    use type_chidg_worker,          only: chidg_worker_t
    use mod_constants,              only: HALF, ONE, TWO
    use type_integral_cache,        only: integral_cache_t
    use mod_functional_operators,   only: integrate_volume
    
    use DNAD_D


    implicit none

    !>  This functional computes the integral of the kientic energy over a volume
    !!
    !!  Ek = 1/2*u**2
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/17/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: kinetic_energy_t
   
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
        procedure, public   :: finalize_functional

    end type kinetic_energy_t
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
        class(kinetic_energy_t), intent(inout)   :: self


        !
        ! Start defining the evaluator information
        !
        call self%set_name("Kinetic Energy")
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
        class(kinetic_energy_t), intent(inout)   :: self

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


        usr_msg_r = "fcl_kinetic_energy: wrong number of reference geometries. The minimum number of reference geoemtries required is 1"
        usr_msg_a = "fcl_kinetic_energy: This functional does not require any auxiliary geometry. Please, use `chidg edit` to remove the input auxiliary geometry"

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check
    !******************************************************************************************






    !>  Compute Kinetic Energy Integral of a single element.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   12/17/2017
    !!
    !!  param[inout]       worker
    !!  param[inout]       cache        ! Storage for integrals computed in each element/face
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(kinetic_energy_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(integral_cache_t),     intent(inout)   :: cache    

        integer(ik)             :: idomain_g, ielement_g, iface, ierr
        
        type(AD_D), allocatable, dimension(:) :: density,                   &
                                                 mom_1, mom_2, mom_3,       &
                                                 u_1, u_2, u_3, norm_vel,   &
                                                 umag, Ek,                  &
                                                 norm_1, norm_2, norm_3,    &
                                                 unorm_1, unorm_2, unorm_3, &
                                                 da, r
        real(rk),   allocatable, dimension(:) :: weights 
        type(AD_D)                            :: energy_integral

        
        ! Get primary fields
        !
        density     = worker%get_field('Density'   ,'value','element')
        mom_1       = worker%get_field('Momentum-1','value','element')
        mom_2       = worker%get_field('Momentum-2','value','element')
        mom_3       = worker%get_field('Momentum-3','value','element')

        !
        ! Account for cylindrical coordinates
        !
        r = worker%coordinate('1','face interior')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom_2 = mom_2/r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,'bad coordinate system')
        end if


        !
        ! Compute velocities
        !
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density
        umag = sqrt(u_1*u_1 + u_2*u_2 + u_3*u_3)

        !
        ! Compute Kinetic Energy
        !
        Ek = HALF * density * umag**TWO

        !
        ! Compute volume integral
        !
        energy_integral = integrate_volume(worker,Ek)
        
        
        !
        ! Store in cache
        !
        call cache%push_back(energy_integral,'kinetic energy')
       

    end subroutine compute_functional
    !******************************************************************************************








    
    !>  Finalize subroutine carries out the final computation before storing the functional
    !!  It comes after the communcation among processore and input argument cache has 
    !!  the over all integrals (in this case only 'kinetic energy') over the reference
    !!  geometry.
    !!
    !!  @auhtor Matteo Ugolotti
    !!  @date   12/17/2017
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function finalize_functional(self,cache) result (res)
        class(kinetic_energy_t),    intent(inout)   :: self
        type(integral_cache_t),     intent(inout)   :: cache    

        type(AD_D)      :: res

        res = cache%get_value('kinetic energy')


    end function finalize_functional
    !*********************************************************************************************









end module fcl_kinetic_energy
