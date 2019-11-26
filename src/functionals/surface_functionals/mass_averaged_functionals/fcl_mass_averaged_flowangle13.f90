module fcl_mass_averaged_flowangle13
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



    !>  This functional computes the mass-averaged flow angle at a given section:
    !!      - reference geom
    !!
    !!
    !!          \--------------- Axis 1
    !!          |\  theta
    !!          | \
    !!          |  \
    !!          |   \ V
    !!          |  
    !!          axis 3
    !!
    !!            sum{int_face_i{theta_gq*rho_gq*v_gq}}
    !!  theta = ----------------------------
    !!              int_geom{n*rho*V}
    !!  
    !!  Where gq refers to the value at a quadrature node, _i is the ith face.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !!
    !------------------------------------------------------------------------------------------
    type,  extends(evaluator_t), public   :: mass_averaged_flowangle13_t
   
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
        procedure, public   :: store_value
        procedure, public   :: store_deriv

    end type mass_averaged_flowangle13_t
    !******************************************************************************************
    


contains




    !>  Initialize the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self)
        class(mass_averaged_flowangle13_t), intent(inout)   :: self

        call self%set_name("MA Flow Angle 13")
        call self%set_eval_type("Functional")
        call self%set_int_type("FACE INTEGRAL")

    end subroutine init
    !******************************************************************************************






    !>  This procedure checks that all the information have been provided by the user to fully 
    !!  compute the functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !------------------------------------------------------------------------------------------
    subroutine check(self)
        class(mass_averaged_flowangle13_t), intent(inout)   :: self

        integer(ik)                     :: aux_geoms, ref_geoms
        character(len=:), allocatable   :: usr_msg_r, usr_msg_a
        logical                         :: ref_geoms_exceed, aux_geoms_exceed

        ref_geoms = self%n_ref_geom()
        aux_geoms = self%n_aux_geom()

        ref_geoms_exceed = (ref_geoms < self%min_ref_geom .or. ref_geoms > self%max_ref_geom)
        aux_geoms_exceed = (aux_geoms < self%min_aux_geom .or. aux_geoms > self%max_aux_geom)

        usr_msg_r = "fcl_mass_average_flowangle13: wrong number of reference geometries. The minimum number of reference geometries is 1."
        usr_msg_a = "fcl_mass_average_flowangle13: wrong number of auxiliary geometries. The minimum number of auxiliary geometries is 1."

        if (ref_geoms_exceed) call chidg_signal(FATAL, usr_msg_r)
        if (aux_geoms_exceed) call chidg_signal(FATAL, usr_msg_a)

    end subroutine check 
    !******************************************************************************************







    !>  Compute 1-2 flow angle flux through an element face.
    !!  The sum of the contributions of the flow angle flux of all the faces belonging to the 
    !!  patch will be divided by the mass flux to achieve the mass-averaged flow angle
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !!  param[inout]       worker
    !!  param[inout]       cache        ! Storage for integrals computed in each element/face
    !!
    !------------------------------------------------------------------------------------------
    subroutine compute_functional(self,worker,cache)
        class(mass_averaged_flowangle13_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(functional_cache_t),       intent(inout)   :: cache

        integer(ik)             :: idomain_g, ielement_g, iface, ierr
        
        type(AD_D), allocatable, dimension(:) :: density,                       &
                                                 mom_1, mom_2, mom_3,           &
                                                 u_1, u_2, u_3, norm_vel,       &
                                                 norm_1, norm_2, norm_3,        &
                                                 unorm_1, unorm_2, unorm_3, r
        type(AD_D)                            :: u1_flux, u3_flux, mass_flux

        
        ! Get primary fields
        density     = worker%get_field('Density'   ,'value','face interior')
        mom_1       = worker%get_field('Momentum-1','value','face interior')
        mom_2       = worker%get_field('Momentum-2','value','face interior')
        mom_3       = worker%get_field('Momentum-3','value','face interior')

        ! Account for cylindrical coordinates
        r = worker%coordinate('1','face interior')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom_2 = mom_2/r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,'bad coordinate system')
        end if

        ! Compute velocities
        u_1 = mom_1/density
        u_2 = mom_2/density
        u_3 = mom_3/density


        ! Compute mass-weighted face integral
        u1_flux = integrate_surface_mass_weighted(worker,u_1)
        u3_flux = integrate_surface_mass_weighted(worker,u_3)

        ! Compute mass flux
        mass_flux =  integrate_surface_mass_weighted(worker)
        
        ! Store in cache 
        call cache%set_value(u1_flux,  'u1 flux',  'reference',worker%function_info)
        call cache%set_value(u3_flux,  'u3 flux',  'reference',worker%function_info)
        call cache%set_value(mass_flux,'mass flux','reference',worker%function_info)

    end subroutine compute_functional
    !******************************************************************************************




    
    
    
    !>  Finalize subroutine carries out the final computation before storing the functional
    !!  It comes after the communcation among processors and input argument 'cache' has 
    !!  the over all integrals (in this case 'flow angle flux' and 'mass flux') over the 
    !!  entire referenc geometry.
    !!
    !!  At this point the final division 'flow angle flux'/'mass flux' can be accomplished
    !! 
    !!
    !!  NOTE: here the integral contribution from each process has been shared and summed.
    !!        Therefore, we are dealing here with integrals over the entire geometry.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !!  param[inout]       cache     cache contains the overall integrals value at this point,
    !!                               since the parallel communication already happened 
    !!
    !---------------------------------------------------------------------------------------------
    subroutine finalize_functional(self,cache) 
        class(mass_averaged_flowangle13_t), intent(inout)   :: self
        type(functional_cache_t),           intent(inout)   :: cache

        type(AD_D)        :: u1_flux, u3_flux, mass_flux, MA_u1, MA_u3, angle_rad, angle_deg
        
        ! Compute MA velocities
        u1_flux   = cache%get_value('u1 flux',  'reference')
        u3_flux   = cache%get_value('u3 flux',  'reference')
        mass_flux = cache%get_value('mass flux','reference')

        MA_u1 = u1_flux/mass_flux
        MA_u3 = u3_flux/mass_flux


        ! Compute flow angle 
        angle_rad = atan(MA_u3/MA_u1)

        ! Convert angle from radians to degrees to be more readable
        angle_deg = (angle_rad*180._rk)/PI

        ! Store MA flow angle
        call cache%set_value(angle_deg,'MA flow angle','reference')

    end subroutine finalize_functional
    !*********************************************************************************************


    
    
    
    
    !>  Store the real value of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_value(self,cache) result(res) 
        class(mass_averaged_flowangle13_t), intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        real(rk)          :: res

        res = cache%ref_cache%get_real('MA flow angle')

    end function store_value
    !*********************************************************************************************



    
    
    !>  Store the derivatives of the actual final functional integral 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   4/4/2019
    !!
    !!  param[inout]       cache     Storage for integrals, this contains the overall 
    !!                               functional (ie after parallel communication).     
    !!
    !---------------------------------------------------------------------------------------------
    function store_deriv(self,cache) result(res) 
        class(mass_averaged_flowangle13_t), intent(inout)   :: self
        type(functional_cache_t),       intent(inout)  :: cache

        type(chidg_vector_t)       :: res
        
        res = cache%ref_cache%get_deriv('MA flow angle')

    end function store_deriv
    !*********************************************************************************************



end module fcl_mass_averaged_flowangle13
