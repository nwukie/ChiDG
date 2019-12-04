module type_integral_cache
#include<messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, NO_ID, NO_DIFF
    use type_mesh,              only: mesh_t
    use type_svector,           only: svector_t
    use mod_string,             only: string_t
    use type_seed,              only: seed_t
    use mod_differentiate,      only: get_nderiv
    use type_chidg_vector,      only: chidg_vector_t
    use type_function_info,     only: function_info_t
    use DNAD_D
    implicit none


    !>  Container to store functional integrals (such as mass-weighted entropy) on a reference geometry.
    !!
    !!  Consider:
    !!      [integral cache] = integrated mass flux
    !!      [integral cache] = integrated mass-weighted entropy
    !!                           |
    !!                           v
    !!      [integral cache] = mass-averaged entropy
    !!
    !!
    !!
    !! 
    !!                                             [integral cache]  (e.g. pressure)
    !!                                            /
    !!                         [geometry cache] -- [integral cache]  (e.g. temperature)
    !!                       /                    \
    !!                      /                      [integral cache]
    !!
    !!  [Functional Cache]
    !!                       
    !!                      \                      [integral cache]  (e.g. T0)
    !!                       \                    /
    !!                         [geometry cache] -- [integral cache]  (e.g. P0)
    !!                                            \
    !!                                             [integral cache]
    !!
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   2/25/2018
    !!
    !!  Restructured
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    type, public :: integral_cache_t
        
        ! Value and derivatives for a particular integral 
        character(:),   allocatable :: name
        real(rk)                    :: integral_value
        type(chidg_vector_t)        :: integral_deriv

        ! Derivative initialization flag
        ! NOTE:
        ! It might be that one processor has no face/element of the auxiliary geometry but it does have
        ! faces/elements of the reference geometry. In this case, the functional has to know about the
        ! overall value of the auxilary integral but it is not needed to store a vector of null derivaties.
        logical                     :: derivatives_initialized
         
    contains
        
        procedure, private :: initialize
        procedure, private :: set_entity_value
        procedure, private :: get_entity_value
        procedure, private :: set_global_value
        procedure, private :: get_global_value
        
        generic,   public :: init       => initialize                         ! Initialize integral
        generic,   public :: get_value  => get_entity_value, get_global_value ! Return the element or global AD_D type
        generic,   public :: set_value  => set_entity_value, set_global_value ! Set element or global AD_D type

    end type integral_cache_t
    !***************************************************************************************************


contains




    !>  Initialize integral 
    !!      - Set name
    !!      - Set real value to zero
    !!      - Set derivatives initialization flag to false
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine initialize(self,i_name)
        class(integral_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: i_name
         
        ! Set name
        self%name = i_name

        ! Initialize real value of the integral
        self%integral_value = ZERO

        ! Set initialized flag
        self%derivatives_initialized = .false.

    end subroutine initialize
    !***************************************************************************************************





    !>  Set local integral. This subroutine is used to store the derivatives of the specific integral
    !!  compute on a single face or element of the input geometry, 
    !!  and accumulate (summation) the real value of the integral.
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_entity_value(self,mesh,integral,vec_model,fcn_info) 
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(AD_D),                 intent(in)      :: integral
        type(chidg_vector_t),       intent(in)      :: vec_model 
        type(function_info_t),      intent(in)      :: fcn_info
         
        integer(ik)     :: itime 
            
        ! Accumulate integral real value, this has to be done for any kind of differentiation
        self%integral_value = self%integral_value + integral%x_ad_
    
        !   - dtype = dQ_DIFF or dX_DIFF or NO_DIFF
        !   - no need to store the chidg vector for NO_DIFF
        if (fcn_info%dtype /= NO_DIFF) then

            ! Initialize chidg_vector of derivatives if not yet done
            if ( .not. self%derivatives_initialized ) then
                itime = 1
                self%integral_deriv = vec_model
                self%derivatives_initialized = .true.
            end if

            associate( int_derivs  => self%integral_deriv%dom(fcn_info%seed%idomain_l)%vecs(fcn_info%seed%ielement_l)%vec )
            
                ! Sanity check of vec size and set derivatives
                if ( size(int_derivs) == size(integral%xp_ad_) ) then
                    int_derivs = int_derivs + integral%xp_ad_
                else
                    call chidg_signal(FATAL,"type_integral_cache.f90: the size of the integral derivatives does not match the worker derivatives size. Implementation Error!")
                end if

            end associate

        end if

    end subroutine set_entity_value
    !***************************************************************************************************








    !>  Get local integral.
    !!  This is called to get the overall integral value over the input geometry with derivatives of the
    !!  given seed.
    !! 
    !!  It returns a AD type.
    !!
    !!  This is, for instance, used when you are computing the integral on a single face/element of the reference
    !!  geometry and you need th overall integral over auxiliary geometry differentiated with respect to given
    !!  seed element.
    !!  (This can be seen for instance in Mass Averaged Entropy, where we use P0_ref computed on the auxiliary 
    !!  geometry)
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_entity_value(self,mesh,fcn_info) result(integral)
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(function_info_t),      intent(in)      :: fcn_info
         
        type(AD_D)      :: integral
        integer(ik)     :: nderivs


        ! Allocate result
        integral = AD_D(get_nderiv(fcn_info))
        
        ! Assign real value
        integral = self%integral_value
        
        ! Assign derivatives, if NO_DIFF then self.derivatives_initialized == .false.
        if (self%derivatives_initialized) then
            integral%xp_ad_ = self%integral_deriv%dom(fcn_info%seed%idomain_l)%vecs(fcn_info%seed%ielement_l)%vec 
        end if

    end function get_entity_value
    !***************************************************************************************************






    !>  Set gloabl integral. This subroutine is used to store the derivatives of the specific integral
    !!  compute on the ENTIRE input geometry, as well as the real value of the overall integral.
    !!
    !!  This is most likely used in the "Finalize" step of functional update. For instance, when you want
    !!  to divide a mass weighted quantity by the overall mass for mass average computation on the auxiliary 
    !!  geometry and store it so that it can be used for the reference geometry calculation.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_global_value(self,mesh,integral,vec_model,dtype) 
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(AD_D),                 intent(in)      :: integral
        type(chidg_vector_t),       intent(in)      :: vec_model
        integer(ik),                intent(in)      :: dtype
        
        integer(ik)             :: istart, iend, idom, ielem
       
        ! Initialize istart and iend indeces
        istart = 0
        iend   = 0
        
        ! Set integral real value
        self%integral_value = integral%x_ad_

         
        if (dtype /= NO_DIFF) then

            ! Initialize chidg_vector of derivatives if not yet done
            if ( .not. self%derivatives_initialized  ) then
                self%integral_deriv = vec_model
                self%derivatives_initialized = .true.
            end if
        
        
            ! Loop thorugh domains and elements
            !do idom = 1,self%integral_deriv%ndomains()
                !do ielem = 1,self%integral_deriv%dom(idom)%nelements()
            do idom = 1,mesh%ndomains()
                do ielem = 1,mesh%domain(idom)%nelements()
                    
                    associate( int_derivs  => self%integral_deriv%dom(idom)%vecs(ielem) )
                        call chidg_signal(FATAL,"integral_cache%set_global_value: needs updated for petsc storage!")

                        ! Find correspondent derivatives in the AD_D vector
                        istart = iend + 1
                        iend   = iend + int_derivs%nentries()
                        
                        ! Set derivatives, check vec size
                        if ( int_derivs%nentries() == size(integral%xp_ad_(istart:iend)) ) then
                            int_derivs%vec = integral%xp_ad_(istart:iend)
                        else
                            call chidg_signal(FATAL,"type_integral_cache.f90: the size of the integral derivatives does not match the worker derivatives size. Implementation Error!")
                        end if

                    end associate
                    
                end do !ielem
            end do !idom 

        end if !NO_DIFF
                    
    end subroutine set_global_value
    !***************************************************************************************************






    !>  Get local integral.
    !!  This is called to get the overall integral value over the input geometry with derivatives of the
    !!  wrt the whole geometry entries.
    !! 
    !!  It returns a AD type.
    !!
    !!  This is, for instance, used when you want to divide the overall integral over the reference geometry
    !!  by the overall integral over the auxiliary geometry.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/6/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_global_value(self,mesh,vec_model,dtype) result(integral)
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(chidg_vector_t),       intent(in)      :: vec_model
        integer(ik),                intent(in)      :: dtype
         
        type(AD_D)              :: integral

        integer(ik)             :: istart, iend, idom, ielem, nderivs


        ! Allocate result
        if (dtype /= NO_DIFF) then
            ! if differentiation is neede, check if the derivatives have been already initialized
            ! if not, initialize it.
            ! Derivatives might not be initialized if the processor does not have any face/element
            ! from the auxiliary/reference geometry, and it only knows the overall integral real value
            ! from the communication process.
            if (.not. self%derivatives_initialized) then
                self%integral_deriv = vec_model
            end if
            nderivs  = self%integral_deriv%nentries()
            integral = AD_D(nderivs)
        else
            ! If no differentiation is needed just allcoate zero derivatives
            integral = AD_D(0)
        end if

        ! Assign real value
        integral%x_ad_  = self%integral_value
       
        ! Assign derivatives
        if (self%derivatives_initialized) then
            
            ! Initialize istart and iend indeces
            istart = 0
            iend   = 0
             
            !do idom = 1,self%integral_deriv%ndomains()
            !    do ielem = 1,self%integral_deriv%dom(idom)%nelements()
            do idom = 1,mesh%ndomains()
                do ielem = 1,mesh%domain(idom)%nelements()
                    
                    associate( int_derivs  => self%integral_deriv%dom(idom)%vecs(ielem) )
                        call chidg_signal(FATAL,"integral_cache%get_global_value: needs updated for petsc storage!")

                        ! Find correspondent derivatives in the AD_D vector
                        istart = iend + 1
                        iend   = iend + int_derivs%nentries()
                        
                        ! Set derivatives, check vec size
                        if ( int_derivs%nentries() == size(integral%xp_ad_(istart:iend)) ) then
                            integral%xp_ad_(istart:iend) = int_derivs%vec 
                        else
                            call chidg_signal(FATAL,"type_integral_cache.f90: the size of the integral derivatives does not match the worker derivatives size. Implementation Error!")
                        end if

                    end associate
                    
                end do !ielem
            end do !idom 
            
        end if            

    end function get_global_value
    !***************************************************************************************************




end module type_integral_cache
