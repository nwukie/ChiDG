module type_integral_cache
#include<messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, NO_ID, NO_DIFF, dX_DIFF, dQ_DIFF
    use mod_io,                 only: backend
    use type_mesh,              only: mesh_t
    use type_svector,           only: svector_t
    use mod_string,             only: string_t
    use type_element_info,      only: element_info_t
    use type_seed,              only: seed_t
    use mod_differentiate,      only: get_nderiv
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
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
        
        procedure, private  :: initialize
        generic,   public   :: init => initialize                       ! Initialize integral
        procedure           :: set_entity_value
        procedure           :: get_entity_value
        procedure           :: set_global_value
        procedure           :: get_global_value
        procedure           :: release

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
    subroutine initialize(self,i_name,derivative_vector_template)
        class(integral_cache_t),    intent(inout)   :: self
        character(*),               intent(in)      :: i_name
        type(chidg_vector_t),       intent(in)      :: derivative_vector_template
         
        ! Set name
        self%name = i_name

        ! Initialize real value of the integral
        self%integral_value = ZERO
        self%integral_deriv = chidg_vector(trim(backend))
        self%integral_deriv = derivative_vector_template
        call self%integral_deriv%assemble()

        ! Set initialized flag
        self%derivatives_initialized = .true.

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
    subroutine set_entity_value(self,mesh,integral,fcn_info) 
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(AD_D),                 intent(in)      :: integral
        type(function_info_t),      intent(in)      :: fcn_info
         
        type(element_info_t)    :: elem_info

        ! Accumulate integral real value, this has to be done for any kind of differentiation
        self%integral_value = self%integral_value + integral%x_ad_
    
        !   - dtype = dQ_DIFF or dX_DIFF or NO_DIFF
        !   - no need to store the chidg vector for NO_DIFF
        if (fcn_info%dtype /= NO_DIFF) then
            elem_info = mesh%get_element_info(fcn_info%seed%idomain_l,fcn_info%seed%ielement_l)
            call self%integral_deriv%add_fields(integral%xp_ad_,elem_info)
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
         
        type(element_info_t)    :: elem_info
        type(AD_D)              :: integral
        integer(ik)             :: nderivs

        call self%integral_deriv%assemble()

        ! Allocate result
        integral = AD_D(get_nderiv(fcn_info))
        
        ! Assign real value
        integral = self%integral_value

        ! Get element_info structure
        elem_info = mesh%get_element_info(fcn_info%seed%idomain_l, fcn_info%seed%ielement_l)
        
        ! Assign derivatives, if NO_DIFF then self.derivatives_initialized == .false.
        if (self%derivatives_initialized) then
            integral%xp_ad_ = self%integral_deriv%get_fields(elem_info)
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
    subroutine set_global_value(self,mesh,integral,dtype) 
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        type(AD_D),                 intent(in)      :: integral
        integer(ik),                intent(in)      :: dtype
        
        integer(ik)             :: istart, iend, idom, ielem
        type(element_info_t)    :: elem_info

        ! Set integral real value
        self%integral_value = integral%x_ad_

        if (dtype /= NO_DIFF) then
        
            ! Initialize istart and iend indeces
            istart = 0
            iend   = 0
        
            ! Loop thorugh domains and elements
            do idom = 1,mesh%ndomains()
                do ielem = 1,mesh%domain(idom)%nelements()

                    ! Find correspondent derivatives in the AD_D vector
                    elem_info = mesh%get_element_info(idom,ielem)

                    if (dtype == dQ_DIFF) then
                        istart = elem_info%dof_local_start
                        iend   = istart + elem_info%nfields*elem_info%nterms_s*elem_info%ntime - 1
                    else if (dtype == dX_DIFF) then
                        istart = iend + 1
                        iend   = istart + 3*elem_info%nterms_c*elem_info%ntime - 1
                    else
                        call chidg_signal_one(FATAL,"integral_cache%set_global_value: differentiation type not implemented.", dtype)
                    end if

                    call self%integral_deriv%set_fields(integral%xp_ad_(istart:iend),elem_info)
                    
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
    function get_global_value(self,mesh,dtype) result(integral)
        class(integral_cache_t),    intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        integer(ik),                intent(in)      :: dtype
         
        type(element_info_t)    :: elem_info
        type(AD_D)              :: integral
        integer(ik)             :: istart, iend, idom, ielem, nderivs


        call self%integral_deriv%assemble()

        ! Allocate result
        if (dtype /= NO_DIFF) then
            !! if differentiation is needed, check if the derivatives have been already initialized
            !! if not, initialize it.
            !! Derivatives might not be initialized if the processor does not have any face/element
            !! from the auxiliary/reference geometry, and it only knows the overall integral real value
            !! from the communication process.
            !if (.not. self%derivatives_initialized) then
            !    !self%integral_deriv = vec_model
            !    call self%integral_deriv%assemble()
            !end if
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
             
            do idom = 1,mesh%ndomains()
                do ielem = 1,mesh%domain(idom)%nelements()
                    
                    elem_info = mesh%get_element_info(idom,ielem)

                    if (dtype == dQ_DIFF) then
                        istart = elem_info%dof_local_start
                        iend   = istart + elem_info%nfields*elem_info%nterms_s*elem_info%ntime - 1
                    else if (dtype == dX_DIFF) then
                        istart = iend + 1
                        iend   = istart + 3*elem_info%nterms_c*elem_info%ntime - 1
                    else
                        call chidg_signal_one(FATAL,"integral_cache%set_global_value: differentiation type not implemented.", dtype)
                    end if

                    integral%xp_ad_(istart:iend) = self%integral_deriv%get_fields(elem_info)

                end do !ielem
            end do !idom 
            
        end if            


    end function get_global_value
    !***************************************************************************************************




    !>  Release memory.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/16/2020
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine release(self)
        class(integral_cache_t),    intent(inout)   :: self

        call self%integral_deriv%release()

    end subroutine release
    !***************************************************************************************************




end module type_integral_cache
