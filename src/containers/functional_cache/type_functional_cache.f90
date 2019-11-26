module type_functional_cache
#include<messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, NO_ID, dQ_DIFF, dX_DIFF, NO_DIFF, dQ_DIFF
    use type_geometry_cache,    only: geometry_cache_t
    use type_function_info,     only: function_info_t
    use type_mesh,              only: mesh_t
    use type_chidg_vector,      only: chidg_vector_t
    use type_svector,           only: svector_t
    use DNAD_D
    implicit none


    !>  Container to store functional integrals for both auxiliary and reference geometry
    !!  It handles high-level operation for storing funtional integrals.
    !!   
    !!  Contains essential access operation used by each single functional
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    type, public :: functional_cache_t
        
        ! Storage of each integral to compute on the geometry 
        type(geometry_cache_t)              :: ref_cache        ! Integrals computed on the reference geom
        type(geometry_cache_t)              :: aux_cache        ! Integrals computed on the auxiliary geom

        ! Info needed for storing derivatives
        type(chidg_vector_t)                :: vector_model
        integer(ik)                         :: dtype
    
    contains
        
        ! Initialization
        procedure,  public  :: init
        
        ! Storage access
        procedure,  public  :: set_value
        procedure,  public  :: get_value
        procedure,  public  :: nentities

        ! Communication
        procedure,  public  :: comm

    end type functional_cache_t
    !***************************************************************************************************


contains

    
    
    
    !>  Initialize geometry cache 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init(self,geometry,mesh,geometries,integral_type,differentiate)
        class(functional_cache_t),  intent(inout)       :: self
        character(*),               intent(in)          :: geometry
        type(mesh_t),               intent(inout)       :: mesh
        type(svector_t),            intent(in)          :: geometries
        character(*),               intent(in)          :: integral_type
        integer(ik),                intent(in)          :: differentiate

        
        ! Initialize spcific geometry cache
        select case (geometry)
            case("reference")
                call self%ref_cache%initialize(mesh,geometries,integral_type,differentiate)

            case("auxiliary")
                call self%aux_cache%initialize(mesh,geometries,integral_type,differentiate)

            case default
                call chidg_signal(FATAL,"functional_cache_t%init: incorrect 'geometry' reference.")
        end select

        ! Create chidg vector model
        call self%vector_model%init(mesh,1,differentiate)

        ! Differentiation type
        self%dtype = differentiate

    end subroutine init
    !***************************************************************************************************




    
    
    !>  Store integral based on the geometry   
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_value(self,integral,int_name,geometry,fcn_info)
        class(functional_cache_t),              intent(inout)   :: self
        type(AD_D),                             intent(in)      :: integral
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry
        type(function_info_t),      optional,   intent(in)      :: fcn_info

        select case (geometry)
            case("reference")
                call self%ref_cache%set_value(int_name,integral,self%vector_model,self%dtype,fcn_info)

            case("auxiliary")
                call self%aux_cache%set_value(int_name,integral,self%vector_model,self%dtype,fcn_info)

            case default
                call chidg_signal(FATAL,"functional_cache_t%set_value: incorrect 'geometry' reference.")
        end select

        
    end subroutine set_value
    !***************************************************************************************************

    
    
    
    
    !>  Get integral based on the geometry   
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_value(self,int_name,geometry,fcn_info) result(integral)
        class(functional_cache_t),              intent(inout)   :: self
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry
        type(function_info_t),      optional,   intent(in)      :: fcn_info

        type(AD_D)      :: integral

        select case (geometry)
            case("reference")
                integral = self%ref_cache%get_value(int_name,self%vector_model,self%dtype,fcn_info)

            case("auxiliary")
                integral = self%aux_cache%get_value(int_name,self%vector_model,self%dtype,fcn_info)

            case default
                call chidg_signal(FATAL,"functional_cache_t%get_value: incorrect 'geometry' reference.")
        end select
        
    end function get_value
    !***************************************************************************************************


    
    
    
    !>  Return the number of entities (faces/elements) in the geometry 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function nentities(self,geometry) result(nentities_)
        class(functional_cache_t),   intent(inout)   :: self
        character(*),                intent(in)      :: geometry

        integer(ik)     :: nentities_

        select case (geometry)
            case("reference")
                nentities_ = self%ref_cache%nentities

            case("auxiliary")
                nentities_ = self%aux_cache%nentities

            case default
                call chidg_signal(FATAL,"functional_cache_t%nentities: incorrect 'geometry' reference.")
        end select
        
    end function nentities
    !***************************************************************************************************


    
    
    
    
    !>  Communication geometry cache
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine comm(self,geometry)
        class(functional_cache_t),  intent(inout)       :: self
        character(*),               intent(in)          :: geometry

        ! Initialize spcific geometry cache
        select case (geometry)
            case("reference")
                call self%ref_cache%comm(self%vector_model,self%dtype)

            case("auxiliary")
                call self%aux_cache%comm(self%vector_model,self%dtype)

            case default
                call chidg_signal(FATAL,"functional_cache_t%comm: incorrect 'geometry' reference.")
        end select


    end subroutine comm
    !***************************************************************************************************


end module type_functional_cache
