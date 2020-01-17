module type_functional_cache
#include<messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, NO_ID, dX_DIFF, NO_DIFF, dD_DIFF
    use mod_io,                 only: backend
    use type_geometry_cache,    only: geometry_cache_t
    use type_function_info,     only: function_info_t
    use type_mesh,              only: mesh_t
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
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
        type(chidg_vector_t)                :: derivative_vector_template
        integer(ik)                         :: dtype
    
    contains
        
        ! Initialization
        procedure,  public  :: init
        
        ! Storage access
        procedure,  public  :: set_entity_value
        procedure,  public  :: set_global_value
        procedure,  public  :: get_entity_value
        procedure,  public  :: get_global_value
        procedure,  public  :: nentities

        ! Communication
        procedure,  public  :: comm

        procedure,  public  :: release

    end type functional_cache_t
    !***************************************************************************************************


contains

    
    
    
    !>  Initialize geometry cache 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init(self,geometry_type,mesh,geometries,integrals,integral_type,differentiate)
        class(functional_cache_t),  intent(inout)       :: self
        character(*),               intent(in)          :: geometry_type
        type(mesh_t),               intent(inout)       :: mesh
        type(svector_t),            intent(in)          :: geometries
        type(svector_t),            intent(in)          :: integrals
        character(*),               intent(in)          :: integral_type
        integer(ik),                intent(in)          :: differentiate

        self%derivative_vector_template = chidg_vector(trim(backend))

        ! Create chidg vector model
        if (differentiate == dX_DIFF) then
            call self%derivative_vector_template%init(mesh,1,dof_type='coordinate')
        else if (differentiate == dD_DIFF) then
            call self%derivative_vector_template%init(mesh,1,dof_type='auxiliary')
        else 
            call self%derivative_vector_template%init(mesh,1,dof_type='primal')
        end if
        call self%derivative_vector_template%assemble()
        
        ! Initialize spcific geometry cache
        select case (geometry_type)
            case("reference")
                call self%ref_cache%initialize(mesh,geometries,integrals,integral_type,differentiate,self%derivative_vector_template)
            case("auxiliary")
                call self%aux_cache%initialize(mesh,geometries,integrals,integral_type,differentiate,self%derivative_vector_template)
            case default
                call chidg_signal(FATAL,"functional_cache_t%init: incorrect input for 'geometry_type'.")
        end select

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
    subroutine set_entity_value(self,mesh,integral,int_name,geometry,fcn_info)
        class(functional_cache_t),              intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh
        type(AD_D),                             intent(in)      :: integral
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry
        type(function_info_t),                  intent(in)      :: fcn_info

        select case (geometry)
            case("reference")
                call self%ref_cache%set_entity_value(mesh,int_name,integral,fcn_info)
            case("auxiliary")
                call self%aux_cache%set_entity_value(mesh,int_name,integral,fcn_info)
            case default
                call chidg_signal(FATAL,"functional_cache_t%set_entity_value: incorrect 'geometry' reference.")
        end select
        
    end subroutine set_entity_value
    !***************************************************************************************************

    
    !>  Store integral based on the geometry   
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine set_global_value(self,mesh,integral,int_name,geometry)
        class(functional_cache_t),              intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh
        type(AD_D),                             intent(in)      :: integral
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry

        select case (geometry)
            case("reference")
                call self%ref_cache%set_global_value(mesh,int_name,integral,self%dtype)
            case("auxiliary")
                call self%aux_cache%set_global_value(mesh,int_name,integral,self%dtype)
            case default
                call chidg_signal(FATAL,"functional_cache_t%set_global_value: incorrect 'geometry' reference.")
        end select
        
    end subroutine set_global_value
    !***************************************************************************************************
    
    
    
    !>  Get integral based on the geometry   
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_entity_value(self,mesh,int_name,geometry,fcn_info) result(integral)
        class(functional_cache_t),              intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry
        type(function_info_t),                  intent(in)      :: fcn_info

        type(AD_D)      :: integral

        select case (geometry)
            case("reference")
                integral = self%ref_cache%get_entity_value(mesh,int_name,fcn_info)
            case("auxiliary")
                integral = self%aux_cache%get_entity_value(mesh,int_name,fcn_info)
            case default
                call chidg_signal(FATAL,"functional_cache_t%get_entity_value: incorrect 'geometry' reference.")
        end select
        
    end function get_entity_value
    !***************************************************************************************************



    !>  Get integral based on the geometry   
    !!
    !!  @author Matteo Ugolotti
    !!  @date   3/7/2019
    !!
    !---------------------------------------------------------------------------------------------------
    function get_global_value(self,mesh,int_name,geometry) result(integral)
        class(functional_cache_t),              intent(inout)   :: self
        type(mesh_t),                           intent(in)      :: mesh
        character(*),                           intent(in)      :: int_name
        character(*),                           intent(in)      :: geometry

        type(AD_D)      :: integral

        select case (geometry)
            case("reference")
                integral = self%ref_cache%get_global_value(mesh,int_name,self%dtype)
            case("auxiliary")
                integral = self%aux_cache%get_global_value(mesh,int_name,self%dtype)
            case default
                call chidg_signal(FATAL,"functional_cache_t%get_global_value: incorrect 'geometry' reference.")
        end select
        
    end function get_global_value
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
                call self%ref_cache%comm(self%dtype)
            case("auxiliary")
                call self%aux_cache%comm(self%dtype)
            case default
                call chidg_signal(FATAL,"functional_cache_t%comm: incorrect 'geometry' reference.")
        end select


    end subroutine comm
    !***************************************************************************************************





    !>  Release memory.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/16/2020
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine release(self)
        class(functional_cache_t),  intent(inout)   :: self

        call self%derivative_vector_template%release()
        call self%ref_cache%release()
        call self%aux_cache%release()

    end subroutine release
    !***************************************************************************************************




end module type_functional_cache
