!-------------------------------------------------------------------------------------
!!
!!                                  A ChiDG Worker
!!
!!  Purpose:
!!  ----------------------------------------
!!  The chidg_worker_t handles the following activities that might occur within
!!  an operator_t:
!!      - interpolate to quadrature nodes. Element and face sets.
!!      - integrate. Volume and Boundaries
!!      - return geometric information such as normals, and coordinates.
!!
!!  The worker knows what element/face is currently being worked on. So, it can then
!!  access that element/face for getting data, performing the correct interpolation,
!!  performing the correct integral. This way, the operator_t flux routines don't
!!  have to worry about where they are getting data from. 
!!
!!  The operator_t's are just concerned with getting information from the worker 
!!  about the element/face, computing a function value, passing that information
!!  back to the worker to be integrated and stored.
!!
!!
!!
!!  @author Nathan A. Wukie
!!  @date   8/22/2016
!!
!!
!-------------------------------------------------------------------------------------
module type_chidg_worker
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: NFACES, ME, NEIGHBOR, BC, ZERO, CHIMERA,  &
                                      ONE, THIRD, TWO, NOT_A_FACE, BOUNDARY,    &
                                      CARTESIAN, CYLINDRICAL, INTERIOR

    use mod_interpolate,        only: interpolate_element_autodiff
    use mod_integrate,          only: integrate_boundary_scalar_flux, &
                                      integrate_volume_vector_flux,   &
                                      integrate_volume_scalar_source

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_element_info,      only: element_info_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none





    !>  The ChiDG worker implementation.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: chidg_worker_t
    
        type(mesh_t),           pointer :: mesh
        type(solverdata_t),     pointer :: solverdata
        type(chidg_cache_t),    pointer :: cache
        !type(properties_t),     pointer :: prop(:)
        type(properties_t), allocatable :: prop(:)

        integer(ik)                 :: iface
        integer(ik)                 :: itime
        type(element_info_t)        :: element_info
        type(function_info_t)       :: function_info
    
        character(:),   allocatable :: interpolation_source
        real(rk)                    :: t    ! Physical time

    contains 
    
        ! Worker state
        procedure   :: init
        procedure   :: set_element          ! Set element_info type
        procedure   :: set_function_info    ! Set function_info type
        procedure   :: set_face             ! Set iface index
        procedure   :: face_info            ! Return a face_info type


        ! Worker get data
        procedure   :: get_primary_field_general
        procedure   :: get_primary_field_face
        procedure   :: get_primary_field_element
        procedure   :: get_model_field_general
        procedure   :: get_model_field_face
        procedure   :: get_model_field_element
        procedure   :: get_auxiliary_field_general
        procedure   :: get_auxiliary_field_face
        procedure   :: get_auxiliary_field_element

        procedure   :: get_field
        !procedure   :: get_field_gradient

        procedure   :: store_bc_state
        procedure   :: store_model_field

        procedure   :: normal
        procedure   :: unit_normal

        procedure   :: coords
        procedure   :: x
        procedure   :: y
        procedure   :: z
        procedure   :: coordinate

        procedure   :: element_size
        procedure   :: solution_order
        procedure   :: quadrature_weights
        procedure   :: inverse_jacobian
        procedure   :: face_area
        procedure   :: coordinate_system

        procedure   :: face_type

        procedure   :: time


        ! Worker process data
        procedure   :: integrate_boundary
        generic     :: integrate_volume => integrate_volume_flux, &
                                           integrate_volume_source
        procedure   :: integrate_volume_flux
        procedure   :: integrate_volume_source
        
        !ALE procedures
        procedure   :: get_grid_velocity_element
        procedure   :: get_grid_velocity_face
        procedure   :: get_jacobian_grid_element
        procedure   :: get_inv_jacobian_grid_element
        procedure   :: get_jacobian_grid_face
        procedure   :: get_inv_jacobian_grid_face
        procedure   :: get_det_jacobian_grid_element
        procedure   :: get_det_jacobian_grid_face


        procedure   :: get_primary_field_value_ale_element
        procedure   :: get_primary_field_grad_ale_element
        procedure   :: get_primary_field_value_ale_face
        procedure   :: get_primary_field_grad_ale_face
        procedure   :: get_primary_field_value_ale_general
        procedure   :: get_primary_field_grad_ale_general

        procedure   :: post_process_volume_advective_flux_ale
        procedure   :: post_process_boundary_advective_flux_ale
        procedure   :: post_process_volume_diffusive_flux_ale
        procedure   :: post_process_boundary_diffusive_flux_ale
        final       :: destructor
    
    end type chidg_worker_t
    !*********************************************************************************






contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self,mesh,prop,solverdata,cache)
        class(chidg_worker_t),  intent(inout)       :: self
        type(mesh_t),       intent(in), target  :: mesh
        type(properties_t),     intent(in), target  :: prop(:)
        type(solverdata_t),     intent(in), target  :: solverdata
        type(chidg_cache_t),    intent(in), target  :: cache

        character(:),   allocatable :: temp_name

        self%mesh       => mesh
        ! having issue with using a pointer here for prop. Theory is that the compiler
        ! creates a temporary array of prop(:) from eqnset(:)%prop when it is passing it in. 
        ! Then after this routine exists, that array ceases to exist and so
        ! points to nothing. For now we will just assign, but probably want this
        ! linked back up in the future.
        self%prop       =  prop
        !self%prop       => prop
        self%solverdata => solverdata
        self%cache      => cache

    end subroutine init
    !**********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine set_element(self,elem_info)
        class(chidg_worker_t),  intent(inout)   :: self
        type(element_info_t),   intent(in)      :: elem_info

        self%element_info = elem_info

    end subroutine set_element
    !**********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine set_function_info(self,fcn_info)
        class(chidg_worker_t),  intent(inout)   :: self
        type(function_info_t),  intent(in)      :: fcn_info

        self%function_info = fcn_info

    end subroutine set_function_info
    !**********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine set_face(self,iface)
        class(chidg_worker_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: iface

        self%iface = iface

    end subroutine set_face
    !***************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function face_info(self) result(face_info_)
        class(chidg_worker_t),  intent(in)  :: self
        
        type(face_info_t)   :: face_info_

        face_info_%idomain_g  = self%element_info%idomain_g
        face_info_%idomain_l  = self%element_info%idomain_l
        face_info_%ielement_g = self%element_info%ielement_g
        face_info_%ielement_l = self%element_info%ielement_l
        face_info_%iface      = self%iface

    end function face_info
    !***************************************************************************************











    !>  Return a primary field evaluated at a quadrature node set. The source here
    !!  is determined by chidg_worker.
    !!
    !!  This routine is specifically for model_t's, because we want them to be evaluated
    !!  on face and element sets the same way. So in a model implementation, we just
    !!  want the model to get some quadrature node set to operate on. The chidg_worker
    !!  handles what node set is currently being returned.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_primary_field_general(self,field,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        type(AD_D), allocatable :: var_gq(:)


        if (self%interpolation_source == 'element') then
            var_gq = self%get_primary_field_element(field,interp_type) 
        else if (self%interpolation_source == 'face interior') then
            var_gq = self%get_primary_field_face(field,interp_type,'face interior')
        else if (self%interpolation_source == 'face exterior') then
            var_gq = self%get_primary_field_face(field,interp_type,'face exterior')
        else if (self%interpolation_source == 'boundary') then
            var_gq = self%get_primary_field_face(field,interp_type,'boundary')
        end if

    end function get_primary_field_general
    !**************************************************************************************








    !>  Return an auxiliary field evaluated at a quadrature node set. The source here
    !!  is determined by chidg_worker.
    !!
    !!  This routine is specifically for model_t's, because we want them to be evaluated
    !!  on face and element sets the same way. So in a model implementation, we just
    !!  want the model to get some quadrature node set to operate on. The chidg_worker
    !!  handles what node set is currently being returned.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/7/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_auxiliary_field_general(self,field,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        !type(AD_D), allocatable :: var_gq(:)
        real(rk), allocatable :: var_gq(:)


        if (self%interpolation_source == 'element') then
            var_gq = self%get_auxiliary_field_element(field,interp_type) 
        else if (self%interpolation_source == 'face interior') then
            var_gq = self%get_auxiliary_field_face(field,interp_type,'face interior')
        else if (self%interpolation_source == 'face exterior') then
            var_gq = self%get_auxiliary_field_face(field,interp_type,'face exterior')
        else if (self%interpolation_source == 'boundary') then
            var_gq = self%get_auxiliary_field_face(field,interp_type,'boundary')
        end if

    end function get_auxiliary_field_general
    !**************************************************************************************















    !>  Return a primary field interpolated to a face quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_primary_field_face(self,field,interp_type,interp_source) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type
        character(*),           intent(in)  :: interp_source

        type(AD_D),     allocatable, dimension(:)   :: var_gq
        character(:),   allocatable                 :: cache_component, cache_type, user_msg
        type(face_info_t)                           :: face_info
        integer(ik)                                 :: idirection, igq



        !
        ! Set cache_component
        !
        if (interp_source == 'face interior') then
            cache_component = 'face interior'
        else if (interp_source == 'face exterior' .or. &
                 interp_source == 'boundary') then
            cache_component = 'face exterior'
        else
            user_msg = "chidg_worker%get_primary_field_face: Invalid value for interpolation source. &
                        Try 'face interior', 'face exterior', or 'boundary'"
            call chidg_signal_two(FATAL,user_msg,trim(field),trim(interp_source))
        end if


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (interp_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (interp_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3


        else if ( (interp_type == 'grad1 + lift') .or. &
                  (interp_type == 'grad1+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 1
        else if ( (interp_type == 'grad2 + lift') .or. &
                  (interp_type == 'grad2+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 2
        else if ( (interp_type == 'grad3 + lift') .or. &
                  (interp_type == 'grad3+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 3
        else
            user_msg = "chidg_worker%get_primary_field_face: Invalid interpolation &
                        type. 'value', 'grad1', 'grad2', 'grad3', 'grad1+lift', 'grad2+lift', 'grad3+lift'"
            call chidg_signal(FATAL,user_msg)
        end if



        !
        ! Retrieve data from cache
        !
        if (cache_type == 'value') then
            var_gq = self%cache%get_data(field,cache_component,'value',idirection,self%function_info%seed,self%iface)

        else if (cache_type == 'gradient') then
            var_gq = self%cache%get_data(field,cache_component,'gradient',idirection,self%function_info%seed,self%iface)

        else if (cache_type == 'gradient + lift') then
            var_gq = self%cache%get_data(field,cache_component,'gradient',idirection,self%function_info%seed,self%iface)

            ! Modify derivative by face lift stabilized by a factor of NFACES
            var_gq = var_gq + real(NFACES,rk)*self%cache%get_data(field,cache_component,'lift face',idirection,self%function_info%seed,self%iface)

        end if


    end function get_primary_field_face
    !***************************************************************************************










    !>  Return a primary field interpolated to an element quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_primary_field_element(self,field,interp_type,Pmin,Pmax) result(var_gq)
        class(chidg_worker_t),  intent(in)              :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: interp_type
        integer(ik),            intent(in), optional    :: Pmin
        integer(ik),            intent(in), optional    :: Pmax

        type(AD_D),     allocatable, dimension(:) :: var_gq, tmp_gq

        type(face_info_t)               :: face_info
        character(:),   allocatable     :: cache_component, cache_type, user_msg
        integer(ik)                     :: idirection, igq, iface, ifield, idomain_l, eqn_ID


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (interp_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (interp_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3


        else if ( (interp_type == 'grad1 + lift') .or. &
                  (interp_type == 'grad1+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 1
        else if ( (interp_type == 'grad2 + lift') .or. &
                  (interp_type == 'grad2+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 2
        else if ( (interp_type == 'grad3 + lift') .or. &
                  (interp_type == 'grad3+lift'  ) ) then
            cache_type = 'gradient + lift'
            idirection = 3
        else
            user_msg = "chidg_worker%get_primary_field_element: Invalid interpolation &
                        type. 'value', 'grad1', 'grad2', 'grad3', 'grad1+lift', 'grad2+lift', 'grad3+lift'"
            call chidg_signal(FATAL,user_msg)
        end if




        !
        ! If we are requesting an interpolation from a subset of the modal expansion, 
        ! then perform a new interpolation rather than using the cache.
        !
        if (present(Pmin) .or. present(Pmax)) then

            if (cache_type == 'value') then
                idomain_l = self%element_info%idomain_l
                eqn_ID    = self%mesh%domain(idomain_l)%eqn_ID
                ifield    = self%prop(eqn_ID)%get_primary_field_index(field)

                var_gq = interpolate_element_autodiff(self%mesh, self%solverdata%q, self%element_info, self%function_info, ifield, self%itime, interp_type, Pmin, Pmax)

            else if ( (cache_type == 'gradient') .or. &
                      (cache_type == 'gradient + lift') ) then
                user_msg = "chidg_worker%get_primary_field_element: On partial field interpolations, &
                            only the 'value' of the field can be interpolated. 'gradient' is not yet implemented."
                call chidg_signal(FATAL,user_msg)
            end if

    

        else


            !
            ! Retrieve data from cache
            !
            if ( cache_type == 'value') then
                var_gq = self%cache%get_data(field,'element','value',idirection,self%function_info%seed)

            else if (cache_type == 'gradient') then
                var_gq = self%cache%get_data(field,'element','gradient',idirection,self%function_info%seed)

            else if (cache_type == 'gradient + lift') then
                var_gq = self%cache%get_data(field,'element','gradient',idirection,self%function_info%seed)

                ! Add lift contributions from each face
                do iface = 1,NFACES
                    tmp_gq = self%cache%get_data(field,'face interior', 'lift element', idirection, self%function_info%seed,iface)
                    var_gq = var_gq + tmp_gq
                end do

            end if


        end if


    end function get_primary_field_element
    !****************************************************************************************










    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/10/2017
    !!
    !---------------------------------------------------------------------------------------
    function get_field(self,field,interp_type,interp_source) result(var_gq)
        class(chidg_worker_t),  intent(in)              :: self
        character(*),           intent(in)              :: field
        character(*),           intent(in)              :: interp_type
        character(*),           intent(in), optional    :: interp_source

        type(AD_D),     allocatable :: var_gq(:), tmp_gq(:)
        character(:),   allocatable :: cache_component, cache_type, lift_source, lift_nodes, user_msg
        integer(ik)                 :: lift_face_min, lift_face_max, idirection, iface
        real(rk)                    :: stabilization
        logical                     :: no_lift


        !
        ! Set cache_component
        !
        if (present(interp_source)) then
            select case(trim(interp_source))
                case('face interior') 
                    cache_component = 'face interior'
                    lift_source     = 'face interior'
                    lift_nodes      = 'lift face'
                    lift_face_min   = self%iface
                    lift_face_max   = self%iface
                    stabilization   = real(NFACES,rk)
                case('face exterior','boundary')
                    cache_component = 'face exterior'
                    lift_source     = 'face exterior'
                    lift_nodes      = 'lift face'
                    lift_face_min   = self%iface
                    lift_face_max   = self%iface
                    stabilization   = real(NFACES,rk)
                case('element')
                    cache_component = 'element'
                    lift_source     = 'face interior'
                    lift_nodes      = 'lift element'
                    lift_face_min   = 1
                    lift_face_max   = NFACES
                    stabilization   = ONE
                case default
                    user_msg = "chidg_worker%get_field: Invalid value for interpolation source. &
                                Try 'face interior', 'face exterior', 'boundary', or 'element'"
                    call chidg_signal_one(FATAL,user_msg,trim(interp_source))
            end select
        else
            cache_component = self%interpolation_source
            if ( (trim(cache_component) /= "face interior") .and. &
                 (trim(cache_component) /= "face exterior") .and. &
                 (trim(cache_component) /= "element") ) then
            user_msg = "chidg_worker%get_field: chidg_worker implicit interpolation source is not valid."
            call chidg_signal(FATAL,user_msg)
            end if
        end if




        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (interp_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (interp_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3
        else
            user_msg = "chidg_worker%get_field: Invalid interpolation &
                        type. 'value', 'grad1', 'grad2', 'grad3'"
            call chidg_signal(FATAL,user_msg)
        end if



        !
        ! Determine when we do not want to lift
        !
        ! Do not lift gradient for a boundary state function. Boundary state functions
        ! interpolate from the 'face interior'. Operators interpolate from 'boundary'.
        ! So, if we are on a BOUNDARY face and interpolating from 'face interior', then
        ! we don't want to lift because there is no lift for the boundary function to use
        ! If we aren't on a face, face_type returns NOT_A_FACE, so this is still valid for 
        ! returning element data.
        no_lift = (self%face_type() == BOUNDARY) .and. (cache_component == 'face interior')


        !
        ! Retrieve data from cache
        !
        if ( cache_type == 'value') then
            var_gq = self%cache%get_data(field,cache_component,'value',idirection,self%function_info%seed,self%iface)

        else if (cache_type == 'gradient') then

            if (self%cache%lift .and. (.not. no_lift)) then
                var_gq = self%cache%get_data(field,cache_component,'gradient',idirection,self%function_info%seed,self%iface)

                ! Add lift contributions from each face
                do iface = lift_face_min,lift_face_max
                    tmp_gq = self%cache%get_data(field,lift_source, lift_nodes, idirection, self%function_info%seed,iface)
                    var_gq = var_gq + stabilization*tmp_gq
                end do

            else
                var_gq = self%cache%get_data(field,cache_component,'gradient',idirection,self%function_info%seed,self%iface)
            end if

        else
            user_msg = "chidg_worker%get_field: invalid cache_type."
            call chidg_signal(FATAL,user_msg)

        end if



    end function get_field
    !***************************************************************************************










    !>  Return a model field evaluated at a quadrature node set. The source here
    !!  is determined by chidg_worker.
    !!
    !!  This routine is specifically for model_t's, because we want them to be evaluated
    !!  on face and element sets the same way. So in a model implementation, we just
    !!  want the model to get some quadrature node set to operate on. The chidg_worker
    !!  handles what node set is currently being returned.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_model_field_general(self,field,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        type(AD_D), allocatable :: var_gq(:)


        if (self%interpolation_source == 'element') then
            var_gq = self%get_model_field_element(field,interp_type) 
        else if (self%interpolation_source == 'face interior') then
            var_gq = self%get_model_field_face(field,interp_type,'face interior')
        else if (self%interpolation_source == 'face exterior') then
            var_gq = self%get_model_field_face(field,interp_type,'face exterior')
        else if (self%interpolation_source == 'boundary') then
            var_gq = self%get_model_field_face(field,interp_type,'boundary')
        end if

    end function get_model_field_general
    !**************************************************************************************










    !>  Return a primary field interpolated to a face quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_model_field_face(self,field,interp_type,interp_source) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type
        character(*),           intent(in)  :: interp_source

        type(AD_D),     allocatable, dimension(:)   :: var_gq
        character(:),   allocatable                 :: cache_component, cache_type, user_msg
        type(face_info_t)                           :: face_info
        integer(ik)                                 :: idirection, igq



        !
        ! Set cache_component
        !
        if (interp_source == 'face interior') then
            cache_component = 'face interior'
        else if (interp_source == 'face exterior' .or. &
                 interp_source == 'boundary') then
            cache_component = 'face exterior'
        else
            user_msg = "chidg_worker%get_model_field_face: Invalid value for interpolation source. &
                        Try 'face interior', 'face exterior', or 'boundary'"
            call chidg_signal_one(FATAL,user_msg,trim(interp_source))
        end if


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if ( (interp_type == 'grad1')          .or. &
                  (interp_type == 'grad2')          .or. &
                  (interp_type == 'grad3')          .or. &
                  (interp_type == 'grad1 + lift')   .or. &
                  (interp_type == 'grad1+lift'  )   .or. &
                  (interp_type == 'grad2 + lift')   .or. &
                  (interp_type == 'grad2+lift'  )   .or. &
                  (interp_type == 'grad3 + lift')   .or. &
                  (interp_type == 'grad3+lift'  ) ) then
                user_msg = 'chidg_worker%get_model_field_face: Computing gradients for model &
                            fields is not yet implemented.'
                call chidg_signal(FATAL,user_msg)
                                    
        end if



        !
        ! Retrieve data from cache
        !
        var_gq = self%cache%get_data(field,cache_component,cache_type,idirection,self%function_info%seed,self%iface)


    end function get_model_field_face
    !***************************************************************************************












    !>  Return a primary field interpolated to an element quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_model_field_element(self,field,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        type(AD_D),     allocatable, dimension(:) :: var_gq

        type(face_info_t)               :: face_info
        character(:),   allocatable     :: cache_component, cache_type, user_msg
        integer(ik)                     :: idirection, igq, iface


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if ( (interp_type == 'grad1')        .or. &
                  (interp_type == 'grad2')        .or. &
                  (interp_type == 'grad3')        .or. &
                  (interp_type == 'grad1+lift'  ) .or. &
                  (interp_type == 'grad2+lift'  ) .or. &
                  (interp_type == 'grad3+lift'  ) .or. &
                  (interp_type == 'grad1 + lift') .or. &
                  (interp_type == 'grad2 + lift') .or. &
                  (interp_type == 'grad3 + lift') ) then
            user_msg = 'chidg_worker%get_model_field_element: Computing gradients for model &
                        fields is not yet implemented.'
            call chidg_signal(FATAL,user_msg)
        end if




        !
        ! Retrieve data from cache
        !
        var_gq = self%cache%get_data(field,'element',cache_type,idirection,self%function_info%seed)



    end function get_model_field_element
    !****************************************************************************************











    !>  Return an auxiliary field interpolated to a face quadrature node set.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   12/7/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_auxiliary_field_face(self,field,interp_type,interp_source) result(var_gq_real)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type
        character(*),           intent(in)  :: interp_source

        type(AD_D),     allocatable, dimension(:)   :: var_gq
        real(rk),       allocatable, dimension(:)   :: var_gq_real
        character(:),   allocatable                 :: cache_component, cache_type, user_msg
        type(face_info_t)                           :: face_info
        integer(ik)                                 :: idirection, igq



        !
        ! Set cache_component
        !
        if (interp_source == 'face interior') then
            cache_component = 'face interior'
        else if (interp_source == 'face exterior' .or. &
                 interp_source == 'boundary') then
            cache_component = 'face exterior'
        else
            user_msg = "chidg_worker%get_auxiliary_field_face: Invalid value for interpolation source. &
                        Try 'face interior', 'face exterior', or 'boundary'"
            call chidg_signal_one(FATAL,user_msg,trim(interp_source))
        end if


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (interp_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (interp_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3

        else if ( (interp_type == 'grad1+lift'  ) .or. &
                  (interp_type == 'grad2+lift'  ) .or. &
                  (interp_type == 'grad3+lift'  ) .or. &
                  (interp_type == 'grad1 + lift') .or. &
                  (interp_type == 'grad2 + lift') .or. &
                  (interp_type == 'grad3 + lift') ) then
                user_msg = 'chidg_worker%get_auxiliary_field_face: Computing lifted derivatives for auxiliary &
                            fields is not supported.'
                call chidg_signal(FATAL,user_msg)
                                    
        end if



        !
        ! Retrieve data from cache
        !
        var_gq = self%cache%get_data(field,cache_component,cache_type,idirection,self%function_info%seed,self%iface)


        !
        ! Return a real array
        !
        var_gq_real = var_gq(:)%x_ad_


    end function get_auxiliary_field_face
    !***************************************************************************************










    !>  Return an auxiliary field interpolated to an element quadrature node set.
    !!
    !!  NOTE: Returns a real array
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   12/7/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_auxiliary_field_element(self,field,interp_type) result(var_gq_real)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        type(AD_D),     allocatable, dimension(:)   :: var_gq
        real(rk),       allocatable, dimension(:)   :: var_gq_real

        type(face_info_t)               :: face_info
        character(:),   allocatable     :: cache_component, cache_type, user_msg
        integer(ik)                     :: idirection, igq, iface


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (interp_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (interp_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3


        else if ( (interp_type == 'grad1+lift'  ) .or. &
                  (interp_type == 'grad2+lift'  ) .or. &
                  (interp_type == 'grad3+lift'  ) .or. &
                  (interp_type == 'grad1 + lift') .or. &
                  (interp_type == 'grad2 + lift') .or. &
                  (interp_type == 'grad3 + lift') ) then

            user_msg = 'chidg_worker%get_auxiliary_field_element: Computing lifted gradients for auxiliary &
                        fields is not supported.'
            call chidg_signal(FATAL,user_msg)

        end if




        !
        ! Retrieve data from cache
        !
        var_gq = self%cache%get_data(field,'element',cache_type,idirection,self%function_info%seed)



        !
        ! Copy real values to be returned.
        !
        var_gq_real = var_gq(:)%x_ad_


    end function get_auxiliary_field_element
    !****************************************************************************************














    !>  Store a primary field being defined from a boundary condition state function
    !!  to the 'face exterior' cache component.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine store_bc_state(self,field,cache_data,data_type)
        class(chidg_worker_t),  intent(inout)   :: self
        character(*),           intent(in)      :: field
        type(AD_D),             intent(in)      :: cache_data(:)
        character(*),           intent(in)      :: data_type

        character(:),   allocatable :: cache_type, user_msg
        integer(ik)                 :: idirection


        !
        ! Set cache_type
        !
        if (data_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (data_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (data_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (data_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3
        else
            user_msg = "chidg_worker%store_bc_state: Invalid data_type specification. &
                        Options are 'value', 'grad1', 'grad2', 'grad3'."
            call chidg_signal_one(FATAL,user_msg,trim(data_type))
        end if



        !
        ! Store bc state in cache, face exterior component
        !
        if (cache_type == 'value') then
            call self%cache%set_data(field,'face exterior',cache_data,'value',0,self%function_info%seed,self%iface)

        else if (cache_type == 'gradient') then
            call self%cache%set_data(field,'face exterior',cache_data,'gradient',idirection,self%function_info%seed,self%iface)

        end if



    end subroutine store_bc_state
    !***************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine store_model_field(self,model_field,data_type,cache_data)
        class(chidg_worker_t),  intent(inout)   :: self
        character(*),           intent(in)      :: model_field
        character(*),           intent(in)      :: data_type
        type(AD_D),             intent(in)      :: cache_data(:)

        type(AD_D),     allocatable, dimension(:)   :: field_current, field_update
        character(:),   allocatable :: cache_type, user_msg
        integer(ik)                 :: idirection


        !
        ! Set cache_type
        !
        if (data_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (data_type == 'grad1') then
            cache_type = 'gradient'
            idirection = 1
        else if (data_type == 'grad2') then
            cache_type = 'gradient'
            idirection = 2
        else if (data_type == 'grad3') then
            cache_type = 'gradient'
            idirection = 3
        else
            user_msg = "chidg_worker%store_model_field: Invalid data_type specification. &
                        Options are 'value', 'grad1', 'grad2', 'grad3'."
            call chidg_signal_one(FATAL,user_msg,trim(data_type))
        end if



        !
        ! Add data to model field cache
        !
        if (cache_type == 'value') then

            field_update = cache_data

            call self%cache%set_data(model_field,self%interpolation_source,field_update,'value',0,self%function_info%seed,self%iface)

        else if (cache_type == 'gradient') then

            field_update = cache_data

            call self%cache%set_data(model_field,self%interpolation_source,field_update,'gradient',idirection,self%function_info%seed,self%iface)

        end if



    end subroutine store_model_field
    !***************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine integrate_boundary(self,primary_field,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        character(*),           intent(in)      :: primary_field
        type(AD_D),             intent(inout)   :: integrand(:)

        integer(ik) :: ifield, idomain_l, eqn_ID

        idomain_l = self%element_info%idomain_l
        eqn_ID    = self%mesh%domain(idomain_l)%eqn_ID
        ifield    = self%prop(eqn_ID)%get_primary_field_index(primary_field)

        call integrate_boundary_scalar_flux(self%mesh,self%solverdata,self%face_info(),self%function_info,ifield,self%itime,integrand)


    end subroutine integrate_boundary
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine integrate_volume_flux(self,primary_field,integrand_x,integrand_y,integrand_z)
        class(chidg_worker_t),  intent(in)      :: self
        character(*),           intent(in)      :: primary_field
        type(AD_D),             intent(inout)   :: integrand_x(:)
        type(AD_D),             intent(inout)   :: integrand_y(:)
        type(AD_D),             intent(inout)   :: integrand_z(:)

        integer(ik) :: ifield, idomain_l, eqn_ID


        idomain_l = self%element_info%idomain_l
        eqn_ID    = self%mesh%domain(idomain_l)%eqn_ID
        ifield    = self%prop(eqn_ID)%get_primary_field_index(primary_field)

        call integrate_volume_vector_flux(self%mesh,self%solverdata,self%element_info,self%function_info,ifield,self%itime,integrand_x,integrand_y,integrand_z)


    end subroutine integrate_volume_flux
    !***************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine integrate_volume_source(self,primary_field,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        character(*),           intent(in)      :: primary_field
        type(AD_D),             intent(inout)   :: integrand(:)

        integer(ik) :: ifield, idomain_l, eqn_ID


        idomain_l = self%element_info%idomain_l
        eqn_ID    = self%mesh%domain(idomain_l)%eqn_ID
        ifield    = self%prop(eqn_ID)%get_primary_field_index(primary_field)

        call integrate_volume_scalar_source(self%mesh,self%solverdata,self%element_info,self%function_info,ifield,self%itime,integrand)


    end subroutine integrate_volume_source
    !***************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function normal(self,direction) result(norm_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: direction

        real(rk), dimension(:), allocatable :: norm_gq

        norm_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%norm(:,direction)

    end function normal
    !***************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function unit_normal(self,direction) result(unorm_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: direction

        real(rk), dimension(:), allocatable :: unorm_gq


        unorm_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%unorm(:,direction)


    end function unit_normal
    !***************************************************************************************






    !>  Return physical coordinates at the support nodes.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function coords(self) result(coords_)
        class(chidg_worker_t),  intent(in)  :: self

        type(point_t), allocatable, dimension(:) :: coords_(:)
        !real(rk), allocatable :: coords_support(:,:)

        !coords_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)
        ! Use constructor to return an array of point_t's from an array of reals
        coords_ = point_t(self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts)

    end function coords
    !***************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function x(self,source) result(x_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: x_gq

        if (source == 'boundary') then
            !x_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c1_
            x_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,1)
        else if (source == 'volume') then
            !x_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c1_
            x_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,1)
        else
            call chidg_signal(FATAL,"chidg_worker%x(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if



    end function x
    !**************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function y(self,source) result(y_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: y_gq


        if (source == 'boundary') then
            !y_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c2_
            y_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,2)
        else if (source == 'volume') then
            !y_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c2_
            y_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,2)
        else
            call chidg_signal(FATAL,"chidg_worker%y(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if



    end function y
    !**************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function z(self,source) result(z_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: z_gq


        if (source == 'boundary') then
            !z_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c3_
            z_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,3)
        else if (source == 'volume') then
            !z_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c3_
            z_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,3)
        else
            call chidg_signal(FATAL,"chidg_worker%z(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if


    end function z
    !**************************************************************************************








    !>  Interface for returning coordinates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/16/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    function coordinate(self,string,user_source) result(coords)
        class(chidg_worker_t),  intent(in)              :: self
        character(*),           intent(in)              :: string
        character(*),           intent(in), optional    :: user_source


        character(:),   allocatable                 :: user_msg, source
        real(rk),       allocatable, dimension(:)   :: gq_1, gq_2, gq_3, coords


        !
        ! Select source
        !
        if ( present(user_source) ) then
            source = user_source
        else
            source = self%interpolation_source
        end if


        !
        ! Get coordinates
        !
        if ( (source == 'boundary') .or. (source == 'face interior') .or. (source == 'face exterior') ) then
            !gq_1 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c1_
            !gq_2 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c2_
            !gq_3 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c3_
            gq_1 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,1)
            gq_2 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,2)
            gq_3 = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:,3)
        else if ( (source == 'volume') .or. (source == 'element') ) then
            !gq_1 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c1_
            !gq_2 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c2_
            !gq_3 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c3_
            gq_1 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,1)
            gq_2 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,2)
            gq_3 = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:,3)
        else
            user_msg = "chidg_worker%coordinate: Invalid source for returning coordinate. Options are 'boundary' and 'volume'."
            call chidg_signal_one(FATAL,user_msg,source)
        end if




        !
        ! Define coordinate to return.
        !
        select case (string)
            case ('1')
                coords = gq_1
            case ('2')
                coords = gq_2
            case ('3')
                coords = gq_3

            
!            case ('x')
!
!            case ('y')
!
!            case ('r')
!
!            case ('theta')
!
!            case ('z')
!
            case default
                call chidg_signal_one(FATAL,"chidg_worker%coordinate: Invalid string for selecting coordinate.",string)
        end select


    end function coordinate
    !**************************************************************************************




    !>  Return the approximate size of an element bounding box.
    !!
    !!  Returns:
    !!      h(3) = [hx, hy, hz]
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function element_size(self,source) result(h)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: source

        integer(ik) :: ineighbor_domain_l, ineighbor_element_l
        real(rk)    :: h(3)
        logical     :: proc_local, chimera_face


        if (source == 'interior') then

            h = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%h

        else if (source == 'exterior') then

            !
            ! If Chimera face, use interior element size. APPROXIMATION
            !
            chimera_face = (self%face_type() == CHIMERA)
            if (chimera_face) then

                h = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%h

            !
            ! If conforming face, check for processor status of neighbor.
            !
            else

                proc_local = (self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_proc  ==  IRANK)
                if (proc_local) then

                    ineighbor_domain_l  = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_domain_l
                    ineighbor_element_l = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_element_l
                    h = self%mesh%domain(ineighbor_domain_l)%elems(ineighbor_element_l)%h

                else

                    h = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%neighbor_h

                end if

            end if


        else
            call chidg_signal(FATAL,"chidg_worker%element_size(source): Invalid value for 'source'. Options are 'interior', 'exterior'")
        end if


    end function element_size
    !**************************************************************************************








    !>  Return the order of the solution polynomial expansion.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function solution_order(self,source) result(order)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: source

        integer(ik) :: ineighbor_domain_l, ineighbor_element_l, nterms_s, order
        logical     :: proc_local, chimera_face


        if (source == 'interior') then

            nterms_s = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%nterms_s


        else if (source == 'exterior') then

            !
            ! If Chimera face, use interior element order. APPROXIMATION
            !
            chimera_face = (self%face_type() == CHIMERA)
            if (chimera_face) then

                nterms_s = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%nterms_s

            !
            ! If conforming face, check for processor status of neighbor.
            !
            else

                proc_local = (self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_proc  ==  IRANK)
                if (proc_local) then

                    ineighbor_domain_l  = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_domain_l
                    ineighbor_element_l = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_element_l
                    nterms_s = self%mesh%domain(ineighbor_domain_l)%elems(ineighbor_element_l)%nterms_s

                else

                    nterms_s = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_nterms_s

                end if

            end if


        else
            call chidg_signal(FATAL,"chidg_worker%solution_order(source): Invalid value for 'source'. Options are 'interior', 'exterior'")
        end if




        !
        ! Compute polynomial order from number of terms in the expansion. 
        !
        ! ASSUMES TENSOR PRODUCT BASIS
        !
        order = nint( real(nterms_s,rk)**(THIRD) - ONE)

    end function solution_order
    !**************************************************************************************



















    !>  Return the quadrature weights for integration.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function quadrature_weights(self,source) result(weights)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: source

        real(rk),   allocatable,    dimension(:)    :: weights

        if (source == 'face') then
            weights = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%basis_s%weights(self%iface)
        else if (source == 'element') then
            weights = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%basis_s%weights()
        else
            call chidg_signal(FATAL,"chidg_worker%quadrature_weights(source): Invalid value for 'source'. Options are 'face', 'element'")
        end if


    end function quadrature_weights
    !**************************************************************************************








    !>  Return the inverse jacobian mapping for integration.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function inverse_jacobian(self,source) result(jinv)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: source

        real(rk),   allocatable,    dimension(:)    :: jinv



        if (source == 'face') then

            jinv = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%jinv

        else if (source == 'element') then

            jinv = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%jinv

        else
            call chidg_signal(FATAL,"chidg_worker%inverse_jacobian(source): Invalid value for 'source'. Options are 'face', 'element'")
        end if


    end function inverse_jacobian
    !**************************************************************************************










    !>  Return the area of the current face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    function face_area(self) result(area)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk)    :: area

        area = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%total_area

    end function face_area
    !**************************************************************************************









    !>  Return the coordinate system of the current geometric object.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   02/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------
    function coordinate_system(self) result(system)
        class(chidg_worker_t),  intent(in)  :: self

        character(:),   allocatable :: system

        if (self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%coordinate_system == CARTESIAN) then
            system = 'Cartesian'
        else if (self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%coordinate_system == CYLINDRICAL) then
            system = 'Cylindrical'
        end if

    end function coordinate_system
    !**************************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function face_type(self) result(face_type_)
        class(chidg_worker_t),  intent(in)  :: self

        integer(ik) :: idom, ielem, iface
        integer(ik) :: face_type_

        idom  = self%element_info%idomain_l
        ielem = self%element_info%ielement_l
        iface = self%iface

        if ( (iface >= 1) .and. (iface <= NFACES) ) then
            face_type_ = self%mesh%domain(idom)%faces(ielem,iface)%ftype
        else
            face_type_ = NOT_A_FACE
        end if


    end function face_type
    !**************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    function time(self) result(solution_time)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk) :: solution_time

        solution_time = self%t

    end function time
    !**************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidg_worker_t),   intent(inout)   :: self

        if (associated(self%mesh))       nullify(self%mesh)
        if (associated(self%solverdata)) nullify(self%solverdata)
        if (associated(self%cache))      nullify(self%cache)

    end subroutine destructor
    !*************************************************************************************


    !
    ! ALE Procedures
    !

    ! Get ALE quantities
    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_grid_velocity_element(self,field) result(grid_vel_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: field

        real(rk), dimension(:), allocatable :: grid_vel_gq

        if (field == 'u_grid') then
            grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%grid_vel(:,1)

        else if (field == 'v_grid') then

            grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%grid_vel(:,2)
        else if (field == 'w_grid') then

            grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%grid_vel(:,3)
        else
            call chidg_signal(FATAL,"chidg_worker%get_grid_velocity_element(field): Invalid value for 'field'.")
        end if
    end function get_grid_velocity_element


    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_grid_velocity_face(self,field,interp_source) result(grid_vel_comp_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: field
        character(*),           intent(in)  :: interp_source


        real(rk),   dimension(:),   allocatable :: grid_vel_comp_gq
        real(rk),   dimension(:,:), allocatable :: grid_vel_gq
        integer(ik)                             :: idomain_l, ielement_l, iface
        logical                                 :: parallel_neighbor

        ! Presumably, the node velocity
        grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%grid_vel
        !if ((interp_source == 'face interior') .or. (interp_source == 'boundary')) then
        !    grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%grid_vel
        !else if (interp_source == 'face exterior') then

        !    if (self%face_type() == INTERIOR) then
        !        parallel_neighbor = (self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_proc /= IRANK) 
        !        if (parallel_neighbor) then
        !            grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%neighbor_grid_vel
        !        else
        !            idomain_l   = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_domain_l
        !            ielement_l  = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_element_l
        !            iface       = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_face
        !            grid_vel_gq = self%mesh%domain(idomain_l)%faces(ielement_l, iface)%grid_vel
        !            !grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%grid_vel
        !        end if


        !    else if (self%face_type() == CHIMERA) then
        !        !ChiID = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ChiID
        !        !grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%grid_vel

        !        ! For Chimera faces, we actually just want to use the interior face velocity
        !        grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%grid_vel

        !    else if (self%face_type() == BOUNDARY) then
        !        grid_vel_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%grid_vel
        !    end if
        !else
        !        call chidg_signal(FATAL,"chidg_worker%get_grid_velocity_face: Invalid value for 'interp_source'.")
        !end if

        if (field == 'u_grid') then
            grid_vel_comp_gq = grid_vel_gq(:,1)
        else if (field == 'v_grid') then
            grid_vel_comp_gq = grid_vel_gq(:,2)
        else if (field == 'w_grid') then
            grid_vel_comp_gq = grid_vel_gq(:,3)
        else
            call chidg_signal(FATAL,"chidg_worker%get_grid_velocity_face: Invalid value for 'field'.")
        end if

    end function get_grid_velocity_face

    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_jacobian_grid_element(self) result(jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk), dimension(:,:,:), allocatable :: jacobian_grid_gq

        jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%jacobian_grid(:,:,:)

    end function get_jacobian_grid_element
    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_inv_jacobian_grid_element(self) result(jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk), dimension(:,:,:), allocatable :: jacobian_grid_gq

        jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%inv_jacobian_grid(:,:,:)

    end function get_inv_jacobian_grid_element


    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_jacobian_grid_face(self) result(jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk), dimension(:,:,:), allocatable :: jacobian_grid_gq

        jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%jacobian_grid(:,:,:)

    end function get_jacobian_grid_face

    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_inv_jacobian_grid_face(self,interp_source) result(jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: interp_source

        real(rk), dimension(:,:,:), allocatable :: jacobian_grid_gq
        integer(ik) :: ChiID, idomain_l, ielement_l, iface
        logical     :: parallel_neighbor


        if ((interp_source == 'face interior') .or. (interp_source == 'boundary')) then
            jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%inv_jacobian_grid
        else if (interp_source == 'face exterior') then

            if (self%face_type() == INTERIOR) then
                parallel_neighbor = (self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_proc /= IRANK) 
                if (parallel_neighbor) then
                    jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%neighbor_inv_jacobian_grid
                else
                    idomain_l  = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_domain_l
                    ielement_l = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_element_l
                    iface      = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_face
                    jacobian_grid_gq = self%mesh%domain(idomain_l)%faces(ielement_l, iface)%inv_jacobian_grid
                end if

            else if (self%face_type() == CHIMERA) then
                ChiID = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ChiID
                jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%inv_jacobian_grid

            else if (self%face_type() == BOUNDARY) then
                 jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%inv_jacobian_grid
            end if

        else
            call chidg_signal(FATAL,"chidg_worker%get_grid_velocity_face: Invalid value for 'interp_source'.")
        end if


    end function get_inv_jacobian_grid_face


    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_det_jacobian_grid_element(self,interp_type) result(det_jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: interp_type


        real(rk), dimension(:), allocatable :: det_jacobian_grid_gq
        real(rk), dimension(:,:), allocatable :: val 

        !det_jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid(:)
        ! Interpolate modes to nodes
        if (interp_type == 'value') then
            det_jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid
        else if (interp_type == 'grad1') then
            det_jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid_grad1
        else if (interp_type == 'grad2') then
            det_jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid_grad2
        else if (interp_type == 'grad3') then
            det_jacobian_grid_gq = self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid_grad3
        end if


!        det_jacobian_grid_gq = &
!        matmul(val,self%mesh%domain(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%det_jacobian_grid_modes)


        !det_jacobian_grid_gq = ONE/det_jacobian_grid_gq
    end function get_det_jacobian_grid_element

    !>
    !!
    !!  @author Eric M. Wolf
    !!  @date 1/9/2017
    !!
    !----------------------------------------------------------------------------------------------------
    function get_det_jacobian_grid_face(self, interp_type,interp_source) result(det_jacobian_grid_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: interp_type
        character(*),           intent(in)  :: interp_source

        real(rk), dimension(:), allocatable :: det_jacobian_grid_gq
        real(rk), dimension(:), allocatable :: val 
        logical                             :: parallel_neighbor
        integer(ik) :: ChiID, idomain_l, ielement_l, iface


        if ( (trim(interp_type) /= 'value') .and. &
             (trim(interp_type) /= 'grad1') .and. &
             (trim(interp_type) /= 'grad2') .and. &
             (trim(interp_type) /= 'grad3') ) call chidg_signal_one(FATAL,"chidg_worker%get_det_jacobian_grid_face: invalid interp_type.",trim(interp_type))

        if ( (trim(interp_source) /= 'face interior') .and. &
             (trim(interp_source) /= 'face exterior') .and. &
             (trim(interp_source) /= 'boundary') ) call chidg_signal_one(FATAL,"chidg_worker%get_det_jacobian_grid_face: invalid interp_type.",trim(interp_type))


        ! Interpolate modes to nodes
        if ((interp_source == 'face interior') .or. (interp_source == 'boundary') )then
            if (interp_type == 'value') then
                val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid
            else if (interp_type == 'grad1') then
                val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad1
            else if (interp_type == 'grad2') then
                val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad2
            else if (interp_type == 'grad3') then
                val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad3
            end if
        else if (interp_source == 'face exterior') then
            if (self%face_type() == INTERIOR) then
                parallel_neighbor = (self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_proc /= IRANK)
                if (parallel_neighbor) then
                    if (interp_type == 'value') then
                        val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid
                    else if (interp_type == 'grad1') then
                        val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad1
                    else if (interp_type == 'grad2') then
                        val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad2
                    else if (interp_type == 'grad3') then
                        val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad3
                    end if
                else
                    idomain_l  = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_domain_l
                    ielement_l = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_element_l
                    iface      = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%ineighbor_face
                    if (interp_type == 'value') then
                        val = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%det_jacobian_grid
                    else if (interp_type == 'grad1') then
                        val = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%det_jacobian_grid_grad1
                    else if (interp_type == 'grad2') then
                        val = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%det_jacobian_grid_grad2
                    else if (interp_type == 'grad3') then
                        val = self%mesh%domain(idomain_l)%faces(ielement_l,iface)%det_jacobian_grid_grad3
                    end if
                end if



            else if (self%face_type() == CHIMERA) then
                ChiID = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ChiID
                if (interp_type == 'value') then
                    val = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%det_jacobian_grid
                else if (interp_type == 'grad1') then
                    val = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%det_jacobian_grid_grad1
                else if (interp_type == 'grad2') then
                    val = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%det_jacobian_grid_grad2
                else if (interp_type == 'grad3') then
                    val = self%mesh%domain(self%element_info%idomain_l)%chimera%recv(ChiID)%det_jacobian_grid_grad3
                end if


            else if (self%face_type() == BOUNDARY) then
                if (interp_type == 'value') then
                    val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid
                else if (interp_type == 'grad1') then
                    val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad1
                else if (interp_type == 'grad2') then
                    val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad2
                else if (interp_type == 'grad3') then
                    val = self%mesh%domain(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad3
                end if

            end if

        else

            call chidg_signal_one(FATAL,"chidg_worker%get_det_jacobian_grid_face: Invalid value for 'interp_source'.", trim(interp_source))
        end if


        det_jacobian_grid_gq = val

    end function get_det_jacobian_grid_face

!    !>
!    !!
!    !!  @author Eric M. Wolf
!    !!  @date 1/9/2017
!    !!
!    !----------------------------------------------------------------------------------------------------
!    function get_det_jacobian_grid_grad_face(self, interp_source) result(det_jacobian_grid_grad_gq)
!        class(chidg_worker_t),  intent(in)  :: self
!        character(*),           intent(in)  :: interp_source
!
!        real(rk), dimension(:,:), allocatable :: det_jacobian_grid_grad_gq
!        real(rk), dimension(:,:), allocatable :: val 
!
!
!
!        allocate(val(self%nnodesf(),3))
!        ! Interpolate modes to nodes
!        if ((interp_source == 'face interior') .or. (interp_source = 'boundary') )then
!                val(:,1) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad1
!
!                val(:,2) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad2
!
!                val(:,3) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad3
!        else if (interp_source == 'face_exterior') then
!            if (self%face_type() == INTERIOR) then
!                val(:,1) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad1
!                val(:,2) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad2
!                val(:,3) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%neighbor_det_jacobian_grid_grad3
!
!
!            else if (self%face_type() == CHIMERA) then
!                ChiID = self%mesh%domain(self%element_info%idomain_l)% &
!                    faces(self%element_info%ielement_l, self%iface)%ChiID
!                val(:,1) = self%mesh%domain(self%element_info%idomain_l)%&
!                chimera%recv(ChiID)%det_jacobian_grid_grad1
!                val(:,2) = self%mesh%domain(self%element_info%idomain_l)%&
!                chimera%recv(ChiID)%det_jacobian_grid_grad2
!                val(:,3) = self%mesh%domain(self%element_info%idomain_l)%&
!                chimera%recv(ChiID)%det_jacobian_grid_grad3
!
!
!            else if (self%face_type() == BOUNDARY) then
!                val(:,1) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad1
!                val(:,2) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad2
!                val(:,3) = &
!                self%mesh%domain(self%element_info%idomain_l)%&
!                faces(self%element_info%ielement_l,self%iface)%det_jacobian_grid_grad3
!
!            end if
!
!        else
!
!            call chidg_signal(FATAL,"chidg_worker%get_det_jacobian_grid_grad_face: Invalid value for 'interp_source'.")
!        end if
!
!
!        det_jacobian_grid_grad_gq = val
!
!    end function get_det_jacobian_grid_face


    !
    ! ALE reference to physical variable conversions
    !

    !>  Return a primary field evaluated at a quadrature node set. The source here
    !!  is determined by chidg_worker.
    !!
    !!  This routine is specifically for model_t's, because we want them to be evaluated
    !!  on face and element sets the same way. So in a model implementation, we just
    !!  want the model to get some quadrature node set to operate on. The chidg_worker
    !!  handles what node set is currently being returned.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_primary_field_value_ale_general(self,field) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field

        type(AD_D), allocatable :: var_gq(:)


        if (self%interpolation_source == 'element') then
            var_gq = self%get_primary_field_value_ale_element(field) 
        else if (self%interpolation_source == 'face interior') then
            var_gq = self%get_primary_field_value_ale_face(field,'face interior')
        else if (self%interpolation_source == 'face exterior') then
            var_gq = self%get_primary_field_value_ale_face(field,'face exterior')
        else if (self%interpolation_source == 'boundary') then
            var_gq = self%get_primary_field_value_ale_face(field,'boundary')
        end if

    end function get_primary_field_value_ale_general
    !**************************************************************************************

    !>  Return a primary field evaluated at a quadrature node set. The source here
    !!  is determined by chidg_worker.
    !!
    !!  This routine is specifically for model_t's, because we want them to be evaluated
    !!  on face and element sets the same way. So in a model implementation, we just
    !!  want the model to get some quadrature node set to operate on. The chidg_worker
    !!  handles what node set is currently being returned.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/30/2016
    !!
    !--------------------------------------------------------------------------------------
    function get_primary_field_grad_ale_general(self,field,gradient_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: gradient_type

        type(AD_D), allocatable :: var_gq(:,:)

        

        if (self%interpolation_source == 'element') then
            var_gq = self%get_primary_field_grad_ale_element(field,gradient_type) 
        else if (self%interpolation_source == 'face interior') then
            var_gq = self%get_primary_field_grad_ale_face(field,gradient_type,'face interior')
        else if (self%interpolation_source == 'face exterior') then
            var_gq = self%get_primary_field_grad_ale_face(field,gradient_type,'face exterior')
        else if (self%interpolation_source == 'boundary') then
            var_gq = self%get_primary_field_grad_ale_face(field,gradient_type,'boundary')
        end if

    end function get_primary_field_grad_ale_general
    !**************************************************************************************



    !>  Returns physical space quantities by tranforming from the reference configuration.
    !!
    !!
    !!  @author Eric M. Wolf
    !!  @date 7/6/2017
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_primary_field_value_ale_element(self, field) result(val_gq)
        class(chidg_worker_t), intent(in)           :: self
        character(*), intent(in)                    :: field

        type(AD_D), allocatable                     :: val_ref(:), val_gq(:)
        real(rk), allocatable                       :: det_jacobian_grid(:)

        val_ref = self%get_primary_field_element(field,'value')
        det_jacobian_grid = self%get_det_jacobian_grid_element('value')

        val_gq = val_ref/det_jacobian_grid
   end function get_primary_field_value_ale_element

    
    !>  Returns physical space quantities by tranforming from the reference configuration.
    !!
    !!
    !!  @author Eric M. Wolf
    !!  @date 7/6/2017
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_primary_field_grad_ale_element(self, field, gradient_type) result(grad_u_gq)
        class(chidg_worker_t), intent(in)           :: self
        character(*), intent(in)                    :: field
        character(*), intent(in)                    :: gradient_type

        type(AD_D), allocatable                     :: u(:), grad1_u(:), grad2_u(:), grad3_u(:), grad_u_gq(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:), &
                                                        det_jacobian_grid_grad1(:), &
                                                        det_jacobian_grid_grad2(:), &
                                                        det_jacobian_grid_grad3(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)


        character(:), allocatable                   :: user_msg

        det_jacobian_grid = self%get_det_jacobian_grid_element('value')
        det_jacobian_grid_grad1 = self%get_det_jacobian_grid_element('grad1')
        det_jacobian_grid_grad2 = self%get_det_jacobian_grid_element('grad2')
        det_jacobian_grid_grad3 = self%get_det_jacobian_grid_element('grad3')

        jacobian_grid = self%get_inv_jacobian_grid_element()
        u = self%get_primary_field_element(field,'value')

        if (gradient_type == 'gradient + lift') then
            grad1_u = self%get_primary_field_element(field,'grad1 + lift')
            grad2_u = self%get_primary_field_element(field,'grad2 + lift')
            grad3_u = self%get_primary_field_element(field,'grad3 + lift')
        elseif (gradient_type == 'gradient') then
            grad1_u = self%get_primary_field_element(field,'grad1')
            grad2_u = self%get_primary_field_element(field,'grad2')
            grad3_u = self%get_primary_field_element(field,'grad3')
        else
            user_msg = "chidg_worker%get_primary_field_grad_ale_element: Invalid interpolation &
                        type. 'gradient' or 'gradient + lift'"
            call chidg_signal(FATAL,user_msg)
        end if


        grad1_u = grad1_u-u/det_jacobian_grid*det_jacobian_grid_grad1
        grad2_u = grad2_u-u/det_jacobian_grid*det_jacobian_grid_grad2
        grad3_u = grad3_u-u/det_jacobian_grid*det_jacobian_grid_grad3

        allocate(grad_u_gq(size(grad1_u,1),3))
        grad_u_gq(:,1) = (jacobian_grid(:,1,1)*grad1_u + jacobian_grid(:,2,1)*grad2_u + jacobian_grid(:,3,1)*grad3_u)/det_jacobian_grid
        grad_u_gq(:,2) = (jacobian_grid(:,1,2)*grad1_u + jacobian_grid(:,2,2)*grad2_u + jacobian_grid(:,3,2)*grad3_u)/det_jacobian_grid
        grad_u_gq(:,3) = (jacobian_grid(:,1,3)*grad1_u + jacobian_grid(:,2,3)*grad2_u + jacobian_grid(:,3,3)*grad3_u)/det_jacobian_grid


   end function get_primary_field_grad_ale_element

    !>  Returns physical space quantities by tranforming from the reference configuration.
    !!
    !!
    !!  @author Eric M. Wolf
    !!  @date 7/6/2017
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_primary_field_value_ale_face(self, field, interp_source) result(val_gq)
        class(chidg_worker_t), intent(in)           :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_source

        type(AD_D), allocatable                     :: val_ref(:), val_gq(:)
        real(rk), allocatable                       :: det_jacobian_grid(:)

        val_ref = self%get_primary_field_face(field,'value', interp_source)
        if (interp_source == 'boundary') then
            ! In this case, the value supplied by the BC is already the physical value!

            val_gq = val_ref
        else
            ! Otherwise, we need to convert the reference configuration value to the physical value.
            det_jacobian_grid = self%get_det_jacobian_grid_face('value', interp_source)

            val_gq = val_ref/det_jacobian_grid

        end if
   end function get_primary_field_value_ale_face

    
    !>  Returns physical space quantities by tranforming from the reference configuration.
    !!
    !!
    !!  @author Eric M. Wolf
    !!  @date 7/6/2017
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    function get_primary_field_grad_ale_face(self, field, gradient_type, interp_source) result(grad_u_gq)
        class(chidg_worker_t), intent(in)           :: self
        character(*),           intent(in)  :: field
        character(*), intent(in)                    :: gradient_type
        character(*),           intent(in)  :: interp_source

        type(AD_D), allocatable                     :: u(:), grad1_u(:), grad2_u(:), grad3_u(:), grad_u_gq(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:), &
                                                        det_jacobian_grid_grad1(:), &
                                                        det_jacobian_grid_grad2(:), &
                                                        det_jacobian_grid_grad3(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)

        character(:), allocatable                   :: user_msg

        if (gradient_type == 'gradient + lift') then
            grad1_u = self%get_primary_field_face(field,'grad1 + lift', interp_source)
            grad2_u = self%get_primary_field_face(field,'grad2 + lift', interp_source)
            grad3_u = self%get_primary_field_face(field,'grad3 + lift', interp_source)
        elseif (gradient_type == 'gradient') then
            grad1_u = self%get_primary_field_face(field,'grad1', interp_source)
            grad2_u = self%get_primary_field_face(field,'grad2', interp_source)
            grad3_u = self%get_primary_field_face(field,'grad3', interp_source)
        else
            user_msg = "chidg_worker%get_primary_field_grad_ale_face: Invalid interpolation &
                        type. 'gradient' or 'gradient + lift'"
            call chidg_signal(FATAL,user_msg)
        end if



        allocate(grad_u_gq(size(grad1_u,1),3))
        if (interp_source == 'boundary') then
            ! In this case, the value supplied by the BC is already the physical value!
            grad_u_gq(:,1) = grad1_u
            grad_u_gq(:,2) = grad2_u
            grad_u_gq(:,3) = grad3_u

        else
            ! Otherwise, we need to convert the reference configuration value to the physical value.
            det_jacobian_grid = self%get_det_jacobian_grid_face('value', interp_source)
            det_jacobian_grid_grad1 = self%get_det_jacobian_grid_face('grad1', interp_source)
            det_jacobian_grid_grad2 = self%get_det_jacobian_grid_face('grad2', interp_source)
            det_jacobian_grid_grad3 = self%get_det_jacobian_grid_face('grad3', interp_source)

            jacobian_grid = self%get_inv_jacobian_grid_face(interp_source)
            u = self%get_primary_field_face(field,'value', interp_source)

            grad1_u = grad1_u-u/det_jacobian_grid*det_jacobian_grid_grad1
            grad2_u = grad2_u-u/det_jacobian_grid*det_jacobian_grid_grad2
            grad3_u = grad3_u-u/det_jacobian_grid*det_jacobian_grid_grad3

            grad_u_gq(:,1) = (jacobian_grid(:,1,1)*grad1_u + jacobian_grid(:,2,1)*grad2_u + jacobian_grid(:,3,1)*grad3_u)/det_jacobian_grid
            grad_u_gq(:,2) = (jacobian_grid(:,1,2)*grad1_u + jacobian_grid(:,2,2)*grad2_u + jacobian_grid(:,3,2)*grad3_u)/det_jacobian_grid
            grad_u_gq(:,3) = (jacobian_grid(:,1,3)*grad1_u + jacobian_grid(:,2,3)*grad2_u + jacobian_grid(:,3,3)*grad3_u)/det_jacobian_grid

        end if

   end function get_primary_field_grad_ale_face



   !
   ! ALE flux post-processing
   !
   
   function post_process_volume_advective_flux_ale(self, flux_1, flux_2, flux_3, advected_quantity) result(flux_ref)
        class(chidg_worker_t),                       intent(in) :: self
        type(AD_D), dimension(:), intent(inout)                                  :: flux_1, flux_2, flux_3
        type(AD_D), dimension(:), intent(in)                                  :: advected_quantity

        type(AD_D), allocatable                                 :: flux_ref(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:), u_grid(:), v_grid(:), w_grid(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)


        u_grid = self%get_grid_velocity_element('u_grid')
        v_grid = self%get_grid_velocity_element('v_grid')
        w_grid = self%get_grid_velocity_element('w_grid')

        det_jacobian_grid = self%get_det_jacobian_grid_element('value')
        jacobian_grid = self%get_inv_jacobian_grid_element()

        flux_1 = flux_1-u_grid*advected_quantity
        flux_2 = flux_2-v_grid*advected_quantity
        flux_3 = flux_3-w_grid*advected_quantity
       
        allocate(flux_ref(size(flux_1,1),3))
        flux_ref(:,1) = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1 + jacobian_grid(:,1,2)*flux_2 + jacobian_grid(:,1,3)*flux_3)
        flux_ref(:,2) = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1 + jacobian_grid(:,2,2)*flux_2 + jacobian_grid(:,2,3)*flux_3)
        flux_ref(:,3) = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1 + jacobian_grid(:,3,2)*flux_2 + jacobian_grid(:,3,3)*flux_3)


   end function post_process_volume_advective_flux_ale


   function post_process_boundary_advective_flux_ale(self, flux_1, flux_2, flux_3, advected_quantity, interp_source) result(flux_ref)
        class(chidg_worker_t),                       intent(in) :: self
        type(AD_D), dimension(:), intent(inout)                                  :: flux_1, flux_2, flux_3
        type(AD_D), dimension(:), intent(in)                                  :: advected_quantity
        character(*),           intent(in)  :: interp_source

        type(AD_D),allocatable, dimension(:)                                  :: flux_1_tmp, flux_2_tmp, flux_3_tmp
        type(AD_D), allocatable                                 :: flux_ref(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:), u_grid(:), v_grid(:), w_grid(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)

        u_grid = self%get_grid_velocity_face('u_grid', interp_source)
        v_grid = self%get_grid_velocity_face('v_grid', interp_source)
        w_grid = self%get_grid_velocity_face('w_grid', interp_source)
!        ! Always choose a consistent velocity
!        u_grid = self%get_grid_velocity_face('u_grid', 'face interior')
!        v_grid = self%get_grid_velocity_face('v_grid', 'face interior')
!        w_grid = self%get_grid_velocity_face('w_grid', 'face interior')


        det_jacobian_grid = self%get_det_jacobian_grid_face('value', interp_source)
        jacobian_grid = self%get_inv_jacobian_grid_face(interp_source)

        flux_1_tmp = flux_1-u_grid*advected_quantity
        flux_2_tmp = flux_2-v_grid*advected_quantity
        flux_3_tmp = flux_3-w_grid*advected_quantity
       
        allocate(flux_ref(size(flux_1,1),3))
        flux_ref(:,1) = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1_tmp + jacobian_grid(:,1,2)*flux_2_tmp + jacobian_grid(:,1,3)*flux_3_tmp)
        flux_ref(:,2) = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1_tmp + jacobian_grid(:,2,2)*flux_2_tmp + jacobian_grid(:,2,3)*flux_3_tmp)
        flux_ref(:,3) = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1_tmp + jacobian_grid(:,3,2)*flux_2_tmp + jacobian_grid(:,3,3)*flux_3_tmp)


   end function post_process_boundary_advective_flux_ale


   function post_process_volume_diffusive_flux_ale(self, flux_1, flux_2, flux_3) result(flux_ref)
        class(chidg_worker_t),                       intent(in) :: self
        type(AD_D), dimension(:), intent(in)                                  :: flux_1, flux_2, flux_3

        type(AD_D), allocatable                                 :: flux_ref(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)


        det_jacobian_grid = self%get_det_jacobian_grid_element('value')
        jacobian_grid = self%get_inv_jacobian_grid_element()

       
        allocate(flux_ref(size(flux_1,1),3))
        flux_ref(:,1) = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1 + jacobian_grid(:,1,2)*flux_2 + jacobian_grid(:,1,3)*flux_3)
        flux_ref(:,2) = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1 + jacobian_grid(:,2,2)*flux_2 + jacobian_grid(:,2,3)*flux_3)
        flux_ref(:,3) = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1 + jacobian_grid(:,3,2)*flux_2 + jacobian_grid(:,3,3)*flux_3)


   end function post_process_volume_diffusive_flux_ale

   function post_process_boundary_diffusive_flux_ale(self, flux_1, flux_2, flux_3, interp_source) result(flux_ref)
        class(chidg_worker_t),                       intent(in) :: self
        type(AD_D), dimension(:), intent(in)                                  :: flux_1, flux_2, flux_3
        character(*),           intent(in)  :: interp_source

        type(AD_D), allocatable                                 :: flux_ref(:,:)
        real(rk), allocatable                       :: det_jacobian_grid(:)
        real(rk), allocatable                       :: jacobian_grid(:,:,:)


        det_jacobian_grid = self%get_det_jacobian_grid_face('value', interp_source)
        jacobian_grid = self%get_inv_jacobian_grid_face(interp_source)

       
        allocate(flux_ref(size(flux_1,1),3))
        flux_ref(:,1) = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1 + jacobian_grid(:,1,2)*flux_2 + jacobian_grid(:,1,3)*flux_3)
        flux_ref(:,2) = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1 + jacobian_grid(:,2,2)*flux_2 + jacobian_grid(:,2,3)*flux_3)
        flux_ref(:,3) = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1 + jacobian_grid(:,3,2)*flux_2 + jacobian_grid(:,3,3)*flux_3)


   end function post_process_boundary_diffusive_flux_ale


end module type_chidg_worker
