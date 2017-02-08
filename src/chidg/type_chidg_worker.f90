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
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: NFACES, ME, NEIGHBOR, BC, ZERO, CHIMERA, ONE, THIRD
    use mod_interpolate,    only: interpolate_element_standard, &
                                  interpolate_element_autodiff, &
                                  interpolate_face_standard,    &
                                  interpolate_face_autodiff

    use mod_integrate,      only: integrate_boundary_scalar_flux, &
                                  integrate_volume_vector_flux,   &
                                  integrate_volume_scalar_source

    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use type_chidg_cache,   only: chidg_cache_t
    use type_properties,    only: properties_t
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
    
        type(mesh_t),           pointer :: mesh(:)
        type(properties_t), allocatable :: prop(:)
        !type(properties_t),     pointer :: prop(:)
        type(solverdata_t),     pointer :: solverdata
        type(chidg_cache_t),    pointer :: cache

        integer(ik)                 :: iface
        integer(ik)                 :: itime
        type(element_info_t)        :: element_info
        type(function_info_t)       :: function_info
    
        character(:),   allocatable :: interpolation_source

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

        procedure   :: store_bc_state
        procedure   :: store_model_field

        procedure   :: normal
        procedure   :: unit_normal

        procedure   :: coords
        procedure   :: x
        procedure   :: y
        procedure   :: z

        procedure   :: element_size
        procedure   :: solution_order
        procedure   :: quadrature_weights
        procedure   :: inverse_jacobian
        procedure   :: face_area

        procedure   :: face_type

        procedure   :: time


        ! Worker process data
        procedure   :: integrate_boundary
        generic     :: integrate_volume => integrate_volume_flux, &
                                           integrate_volume_source
        procedure   :: integrate_volume_flux
        procedure   :: integrate_volume_source


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
        type(mesh_t),           intent(in), target  :: mesh(:)
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
            call chidg_signal_one(FATAL,user_msg,trim(interp_source))
        end if


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (interp_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (interp_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3


        else if ( (interp_type == 'ddx + lift') .or. &
                  (interp_type == 'ddx+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 1
        else if ( (interp_type == 'ddy + lift') .or. &
                  (interp_type == 'ddy+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 2
        else if ( (interp_type == 'ddz + lift') .or. &
                  (interp_type == 'ddz+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 3
        end if



        !
        ! Retrieve data from cache
        !
        if (cache_type == 'value') then
            var_gq = self%cache%get_data(field,cache_component,'value',idirection,self%function_info%seed,self%iface)

        else if (cache_type == 'derivative') then
            var_gq = self%cache%get_data(field,cache_component,'derivative',idirection,self%function_info%seed,self%iface)

        else if (cache_type == 'derivative + lift') then
            var_gq = self%cache%get_data(field,cache_component,'derivative',idirection,self%function_info%seed,self%iface)

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
        integer(ik)                     :: idirection, igq, iface, ifield, idomain_l


        !
        ! Set cache_type
        !
        if (interp_type == 'value') then
            cache_type = 'value'
            idirection = 0
        else if (interp_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (interp_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (interp_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3


        else if ( (interp_type == 'ddx + lift') .or. &
                  (interp_type == 'ddx+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 1
        else if ( (interp_type == 'ddy + lift') .or. &
                  (interp_type == 'ddy+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 2
        else if ( (interp_type == 'ddz + lift') .or. &
                  (interp_type == 'ddz+lift'  ) ) then
            cache_type = 'derivative + lift'
            idirection = 3
        end if




        !
        ! If we are requesting an interpolation from a subset of the modal expansion, 
        ! then perform a new interpolation rather than using the cache.
        !
        if (present(Pmin) .or. present(Pmax)) then

            if (cache_type == 'value') then
                idomain_l = self%element_info%idomain_l
                ifield    = self%prop(idomain_l)%get_primary_field_index(field)

                var_gq = interpolate_element_autodiff(self%mesh, self%solverdata%q, self%element_info, self%function_info, ifield, self%itime, interp_type, Pmin, Pmax)

            else if ( (cache_type == 'derivative') .or. &
                      (cache_type == 'derivative + lift') ) then
                user_msg = "chidg_worker%get_primary_field_element: On partial field interpolations, &
                            only the 'value' of the field can be interpolated. 'derivative' is not yet implemented."
                call chidg_signal(FATAL,user_msg)
            end if

    

        else


            !
            ! Retrieve data from cache
            !
            if ( cache_type == 'value') then
                var_gq = self%cache%get_data(field,'element','value',idirection,self%function_info%seed)

            else if (cache_type == 'derivative') then
                var_gq = self%cache%get_data(field,'element','derivative',idirection,self%function_info%seed)

            else if (cache_type == 'derivative + lift') then
                var_gq = self%cache%get_data(field,'element','derivative',idirection,self%function_info%seed)

                ! Add lift contributions from each face
                do iface = 1,NFACES
                    tmp_gq = self%cache%get_data(field,'face interior', 'lift element', idirection, self%function_info%seed,iface)
                    var_gq = var_gq + tmp_gq
                end do

            end if


        end if


    end function get_primary_field_element
    !****************************************************************************************






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
        else if ( (interp_type == 'ddx')          .or. &
                  (interp_type == 'ddy')          .or. &
                  (interp_type == 'ddz')          .or. &
                  (interp_type == 'ddx + lift')   .or. &
                  (interp_type == 'ddx+lift'  )   .or. &
                  (interp_type == 'ddy + lift')   .or. &
                  (interp_type == 'ddy+lift'  )   .or. &
                  (interp_type == 'ddz + lift')   .or. &
                  (interp_type == 'ddz+lift'  ) ) then
                user_msg = 'chidg_worker%get_model_field_face: Computing derivatives for model &
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
        else if ( (interp_type == 'ddx')        .or. &
                  (interp_type == 'ddy')        .or. &
                  (interp_type == 'ddz')        .or. &
                  (interp_type == 'ddx+lift'  ) .or. &
                  (interp_type == 'ddy+lift'  ) .or. &
                  (interp_type == 'ddz+lift'  ) .or. &
                  (interp_type == 'ddx + lift') .or. &
                  (interp_type == 'ddy + lift') .or. &
                  (interp_type == 'ddz + lift') ) then
            user_msg = 'chidg_worker%get_model_field_element: Computing derivatives for model &
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
        else if (interp_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (interp_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (interp_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3

        else if ( (interp_type == 'ddx+lift'  ) .or. &
                  (interp_type == 'ddy+lift'  ) .or. &
                  (interp_type == 'ddz+lift'  ) .or. &
                  (interp_type == 'ddx + lift') .or. &
                  (interp_type == 'ddy + lift') .or. &
                  (interp_type == 'ddz + lift') ) then
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
        else if (interp_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (interp_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (interp_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3


        else if ( (interp_type == 'ddx+lift'  ) .or. &
                  (interp_type == 'ddy+lift'  ) .or. &
                  (interp_type == 'ddz+lift'  ) .or. &
                  (interp_type == 'ddx + lift') .or. &
                  (interp_type == 'ddy + lift') .or. &
                  (interp_type == 'ddz + lift') ) then

            user_msg = 'chidg_worker%get_auxiliary_field_element: Computing lifted derivatives for auxiliary &
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
        else if (data_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (data_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (data_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3
        else
            user_msg = "chidg_worker%store_bc_state: Invalid data_type specification. &
                        Options are 'value', 'ddx', 'ddy', 'ddz'."
            call chidg_signal_one(FATAL,user_msg,trim(data_type))
        end if



        !
        ! Store bc state in cache, face exterior component
        !
        if (cache_type == 'value') then
            call self%cache%set_data(field,'face exterior',cache_data,'value',0,self%function_info%seed,self%iface)

        else if (cache_type == 'derivative') then
            call self%cache%set_data(field,'face exterior',cache_data,'derivative',idirection,self%function_info%seed,self%iface)

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
        else if (data_type == 'ddx') then
            cache_type = 'derivative'
            idirection = 1
        else if (data_type == 'ddy') then
            cache_type = 'derivative'
            idirection = 2
        else if (data_type == 'ddz') then
            cache_type = 'derivative'
            idirection = 3
        else
            user_msg = "chidg_worker%store_model_field: Invalid data_type specification. &
                        Options are 'value', 'ddx', 'ddy', 'ddz'."
            call chidg_signal_one(FATAL,user_msg,trim(data_type))
        end if



        !
        ! Add data to model field cache
        !
        if (cache_type == 'value') then

!            field_current = self%cache%get_data(model_field,self%interpolation_source,'value',0,self%function_info%seed,self%iface)
!            field_update = field_current + cache_data
            field_update = cache_data

            call self%cache%set_data(model_field,self%interpolation_source,field_update,'value',0,self%function_info%seed,self%iface)

        else if (cache_type == 'derivative') then

!            field_current = self%cache%get_data(model_field,self%interpolation_source,'derivative',idirection,self%function_info%seed,self%iface)
!            field_update = field_current + cache_data
            field_update = cache_data

            call self%cache%set_data(model_field,self%interpolation_source,field_update,'derivative',idirection,self%function_info%seed,self%iface)

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

        integer(ik) :: ifield, idomain_l

        idomain_l = self%element_info%idomain_l
        ifield    = self%prop(idomain_l)%get_primary_field_index(primary_field)

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

        integer(ik) :: ifield, idomain_l


        idomain_l = self%element_info%idomain_l
        ifield    = self%prop(idomain_l)%get_primary_field_index(primary_field)

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

        integer(ik) :: ifield, idomain_l


        idomain_l = self%element_info%idomain_l
        ifield    = self%prop(idomain_l)%get_primary_field_index(primary_field)

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


        norm_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%norm(:,direction)


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


        unorm_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%unorm(:,direction)


    end function unit_normal
    !***************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function coords(self) result(coords_gq)
        class(chidg_worker_t),  intent(in)  :: self

        type(point_t), allocatable, dimension(:) :: coords_gq

        coords_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)

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
            x_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c1_
        else if (source == 'volume') then
            x_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c1_
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
            y_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c2_
        else if (source == 'volume') then
            y_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c2_
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
            z_gq = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l,self%iface)%quad_pts(:)%c3_
        else if (source == 'volume') then
            z_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c3_
        else
            call chidg_signal(FATAL,"chidg_worker%z(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if


    end function z
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

            h = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%h

        else if (source == 'exterior') then

            !
            ! If Chimera face, use interior element size. APPROXIMATION
            !
            chimera_face = (self%face_type() == CHIMERA)
            if (chimera_face) then

                h = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%h

            !
            ! If conforming face, check for processor status of neighbor.
            !
            else

                proc_local = (self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_proc  ==  IRANK)
                if (proc_local) then

                    ineighbor_domain_l  = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_domain_l
                    ineighbor_element_l = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_element_l
                    h = self%mesh(ineighbor_domain_l)%elems(ineighbor_element_l)%h

                else

                    h = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%neighbor_h

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

            nterms_s = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%nterms_s


        else if (source == 'exterior') then

            !
            ! If Chimera face, use interior element order. APPROXIMATION
            !
            chimera_face = (self%face_type() == CHIMERA)
            if (chimera_face) then

                nterms_s = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%nterms_s

            !
            ! If conforming face, check for processor status of neighbor.
            !
            else

                proc_local = (self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_proc  ==  IRANK)
                if (proc_local) then

                    ineighbor_domain_l  = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_domain_l
                    ineighbor_element_l = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_element_l
                    nterms_s = self%mesh(ineighbor_domain_l)%elems(ineighbor_element_l)%nterms_s

                else

                    nterms_s = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%ineighbor_nterms_s

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

            weights = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%gq%face%weights(:,self%iface)

        else if (source == 'element') then

            weights = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%gq%vol%weights

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

            jinv = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%jinv

        else if (source == 'element') then

            jinv = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%jinv

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

        area = self%mesh(self%element_info%idomain_l)%faces(self%element_info%ielement_l, self%iface)%area

    end function face_area
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


        face_type_ = self%mesh(idom)%faces(ielem,iface)%ftype


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

        solution_time = self%solverdata%t

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

end module type_chidg_worker
