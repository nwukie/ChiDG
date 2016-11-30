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
    use mod_constants,      only: NFACES, ME, NEIGHBOR, BC, ZERO
    use mod_interpolate,    only: interpolate_element_standard, &
                                  interpolate_face_standard

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
        type(solverdata_t),     pointer :: solverdata
        type(chidg_cache_t),    pointer :: cache

        integer(ik)                 :: iface
        type(element_info_t)        :: element_info
        type(function_info_t)       :: function_info
    

    contains 
    
        ! Worker state
        procedure   :: init
        procedure   :: set_element          ! Set element_info type
        procedure   :: set_function_info    ! Set function_info type
        procedure   :: set_face             ! Set iface index
        procedure   :: face_info            ! Return a face_info type


        ! Worker get data
        procedure   :: get_face_variable
        procedure   :: get_element_variable
        procedure   :: get_element_auxiliary_field

        procedure   :: store_bc_state

        procedure   :: normal
        procedure   :: unit_normal

        procedure   :: coords
        procedure   :: x
        procedure   :: y
        procedure   :: z

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
    subroutine init(self,mesh,solverdata,cache)
        class(chidg_worker_t),  intent(inout)       :: self
        type(mesh_t),           intent(in), target  :: mesh(:)
        type(solverdata_t),     intent(in), target  :: solverdata
        type(chidg_cache_t),    intent(in), target  :: cache


        self%mesh       => mesh
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










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function get_face_variable(self,ieqn,interp_type,interp_source) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: ieqn
        character(len=*),       intent(in)  :: interp_type
        integer(ik),            intent(in)  :: interp_source

        type(AD_D), allocatable, dimension(:) :: &
            var_gq

        type(face_info_t)               :: face_info
        character(len=:), allocatable   :: cache_component, cache_type
        integer(ik)                     :: idirection, igq
        logical                         :: keep_linearization



        !
        ! Set cache_component
        !
        if (interp_source == ME) then
            cache_component = 'face interior'
        else if (interp_source == NEIGHBOR .or. &
                 interp_source == BC) then
            cache_component = 'face exterior'
        else
            call chidg_signal(FATAL,"chidg_worker%get_face_variable: Invalid value for interp_source. Try ME, NEIGHBOR, or BC.")
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
            var_gq = self%cache%get_data(cache_component,cache_type,idirection,self%function_info%seed,ieqn,self%iface)

        else if (cache_type == 'derivative') then

            ! Get DG derivative on face
            var_gq = self%cache%get_data(cache_component,'derivative',idirection,self%function_info%seed,ieqn,self%iface)


        else if (cache_type == 'derivative + lift') then

            ! Get DG derivative on face
            var_gq = self%cache%get_data(cache_component,'derivative',idirection,self%function_info%seed,ieqn,self%iface)

            ! Modify derivative by face lift stabilized by a factor of NFACES
            var_gq = var_gq + real(NFACES,rk)*self%cache%get_data(cache_component,'lift face',idirection,self%function_info%seed,ieqn,self%iface)

        end if


    end function get_face_variable
    !***************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_element_variable(self,ieqn,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: ieqn
        character(len=*),       intent(in)  :: interp_type

        type(AD_D), allocatable, dimension(:) :: &
            var_gq, tmp_gq

        type(face_info_t)               :: face_info
        character(len=:), allocatable   :: cache_component, cache_type
        integer(ik)                     :: idirection, igq, iface
        logical                         :: keep_linearization




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
        if ( cache_type == 'value') then
            var_gq = self%cache%get_data('element',cache_type,idirection,self%function_info%seed,ieqn)



        else if (cache_type == 'derivative') then
            ! Get DG derivative
            var_gq = self%cache%get_data('element','derivative',idirection,self%function_info%seed,ieqn)



        else if (cache_type == 'derivative + lift') then
            ! Get DG derivative
            var_gq = self%cache%get_data('element','derivative',idirection,self%function_info%seed,ieqn)

            ! Add lift contributions from each face
            do iface = 1,NFACES
                tmp_gq = self%cache%get_data('face interior', 'lift element', idirection, self%function_info%seed,ieqn,iface)
                var_gq = var_gq + tmp_gq
            end do

        end if






    end function get_element_variable
    !****************************************************************************************













    !>  Return an interpolation of an auxiliary chidgVector on the element quadrature node set.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/4/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_element_auxiliary_field(self,field,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(*),           intent(in)  :: field
        character(*),           intent(in)  :: interp_type

        character(:),   allocatable :: user_msg
        real(rk),       allocatable :: var_gq(:)

        integer(ik) :: ifield


        !
        ! Get index of the auxiliary field vector.
        !
        ifield = self%solverdata%get_auxiliary_field_index(field)


        user_msg = "chidg_worker%get_element_auxiliary_field: There was no field data found for the &
                    specified field string."
        if (ifield == 0) call chidg_signal_one(FATAL,user_msg,trim(field))


        if ( (trim(interp_type) == 'value') .or. &
             (trim(interp_type) == 'ddx'  ) .or. &
             (trim(interp_type) == 'ddy'  ) .or. &
             (trim(interp_type) == 'ddz'  ) ) then

            ! Here, we implicitly assume that all auxiliary field vectors contain only
            ! one variable expansion. Hence, ieqn = 1.
            var_gq = interpolate_element_standard(self%mesh,self%solverdata%auxiliary_field(ifield),    &
                                                            idomain_l  = self%element_info%idomain_l,   &
                                                            ielement_l = self%element_info%ielement_l,  &
                                                            ieqn = 1,                                   &
                                                            interpolation_type = interp_type)

        elseif ( (trim(interp_type) == 'ddx+lift')   .or. &
                 (trim(interp_type) == 'ddy+lift')   .or. &
                 (trim(interp_type) == 'ddz+lift')   .or. &
                 (trim(interp_type) == 'ddx + lift') .or. &
                 (trim(interp_type) == 'ddy + lift') .or. &
                 (trim(interp_type) == 'ddz + lift') ) then

            user_msg = "chidg_worker%get_element_auxiliary_field: Lifted derivatives('ddx + lift') are&
                        not supported. Only the standard derivatives('ddx','ddy','ddz') are supported."
            call chidg_signal_one(FATAL,user_msg,trim(interp_type))
        end if






    end function get_element_auxiliary_field
    !****************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine store_bc_state(self,ieqn,cache_data,data_type)
        class(chidg_worker_t),  intent(inout)   :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(in)      :: cache_data(:)
        character(len=*),       intent(in)      :: data_type

        character(len=:), allocatable   :: cache_type
        integer(ik)                     :: idirection


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
            call chidg_signal(FATAL,"worker%store_bc_state: Invalid data_type specification. Options are 'value', 'ddx', 'ddy', 'ddz'.")
        end if







        !
        ! Store bc state in cache, face exterior component
        !
        if (cache_type == 'value') then
            call self%cache%set_data('face exterior',cache_data,'value',0,self%function_info%seed,ieqn,self%iface)

        else if (cache_type == 'derivative') then
            call self%cache%set_data('face exterior',cache_data,'derivative',idirection,self%function_info%seed,ieqn,self%iface)

        end if



    end subroutine store_bc_state
    !***************************************************************************************















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine integrate_boundary(self,ieqn,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand(:)


        call integrate_boundary_scalar_flux(self%mesh,self%solverdata,self%face_info(),self%function_info,ieqn,integrand)


    end subroutine integrate_boundary
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine integrate_volume_flux(self,ieqn,integrand_x,integrand_y,integrand_z)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand_x(:)
        type(AD_D),             intent(inout)   :: integrand_y(:)
        type(AD_D),             intent(inout)   :: integrand_z(:)


        call integrate_volume_vector_flux(self%mesh,self%solverdata,self%element_info,self%function_info,ieqn,integrand_x,integrand_y,integrand_z)


    end subroutine integrate_volume_flux
    !***************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine integrate_volume_source(self,ieqn,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand(:)


        call integrate_volume_scalar_source(self%mesh,self%solverdata,self%element_info,self%function_info,ieqn,integrand)


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
