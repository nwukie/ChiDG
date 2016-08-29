module type_chidg_worker
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: NFACES
    use mod_interpolate,    only: interpolate_element_autodiff, &
                                  interpolate_face_autodiff
    use mod_integrate,      only: integrate_boundary_scalar_flux, &
                                  integrate_volume_vector_flux,   &
                                  integrate_volume_scalar_source

    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use type_solverdata,    only: solverdata_t
    use type_BR2,           only: BR2_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use DNAD_D
    implicit none





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: chidg_worker_t
    
        type(mesh_t),           pointer :: mesh(:)
        type(solverdata_t),     pointer :: solverdata
        type(BR2_t),            pointer :: BR2
    

        type(element_info_t)        :: element_info
        type(face_info_t)           :: face_info
        type(function_info_t)       :: function_info
    
    contains 
    
        ! Worker state
        procedure   :: init
        procedure   :: set_element_info
        procedure   :: set_face_info
        procedure   :: set_function_info


        ! Worker get data
        procedure   :: interpolate_face
        procedure   :: interpolate_element
        generic     :: interpolate => interpolate_face, &
                                      interpolate_element

        procedure   :: normal
        procedure   :: unit_normal

        procedure   :: coords
        procedure   :: x
        procedure   :: y
        procedure   :: z

        procedure   :: time


        ! Worker process data
        procedure   :: integrate_boundary

        procedure   :: integrate_volume_flux
        procedure   :: integrate_volume_source
        generic     :: integrate_volume => integrate_volume_flux, &
                                           integrate_volume_source


        final       :: destructor
    
    end type chidg_worker_t
    !*****************************************************************************    






contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !-------------------------------------------------------------------------------
    subroutine init(self,mesh,solverdata,BR2)
        class(chidg_worker_t),  intent(inout)       :: self
        type(mesh_t),           intent(in), target  :: mesh(:)
        type(solverdata_t),     intent(in), target  :: solverdata
        type(BR2_t),            intent(in), target, optional    :: BR2

        self%mesh       => mesh
        self%solverdata => solverdata

        if (present(BR2)) then
            self%BR2 => BR2
        end if

    end subroutine init
    !********************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine set_element_info(self,elem_info)
        class(chidg_worker_t),  intent(inout)   :: self
        type(element_info_t),   intent(in)      :: elem_info

        self%element_info = elem_info

    end subroutine set_element_info
    !**********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine set_face_info(self,face_info)
        class(chidg_worker_t),  intent(inout)   :: self
        type(face_info_t),      intent(in)      :: face_info

        self%face_info = face_info

    end subroutine set_face_info
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
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function interpolate_face(self,ieqn,interp_type,interp_source) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: ieqn
        character(len=*),       intent(in)  :: interp_type
        integer(ik),            intent(in)  :: interp_source

        !type(AD_D), dimension(self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%gq%face%nnodes) :: &
        type(AD_D), allocatable, dimension(:) :: &
            var_gq, deriv, lift


        if (interp_type == 'value') then
            var_gq = interpolate_face_autodiff(self%mesh,self%solverdata%q,self%face_info,self%function_info,ieqn,interp_type,interp_source)

        elseif ((interp_type == 'ddx') .or. &
                (interp_type == 'ddy') .or. &
                (interp_type == 'ddz') ) then

            deriv = interpolate_face_autodiff(self%mesh,self%solverdata%q,self%face_info,self%function_info,ieqn,interp_type,interp_source)

!            lift = self%solverdata%BR2%interpolate_lift_face(self%mesh,self%face_info,self%function_info,ieqn,interp_type,interp_source)
            lift = self%BR2%interpolate_lift_face(self%mesh,self%face_info,self%function_info,ieqn,interp_type,interp_source)

            var_gq = deriv + real(NFACES,rk)*lift

        end if

    end function interpolate_face
    !**********************************************************************************************



    !>
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function interpolate_element(self,ieqn,interp_type) result(var_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: ieqn
        character(len=*),       intent(in)  :: interp_type

        !type(AD_D), dimension(self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%gq%vol%nnodes) :: &
        type(AD_D), allocatable, dimension(:) :: &
            var_gq, deriv, lift

        if (interp_type == 'value') then
            var_gq = interpolate_element_autodiff(self%mesh,self%solverdata%q,self%element_info,self%function_info,ieqn,interp_type)

        elseif ((interp_type == 'ddx') .or. &
                (interp_type == 'ddy') .or. &
                (interp_type == 'ddz') ) then

            deriv = interpolate_element_autodiff(self%mesh,self%solverdata%q,self%element_info,self%function_info,ieqn,interp_type)

!            lift  = self%solverdata%BR2%interpolate_lift_element(self%mesh,self%element_info,self%function_info,ieqn,interp_type)
            lift  = self%BR2%interpolate_lift_element(self%mesh,self%element_info,self%function_info,ieqn,interp_type)

            var_gq = deriv + lift

        end if

    end function interpolate_element
    !**********************************************************************************************












    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine integrate_boundary(self,ieqn,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand(:)


        call integrate_boundary_scalar_flux(self%mesh,self%solverdata,self%face_info,self%function_info,ieqn,integrand)


    end subroutine integrate_boundary
    !*********************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine integrate_volume_flux(self,ieqn,integrand_x,integrand_y,integrand_z)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand_x(:)
        type(AD_D),             intent(inout)   :: integrand_y(:)
        type(AD_D),             intent(inout)   :: integrand_z(:)


        call integrate_volume_vector_flux(self%mesh,self%solverdata,self%element_info,self%function_info,ieqn,integrand_x,integrand_y,integrand_z)


    end subroutine integrate_volume_flux
    !*********************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine integrate_volume_source(self,ieqn,integrand)
        class(chidg_worker_t),  intent(in)      :: self
        integer(ik),            intent(in)      :: ieqn
        type(AD_D),             intent(inout)   :: integrand(:)


        call integrate_volume_scalar_source(self%mesh,self%solverdata,self%element_info,self%function_info,ieqn,integrand)


    end subroutine integrate_volume_source
    !*********************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function normal(self,direction) result(norm_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: direction

        real(rk), dimension(:), allocatable :: norm_gq


        norm_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%norm(:,direction)


    end function normal
    !**********************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function unit_normal(self,direction) result(unorm_gq)
        class(chidg_worker_t),  intent(in)  :: self
        integer(ik),            intent(in)  :: direction

        real(rk), dimension(:), allocatable :: unorm_gq


        unorm_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%unorm(:,direction)


    end function unit_normal
    !**********************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function coords(self) result(coords_gq)
        class(chidg_worker_t),  intent(in)  :: self

        type(point_t), allocatable, dimension(:) :: coords_gq

        coords_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%quad_pts(:)

    end function coords
    !**********************************************************************************************













    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function x(self,source) result(x_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: x_gq

        if (source == 'boundary') then
            x_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%quad_pts(:)%c1_
        else if (source == 'volume') then
            x_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c1_
        else
            call chidg_signal(FATAL,"chidg_worker%x(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if



    end function x
    !**********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function y(self,source) result(y_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: y_gq


        if (source == 'boundary') then
            y_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%quad_pts(:)%c2_
        else if (source == 'volume') then
            y_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c2_
        else
            call chidg_signal(FATAL,"chidg_worker%y(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if



    end function y
    !**********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function z(self,source) result(z_gq)
        class(chidg_worker_t),  intent(in)  :: self
        character(len=*),       intent(in)  :: source

        real(rk), dimension(:), allocatable :: z_gq


        if (source == 'boundary') then
            z_gq = self%mesh(self%face_info%idomain_l)%faces(self%face_info%ielement_l,self%face_info%iface)%quad_pts(:)%c3_
        else if (source == 'volume') then
            z_gq = self%mesh(self%element_info%idomain_l)%elems(self%element_info%ielement_l)%quad_pts(:)%c3_
        else
            call chidg_signal(FATAL,"chidg_worker%z(source): Invalid value for 'source'. Options are 'boundary', 'volume'")
        end if


    end function z
    !**********************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/22/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function time(self) result(solution_time)
        class(chidg_worker_t),  intent(in)  :: self

        real(rk) :: solution_time

        solution_time = self%solverdata%t

    end function time
    !**********************************************************************************************








    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine destructor(self)
        type(chidg_worker_t),   intent(inout)   :: self

        if (associated(self%mesh))       nullify(self%mesh)
        if (associated(self%solverdata)) nullify(self%solverdata)

    end subroutine destructor
    !**********************************************************************************

end module type_chidg_worker
