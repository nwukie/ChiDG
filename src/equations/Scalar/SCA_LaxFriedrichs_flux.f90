module SCA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF,ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use DNAD_D

    use SCA_properties,          only: SCA_properties_t
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: SCA_LaxFriedrichs_flux_t


    contains
        procedure   :: compute

    end type SCA_LaxFriedrichs_flux_t
    !*****************************************************************************************************

contains




    !>
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(SCA_LaxFriedrichs_flux_t),    intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        type(face_info_t),                  intent(in)      :: face_info
        type(function_info_t),              intent(in)      :: function_info


        integer(ik)                 :: idom, ielem, iface
        integer(ik)                 :: ifcn, idomor, iblk

        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, ierr, nnodes

        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)  :: &
                            u_l, u_r, flux_x, flux_y, flux_z, integrand


        !
        ! Get variable index
        !
        iu = prop%get_eqn_index("u")

        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface

        



        associate ( norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm )

        !
        ! Get equation set property data
        !
        select type(prop)
            type is (SCA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select


        !
        ! Interpolate solution to quadrature nodes
        !
        u_r = interpolate(mesh,sdata,face_info,function_info,iu, 'value', ME)
        u_l = interpolate(mesh,sdata,face_info,function_info,iu, 'value', NEIGHBOR)

        !
        ! Compute boundary upwind flux
        !
        flux_x = (cx * (u_l - u_r)/TWO )  *  norms(:,1) * unorms(:,1)
        flux_y = (cy * (u_l - u_r)/TWO )  *  norms(:,2) * unorms(:,2)
        flux_z = (cz * (u_l - u_r)/TWO )  *  norms(:,3) * unorms(:,3)

        integrand = flux_x + flux_y + flux_z

        !
        ! Integrate flux
        !
        call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu,integrand)

        end associate
    end subroutine compute
    !*********************************************************************************************************









end module SCA_LaxFriedrichs_flux
