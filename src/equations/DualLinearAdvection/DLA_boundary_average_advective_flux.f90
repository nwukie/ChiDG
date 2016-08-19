module DLA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                      ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t


    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use DNAD_D

    use type_properties,        only: properties_t
    use DLA_properties,         only: DLA_properties_t
    implicit none

    private



    !> This equation set exists really just to test equationsets with more than one equation. 
    !! The idea is just to compute the linear advection solution twice at the same time. 
    !! The equations are independent of each other. So, we can verify, for example, the volume 
    !! flux jacobians for each equation. They should be the same as for the single 
    !! LinearAdvection equation set.
    !!
    !!  @author Nathan A. Wukie
    !!
    !---------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: DLA_boundary_average_advective_flux_t

    contains
        procedure   :: compute

    end type DLA_boundary_average_advective_flux_t
    !*********************************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(DLA_boundary_average_advective_flux_t),   intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        type(face_info_t),                              intent(in)      :: face_info
        type(function_info_t),                          intent(in)      :: function_info


        integer(ik)             :: idom, ielem, iface
        integer(ik)             :: ifcn, idonor, iblk

        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu_a, iu_b
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes) :: &
                        ua_l, ua_r, ub_l, ub_r, flux_x, flux_y, flux_z, integrand


        !
        ! Get integer data
        !
        iu_a = prop%get_eqn_index("u_a")
        iu_b = prop%get_eqn_index("u_b")


        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface


        associate ( norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm )

        !
        ! Get equation set properties
        !
        select type(prop)
            type is (DLA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select



        !
        ! Interpolate solution to quadrature nodes
        !
        ua_r = interpolate(mesh,sdata,face_info,function_info,iu_a, 'value', ME)
        ua_l = interpolate(mesh,sdata,face_info,function_info,iu_a, 'value', NEIGHBOR)

        ub_r = interpolate(mesh,sdata,face_info,function_info,iu_b, 'value', ME)
        ub_l = interpolate(mesh,sdata,face_info,function_info,iu_b, 'value', NEIGHBOR)


        !
        ! Compute boundary average flux for u_a
        !
        flux_x = ((cx*ua_l + cx*ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = ((cy*ua_l + cy*ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = ((cz*ua_l + cz*ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)

        integrand = flux_x + flux_y + flux_z
        call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu_a,integrand)



        !
        ! Compute boundary average flux for u_b
        !
        flux_x = ((cx*ub_l + cx*ub_r)/TWO )  *  norms(:,1)
        flux_y = ((cy*ub_l + cy*ub_r)/TWO )  *  norms(:,2)
        flux_z = ((cz*ub_l + cz*ub_r)/TWO )  *  norms(:,3)

        integrand = flux_x + flux_y + flux_z
        call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu_b,integrand)



        end associate
    end subroutine compute
    !************************************************************************************







end module DLA_boundary_average_advective_flux
