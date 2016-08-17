module BU_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                      LOCAL, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
    use DNAD_D

    use BU_properties,          only: BU_properties_t
    implicit none

    private



    !--------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: BU_LaxFriedrichs_flux_t


    contains
        procedure   :: compute

    end type BU_LaxFriedrichs_flux_t
    !---------------------------------------------------------------------------

contains




    !
    !
    !
    !
    !
    !---------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(BU_LaxFriedrichs_flux_t),     intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        type(face_info_t),                  intent(in)      :: face_info
        type(function_info_t),              intent(in)      :: function_info

        integer(ik)              :: idom, ielem, iface
        integer(ik)              :: iblk, idonor, ifcn

        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu, ierr, nnodes
        type(seed_t)             :: seed
        type(AD_D), dimension(mesh(face_info%idomain)%faces(face_info%ielement,face_info%iface)%gq%face%nnodes) :: &
                            u_l, u_r, flux_x, flux_y, flux_z, integrand


        !
        ! Get integer data
        !
        iu = prop%get_eqn_index("u")



        idom  = face_info%idomain
        ielem = face_info%ielement
        iface = face_info%iface



        associate ( norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm )

        !
        ! Get equation set property data
        !
        select type(prop)
            type is (BU_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select



        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_face(mesh,face_info,sdata%q,iu, u_r, LOCAL)
        call interpolate_face(mesh,face_info,sdata%q,iu, u_l, NEIGHBOR)
 

        !
        ! Compute boundary upwind flux
        !
        flux_x = cx * max(abs(u_l),abs(u_r))*( (u_l - u_r)/TWO )  *  norms(:,1) * unorms(:,1)
        flux_y = cy * max(abs(u_l),abs(u_r))*( (u_l - u_r)/TWO )  *  norms(:,2) * unorms(:,2)
        flux_z = cz * max(abs(u_l),abs(u_r))*( (u_l - u_r)/TWO )  *  norms(:,3) * unorms(:,3)

        integrand = flux_x + flux_y + flux_z



        !
        ! Integrate flux
        !
        call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu,integrand)


        end associate


    end subroutine









end module BU_LaxFriedrichs_flux
