module LA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NFACES,ZERO,ONE,TWO,HALF, &
                                          XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                          LOCAL, NEIGHBOR

    use atype_boundary_flux,        only: boundary_flux_t
    use type_mesh,                  only: mesh_t
    use type_solverdata,            only: solverdata_t
    use type_properties,            only: properties_t
    use type_seed,                  only: seed_t
    use type_face_info,             only: face_info_t
    use type_function_info,         only: function_info_t


    use mod_interpolate,            only: interpolate_face
    use mod_integrate,              only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
    use DNAD_D

    use LA_properties,              only: LA_properties_t
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LA_boundary_average_advective_flux_t


    contains
        procedure   :: compute

    end type LA_boundary_average_advective_flux_t
    !********************************************************************************



contains

    !> Compute the average advective boundary flux for scalar linear advection
    !!
    !!   @author Nathan A. Wukie
    !!
    !!   @param[in]      mesh    Mesh data
    !!   @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !!   @param[in]      ielem   Element index
    !!   @param[in]      iface   Face index
    !!   @param[in]      iblk    Block index indicating the linearization direction
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(LA_boundary_average_advective_flux_t),    intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        type(face_info_t),                              intent(in)      :: face_info
        type(function_info_t),                          intent(in)      :: function_info


        integer(ik)                 :: idom, ielem, iface
        integer(ik)                 :: iblk, ifcn, idonor

        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, ierr, nnodes, i
        type(seed_t)                :: seed
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                                        u_l, u_r, flux_x, flux_y, flux_z, integrand


        !
        ! Get variable index
        !
        iu = prop%get_eqn_index("u")



        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface


        associate ( norms => mesh(idom)%faces(ielem,iface)%norm )


        !
        ! Get equation set properties
        !
        select type(prop)
            type is (LA_properties_t)
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
        ! Compute boundary average flux
        !
        flux_x = ((cx*u_r + cx*u_l)/TWO ) * norms(:,1)
        flux_y = ((cy*u_r + cy*u_l)/TWO ) * norms(:,2)
        flux_z = ((cz*u_r + cz*u_l)/TWO ) * norms(:,3)

        integrand = flux_x + flux_y + flux_z


        !
        ! Integrate flux
        !
        call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu,integrand)


        end associate
    end subroutine compute
    !**************************************************************************************************




end module LA_boundary_average_advective_flux
