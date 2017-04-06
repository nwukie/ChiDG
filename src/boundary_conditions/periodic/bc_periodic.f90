module bc_periodic
#include <messenger.h>
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use type_mesh,          only: mesh_t
    use type_point,         only: point_t
    implicit none





    !> Periodic boundary condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(bc_state_t) :: periodic_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type periodic_t
    !********************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(periodic_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Periodic")
        call self%set_family("Periodic")




        !
        ! Add parameters
        !
        call self%bcproperties%add('Offset-1', 'Required')
        call self%bcproperties%add('Offset-2', 'Required')
        call self%bcproperties%add('Offset-3', 'Required')


        !
        ! Set default values
        !
        call self%set_fcn_option('Offset-1', 'val', 0._rk)
        call self%set_fcn_option('Offset-2', 'val', 0._rk)
        call self%set_fcn_option('Offset-3', 'val', 0._rk)


    end subroutine init
    !********************************************************************************








    !>  Set periodic offset data on the mesh.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/6/2017
    !!
    !---------------------------------------------------------------------------------
    subroutine init_bc_specialized(self,mesh)
        class(periodic_t),  intent(in)  :: self
        type(mesh_t),       intent(in)  :: mesh

        integer(ik)     :: ipatch, ipatch_face, idom, ielem, iface
        real(rk)        :: time
        type(point_t)   :: pnt

!        !
!        ! Loop over patches and set periodic offset parameters.
!        !
!        do ipatch = 1,mesh%bc_patch_group(self%group_ID)%npatches()
!            do ipatch_face = 1,mesh%bc_patch_group(self%group_ID)%patch(ipatch)%nfaces()
!
!                !
!                ! Get geometry locator
!                !
!                idom  = mesh%bc_patch_group(self%group_ID)%patch(ipatch)%idomain_l(ipatch_face)
!                ielem = mesh%bc_patch_group(self%group_ID)%patch(ipatch)%ielement_l(ipatch_face)
!                iface = mesh%bc_patch_group(self%group_ID)%patch(ipatch)%iface(ipatch_face)
!
!
!                !
!                ! Set face type to ORPHAN face so it will be recognized as chimera
!                !
!                mesh%domain(idom)%faces(ielem,iface)%ftype = ORPHAN
!
!
!                !
!                ! time, pnt do nothing here, but interface for function requires them.
!                !
!                mesh%domain(idom)%faces(ielem,iface)%periodic_offset  = .true.
!                mesh%domain(idom)%faces(ielem,iface)%chimera_offset_1 = self%state%bcproperties%compute('Offset-1',time,pnt)
!                mesh%domain(idom)%faces(ielem,iface)%chimera_offset_2 = self%state%bcproperties%compute('Offset-2',time,pnt)
!                mesh%domain(idom)%faces(ielem,iface)%chimera_offset_3 = self%state%bcproperties%compute('Offset-3',time,pnt)
!
!
!            end do !ipatch_face
!        end do !ipatch



    end subroutine init_bc_specialized
    !**********************************************************************************










    !> Boundary condition compute routine called by spatial scheme
    !!      - Matching periodic boundary condition, so the interior scheme 
    !!        is just adjusted to connect elements and no extra calculation
    !!        routine is needed here. Hence the empty routine below.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(periodic_t),              intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! DO NOTHING IN PERIODIC BOUNDARY CONDITION

    end subroutine compute_bc_state
    !***********************************************************************************













end module bc_periodic
