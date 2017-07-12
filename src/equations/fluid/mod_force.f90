!>  Compute
!!
!!  @author Nathan A. Wukie (AFRL)
!!  @date   7/12/2017
!!  @note   Modified directly from 'chidg airfoil' action
!!
!!
!---------------------------------------------------------------------------------------------
module mod_force
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, TWO, NO_ID
    use mod_chidg_mpi,          only: ChiDG_COMM
    use type_chidg_data,        only: chidg_data_t
    use type_element_info,      only: element_info_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use mpi_f08,                only: MPI_AllReduce, MPI_REAL8, MPI_SUM
    use DNAD_D
    implicit none







contains



    !>  Compute force integrated over a specified patch group.
    !!
    !!
    !!  F = int[ (tau-p) dot n ] dPatch
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/12/2017
    !!  @note   Modified directly from 'chidg airfoil' action
    !!
    !!
    !!  @param[in]      data            chidg_data instance
    !!  @param[in]      patch_group     Name of patch group over which the force will be integrated.
    !!  @result[out]    force           Integrated force vector: force = [f1, f2, f3]
    !!
    !-----------------------------------------------------------------------------------
    function compute_force(data,patch_group) result(force_reduced)
        type(chidg_data_t), intent(inout)   :: data
        character(*),       intent(in)      :: patch_group
    
        integer(ik)                 :: group_ID, patch_ID, face_ID, &
                                       idomain_g,  idomain_l,        &
                                       ielement_g, ielement_l, iface

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler
        type(element_info_t)        :: elem_info


        real(rk)                                ::  &
            force(3), force_reduced(3)


        real(rk),   allocatable, dimension(:)   ::  &
            norm_1,  norm_2,  norm_3,               &
            unorm_1, unorm_2, unorm_3,              &
            weights, areas


        type(AD_D), allocatable, dimension(:)   ::  &
            tau_11,     tau_12,     tau_13,         &
            tau_21,     tau_22,     tau_23,         &
            tau_31,     tau_32,     tau_33,         &
            stress_x,   stress_y,   stress_z,       &
            pressure,   normal_stress

        integer(ik) :: ierr



        !
        ! Initialize Chidg Worker references
        !
        call worker%init(data%mesh, data%eqnset(:)%prop, data%sdata, cache)


        !
        ! Get patch_group boundary group ID
        !
        group_ID = data%mesh%get_bc_patch_group_id(trim(patch_group))


        !
        ! Check if a group matching "patch_group" was found
        !
        if (group_ID == NO_ID) call chidg_signal(FATAL,"chidg airfoil: No airfoil boundary was found.")


        !
        ! Loop over domains/elements/faces for "patch_group" 
        !
        force = ZERO
        do patch_ID = 1,size(data%mesh%bc_patch_group(group_ID)%patch)

            !
            ! Loop over faces in the patch
            !
            do face_ID = 1,data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                idomain_g  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                idomain_l  = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                ielement_g = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                ielement_l = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                iface      = data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                !
                ! Initialize element location object
                ! 
                elem_info%idomain_g  = idomain_g
                elem_info%idomain_l  = idomain_l
                elem_info%ielement_g = ielement_g
                elem_info%ielement_l = ielement_l
                call worker%set_element(elem_info)
                worker%itime = 1


                !
                ! Update the element cache and all models so they are available
                !
                call cache_handler%update(worker,data%eqnset,data%bc_state_group, components    = 'all',   &
                                                                                  face          = NO_ID,   &
                                                                                  differentiate = .false., &
                                                                                  lift          = .true.)



                call worker%set_face(iface)


                !
                ! Get pressure
                !
                pressure = worker%get_model_field_face('Pressure', 'value', 'face interior')


                
                !
                ! Get shear stress tensor
                !
                tau_11 = worker%get_model_field_face('Shear-11', 'value', 'face interior')
                tau_22 = worker%get_model_field_face('Shear-22', 'value', 'face interior')
                tau_33 = worker%get_model_field_face('Shear-33', 'value', 'face interior')
                tau_12 = worker%get_model_field_face('Shear-12', 'value', 'face interior')
                tau_13 = worker%get_model_field_face('Shear-13', 'value', 'face interior')
                tau_23 = worker%get_model_field_face('Shear-23', 'value', 'face interior')

                ! From symmetry
                tau_21 = tau_12
                tau_31 = tau_13
                tau_32 = tau_23

                
                !
                ! Add pressure component
                !
                tau_11 = tau_11 - pressure
                tau_22 = tau_22 - pressure
                tau_33 = tau_33 - pressure


                !
                ! Get normal vectors and reverse, because we want outward-facing vector from
                ! the geometry.
                !
                norm_1  = -worker%normal(1)
                norm_2  = -worker%normal(2)
                norm_3  = -worker%normal(3)

                unorm_1 = -worker%unit_normal(1)
                unorm_2 = -worker%unit_normal(2)
                unorm_3 = -worker%unit_normal(3)
                

                !
                ! Compute \vector{n} dot \tensor{tau}
                !   : These should produce the same result since the tensor is 
                !   : symmetric. Not sure which is more correct.
                !
                !stress_x = unorm_1*tau_11 + unorm_2*tau_21 + unorm_3*tau_31
                !stress_y = unorm_1*tau_12 + unorm_2*tau_22 + unorm_3*tau_32
                !stress_z = unorm_1*tau_13 + unorm_2*tau_23 + unorm_3*tau_33
                stress_x = tau_11*unorm_1 + tau_12*unorm_2 + tau_13*unorm_3
                stress_y = tau_21*unorm_1 + tau_22*unorm_2 + tau_23*unorm_3
                stress_z = tau_31*unorm_1 + tau_32*unorm_2 + tau_33*unorm_3



                !
                ! Integrate
                !
                weights = worker%mesh%domain(idomain_l)%faces(ielement_l,iface)%basis_s%weights(iface)
                areas   = sqrt(norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)

                force(1) = force(1) + sum( stress_x(:)%x_ad_ * weights * areas)
                force(2) = force(2) + sum( stress_y(:)%x_ad_ * weights * areas)
                force(3) = force(3) + sum( stress_z(:)%x_ad_ * weights * areas)

            end do !iface

        end do !ipatch



        !
        ! Reduce result across processors
        !
        call MPI_AllReduce(force,force_reduced,3,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)


    end function compute_force
    !******************************************************************************************


    



end module mod_force
