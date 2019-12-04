!>  asdf
!!
!!  @author Nathan A. Wukie
!!  @date   3/8/2017
!!
!!
!! Usage:   chidg airfoil 'chidgfile'
!!
!!
!---------------------------------------------------------------------------------------------
module mod_chidg_forces
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, TWO, NO_ID, NO_DIFF
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf
    use type_element_info,      only: element_info_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_cache_handler,     only: cache_handler_t
    use mod_io
    use DNAD_D
    implicit none







contains



    !>  Post-processing tool for computing airfoil relevant quantities.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/8/2017
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine chidg_forces(filename,patch_group)
        character(*)    :: filename
        character(*)    :: patch_group
    
        character(:),   allocatable :: user_msg
        type(chidg_t)               :: chidg
        type(file_properties_t)     :: file_props
        integer(ik)                 :: nterms_s, solution_order, group_ID, &
                                       ibc, patch_ID, face_ID, idomain_g, &
                                       ielement_g, iface, itime, myunit, idomain_l, ielement_l
        logical                     :: exists

        type(chidg_worker_t)        :: worker
        type(chidg_cache_t)         :: cache
        type(cache_handler_t)       :: cache_handler

        type(element_info_t)        :: elem_info

        real(rk),   allocatable, dimension(:)   ::  &
            norm_1,  norm_2,  norm_3,               &
            unorm_1, unorm_2, unorm_3,              &
            weights, areas

        type(AD_D), allocatable, dimension(:)   ::  &
            tau_11,     tau_12,     tau_13,         &
            tau_21,     tau_22,     tau_23,         &
            tau_31,     tau_32,     tau_33,         &
            stress_x,   stress_y,   stress_z,       &
            pressure, normal_stress

        type(AD_D)  :: lift, drag, force(3)


        gq_rule = 3

        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('core')

        user_msg = "NOTE: if using a constant-viscosity model, make sure the models.nml &
                    file with the appropriate viscosity constant is available in the &
                    working directory. Otherwise the default viscosity will be used, &
                    which will cause forces to be computed incorrectly."
        call chidg_signal(WARN,user_msg)

        ! Get nterms_s from file
        file_props = get_properties_hdf(filename)
        nterms_s   = file_props%nterms_s(1)
        solution_order = 0
        do while (solution_order*solution_order*solution_order < nterms_s)
            solution_order = solution_order + 1
        end do


        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%set('Time Integrator' , algorithm=time_integrator)
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions)
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions)
        call chidg%set('Preconditioner'  , algorithm=preconditioner                    )



        ! Initialize solution data storage
        ! Read grid data from file
        gridfile = filename
        call chidg%read_mesh(filename)

        ! Read solution modes from HDF5
        call chidg%read_fields(filename)
        
        ! Process for getting wall distance
        call chidg%process()

        ! Initialize time integrator state
        call chidg%time_integrator%initialize_state(chidg%data)

        ! Initialize Chidg Worker references
        call worker%init(chidg%data%mesh, chidg%data%eqnset(:)%prop, chidg%data%sdata, chidg%data%time_manager, cache)



        ! Get patch group ID and check to make sure it was found.
        group_ID = chidg%data%mesh%get_bc_patch_group_id(trim(patch_group))
        if (group_ID == NO_ID) call chidg_signal_one(FATAL,"chidg forces: Patch was not found.", trim(patch_group))


        ! Loop over domains/elements/faces for 'Airfoil' patches
        force(1:3) = AD_D(1)
        force(1:3) = ZERO
        do patch_ID = 1,chidg%data%mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                idomain_g  = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g()
                idomain_l  = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                ielement_g = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                ielement_l = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                iface      = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)
                call write_line(trim(patch_group), idomain_g, ielement_g, iface)


                ! Initialize element location with worker
                elem_info = chidg%data%mesh%get_element_info(idomain_l,ielement_l)
                call worker%set_element(elem_info)
                worker%itime = 1

                ! Update the element cache and all models so they are available
                call cache_handler%update(worker,chidg%data%eqnset,chidg%data%bc_state_group,   &
                                                                   components    = 'all',       &
                                                                   face          = NO_ID,       &
                                                                   differentiate = NO_DIFF,     &
                                                                   lift          = .true.)     

                call worker%set_face(iface)

                ! Get pressure
                pressure = worker%get_field('Pressure', 'value', 'face interior')

                ! Get shear stress tensor
                tau_11 = worker%get_field('Shear-11', 'value', 'face interior')
                tau_22 = worker%get_field('Shear-22', 'value', 'face interior')
                tau_33 = worker%get_field('Shear-33', 'value', 'face interior')
                tau_12 = worker%get_field('Shear-12', 'value', 'face interior')
                tau_13 = worker%get_field('Shear-13', 'value', 'face interior')
                tau_23 = worker%get_field('Shear-23', 'value', 'face interior')

                ! From symmetry
                tau_21 = tau_12
                tau_31 = tau_13
                tau_32 = tau_23

                ! Add pressure component
                tau_11 = tau_11 - pressure
                tau_22 = tau_22 - pressure
                tau_33 = tau_33 - pressure


                ! Get normal vectors and reverse, because we want outward-facing vector from
                ! the geometry.
                norm_1  = -worker%normal(1)
                norm_2  = -worker%normal(2)
                norm_3  = -worker%normal(3)

                unorm_1 = -worker%unit_normal(1)
                unorm_2 = -worker%unit_normal(2)
                unorm_3 = -worker%unit_normal(3)
                

                ! Compute \vector{n} dot \tensor{tau}
                !   These should produce the same result since the tensor is 
                !   symmetric. Not sure which is more correct.
                !stress_x = unorm_1*tau_11 + unorm_2*tau_21 + unorm_3*tau_31
                !stress_y = unorm_1*tau_12 + unorm_2*tau_22 + unorm_3*tau_32
                !stress_z = unorm_1*tau_13 + unorm_2*tau_23 + unorm_3*tau_33
                stress_x = tau_11*unorm_1 + tau_12*unorm_2 + tau_13*unorm_3
                stress_y = tau_21*unorm_1 + tau_22*unorm_2 + tau_23*unorm_3
                stress_z = tau_31*unorm_1 + tau_32*unorm_2 + tau_33*unorm_3


                ! Integrate
                weights = worker%mesh%domain(idomain_g)%faces(ielement_g,iface)%basis_s%weights_face(iface)
                areas   = sqrt(norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)

                force(1) = force(1) + sum( stress_x * weights * areas)
                force(2) = force(2) + sum( stress_y * weights * areas)
                force(3) = force(3) + sum( stress_z * weights * areas)


            end do !iface
        end do !ipatch


        call write_line('Force: ', force(1)%x_ad_, force(2)%x_ad_, force(3)%x_ad_)


        if (IRANK == GLOBAL_MASTER) then
            inquire(file="forces.txt", exist=exists)
            if (exists) then
                open(newunit=myunit, file="forces.txt", status="old", position="append",action="write")
            else
                open(newunit=myunit, file="forces.txt", status="new",action="write")
                write(myunit,*) 'force-1', 'force-2', 'force-3'
            end if
            write(myunit,*) force(1)%x_ad_, force(2)%x_ad_, force(3)%x_ad_
            close(myunit)
        end if


        ! Close ChiDG
        call chidg%shut_down('core')


    end subroutine chidg_forces
    !******************************************************************************************


    



end module mod_chidg_forces
