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
module mod_chidg_airfoil
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, TWO
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf
    use type_element_info,      only: element_info_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_chidg_cache,       only: chidg_cache_t
    use type_chidg_manager,     only: chidg_manager_t
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
    subroutine chidg_airfoil(filename)
        character(*)    :: filename
    
    type(chidg_manager_t)                       :: manager
        type(chidg_t)               :: chidg
        type(file_properties_t)     :: file_props
        integer(ik)                 :: nterms_s, spacedim, solution_order, group_ID, &
                                       ibc, patch_ID, face_ID, idomain_g, &
                                       ielement_g, iface, itime
        logical                     :: found_airfoil

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

        type(AD_D)  :: lift, drag

        gq_rule = 3

        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')


        !
        ! Get nterms_s from file
        !
        file_props = get_properties_hdf(filename)
        nterms_s   = file_props%nterms_s(1)
        spacedim   = file_props%spacedim(1)

        solution_order = 0
        do while (solution_order*solution_order*solution_order < nterms_s)
            solution_order = solution_order + 1
        end do





        ! Set linear solver options to pass during initialization
        call loptions%set('tol',1.e-9_rk)

        ! Set nonlinear solver options
        call noptions%set('tol',3.e-5_rk)
        call noptions%set('cfl0',1.0_rk)
        call noptions%set('nsteps',100)




        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%set('Time Integrator' , algorithm=time_integrator)
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions)
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions)
        call chidg%set('Preconditioner'  , algorithm=preconditioner                    )




        !
        ! Initialize solution data storage
        !
        !
        ! Read grid data from file
        !
        gridfile = filename
        call chidg%read_grid(filename,spacedim)

        call manager%process(chidg)



        !
        ! Read solution modes from HDF5
        !
        call chidg%read_solution(filename)
        chidg%data%sdata%q = chidg%data%sdata%q_in



        !
        ! Initialize Chidg Worker references
        !
        call worker%init(chidg%data%mesh, chidg%data%eqnset(:)%prop, chidg%data%sdata, cache)




        !
        ! Get 'Airfoil' boundary group ID
        !
        group_ID = chidg%data%mesh%get_bc_patch_group_id('Airfoil')


        !
        ! Check if an 'Airfoil' boundary was found
        !
        if (group_ID == 0) call chidg_signal(FATAL,"chidg airfoil: No airfoil boundary was found.")




        !
        ! Loop over domains/elements/faces for 'Airfoil' patches
        !
        lift = AD_D(1)
        drag = AD_D(1)

        lift = ZERO
        drag = ZERO
        do patch_ID = 1,size(chidg%data%mesh%bc_patch_group(group_ID)%patch)

            !
            ! Loop over faces in the patch
            !
            do face_ID = 1,chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                idomain_g  = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_g(face_ID)
                ielement_g = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_g(face_ID)
                iface      = chidg%data%mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                call write_line('Airfoil: ', idomain_g, ielement_g, iface)


                !
                ! Initialize element location object
                ! 
                elem_info%idomain_g  = idomain_g
                elem_info%idomain_l  = idomain_g
                elem_info%ielement_g = ielement_g
                elem_info%ielement_l = ielement_g
                call worker%set_element(elem_info)
                worker%itime = 1


                !
                ! Update the element cache and all models so they are available
                !
                call cache_handler%update(worker,chidg%data%eqnset,chidg%data%bc_state_group, differentiate=.false.)



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
                weights = worker%mesh%domain(idomain_g)%faces(ielement_g,iface)%gq%face%weights(:,iface)
                areas   = sqrt(norm_1*norm_1 + norm_2*norm_2 + norm_3*norm_3)

                lift = lift + sum( stress_y * weights * areas)
                drag = drag + sum( stress_x * weights * areas)






            end do !iface

        end do !ipatch



        call write_line('Lift: ', lift%x_ad_)
        call write_line('Drag: ', drag%x_ad_)



        !
        ! Close ChiDG
        !
        call chidg%shut_down('core')



    end subroutine chidg_airfoil
    !******************************************************************************************


    



end module mod_chidg_airfoil
