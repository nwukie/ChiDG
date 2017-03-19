!>  ChiDG HDF5 to vtk conversion utility
!!
!!  @author Mayank Sharma
!!  @date 10/31/2016
!!
!!
!!  Usage: chidg post 'chidgfile'
!!  
!
!----------------------------------------------------------------------------------------------
module mod_chidg_post_hdf2vtk
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_vtk_calc_func,      only: get_Fourier_coeff_vector, get_Fourier_coeffs, get_interp_time_levels, &
                                      get_interp_solution, get_interp_solution_vector
    use mod_vtkio,              only: write_vtk_file
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf
    use type_chidg_vector,      only: chidg_vector_t
    use mod_HB_matrices,        only: calc_E

    implicit none



contains



    !>  Post processing tool for writing a tecplot file of sampled modal data
    !!
    !!  @author Mayank Sharma
    !!  @date 10/31/2016
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_post_hdf2vtk(filename)

        character(*)                    ::  filename


        type(chidg_t)                       ::  chidg
        type(file_properties_t)             ::  file_props
        character(:),           allocatable ::  eqnset
        character(:),           allocatable ::  time_string
        real(rk),               allocatable ::  freq(:), time_lev(:)
        integer(ik)                         ::  nterms_s,spacedim,solution_order,nfreq,ntime,ierr
        type(chidg_vector_t)                ::  q_HB               ! Storage vector for original HB solution
        real(rk),               allocatable ::  E(:,:)
        type(chidg_vector_t)                ::  q_coeff_vector     ! Vector containing coefficients of solution &
                                                                   ! Fourier expansion in time
        type(chidg_vector_t),   allocatable ::  q_coeffs(:)
        type(chidg_vector_t),   allocatable ::  q_interp(:)
        integer(ik)                         ::  ntime_interp, itime_interp
        real(rk),               allocatable ::  time_interp(:)
        character(:),           allocatable ::  flag
        character(len = 100)                ::  new_dir_path_1, new_dir_path_2, &
                                                pvd_filename_1, pvd_filename_2



        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')

        !
        ! Get nterms_s and eqnset
        !
        file_props  = get_properties_hdf(filename)

        nterms_s    = file_props%nterms_s(1)     ! Global variable from mod_io 
        eqnset      = file_props%eqnset(1)       ! Global variable from mod_io
        spacedim    = file_props%spacedim(1)     ! Global variable from mod_io
        time_string = file_props%time_integrator
        freq        = file_props%HB_frequencies; nfreq = size(freq)
        time_lev    = file_props%HB_time_lev;    ntime = size(time_lev)


        !
        ! Read grid data from file
        !
        call chidg%read_grid(filename,spacedim)


        solution_order = 0
        do while (solution_order*solution_order*solution_order < nterms_s)
            solution_order = solution_order + 1
        end do

        !
        ! Initialize solution data storage
        !
        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%init('domains')
        call chidg%init('communication')
        call chidg%init('solvers')


        !
        ! Read solution modes from HDF5
        !
        call chidg%read_solution(filename)
        
        
        !
        ! If harmonic balance is being used, interpolate the solution at given times and write output
        !
        flag = trim(time_string)

        if (flag == 'Harmonic Balance' .or. flag == 'Harmonic_Balance' .or. flag == 'harmonic balance' .or. &
            flag == 'harmonic_balance' .or. flag == 'HB') then 

            !
            ! Write original solution in vtk format
            !
            new_dir_path_1 = 'ChiDG_HB_results'
            pvd_filename_1 = 'chidg_HB_results.pvd'
            call write_vtk_file(chidg%data,new_dir_path_1,pvd_filename_1)

            
            !
            ! Store original solution in a separate vector
            !
            call q_HB%init(chidg%data%mesh,chidg%data%sdata%q_in%get_ntime())
            call q_HB%clear()
            q_HB = chidg%data%sdata%q_in


            !
            ! Compute Fourier transform matrix
            !
            if (allocated(E)) deallocate(E)
            allocate(E(ntime,ntime), stat=ierr)
            if (ierr /= 0) call AllocationError
            
            call calc_E(nfreq,ntime,freq,time_lev,E)


            !
            ! Compute coefficients of Fourier expansion of solution in time
            !
            call get_Fourier_coeff_vector(chidg%data,nterms_s,E,q_coeff_vector)
            call get_Fourier_coeffs(chidg%data,nterms_s,q_coeff_vector,q_coeffs)


            !
            ! Compute interpolation times
            !
            call get_interp_time_levels(ntime_interp,time_interp)
            

            !
            ! Compute interpolated solutions
            !
            if (allocated(q_interp)) deallocate(q_interp)
            allocate(q_interp(ntime_interp), stat=ierr)
            if (ierr /= 0) call AllocationError

            do itime_interp = 1,ntime_interp

                call get_interp_solution(chidg%data,nterms_s,time_interp(itime_interp), &
                                         freq,q_coeffs,q_interp(itime_interp))

            end do
           
           
            !
            ! Store the interpolated solutions in resized vector, q_in
            ! 
            call get_interp_solution_vector(chidg%data,nterms_s,q_interp)


            !
            ! Write interpolated solution in vtk format
            !
            new_dir_path_2 = 'ChiDG_interpolated_results'
            pvd_filename_2 = 'chidg_interp_results.pvd'
            call write_vtk_file(chidg%data,new_dir_path_2,pvd_filename_2)
        
        else
            
            new_dir_path_1 = 'ChiDG_results'
            pvd_filename_1 = 'chidg_results.pvd'
            call write_vtk_file(chidg%data, new_dir_path_1,pvd_filename_1)

        end if


        !
        ! Shut down ChiDG
        !
        call chidg%shut_down('core')


    end subroutine chidg_post_hdf2vtk 
    !******************************************************************************************




















end module mod_chidg_post_hdf2vtk
