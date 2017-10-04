module mod_HB_post
#include<messenger.h>
    use mod_kinds,          only: rk,ik,rdouble
    use mod_constants,      only: ONE,HALF,TWO,OUTPUT_RES
    use type_chidg_data,    only: chidg_data_t
    use type_chidg_vector,  only: chidg_vector_t
    use mod_HB_matrices,    only: calc_E

    implicit none



contains



    !>  Compute vector of coefficients of Fourier expansion of solution in time
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_Fourier_coeff_vector(data,E,q_coeff_vector)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk)                                :: E(:,:)
        type(chidg_vector_t),   intent(inout)   :: q_coeff_vector

        real(rk),   allocatable     :: temp_1(:),temp_2(:)
        integer(ik)                 :: otime, itime, ielem, idom, ierr, ivar, eqn_ID

        
        associate ( q_in => data%sdata%q_in)

         
            !
            ! Initialize Fourier coefficient chidg_vector
            !
            call q_coeff_vector%init(data%mesh,q_in%get_ntime())
            call q_coeff_vector%set_ntime(q_in%get_ntime())
            call q_coeff_vector%clear()


            !
            ! Generate Fourier coefficients for solution expansion in time
            ! Reference: Knapke, R.D., "High-Order Unsteady Heat Transfer with the Harmonic
            !            Balance Method", Ph.D. Dissertation, University of Cincinnati, pg. 40-41
            !
            do otime = 1,q_in%get_ntime()
                do idom = 1,data%mesh%ndomains()
                    
                    eqn_ID = data%mesh%domain(idom)%eqn_ID

                    associate ( domain   => data%mesh%domain(idom),          &
                                nelem    => data%mesh%domain(idom)%nelem,    &
                                nterms_s => data%mesh%domain(idom)%nterms_s, &
                                nvars    => data%eqnset(eqn_ID)%prop%nprimary_fields()) 
                    
                    if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                    allocate(temp_1(nterms_s), temp_2(nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    do ielem = 1,nelem
                        do itime = 1,q_in%get_ntime()
                            do ivar = 1,nvars

                            temp_1 = E(otime,itime)*q_in%dom(idom)%vecs(ielem)%getvar(ivar,itime)
                            temp_2 = q_coeff_vector%dom(idom)%vecs(ielem)%getvar(ivar,otime) + temp_1
                            call q_coeff_vector%dom(idom)%vecs(ielem)%setvar(ivar,otime,temp_2)

                            end do  ! ivar
                        end do  ! itime
                    end do  ! ielem

                    end associate  

                end do  ! idom
            end do  ! otime

        end associate


    end subroutine get_Fourier_coeff_vector
    !*******************************************************************************************



    !>  Compute the individual coefficients of Fourier expansion of solution in time
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_Fourier_coeffs(data,q_coeff_vector,q_coeffs)
        type(chidg_data_t),                     intent(inout)   :: data
        type(chidg_vector_t),                   intent(inout)   :: q_coeff_vector
        type(chidg_vector_t),   allocatable,    intent(inout)   :: q_coeffs(:)

        real(rk),   allocatable :: temp(:)
        integer(ik)             :: ntime, itime, idom, ielem, ivar, ierr, eqn_ID


        !
        ! Get number of time levels
        !
        ntime = q_coeff_vector%get_ntime()
        
        if (allocated(q_coeffs)) deallocate(q_coeffs)
        allocate(q_coeffs(ntime), stat = ierr)
        if (ierr /= 0) call AllocationError


        do itime = 1,ntime

            !
            ! Initialize coefficient vector at itime
            !
            call q_coeffs(itime)%init(data%mesh,1)
            call q_coeffs(itime)%clear()

            do idom = 1,data%mesh%ndomains()

                eqn_ID = data%mesh%domain(idom)%eqn_ID

                associate ( nelem    => data%mesh%domain(idom)%nelem,    &
                            nterms_s => data%mesh%domain(idom)%nterms_s, &
                            nvars    => data%eqnset(eqn_ID)%prop%nprimary_fields() )

                if (allocated(temp)) deallocate(temp)
                allocate(temp(nterms_s), stat = ierr)
                if (ierr /= 0) call AllocationError
            
                do ielem = 1,nelem
                    do ivar = 1,nvars
                        
                        !
                        ! Set the values in q_coeffs(itime) = values in itime component
                        ! of the main coefficient vector
                        !
                        temp = q_coeff_vector%dom(idom)%vecs(ielem)%getvar(ivar,itime)
                        call q_coeffs(itime)%dom(idom)%vecs(ielem)%setvar(ivar,1,temp)

                    end do  ! ivar
                end do  ! ielem

                end associate

            end do  ! idom
        end do  ! itime


    end subroutine get_Fourier_coeffs
    !*******************************************************************************************



    !>  Compute interpolation times
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_interp_time_levels(data,ntime_interp,time_interp)
        type(chidg_data_t),         intent(inout)   :: data
        integer(ik),                intent(inout)   :: ntime_interp
        real(rk),   allocatable,    intent(inout)   :: time_interp(:)
        
        integer(ik)     :: itime_interp, ierr
        

        ntime_interp = data%time_manager%ntime
        !ntime_interp = 21

        if (allocated(time_interp)) deallocate(time_interp)
        allocate(time_interp(ntime_interp), stat=ierr) 
        if (ierr /= 0) call AllocationError

        !do itime_interp = 1,ntime_interp

        !    !
        !    ! TODO: Change post processing from 0 - 1 to 0 - t_max
        !    ! t_max is the maximum time level obtained from the HB frequencies (2*PI/omega_min)
        !    ! TODO: Can also be changed to 0 - t_desired
        !    !
        !do itime_interp = 1,ntime_interp
        !    time_interp(itime_interp) = real(itime_interp - 1,rk)*(1.0_rk&
        !                                /real(ntime_interp - 1,rk))

        !end do

        time_interp = data%time_manager%times

    end subroutine get_interp_time_levels
    !*******************************************************************************************



    !>  Compute interpolated solution at a given time
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_interp_solution(data,time,q_coeffs,q_time)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk),               intent(in)      :: time
        type(chidg_vector_t),   intent(inout)   :: q_coeffs(:)
        type(chidg_vector_t),   intent(inout)   :: q_time

        integer(ik)                 :: ntime, nfreq, itime, idom, ielem, ivar, ierr, eqn_ID
        real(rk),   allocatable     :: temp_1(:), temp_2(:), freq_data(:)
        character(:),   allocatable :: user_msg, dev_msg

        
        !
        ! Get ntime and nfreq
        !
        !ntime = size(q_coeffs)
        !nfreq = data%time_manager%freq_data%size()
        ntime = size(data%time_manager%times)
        nfreq = size(data%time_manager%freqs)
        

        !
        ! Get frequency data
        !
        !if (allocated(freq_data)) deallocate(freq_data)
        !allocate(freq_data(nfreq), stat=ierr)
        !if (ierr /= 0) call AllocationError
        !freq_data = data%time_manager%freq_data%data()
        freq_data = data%time_manager%freqs

        
        !
        ! Initialize the interpolated solution vector
        !
        call q_time%init(data%mesh,1)
        call q_time%clear()


        do itime = 1,ntime
            do idom = 1,data%mesh%ndomains()

                eqn_ID = data%mesh%domain(idom)%eqn_ID

                associate ( nelem    => data%mesh%domain(idom)%nelem,    &
                            nterms_s => data%mesh%domain(idom)%nterms_s, &
                            nvars    => data%eqnset(eqn_ID)%prop%nprimary_fields() )
                    
                    ! Allocate temporary arrays
                    if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                    allocate(temp_1(nterms_s), temp_2(nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    do ielem = 1,nelem
                        do ivar = 1,nvars

                            !
                            ! Compute interpolated solution
                            ! Reference: Knapke, R.D., "High-Order Unsteady Heat Transfer with the Harmonic
                            !            Balance Method", Ph.D. Dissertation, University of Cincinnati, pg. 40-41
                            !
                            if (itime == 1) then

                                temp_1 = q_coeffs(itime)%dom(idom)%vecs(ielem)%getvar(ivar,1)
                                call q_time%dom(idom)%vecs(ielem)%setvar(ivar,1,temp_1)

                            else if (itime .ge. 2 .and. itime .le. nfreq + 1) then

                                temp_1 = q_coeffs(itime)%dom(idom)%vecs(ielem)%getvar(ivar,1)*&
                                         sin(freq_data(itime - 1)*time)
                                temp_2 = q_time%dom(idom)%vecs(ielem)%getvar(ivar,1) + temp_1
                                call q_time%dom(idom)%vecs(ielem)%setvar(ivar,1,temp_2)

                            else if (itime .ge. nfreq + 2 .and. itime .le. ntime) then

                                temp_1 = q_coeffs(itime)%dom(idom)%vecs(ielem)%getvar(ivar,1)*&
                                         cos(freq_data(itime - nfreq - 1)*time)
                                temp_2 = q_time%dom(idom)%vecs(ielem)%getvar(ivar,1) + temp_1
                                call q_time%dom(idom)%vecs(ielem)%setvar(ivar,1,temp_2)

                            else

                                user_msg = 'There is a discrepancy between the number of HB frequencies and &
                                            number of HB time levels'
                                dev_msg  = 'q_in might be initialized incorrectly in read_solution'
                                call chidg_signal_two(FATAL, user_msg, ntime, dev_msg)

                            end if

                        end do  ! ivar
                    end do  ! ielem

                end associate
                
            end do  ! idom
        end do  ! itime


    end subroutine get_interp_solution
    !*******************************************************************************************



    !>  Store all interpolated solutions in a resized chidg_vector, q_in
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_interp_solution_vector(data,q_interp)
        type(chidg_data_t),     intent(inout)   :: data
        type(chidg_vector_t),   intent(in)      :: q_interp(:)

        real(rk),   allocatable:: temp(:)
        integer(ik)            :: ntime_interp, itime_interp, idom, ielem, &
                                  ivar, ierr, eqn_ID

        
        !
        ! Get number of interpolated solutions
        ! 
        ntime_interp = size(q_interp)

        
        !
        ! Reinitialize q_out
        !
        call data%sdata%q_out%init(data%mesh,ntime_interp)
        call data%sdata%q_out%set_ntime(ntime_interp)
        call data%sdata%q_out%clear()


        !
        ! Generate new q_in where itime component of q_in = q_interp(itime)
        !
        do itime_interp = 1,ntime_interp
            do idom = 1,data%mesh%ndomains()

                eqn_ID = data%mesh%domain(idom)%eqn_ID

                associate ( nelem    => data%mesh%domain(idom)%nelem,    &
                            nterms_s => data%mesh%domain(idom)%nterms_s, &
                            nvars    => data%eqnset(eqn_ID)%prop%nprimary_fields() )
                
                    if (allocated(temp)) deallocate(temp)
                    allocate(temp(nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    do ielem = 1,nelem
                        do ivar = 1,nvars
                            
                            temp = q_interp(itime_interp)%dom(idom)%vecs(ielem)%getvar(ivar,1)
                            call data%sdata%q_out%dom(idom)%vecs(ielem)%setvar(ivar,itime_interp,temp)

                        end do
                    end do

                end associate

            end do
        end do


    end subroutine get_interp_solution_vector
    !*******************************************************************************************




    !>  
    !!
    !!  @author Mayank Sharma
    !!  @date   3/20/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine get_post_processing_data(data)
        type(chidg_data_t),     intent(inout)   :: data

        real(rk),   allocatable                 :: E(:,:), freq(:), time_lev(:)
        type(chidg_vector_t)                    :: q_coeff_vector
        type(chidg_vector_t),   allocatable     :: q_coeffs(:), q_interp(:)
        real(rk),               allocatable     :: time_interp(:)
        integer(ik)                             :: ntime, nfreq, ierr, ntime_interp, &
                                                   itime_interp    
    
        
        !
        ! Set no. of frequencies and no. of time levels
        ! Also set frequency and time level data
        !
        !nfreq    = data%time_manager%freq_data%size()
        !ntime    = data%time_manager%time_lev%size()
        !freq     = data%time_manager%freq_data%data()
        !time_lev = data%time_manager%time_lev%data()

        nfreq    = size(data%time_manager%freqs)
        ntime    = size(data%time_manager%times)
        freq     = data%time_manager%freqs
        time_lev = data%time_manager%times

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
        call get_Fourier_coeff_vector(data,E,q_coeff_vector)
        call get_Fourier_coeffs(data,q_coeff_vector,q_coeffs)


        !
        ! Compute interpolation times
        !
        call get_interp_time_levels(data,ntime_interp,time_interp)


        !
        ! Compute interpolated solutions
        !
        if (allocated(q_interp)) deallocate(q_interp)
        allocate(q_interp(ntime_interp), stat=ierr)
        if (ierr /= 0) call AllocationError
    
        
        do itime_interp = 1,ntime_interp

            call get_interp_solution(data,time_interp(itime_interp), &
                                     q_coeffs,q_interp(itime_interp))

        end do


        !
        ! Store the interpolated solutions in resized vector, q_in
        !
        call get_interp_solution_vector(data,q_interp)


    end subroutine get_post_processing_data
    !*******************************************************************************************




















end module mod_HB_post
