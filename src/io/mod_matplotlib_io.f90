module mod_matplotlib_io
#include<messenger.h>
    use mod_kinds,          only: rk,ik,rdouble
    use mod_constants,      only: ZERO,ONE,HALF,TWO,OUTPUT_RES
    use type_point,         only: point_t
    use type_chidg_data,    only: chidg_data_t
    use type_chidg_vector,  only: chidg_vector_t
    use mod_chimera,        only: find_gq_donor
    use type_face_info,     only: face_info_t, face_info
    use type_element_info,  only: element_info_t
    use mod_HB_matrices,    only: calc_E
    use mod_HB_post,        only: get_Fourier_coeff_vector



contains



    !>  Get points along the domain for matplotlib computations
    !!
    !!  @author Mayank Sharma   
    !!  @date   3/24/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_points(points)
        type(point_t),  allocatable,    intent(inout)     :: points(:)

        integer(ik),    parameter       :: npts = 81
        real(rk),       allocatable     :: xpt(:), ypt(:), zpt(:)
        integer(ik)                     :: ipt, ierr


        !
        ! Allocate points array
        !
        if (allocated(xpt) .and. allocated(ypt) .and. allocated(zpt) .and. allocated(points)) deallocate(xpt,ypt,zpt,points)
        allocate(xpt(npts), ypt(npts), zpt(npts), points(npts), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Generate points
        ! TODO: Generalize to read from file
        !
        do ipt = 1,npts

            xpt(ipt) = real(ipt - 1)/real(npts - 1)
            ypt(ipt) = ZERO
            zpt(ipt) = ZERO

            ! Read points into type_point array
            call points(ipt)%set(xpt(ipt),ypt(ipt),zpt(ipt)) 

        end do


    end subroutine get_points
    !***********************************************************************************
    


    !>  Get donor point coordinates (xi, eta, zeta) for all points read in
    !!
    !!  @author Mayank Sharma
    !!  @date   3/24/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_donors(data,donors,donor_coords)
        type(chidg_data_t),                     intent(inout)   :: data
        type(element_info_t),   allocatable,    intent(inout)   :: donors(:)
        type(point_t),          allocatable,    intent(inout)   :: donor_coords(:)

        type(point_t),      allocatable     :: gq_nodes(:)
        type(face_info_t),  allocatable     :: receivers(:)
        integer(ik)                         :: npts, ipt, ierr
        real(rk)                            :: donor_coord(3)
        logical                             :: donor_found


        !
        ! Get gq_nodes and no. of gq_nodes
        !
        call get_points(gq_nodes)
        npts = size(gq_nodes)


        !
        ! Allocate donor coordinate arrays
        ! Alos allocate receiver_face array
        !
        if (allocated(donors) .and. allocated(donor_coords) .and. allocated(receivers)) deallocate(donors,donor_coords,receivers)
        allocate(donors(npts), donor_coords(npts), receivers(npts), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Generate donor point coordinates
        !
        do ipt = 1,npts
            receivers(ipt) = face_info(0_ik,0_ik,0_ik,0_ik,0_ik)
            !call find_gq_donor(data%mesh,gq_nodes(ipt),point_t(0._rk,0._rk,0._rk),receivers(ipt),donors(ipt),donor_coord, donor_found)
            call find_gq_donor(data%mesh,                                               &
                               [gq_nodes(ipt)%c1_,gq_nodes(ipt)%c2_,gq_nodes(ipt)%c3_], &
                               [0._rk,0._rk,0._rk],                                     &
                               receivers(ipt),                                          &
                               donors(ipt),                                             &
                               donor_coord,                                             &
                               donor_found)
            donor_coords(ipt)%c1_ = donor_coord(1)
            donor_coords(ipt)%c2_ = donor_coord(2)
            donor_coords(ipt)%c3_ = donor_coord(3)
            if (.not. donor_found) call chidg_signal(FATAL,"matplotlib_io::get_donors: interpolation donor not found.")
        end do


    end subroutine get_donors
    !***********************************************************************************



    !>  Get original solution values at donor points
    !!
    !!  @author Mayank Sharma
    !!  @date   3/27/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_original_solution_at_donor_points(data,ieqn,solution_values)
        type(chidg_data_t),         intent(inout)   :: data
        integer(ik),                intent(in)      :: ieqn
        real(rk),   allocatable,    intent(inout)   :: solution_values(:,:)

        type(element_info_t),   allocatable     :: donors(:)
        type(point_t),          allocatable     :: donor_coords(:)
        real(rdouble)                           :: val(1)                     
        integer(ik)                             :: ntime, npts, ipt, itime, ierr


        !
        ! Get donor points
        !
        call get_donors(data,donors,donor_coords)

        
        !
        ! Compute no. of points and no. of time levels and allocate solution values array
        !
        npts  = size(donors)
        ntime = data%sdata%q_in%get_ntime()

        if (allocated(solution_values)) deallocate(solution_values)
        allocate(solution_values(ntime,npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        
        !
        ! Compute solutions at donor points
        ! 
        do itime = 1,ntime
            do ipt = 1,npts

                associate( idom  => donors(ipt)%idomain_g,  &
                           ielem => donors(ipt)%ielement_g, &
                           xi    => donor_coords(ipt)%c1_,  &
                           eta   => donor_coords(ipt)%c2_,  &
                           zeta  => donor_coords(ipt)%c3_)

                    val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_in%dom(idom)%&
                                   vecs(ielem),ieqn,itime,xi,eta,zeta),rdouble)
                    solution_values(itime,ipt) = val(1)

                end associate

            end do  ! ipt
        end do  ! itime


    end subroutine get_original_solution_at_donor_points
    !***********************************************************************************



    !>  Get Fourier coefficient values at donor points
    !!
    !!  @author Mayank Sharma
    !!  @date   3/27/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_Fourier_coeffs_at_donor_points(data,ieqn,solution_values)
        type(chidg_data_t),         intent(inout)   :: data
        integer(ik),                intent(in)      :: ieqn
        real(rk),   allocatable,    intent(inout)   :: solution_values(:,:)

        type(chidg_vector_t)                    :: q_coeff_vector
        real(rk),               allocatable     :: E(:,:)
        type(element_info_t),   allocatable     :: donors(:)
        type(point_t),          allocatable     :: donor_coords(:)
        real(rdouble)                           :: val(1)                     
        integer(ik)                             :: ntime, npts, ipt, itime, icoeff, ierr


        !
        ! Get donor points
        !
        call get_donors(data,donors,donor_coords)


        !
        ! Compute Fourier coefficient modes
        !
        associate ( ntimes   => size(data%time_manager%times),  &
                    nfreq    => size(data%time_manager%freqs),  &
                    time_lev => data%time_manager%times,        &
                    freq     => data%time_manager%freqs)
          
            if (allocated(E)) deallocate(E)
            allocate(E(ntimes,ntimes), stat=ierr)
            if (ierr /= 0) call AllocationError

            call calc_E(nfreq,ntimes,freq,time_lev,E)

            call get_Fourier_coeff_vector(data,E,q_coeff_vector)
        
        end associate            


        !
        ! Compute no. of points and no. of time levels and allocate solution values array
        !
        npts  = size(donors)
        ntime = q_coeff_vector%get_ntime()

        if (allocated(solution_values)) deallocate(solution_values)
        allocate(solution_values(ntime,npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        
        !
        ! Compute Fourier coefficients at donor points
        ! 
        do itime = 1,ntime
            do ipt = 1,npts

                associate( idom  => donors(ipt)%idomain_g,  &
                           ielem => donors(ipt)%ielement_g, &
                           xi    => donor_coords(ipt)%c1_,  &
                           eta   => donor_coords(ipt)%c2_,  &
                           zeta  => donor_coords(ipt)%c3_)

                    val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(q_coeff_vector%dom(idom)%&
                                   vecs(ielem),ieqn,itime,xi,eta,zeta),rdouble)
                    solution_values(itime,ipt) = val(1)

                end associate

            end do  ! ipt
        end do  ! itime


    end subroutine get_Fourier_coeffs_at_donor_points
    !***********************************************************************************



    !>  Get post processed solution values at donor points
    !!
    !!  @author Mayank Sharma
    !!  @date   3/24/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_solution_at_donor_points(data,ieqn,solution_values)
        type(chidg_data_t),         intent(inout)   :: data
        integer(ik),                intent(in)      :: ieqn
        real(rk),   allocatable,    intent(inout)   :: solution_values(:,:)

        type(element_info_t),   allocatable     :: donors(:)
        type(point_t),          allocatable     :: donor_coords(:)
        real(rdouble)                           :: val(1)                     
        integer(ik)                             :: ntime, npts, ipt, itime, ierr


        !
        ! Get donor points
        !
        call get_donors(data,donors,donor_coords)

        
        !
        ! Compute no. of points and no. of time levels and allocate solution values array
        !
        npts  = size(donors)
        ntime = data%sdata%q_out%get_ntime()

        if (allocated(solution_values)) deallocate(solution_values)
        allocate(solution_values(ntime,npts), stat=ierr)
        if (ierr /= 0) call AllocationError

        
        !
        ! Compute solutions at donor points
        ! 
        do itime = 1,ntime
            do ipt = 1,npts

                associate( idom  => donors(ipt)%idomain_g,  &
                           ielem => donors(ipt)%ielement_g, &
                           xi    => donor_coords(ipt)%c1_,  &
                           eta   => donor_coords(ipt)%c2_,  &
                           zeta  => donor_coords(ipt)%c3_)

                    val = real(data%mesh%domain(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%&
                                   vecs(ielem),ieqn,itime,xi,eta,zeta),rdouble)
                    solution_values(itime,ipt) = val(1)

                end associate

            end do  ! ipt
        end do  ! itime


    end subroutine get_solution_at_donor_points
    !***********************************************************************************



    !>  Write matplotlib .dat file
    !!
    !!  @author Mayank Sharma
    !!  @date   3/24/2017
    !!
    !!  @param[in]  ieqn     - equation variable which needs to be plotted (input parameter 
    !!                         because plot wrt only one variable might be needed)
    !!  @param[in]  coord    - coordinate wrt to which plot is needed
    !!  @param[in]  filename - output file name
    !!
    !-----------------------------------------------------------------------------------
    subroutine write_matplotlib_file(data,ieqn,coord,filename_1,filename_2,filename_3)
        type(chidg_data_t),             intent(inout)           :: data
        integer(ik),                    intent(in)              :: ieqn
        character(:),   allocatable,    intent(in)              :: coord
        character(:),   allocatable,    intent(in)              :: filename_1
        character(:),   allocatable,    intent(in), optional    :: filename_2
        character(:),   allocatable,    intent(in), optional    :: filename_3

        type(point_t),  allocatable     :: points(:)
        real(rk),       allocatable     :: solution_values(:,:), orig_sol_val(:,:), Fourier_coeffs(:,:)
        integer(ik)                     :: npts, ipt, funit = 10
        logical                         :: exist


        !
        ! Get solution points and interpolated solutions
        !
        call get_points(points)
        call get_solution_at_donor_points(data,ieqn,solution_values)
        if (present(filename_2)) call get_original_solution_at_donor_points(data,ieqn,orig_sol_val)


        !
        ! Get number of points
        ! 
        npts = size(solution_values,2)


        !
        ! Write output file for post processed data
        !
        inquire(file = trim(filename_1), exist=exist)
        if (exist) then
            open(funit, file = trim(filename_1), status = 'old', action = 'write')
        else
            open(funit, file = trim(filename_1), status = 'new', action = 'write')
        end if

            select case(coord)
                case('x')
                    do ipt = 1,npts

                        write(funit,*) points(ipt)%c1_, solution_values(:,ipt)

                    end do
                case('y')
                    do ipt = 1,npts

                        write(funit,*) points(ipt)%c2_, solution_values(:,ipt)

                    end do
                case('z')
                    do ipt = 1,npts

                        write(funit,*) points(ipt)%c3_, solution_values(:,ipt)

                    end do
            end select

        close(funit)


        !
        ! Write output file for original data
        !
        if ( present(filename_2) ) then

            call get_original_solution_at_donor_points(data,ieqn,orig_sol_val)

            inquire(file = trim(filename_2), exist=exist)
            if (exist) then
                open(funit, file = trim(filename_2), status = 'old', action = 'write')
            else
                open(funit, file = trim(filename_2), status = 'new', action = 'write')
            end if

                select case(coord)
                    case('x')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c1_, orig_sol_val(:,ipt)

                        end do
                    case('y')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c2_, orig_sol_val(:,ipt)

                        end do
                    case('z')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c3_, orig_sol_val(:,ipt)

                        end do
                end select

            close(funit)  
            
        end if
        

        !
        ! Write output file for Fourier coefficients
        !
        if ( present(filename_3) ) then
    
            call get_Fourier_coeffs_at_donor_points(data,ieqn,Fourier_coeffs)
            

            inquire(file = trim(filename_3), exist=exist)
            if (exist) then
                open(funit, file = trim(filename_3), status = 'old', action = 'write')
            else
                open(funit, file = trim(filename_3), status = 'new', action = 'write')
            end if

                select case(coord)
                    case('x')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c1_, Fourier_coeffs(:,ipt)

                        end do
                    case('y')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c2_, Fourier_coeffs(:,ipt)

                        end do
                    case('z')
                        do ipt = 1,npts

                            write(funit,*) points(ipt)%c3_, Fourier_coeffs(:,ipt)

                        end do
                end select

            close(funit)  
            
        end if
          

    end subroutine write_matplotlib_file
    !***********************************************************************************




















end module mod_matplotlib_io
