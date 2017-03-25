module mod_matplotlib_io
#include<messenger.h>
    use mod_kinds,          only: rk,ik,rdouble
    use mod_constants,      only: ZERO,ONE,HALF,TWO,OUTPUT_RES
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_blockvector,   only: blockvector_t   
    use type_solverdata,    only: solverdata_t
    use type_chidg_data,    only: chidg_data_t
    use type_chidg_vector,  only: chidg_vector_t
    use mod_chimera,        only: find_gq_donor
    use type_face_info,     only: face_info_t
    use type_element_info,  only: element_info_t



contains



    !>  Get points along the domain for matplotlib computations
    !!
    !!  @author Mayank Sharma   
    !!  @date   3/24/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine get_points(points)
        type(point_t),  allocatable,    intent(inout)     :: points(:)

        integer(ik),    parameter       :: npts = 101
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

            ! Set receiver face indices to zero
            call receivers(ipt)%init(0_ik,0_ik,0_ik,0_ik,0_ik)
            call find_gq_donor(data%mesh,gq_nodes(ipt),receivers(ipt),donors(ipt),donor_coords(ipt))

        end do


    end subroutine get_donors
    !***********************************************************************************



    !>  Get solution values at donor points
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

                    val = real(data%mesh(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%&
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
    subroutine write_matplotlib_file(data,ieqn,coord,filename)
        type(chidg_data_t),             intent(inout)   :: data
        integer(ik),                    intent(in)      :: ieqn
        character(:),   allocatable,    intent(in)      :: coord
        character(:),   allocatable,    intent(in)      :: filename

        type(point_t),  allocatable     :: points(:)
        real(rk),       allocatable     :: solution_values(:,:)
        integer(ik)                     :: npts, ipt, funit = 10
        logical                         :: exist


        !
        ! Get solution points and interpolated solutions
        !
        call get_points(points)
        call get_solution_at_donor_points(data,ieqn,solution_values)


        !
        ! Get number of points
        ! 
        npts = size(solution_values,2)


        !
        ! Write output file
        !
        inquire(file = trim(filename), exist=exist)
        if (exist) then
            open(funit, file = trim(filename), status = 'old', action = 'write')
        else
            open(funit, file = trim(filename), status = 'new', action = 'write')
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


    end subroutine write_matplotlib_file
    !***********************************************************************************




















end module mod_matplotlib_io
