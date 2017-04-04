module mod_vtk_calc_func
#include <messenger.h>
    
! Module containing functions and subroutines for arranging a ChiDG instance data for a vtk file

    use mod_kinds,           only: rk,ik,rdouble
    use mod_constants,       only: ONE,HALF,TWO,OUTPUT_RES

    use type_element,        only: element_t
    use type_blockvector,    only: blockvector_t
    use type_solverdata,     only: solverdata_t
    use type_chidg_data,     only: chidg_data_t
    use type_chidg_vector,   only: chidg_vector_t

    implicit none



contains




    !>
    !!
    !!
    !! @author Mayank Sharma
    !! @date 8/30/2016
    !!
    !!
    !!
    ! Get equation variable names
    !-------------------------------------------------------------------------------------------
    subroutine get_cons_var(data,idom,noeq,cons_var)
        
        type(chidg_data_t),                           intent(inout)        :: data
        integer(ik),                                  intent(in)           :: idom
        integer(ik),                                  intent(in)           :: noeq
        character(len = 100),dimension(:),allocatable,intent(inout)        :: cons_var


        integer(ik)                                                        :: ieq,ierr
 
  
        ! 
        ! Check if the array storing the conservative variable names is allocated or not and deallocate accordingly
        !  
        if (allocated(cons_var)) deallocate(cons_var)
        allocate(cons_var(noeq), stat = ierr)    
        if (ierr /= 0) call AllocationError
        
        ieq = 1

        !
        ! Store conservative variable names in an array of strings (to be used when writing out point data for a given piece)
        !
        do while (ieq .le. noeq)
            !cons_var(ieq) = trim(data%eqnset(1)%prop%eqns(ieq)%name)
            cons_var(ieq) = trim(data%eqnset(1)%prop%get_primary_field_name(ieq))
            ieq  = ieq + 1
        end do

 

    end subroutine get_cons_var
    !*******************************************************************************************



    !> Get the number of points and number of cells for each block (piece)
    !!
    !!
    !! @author Mayank Sharma
    !! @date 10/28/2016
    !!
    !!
    !!
    ! 
    !-------------------------------------------------------------------------------------------
    subroutine get_piece_nums(data,nelem,num_pts,num_cells)

        type(chidg_data_t),intent(inout)        :: data
        integer(ik),       intent(in)           :: nelem
        integer(ik),       intent(inout)        :: num_pts, num_cells


        num_pts   = nelem*((OUTPUT_RES + 1)**3)

        num_cells = nelem*((OUTPUT_RES)**3)
        
    
    end subroutine get_piece_nums
    !*******************************************************************************************

    
    
    
    !>
    !!
    !!
    !! @author Mayank Sharma
    !! @date 8/30/2016
    !!
    !!
    !!
    ! Get point coordinates for a structured grid piece
    ! In an individual .vts file, these are the coordinates for each block which consists of only one piece
    !-------------------------------------------------------------------------------------------
    subroutine get_piece_coord(data,idom,nelem,num_pts,X,Y,Z)
        
        type(chidg_data_t),               intent(inout)       :: data
        integer(ik),                      intent(in)          :: idom
        integer(ik),                      intent(in)          :: nelem      ! No. of elements in a block
        integer(ik),                      intent(in)          :: num_pts
        real(rk),dimension(:),allocatable,intent(inout)       :: X,Y,Z


        integer(ik)                                           :: npts_xi, npts_eta, npts_zeta
        integer(ik)                                           :: ipt_xi, ipt_eta, ipt_zeta
        integer(ik)                                           :: xilim, etalim, zetalim
        integer(ik)                                           :: npts, icoord
        real(rdouble)                                         :: val(1), r_coord, theta_coord, z_coord
        real(rk)                                              :: xi, eta, zeta
        integer(ik)                                           :: ival, ielem, ierr


        npts = OUTPUT_RES + 1

 
        xilim   = npts
        etalim  = npts
        zetalim = npts

    

        ! 
        ! Check if the arrays storing the x,y,z coordinates are allocated or not and deallocate accordingly
        !
        if (allocated(X) .and. allocated(Y) .and. allocated(Z)) deallocate(X,Y,Z)
        allocate(X(num_pts),Y(num_pts),Z(num_pts), stat = ierr)
        if (ierr /= 0) call AllocationError

        ival = 0

        !
        ! For each coordinate, compute pointwise values and store in respective arrays
        !
         do icoord = 1,3

            !
            ! For each actual element, create a sub-sampling of elements to resolve solution variation
            !
            do ielem = 1,nelem
            
                !Write sampling for current element  
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO


                                    ! Get coordinate value at point
                                    if ( data%mesh(idom)%elems(ielem)%coordinate_system == 'Cylindrical' ) then
                                        r_coord     = real(data%mesh(idom)%elems(ielem)%grid_point(1,xi,eta,zeta),rdouble)
                                        theta_coord = real(data%mesh(idom)%elems(ielem)%grid_point(2,xi,eta,zeta),rdouble)
                                        z_coord     = real(data%mesh(idom)%elems(ielem)%grid_point(3,xi,eta,zeta),rdouble)

                                        if (icoord == 1) then
                                            val = r_coord*cos(theta_coord)
                                        else if (icoord == 2) then
                                            val = r_coord*sin(theta_coord)
                                        else if (icoord == 3) then
                                            val = z_coord
                                        end if

                                    else

                                        val = real(data%mesh(idom)%elems(ielem)%grid_point(icoord,xi,eta,zeta),rdouble)

                                    end if


                                    ival = ival + 1           ! Counter for the coordinate arrays
                                    select case(icoord)
                                        case(1)
                                            X(ival) = val(1)
                                        case(2)               
                                            Y(ival) = val(1)  ! Store X,Y,Z pointwise values in their respective arrays                              
                                        case(3)
                                            Z(ival) = val(1)
                                    end select

                                end do  ! ipt_xi
                        end do  ! ipt_eta
                end do  ! ipt_zeta

            end do  ! ielem

            ival = 0  ! Reset counter to zero before starting next icoord iteration

        end do  ! coords              



    end subroutine get_piece_coord
    !*******************************************************************************************



    !>
    !!
    !!
    !! @author Mayank Sharma
    !! @date 8/30/2016
    !!
    !!
    !!
    ! Get point data for conservative variables at all points for a structured grid piece
    ! In an individual .vts file, these are the conservative variable values for a block, which consists of only piece
    !-------------------------------------------------------------------------------------------
    subroutine get_piece_data(data,idom,itime,nelem,num_pts,noeq,cons_var_val)

        type(chidg_data_t),                 intent(inout)       :: data
        integer(ik),                        intent(in)          :: idom
        integer(ik),                        intent(in)          :: itime
        integer(ik),                        intent(in)          :: nelem
        integer(ik),                        intent(in)          :: num_pts
        integer(ik),                        intent(in)          :: noeq
        real(rk),dimension(:,:),allocatable,intent(inout)       :: cons_var_val


        integer(ik)                                             :: npts_xi, npts_eta, npts_zeta
        integer(ik)                                             :: ipt_xi, ipt_eta, ipt_zeta
        integer(ik)                                             :: xilim, etalim, zetalim
        integer(ik)                                             :: npts
        real(rdouble)                                           :: val(1)
        real(rk)                                                :: xi, eta, zeta
        integer(ik)                                             :: ieq, ivar, ival, ielem, ierr, eqn_ID



        npts = OUTPUT_RES + 1


        xilim   = npts
        etalim  = npts
        zetalim = npts


        ! 
        ! Check if the array storing the conservative variable values is allocated or not and deallocate accordingly
        !
        if (allocated(cons_var_val)) deallocate(cons_var_val)
        allocate(cons_var_val(noeq,num_pts), stat = ierr)
        if (ierr /= 0) call AllocationError

        ival = 0


        !
        ! For each conservative variable in equation set, compute values pointwise and save in the conservative variable array
        !
        !do ivar = 1,data%eqnset(idom)%prop%nprimary_fields()
        eqn_ID = data%mesh(idom)%eqn_ID
        do ivar = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

            !
            ! For each actual element, create a sub-sampling of elements to resolve solution variation
            !
            do ielem = 1, nelem

                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO
                    do ipt_eta = 1,etalim
                        eta = (((real(ipt_eta,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO
                        do ipt_xi = 1,xilim
                            xi = (((real(ipt_xi,rk) - ONE)/(real(npts,rk) - ONE)) - HALF)*TWO

                            ! Get solution value at a point 
                            val = real(data%mesh(idom)%elems(ielem)%solution_point(data%sdata%q_out%dom(idom)%vecs(ielem),ivar,itime,xi,eta,zeta),rdouble)
                            ival = ival + 1                   ! Counter for conservative variable array
                            cons_var_val(ivar,ival) = val(1)  ! Store values in the array
                                                              ! Each row of the array contains the values for one conservative variable

                        end do  ! ipt_xi
                    end do  ! ipt_eta
                end do  ! ipt_zeta

            end do  ! ielem

            ival = 0  ! Reset array counter before starting another ivar iteration

        end do  ! ivar


    end subroutine get_piece_data
    !*******************************************************************************************



    !> Get connectivity, offsets and element type for each element of a block
    !!
    !!
    !! @author Mayank Sharma
    !! @date 10/28/2016
    !!
    !!
    !!
    !  
    !-------------------------------------------------------------------------------------------
    subroutine get_piece_connectivity_data(data,idom,nelem,num_pts,num_cells,connectivity,offsets,types)

        type(chidg_data_t),                    intent(inout)       :: data
        integer(ik),                           intent(in)          :: idom
        integer(ik),                           intent(in)          :: nelem,num_pts,num_cells
        integer(ik),dimension(:,:),allocatable,intent(inout)       :: connectivity
        integer(ik),dimension(:),allocatable  ,intent(inout)       :: offsets,types

        integer(ik),dimension(:,:),allocatable                     :: connectivity_temp     ! Temporary connectivity array which will transposed to
                                                                                            ! get final connectivity values

        integer(ik)                                                :: nelem_xi, nelem_eta, nelem_zeta
        integer(ik)                                                :: ielem_xi, ielem_eta, ielem_zeta
        integer(ik)                                                :: ielem_global, ielem_offset, istart, ielem_start, ierr
        integer(ik)                                                :: npts_element, nsub_per_element, nsub_elements
        integer(ik)                                                :: icell_var,ielem, iconn

        logical                                                    :: exist
       
        
        npts_element     = (OUTPUT_RES + 1)**3
        nsub_per_element = (OUTPUT_RES)**3
        nsub_elements    = nsub_per_element*nelem

        nelem_xi   = OUTPUT_RES
        nelem_eta  = OUTPUT_RES
        nelem_zeta = OUTPUT_RES

        
        !
        ! Check if the temporary array is allocated or not
        ! If it is allocated, deallocate it and allocate it again
        !
        if (allocated(connectivity_temp)) deallocate(connectivity_temp)
        allocate(connectivity_temp(8,num_cells), stat = ierr)
        if (ierr /= 0) call AllocationError

        do ielem = 1,nelem
            
            !
            ! Within each element, sample into sub-elements to account for variation within element
            !
            do ielem_zeta = 1, nelem_zeta

                do ielem_eta = 1, nelem_eta

                    do ielem_xi = 1, nelem_xi

                        ielem_start  = ielem_xi + (ielem_eta - 1)*(nelem_xi) + (ielem_zeta - 1)*(nelem_xi*nelem_eta)
                        ielem_global = ielem_start + (ielem - 1)*(nelem_zeta*nelem_eta*nelem_xi)
                        ielem_offset = (npts_element)*(ielem - 1)

                        istart = ielem_xi + (ielem_eta  - 1)*(nelem_xi + 1) + (ielem_zeta - 1)*((nelem_xi + 1)*(nelem_eta + 1))

                        connectivity_temp(1,ielem_global) = istart                                              + ielem_offset
                        connectivity_temp(2,ielem_global) = istart + 1                                          + ielem_offset
                        connectivity_temp(3,ielem_global) = istart + (OUTPUT_RES + 1) + 1                       + ielem_offset
                        connectivity_temp(4,ielem_global) = istart + (OUTPUT_RES + 1)                           + ielem_offset
                        connectivity_temp(5,ielem_global) = istart                        + (OUTPUT_RES + 1)**2 + ielem_offset
                        connectivity_temp(6,ielem_global) = istart + 1                    + (OUTPUT_RES + 1)**2 + ielem_offset
                        connectivity_temp(7,ielem_global) = istart + (OUTPUT_RES + 1) + 1 + (OUTPUT_RES + 1)**2 + ielem_offset
                        connectivity_temp(8,ielem_global) = istart + (OUTPUT_RES + 1)     + (OUTPUT_RES + 1)**2 + ielem_offset

                    end do  ! ielem_xi

                end do  ! ielem_eta

            end do  ! ielem_zeta

        end do  ! ielem


        !
        ! Check if the connectivity array is allocaated or not
        ! If it is allocated, deallocate and allocate it again
        !

        if (allocated(connectivity)) deallocate(connectivity)
        allocate(connectivity(num_cells,8), stat = ierr)
        if (ierr /= 0) call AllocationError

        !
        ! Transposed to fit vtk connectivity output requirements
        !
        connectivity = transpose(connectivity_temp) - 1.0

        !
        ! Check if offsets and types arrays are allocated or not
        ! If they are allocated, deallocate and allocate them again
        !

        if (allocated(offsets) .and. allocated(types)) deallocate(offsets,types)
        allocate(offsets(num_cells),types(num_cells), stat = ierr)
        if (ierr /= 0) call AllocationError

        !
        ! Offsets and types fields required in the vtk format
        !
        do icell_var = 1,num_cells

            offsets(icell_var) = 8*icell_var  ! 8 vertices in each sub-element
 
            types(icell_var)   = 12           ! VTK_HEXAHEDRON

        end do


    end subroutine get_piece_connectivity_data
    !*******************************************************************************************





















end module mod_vtk_calc_func
