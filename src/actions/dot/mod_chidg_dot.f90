!>  ChiDG DOT: this action allow to compute the final derivatives of an objective function
!!  wrt to a design parameter. The grid-node senstivities for a given objective function are
!!  multiplied (dot product) with the mesh sensitivities (dX/dY), which are computed by means 
!!  of Finite-Difference.
!! 
!!
!!  @author Matteo Ugolotti
!!  @date 12/5/2018
!!
!!
!!  Usage: chidg dot
!!  
!
!----------------------------------------------------------------------------------------------
module mod_chidg_dot
#include <messenger.h>

    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: ZERO, TWO, IO_DESTINATION
    use type_chidg,                 only: chidg_t
    use mod_chidg_mpi,              only: GLOBAL_MASTER, ChiDG_COMM, IRANK, NRANK
    use mod_string,                 only: string_t
    use type_file_properties,       only: file_properties_t
    use mod_file_utilities,         only: delete_file


    implicit none





contains



    !>  Driver for computing the sensitivities of an objective function wrt design parameters
    !!
    !!  @author Matteo Ugolotti
    !!  @date 10/29/2019
    !!
    !!  Required computation
    !!
    !!      DJ/DY = dot(dJ/dX,dX/dY) 
    !!
    !!      where "Y" is the design parameter 
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_dot(func_sens,mesh_sens)
        character(*),   intent(in)     :: func_sens
        character(*),   intent(in)     :: mesh_sens
         

        ! ChiDG environment
        type(chidg_t)                               :: chidg


        integer(ik)                                 :: fileunit1, fileunit2,      &
                                                       nblocks_1, nblocks_2,      &
                                                       igrid, ierr, npts, npt_i,  &
                                                       npt_j, npt_k, i, j, k,     &
                                                       unit_file
        logical                                     :: file1_exists, file2_exists
        real(rk)                                    :: dJdY
        integer(ik),                    allocatable :: blkdims_1(:,:), blkdims_2(:,:), blkdims_3(:,:)
        real(rk),                       allocatable :: msens_1(:,:,:), msens_2(:,:,:), msens_3(:,:,:), &
                                                       fsens_1(:,:,:), fsens_2(:,:,:), fsens_3(:,:,:), &
                                                       xsens(:,:,:), ysens(:,:,:), zsens(:,:,:)
        type(string_t)                              :: output_file_str
        character(len=18)                           :: FMT = '(ES20.14)'
        logical                                     :: stat
        integer(ik)                                 :: f_vars, m_vars, ivar
        character(50)                               :: ivar_str
        

        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('core')


        ! Write out input files
        call write_line("Functional sensitivities :", func_sens      ,delimiter=" ")
        call write_line("Mesh sensitivities       :", mesh_sens      ,delimiter=" ")

        ! Send output to screen
        IO_DESTINATION = 'screen'


        ! Check if input files exist.
        inquire(file=func_sens, exist=file1_exists)
        inquire(file=mesh_sens, exist=file2_exists)
        if (.not. file1_exists) call chidg_signal(FATAL,"Functional sensitivities file is missing. Check if the file exists or if has been missplelled.")
        if (.not. file2_exists) call chidg_signal(FATAL,"Mesh sensitivities file is missing. Check if the file exists or if has been missplelled.")


        ! Start reading files
        open(newunit=fileunit1, file=trim(func_sens), form='unformatted')
        open(newunit=fileunit2, file=trim(mesh_sens), form='unformatted')


        ! Read number of grid blocks
        call write_line("Reading number of block...")
        read(fileunit1) nblocks_1
        read(fileunit2) nblocks_2
        
        if ( (nblocks_1 /= nblocks_2) ) call chidg_signal(FATAL,"Apparently there is mismatch in the number of grid blocks between the file provided")


        ! Read block dimensions
        call write_line("Reading block dimensions...")
        allocate(blkdims_1(4,nblocks_1),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit1) (blkdims_1(1,igrid), blkdims_1(2,igrid), blkdims_1(3,igrid), blkdims_1(4,igrid), igrid=1,nblocks_1)

        
        allocate(blkdims_2(4,nblocks_2),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit2) (blkdims_2(1,igrid), blkdims_2(2,igrid), blkdims_2(3,igrid), blkdims_2(4,igrid), igrid=1,nblocks_2)
        
        do igrid = 1,nblocks_1
            if ( (blkdims_1(1,igrid) /= blkdims_2(1,igrid)) .or. (blkdims_1(2,igrid) /= blkdims_2(2,igrid)) ) then
                call chidg_signal_one(FATAL,"Apparently, there is mismatch in the number of points in grid: ",igrid) 
            end if 
        end do 
        

        ! Check the number of variables in the func sensitivities file is 3
        f_vars = blkdims_1(4,1)
        if ( f_vars /= 3 ) call chidg_signal(FATAL,"The functional sensitivities file should contain 3 variables. Implementation error.")
        

        ! Check the number of variables in the mesh sensitivities file
        m_vars = blkdims_2(4,1)/3


        ! The mesh sensitivities can contain more than one design parameter.
        ! Each design parameter will have x, y, z sensitivities. So 3*nparams will be the number of mesh sensitivities
        ! 
        ! Here we need to close the functional sensitivities file to read it again and again for each design paramter
        close(fileunit1)


        ! Now we loop thorugh all the design variables
        do ivar = 1,m_vars 
            
            ! Convert ivar into string for io
            if (ivar .lt. 10) then 
                write(ivar_str, "(I1)") ivar
            else if ( (ivar .ge. 10) .and. (ivar .lt. 100) ) then
                write(ivar_str, "(I2)") ivar
            else
                write(ivar_str, "(I3)") ivar
            end if

            ! Reopen the functional sensitivities file and read headers.
            open(newunit=fileunit1, file=trim(func_sens), form='unformatted')
            read(fileunit1) nblocks_1
            if (allocated(blkdims_1)) deallocate(blkdims_1)
            allocate(blkdims_1(4,nblocks_1),stat=ierr)
            if (ierr /= 0) call AllocationError
            read(fileunit1) (blkdims_1(1,igrid), blkdims_1(2,igrid), blkdims_1(3,igrid), blkdims_1(4,igrid), igrid=1,nblocks_1)

            ! Reading node coordiantes
            call write_line("Computing sensitivities for design parameter ",trim(ivar_str))
            
            
            ! Initailize dJdY = 0
            dJdY = ZERO
            
            do igrid = 1,nblocks_1

                ! Dimensions for reading plot3d grid
                npt_i = blkdims_1(1,igrid)
                npt_j = blkdims_1(2,igrid)
                npt_k = blkdims_1(3,igrid)
                npts  = blkdims_1(1,igrid) * blkdims_1(2,igrid) * blkdims_1(3,igrid)


                if (allocated(fsens_1)) deallocate(fsens_1,fsens_2,fsens_3)
                allocate(fsens_1(npt_i,npt_j,npt_k),fsens_2(npt_i,npt_j,npt_k),fsens_3(npt_i,npt_j,npt_k), stat=ierr)
                if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

                if (allocated(msens_1)) deallocate(msens_1,msens_2,msens_3)
                allocate(msens_1(npt_i,npt_j,npt_k),msens_2(npt_i,npt_j,npt_k),msens_3(npt_i,npt_j,npt_k), stat=ierr)
                if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"
                 
                
                ! Read func sensitivities of the function file .q
                read(fileunit1) ((( fsens_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                                ((( fsens_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                                ((( fsens_3(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

                ! Read mesh sensitivities of the function file .q
                read(fileunit2) ((( msens_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                                ((( msens_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                                ((( msens_3(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

                ! Multiply each dX_i/dY for the correspondent dJ/dX_i
                xsens = fsens_1*msens_1
                ysens = fsens_2*msens_2
                zsens = fsens_3*msens_3

                ! Sum all derivatives contributions
                dJdY = dJdY + sum(xsens) + sum(ysens) + sum(zsens)

            end do
            

            ! Print total sensititivites to a file 
            call write_line("Total derivatives for design parameter ",trim(ivar_str),":",dJdY)
            
            
            ! Write total sensititivites to a file 
            
            ! Here we might have the mesh sens file that has more paramters
            ! Replace dJdX in the mesh_sens file with dJd1
            ! dJdX_drag.q -> dJd1_drag.txt
            call output_file_str%set(func_sens)
            output_file_str = output_file_str%replace('.q','.txt')
            output_file_str = output_file_str%replace('dX','d'//trim(ivar_str))
            

            ! Delete old file if exists
            inquire(file=output_file_str%get(), exist=stat)
            if (stat) call delete_file(output_file_str%get())

            ! Create and open a new one
            open(newunit=unit_file, file=trim(output_file_str%get()), status='new', action='write', iostat=ierr)
            if ( ierr /= 0 ) call chidg_signal(FATAL,"mod_chidg_dot: error opening target file.")


            ! Write data
            write(unit_file, *) dJdY

            ! Close files
            close(unit_file)
            close(fileunit1)

        end do !ivar

        ! Close mesh sensitivities file
        close(fileunit2)
        
        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')

    end subroutine chidg_dot
    !******************************************************************************************






    !>  Driver for computing the sensitivities of an objective function wrt design parameters
    !!
    !!  @author Matteo Ugolotti
    !!  @date 12/5/2018
    !!
    !!  Required computation
    !!
    !!      DJ/DY = dot(dJ/dX,dX/dY) 
    !!
    !!      where "Y" is the design parameter 
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_dot_fd(original_grid,perturbed_grid,mesh_sens,FD_delta)
        character(*),   intent(in)     :: original_grid
        character(*),   intent(in)     :: perturbed_grid
        character(*),   intent(in)     :: mesh_sens
        character(*),   intent(in)     :: FD_delta

        ! ChiDG environment
        type(chidg_t)                               :: chidg


        integer(ik)                                 :: fileunit1, fileunit2, fileunit3,         &
                                                       nblocks_1, nblocks_2, nblocks_3,         &
                                                       igrid, ierr, npts, npt_i, npt_j, npt_k,  &
                                                       unit_file, i, j ,k
        logical                                     :: file1_exists, file2_exists, file3_exists
        real(rk)                                    :: delta, dJdY
        integer(ik),                    allocatable :: blkdims_1(:,:), blkdims_2(:,:), blkdims_3(:,:)
        real(rk),                       allocatable :: xcoords_1(:,:,:), xcoords_2(:,:,:), ycoords_1(:,:,:), ycoords_2(:,:,:), &
                                                       zcoords_1(:,:,:), zcoords_2(:,:,:), xsens(:,:,:), ysens(:,:,:),     &
                                                       zsens(:,:,:), dxcoords(:,:,:), dycoords(:,:,:), dzcoords(:,:,:)
        type(string_t)                              :: output_file_str
        character(len=18)                           :: FMT = '(ES20.14)'
        logical                                     :: stat
        

        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('core')

        ! Write out input files
        call write_line("Original grid file :", original_grid  ,delimiter=" ")
        call write_line("Perturbed grid file:", perturbed_grid ,delimiter=" ")
        call write_line("Mesh sensitivities :", mesh_sens      ,delimiter=" ")
        call write_line("Delta perturbation :", FD_delta       ,delimiter=" ")

        ! Send output to screen
        IO_DESTINATION = 'screen'

        ! Read Finte Difference delta
        read(FD_delta,'(f10.0)') delta

        ! Check if input files exist.
        inquire(file=original_grid, exist=file1_exists)
        inquire(file=perturbed_grid,exist=file2_exists)
        inquire(file=mesh_sens,     exist=file3_exists)
        if (.not. file1_exists) call chidg_signal(FATAL,"Original grid file is missing. Check if the grid file exists or if has been missplelled.")
        if (.not. file2_exists) call chidg_signal(FATAL,"Perturbed grid file is missing. Check if the grid file exists or if has been missplelled.")
        if (.not. file3_exists) call chidg_signal(FATAL,"Mesh sensitivities file is missing. Check if the grid file exists or if has been missplelled.")


        ! Start reading files
        open(newunit=fileunit1, file=trim(original_grid),  form='unformatted')
        open(newunit=fileunit2, file=trim(perturbed_grid), form='unformatted')
        open(newunit=fileunit3, file=trim(mesh_sens),      form='unformatted')


        ! Read number of grid blocks
        call write_line("Reading number of block...")
        read(fileunit1) nblocks_1
        read(fileunit2) nblocks_2
        read(fileunit3) nblocks_3
        
        if ( (nblocks_1 /= nblocks_2) .or. (nblocks_1 /= nblocks_3) ) call chidg_signal(FATAL,"Apparently there is mismatch in the number of grid blocks between the file provided")


        ! Read block dimensions
        call write_line("Reading block dimensions...")
        allocate(blkdims_1(3,nblocks_1),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit1) (blkdims_1(1,igrid), blkdims_1(2,igrid), blkdims_1(3,igrid), igrid=1,nblocks_1)
        
        allocate(blkdims_2(3,nblocks_2),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit2) (blkdims_2(1,igrid), blkdims_2(2,igrid), blkdims_2(3,igrid), igrid=1,nblocks_2)
        
        allocate(blkdims_3(3,nblocks_3),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit3) (blkdims_3(1,igrid), blkdims_3(2,igrid), blkdims_3(3,igrid), igrid=1,nblocks_3)

        do igrid = 1,nblocks_1
            if ( (blkdims_1(1,igrid) /= blkdims_2(1,igrid)) .or. (blkdims_1(2,igrid) /= blkdims_2(2,igrid)) .or. (blkdims_1(3,igrid) /= blkdims_2(3,igrid)) ) then
                call chidg_signal_one(FATAL,"Apparently, there is mismatch in the number of points in grid: ",igrid) 
            end if 
        end do 
        
        
        ! Reading node coordiantes
        call write_line("Reading node coordiantes and sensitivities, and accumulating derivatives...")
        
        
        ! Initailize dJdY = 0
        dJdY = ZERO
        
        do igrid = 1,nblocks_1

            ! Dimensions for reading plot3d grid
            npt_i = blkdims_1(1,igrid)
            npt_j = blkdims_1(2,igrid)
            npt_k = blkdims_1(3,igrid)
            npts  = blkdims_1(1,igrid) * blkdims_1(2,igrid) * blkdims_1(3,igrid)


            ! Read block coordinates
            if (allocated(xcoords_1)) deallocate(xcoords_1,ycoords_1,zcoords_1)
            allocate(xcoords_1(npt_i,npt_j,npt_k),ycoords_1(npt_i,npt_j,npt_k),zcoords_1(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            if (allocated(xcoords_2)) deallocate(xcoords_2,ycoords_2,zcoords_2)
            allocate(xcoords_2(npt_i,npt_j,npt_k),ycoords_2(npt_i,npt_j,npt_k),zcoords_2(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            if (allocated(xsens)) deallocate(xsens,ysens,zsens)
            allocate(xsens(npt_i,npt_j,npt_k),ysens(npt_i,npt_j,npt_k),zsens(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"
            
            if (allocated(dxcoords)) deallocate(dxcoords,dycoords,dzcoords)
            allocate(dxcoords(npt_i,npt_j,npt_k),dycoords(npt_i,npt_j,npt_k),dzcoords(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            ! Read coordinates of the original/unperturbed grid
            read(fileunit1) ((( xcoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ycoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zcoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

            ! Read coordinates of the unperturbed grid
            read(fileunit2) ((( xcoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ycoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zcoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)
            
            ! Read sensitivities of the function file .q
            read(fileunit3) ((( xsens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ysens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zsens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

            ! Compute finite difference
            dxcoords = (xcoords_2 - xcoords_1)/delta
            dycoords = (ycoords_2 - ycoords_1)/delta
            dzcoords = (zcoords_2 - zcoords_1)/delta

            
            ! Multiply each dX_i/dY for the correspondent dJ/dX_i
            xsens = xsens*dxcoords
            ysens = ysens*dycoords
            zsens = zsens*dzcoords


            ! Sum all derivatives contributions
            dJdY = dJdY + sum(xsens) + sum(ysens) + sum(zsens)

        end do
        

        ! Print total sensititivites to a file 
        print*, "Total sensitivities: ", dJdY 
        
        
        ! Write total sensititivites to a file 
        
        ! Replace dJdX in the mesh_sens file with dJdY
        ! dJdX_drag.q -> dJdY_drag.txt
        call output_file_str%set(mesh_sens)
        output_file_str = output_file_str%replace('.q','.txt')
        output_file_str = output_file_str%replace('X','Y')

        ! Delete old file if exists
        inquire(file=output_file_str%get(), exist=stat)
        if (stat) call delete_file(output_file_str%get())

        ! Create and open a new one
        open(newunit=unit_file, file=trim(output_file_str%get()), status='new', action='write', iostat=ierr)
        if ( ierr /= 0 ) call chidg_signal(FATAL,"mod_chidg_dot: error opening target file.")


        ! Write data
        write(unit_file, *) dJdY


        ! Close files
        close(unit_file)
        close(fileunit1)
        close(fileunit2)
        close(fileunit3)


        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')


    end subroutine chidg_dot_fd
    !******************************************************************************************






    !>  Driver for computing the sensitivities of an objective function wrt design parameters
    !!
    !!  @author Matteo Ugolotti
    !!  @date 12/5/2018
    !!
    !!  Required computation
    !!
    !!      DJ/DY = dot(dJ/dX,dX/dY) 
    !!
    !!      where "Y" is the design parameter 
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_dot_cd(neg_perturbed_grid,pos_perturbed_grid,mesh_sens,FD_delta)
        character(*),   intent(in)     :: neg_perturbed_grid
        character(*),   intent(in)     :: pos_perturbed_grid
        character(*),   intent(in)     :: mesh_sens
        character(*),   intent(in)     :: FD_delta
         
        ! ChiDG environment
        type(chidg_t)                               :: chidg


        integer(ik)                                 :: fileunit1, fileunit2, fileunit3,         &
                                                       nblocks_1, nblocks_2, nblocks_3,         &
                                                       igrid, ierr, npts, npt_i, npt_j, npt_k,  &
                                                       unit_file, i, j ,k
        logical                                     :: file1_exists, file2_exists, file3_exists
        real(rk)                                    :: delta, dJdY
        integer(ik),                    allocatable :: blkdims_1(:,:), blkdims_2(:,:), blkdims_3(:,:)
        real(rk),                       allocatable :: xcoords_1(:,:,:), xcoords_2(:,:,:), ycoords_1(:,:,:), ycoords_2(:,:,:), &
                                                       zcoords_1(:,:,:), zcoords_2(:,:,:), xsens(:,:,:), ysens(:,:,:),     &
                                                       zsens(:,:,:), dxcoords(:,:,:), dycoords(:,:,:), dzcoords(:,:,:)
        type(string_t)                              :: output_file_str
        character(len=18)                           :: FMT = '(ES20.14)'
        logical                                     :: stat
        

        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('core')


        ! Write out input files
        call write_line("Negatively perturbed grid file :", neg_perturbed_grid ,delimiter=" ")
        call write_line("Positively perturbed grid file :", pos_perturbed_grid ,delimiter=" ")
        call write_line("Mesh sensitivities :"            , mesh_sens          ,delimiter=" ")
        call write_line("Delta perturbation :"            , FD_delta           ,delimiter=" ")

        ! Send output to screen
        IO_DESTINATION = 'screen'


        ! Read Finte Difference delta
        read(FD_delta,'(f10.0)') delta


        ! Check if input files exist.
        inquire(file=neg_perturbed_grid, exist=file1_exists)
        inquire(file=pos_perturbed_grid, exist=file2_exists)
        inquire(file=mesh_sens,          exist=file3_exists)
        if (.not. file1_exists) call chidg_signal(FATAL,"Negatively perturbed grid file is missing. Check if the grid file exists or if has been missplelled.")
        if (.not. file2_exists) call chidg_signal(FATAL,"Positively perturbed grid file is missing. Check if the grid file exists or if has been missplelled.")
        if (.not. file3_exists) call chidg_signal(FATAL,"Mesh sensitivities file is missing. Check if the grid file exists or if has been missplelled.")


        ! Start reading files
        open(newunit=fileunit1, file=trim(neg_perturbed_grid), form='unformatted')
        open(newunit=fileunit2, file=trim(pos_perturbed_grid), form='unformatted')
        open(newunit=fileunit3, file=trim(mesh_sens),          form='unformatted')


        ! Read number of grid blocks
        call write_line("Reading number of block...")
        read(fileunit1) nblocks_1
        read(fileunit2) nblocks_2
        read(fileunit3) nblocks_3
        
        if ( (nblocks_1 /= nblocks_2) .or. (nblocks_1 /= nblocks_3) ) call chidg_signal(FATAL,"Apparently there is mismatch in the number of grid blocks between the file provided")


        ! Read block dimensions
        call write_line("Reading block dimensions...")
        allocate(blkdims_1(3,nblocks_1),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit1) (blkdims_1(1,igrid), blkdims_1(2,igrid), blkdims_1(3,igrid), igrid=1,nblocks_1)
        
        allocate(blkdims_2(3,nblocks_2),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit2) (blkdims_2(1,igrid), blkdims_2(2,igrid), blkdims_2(3,igrid), igrid=1,nblocks_2)
        
        allocate(blkdims_3(3,nblocks_3),stat=ierr)
        if (ierr /= 0) call AllocationError
        read(fileunit3) (blkdims_3(1,igrid), blkdims_3(2,igrid), blkdims_3(3,igrid), igrid=1,nblocks_3)

        do igrid = 1,nblocks_1
            if ( (blkdims_1(1,igrid) /= blkdims_2(1,igrid)) .or. (blkdims_1(2,igrid) /= blkdims_2(2,igrid)) .or. (blkdims_1(3,igrid) /= blkdims_2(3,igrid)) ) then
                call chidg_signal_one(FATAL,"Apparently, there is mismatch in the number of points in grid: ",igrid) 
            end if 
        end do 
        
        
        ! Reading node coordiantes
        call write_line("Reading node coordiantes and sensitivities, and accumulating derivatives...")
        
        
        ! Initailize dJdY = 0
        dJdY = ZERO
        
        do igrid = 1,nblocks_1

            ! Dimensions for reading plot3d grid
            npt_i = blkdims_1(1,igrid)
            npt_j = blkdims_1(2,igrid)
            npt_k = blkdims_1(3,igrid)
            npts  = blkdims_1(1,igrid) * blkdims_1(2,igrid) * blkdims_1(3,igrid)


            ! Read block coordinates
            if (allocated(xcoords_1)) deallocate(xcoords_1,ycoords_1,zcoords_1)
            allocate(xcoords_1(npt_i,npt_j,npt_k),ycoords_1(npt_i,npt_j,npt_k),zcoords_1(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            if (allocated(xcoords_2)) deallocate(xcoords_2,ycoords_2,zcoords_2)
            allocate(xcoords_2(npt_i,npt_j,npt_k),ycoords_2(npt_i,npt_j,npt_k),zcoords_2(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            if (allocated(xsens)) deallocate(xsens,ysens,zsens)
            allocate(xsens(npt_i,npt_j,npt_k),ysens(npt_i,npt_j,npt_k),zsens(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"
            
            if (allocated(dxcoords)) deallocate(dxcoords,dycoords,dzcoords)
            allocate(dxcoords(npt_i,npt_j,npt_k),dycoords(npt_i,npt_j,npt_k),dzcoords(npt_i,npt_j,npt_k), stat=ierr)
            if (ierr /= 0) stop "memory allocation error: mod_chidg_dot"

            ! Read coordinates of the original/unperturbed grid
            read(fileunit1) ((( xcoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ycoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zcoords_1(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)

            ! Read coordinates of the unperturbed grid
            read(fileunit2) ((( xcoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ycoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zcoords_2(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)
            
            ! Read sensitivities of the function file .q
            read(fileunit3) ((( xsens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( ysens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k), &
                            ((( zsens(i,j,k), i=1,npt_i), j=1,npt_j), k=1,npt_k)


            ! Compute finite difference
            dxcoords = (xcoords_2 - xcoords_1)/(TWO*delta)
            dycoords = (ycoords_2 - ycoords_1)/(TWO*delta)
            dzcoords = (zcoords_2 - zcoords_1)/(TWO*delta)

            
            ! Multiply each dX_i/dY for the correspondent dJ/dX_i
            xsens = xsens*dxcoords
            ysens = ysens*dycoords
            zsens = zsens*dzcoords


            ! Sum all derivatives contributions
            dJdY = dJdY + sum(xsens) + sum(ysens) + sum(zsens)

        end do
        

        ! Print total sensititivites to a file 
        print*, "Total sensitivities: ", dJdY 
        
        
        ! Write total sensititivites to a file 
        
        ! Replace dJdX in the mesh_sens file with dJdY
        ! dJdX_drag.q -> dJdY_drag.txt
        call output_file_str%set(mesh_sens)
        output_file_str = output_file_str%replace('.q','.txt')
        output_file_str = output_file_str%replace('X','Y')

        ! Delete old file if exists
        inquire(file=output_file_str%get(), exist=stat)
        if (stat) call delete_file(output_file_str%get())

        ! Create and open a new one
        open(newunit=unit_file, file=trim(output_file_str%get()), status='new', action='write', iostat=ierr)
        if ( ierr /= 0 ) call chidg_signal(FATAL,"mod_chidg_dot: error opening target file.")


        ! Write data
        write(unit_file, *) dJdY


        ! Close files
        close(unit_file)
        close(fileunit1)
        close(fileunit2)
        close(fileunit3)


        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')

    end subroutine chidg_dot_cd
    !******************************************************************************************





end module mod_chidg_dot
