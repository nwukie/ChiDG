module type_chidg_adjoint
#include<messenger.h>
    use type_chidg_vector,      only: chidg_vector_t, chidg_vector
    use type_chidg_matrix,      only: chidg_matrix_t, chidg_matrix
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: ZERO, dD_DIFF 
    use mod_io,                 only: backend
    use type_mesh,              only: mesh_t
    use type_svector,           only: svector_t
    use type_storage_flags,     only: storage_flags_t
    use mod_string

    implicit none


    !>  Storage for adjoint variables
    !!
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/14/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    type,   public  :: chidg_adjoint_t

        ! PRIMARY ADJOINT VARIABLES
        ! Adjoitn variables of primary fields
        type(chidg_vector_t),   allocatable     :: v(:,:)           ! Primary adjoint variables (nfunc,nstep)
        type(chidg_vector_t),   allocatable     :: v_in(:)          ! Primary adjoint variables input from hdf(nfunc)
        type(chidg_vector_t),   allocatable     :: q_time(:)        ! primal solutions vectors (istep)
        type(chidg_vector_t),   allocatable     :: Jq(:)            ! Functional derivatives (nfunc)

        ! AUXILIARY ADJOINT VARIABLES
        ! Adjoint variables of auxiliary fields: wall_distance
        type(chidg_matrix_t),   allocatable     :: Rd(:)            ! Flow Jacobian distance field linearization (naux) 
        type(chidg_vector_t),   allocatable     :: vRd(:,:)         ! Dot product of primary adjoint variables
                                                                    ! and auxiliary Jacobian (naux,nfunc)
        ! Storage for solver report info (matrix_time, niter)
        integer(ik),            allocatable     :: solver_iter(:,:) ! Linear solver iterations (ifunc,istep) 
        real(rk),               allocatable     :: solver_time(:,:) ! Linear solver time (ifunc,istep)
        real(rk),               allocatable     :: solver_err(:,:)  ! Linear solver error (ifunc,istep)

        ! Adjoint initialization completed
        logical         :: primary_adjoint_initialized = .false.
        logical         :: auxiliary_adjoint_initialized = .false.

    contains

        procedure       :: init
        procedure       :: init_vector
!        procedure       :: check_adjoint_stored
        procedure       :: process_adjoint_solution
        procedure       :: process_primal_solution
        procedure       :: store_solver_info
        procedure       :: release
        procedure       :: report

    end type chidg_adjoint_t
    !************************************************************************************


contains


    !>  Allocate storage for adjoint variables and solution vector at nstep
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/14/2017
    !!
    !!  param[in]   nstep       Number of steps of the solution (nstep = 1 for steady and HB)
    !!  param[in]   nfunc       Number of functionals
    !!
    !------------------------------------------------------------------------------------
    subroutine init(self,nfunc,nstep,sflags)
        class(chidg_adjoint_t), intent(inout)   :: self
        integer(ik),            intent(in)      :: nstep
        integer(ik),            intent(in)      :: nfunc
        type(storage_flags_t),  intent(in)      :: sflags

        integer(ik)     :: ierr

        ! Deallocate any data previously allocated
        call self%release()

        ! Default, no error
        ierr = 0

        if (sflags%v) allocate( self%v(nfunc,nstep), stat=ierr )
        if (ierr/=0) call AllocationError

        if (sflags%q_time) allocate( self%q_time(nstep), stat=ierr )
        if (ierr/=0) call AllocationError

        if (sflags%Jq) allocate( self%Jq(nfunc), stat=ierr )
        if (ierr/=0) call AllocationError

        if (sflags%solver_iter) allocate( self%solver_iter(nfunc,nstep), stat=ierr )
        if (ierr/=0) call AllocationError

        if (sflags%solver_time) allocate( self%solver_time(nfunc,nstep), stat=ierr )
        if (ierr/=0) call AllocationError

        if (sflags%solver_err) allocate( self%solver_err(nfunc,nstep), stat=ierr )
        if (ierr/=0) call AllocationError
        
        
        self%primary_adjoint_initialized = .true.
    
        
        ! Auxiliary adjoint storage
        !   assumption: only one auxiliary field (wall_distance)
        if (sflags%vRd) then
            allocate( self%vRd(1,nfunc), stat=ierr )
            if (ierr/=0) call AllocationError
            self%auxiliary_adjoint_initialized = .true.
        end if
        
        if (sflags%Rd) then
            allocate( self%Rd(1), stat=ierr )
            if (ierr/=0) call AllocationError
        end if

    end subroutine init
    !************************************************************************************






    !>  Allocate storage necessary for adjoint chidg_vectors
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/14/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine init_vector(self,mesh,ntime,sflags)
        class(chidg_adjoint_t),     intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(inout)   :: ntime  
        type(storage_flags_t),      intent(in)      :: sflags
             
        integer(ik)     :: ifunc,istep,iaux !specialization

        if (sflags%v) then
            do ifunc = 1,size(self%v,1)
                do istep = 1,size(self%v,2)
                    self%v(ifunc,istep) = chidg_vector(trim(backend)) 
                    call self%v(ifunc,istep)%init(mesh,ntime,'primal')
                end do !istep
            end do !ifunc
        end if
        
        if (sflags%q_time) then
            do istep = 1,size(self%q_time,1)
                self%q_time(istep) = chidg_vector(trim(backend)) 
                call self%q_time(istep)%init(mesh,ntime,'primal')
            end do !istep
        end if

        if (sflags%Jq) then
            do ifunc = 1,size(self%Jq,1)
                self%Jq(ifunc) = chidg_vector(trim(backend)) 
                call self%Jq(ifunc)%init(mesh,ntime,'primal')
            end do !ifunc
        end if
 
        ! Allocate storage for auxiliary
!        specialization = dD_DIFF
        
        if (sflags%vRd) then
            do iaux = 1,size(self%vRd,1)
                do ifunc = 1,size(self%vRd,2)
                     !call self%vRd(iaux,ifunc)%init(mesh,ntime,specialization)
                    self%vRd(iaux,ifunc) = chidg_vector(trim(backend)) 
                    call self%vRd(iaux,ifunc)%init(mesh,ntime,'auxiliary')
                end do !istep
            end do !ifunc
        end if
    
        if (sflags%Rd) then
            do iaux = 1,size(self%Rd)
                self%Rd(iaux) = chidg_matrix(trim(backend)) 
                call self%Rd(iaux)%init(mesh,storage_config='dD',dof_type='auxiliary')
                call self%Rd(iaux)%init_recv(self%vRd(1,1))
                if (sflags%Rd_trans) self%Rd(iaux)%transposed = .true.
            end do
        end if


        ! Set Rd as transpose matrix
    


    end subroutine init_vector
    !************************************************************************************





!    !>  Check if at least one chidg_vector of adjoint variables 
!    !!  Return .true. if at least one adjoint variable vector is stored
!    !!
!    !!  @author Matteo Ugolotti
!    !!  @date   7/14/2017
!    !!
!    !!
!    !------------------------------------------------------------------------------------
!    function check_adjoint_stored(self) result(stat)
!        class(chidg_adjoint_t),     intent(in)   :: self
!        
!        logical     :: stat,not_zero
!        integer(ik) :: idom,ielem,iterm
!
!        stat = .false.
!
!        if ( self%primary_adjoint_initialized ) then
!
!            ! check if v(1,1) is ZERO, if the chidg_vector is ZERO then nothing 
!            ! has been written in any of the chidg_vectors belonging to v
!            do idom = 1, size(self%v(1,1)%dom)
!                do ielem = 1, size(self%v(1,1)%dom(idom)%vecs)
!                   do iterm = 1,size(self%v(1,1)%dom(idom)%vecs(ielem)%vec)
!                        
!                        not_zero = ( self%v(1,1)%dom(idom)%vecs(ielem)%vec(iterm) /= ZERO ) 
!                        if (not_zero) then
!                            stat = .true.
!                            exit
!                        end if
!
!                    end do !iterm 
!                end do !ielem
!            end do !idom
!
!        end if
!
!
!    end function check_adjoint_stored 
!    !************************************************************************************







    !>  This procedure copies an input primal solution stored in sdata%q_in (this is done
    !!  by read_primary_fields_hdf) into adjoint%q_time
    !!
    !!  For unsteady mesh-senstivities, we will read in each primal solution at the time (istep)
    !!  and the q_time vector matrix will be already initialized based on the number of
    !!  of time steps (nsteps) defined in the namelist. This means we need to pass in the istep
    !!  index so that we can locate the istep-th solution read in in the right spot in
    !!  q_time vector.
    !!
    !!  NOTE1: This subroutine is executed if we are actaully reading in primal solutions
    !!         for adjoint and mesh-sensitivites calculations, where we do need solution at
    !!         every time step.
    !!  NOTE2: q_in is re/initialized in read_primary_fields_hdf
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/9/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine process_primal_solution(self,q_in,istep)
        class(chidg_adjoint_t),     intent(inout)   :: self
        type(chidg_vector_t),       intent(in)      :: q_in
        integer(ik),                intent(in)      :: istep
        
        logical :: storage_available

        ! Check if there is room to store the new input primal solution
        storage_available = ( istep <= size(self%q_time) )
        if (.not. storage_available) then
            call chidg_signal(FATAL,'adjoint%process_primal_solution: the last primal    &
                                     solution read in cannot be stored in q_time vector. &
                                     Apparently, the q_time vector has not been properly &
                                     allocated. Check if the number of solutions read in &
                                     is euqal to the time_steps input in the chidg.nml')

        end if

        ! Copy q_in over to q_time(istep)
        self%q_time(istep) = q_in
        call self%q_time(istep)%assemble()
        
    end subroutine process_primal_solution
    !************************************************************************************







    !>  This procedure copies an input adjoint solution stored in v-in (this is done
    !!  by read_adjoint_fields_hdf) and copies it into v.
    !!
    !!  For post-process purposes, not input argument is necessary since only one step is written
    !!  out. This means that the v vector will be of size (nfuncs,1). In this case
    !!  of post-processing chidg%init('adjoint_storage') is NOT called, therefore we need to 
    !!  initialize the containers now.
    !!
    !!  For unsteady mesh-senstivities, we will read each adjoint solution at the time (istep)
    !!  and the v matrix will be already initialized. This means we need to pass in the istep
    !!  index so that we can locate the istep-th solution read in in the right spot in
    !!  v matrix.
    !!
    !!  NOTE1: This subroutine is executed if we are actaully reading in an adjoint solution.
    !!         This routine is called every time we are reading a solution file, however
    !!         we need to consider that the solution file might not contain an adjoint solution.
    !!  NOTE2: v_in is allocated in read_adjoint_fields_hdf based on the number of functionals
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/9/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine process_adjoint_solution(self,istep)
        class(chidg_adjoint_t),     intent(inout)   :: self
        integer(ik),    optional,   intent(in)      :: istep
        
        integer(ik)     :: ifunc,ierr
        logical         :: process_solution, nfuncs_mismatch


        ! Check if we have read an adjoint solution in
        process_solution = .false.
        if (allocated(self%v_in)) process_solution = .true.

        ! If so, process it
        if (process_solution) then
            
            if (present(istep)) then
                
                ! Check that v and v_in have the same number of functionals
                nfuncs_mismatch = ( size(self%v,1) /= size(self%v_in) )
                if (nfuncs_mismatch) then
                    call chidg_signal(FATAL,'adjoint%process_adjoint_solution: the number of &
                                    functionals in the input file provided does not match    &
                                    the number of functionals used for initializing the      &
                                    the chidg containers.                                    &
                                    Check that the first solution file read in by AdjointX   &
                                    contains the right number of functionals.')
                end if
                         
                ! For reading in solutions for mesh-sensitivity calculation
                do ifunc = 1,size(self%v,1)
                    self%v(ifunc,istep) = self%v_in(ifunc)
                end do

            else 

                ! For reading in solution for post-process (istep = 1)
                if (allocated(self%v)) deallocate (self%v)
                allocate( self%v(size(self%v_in,1),1), stat=ierr )
                if (ierr/=0) call AllocationError
                
                do ifunc = 1,size(self%v,1)
                    self%v(ifunc,1) = self%v_in(ifunc)
                end do

            end if

        end if !process-solution

    end subroutine process_adjoint_solution
    !************************************************************************************




    
    
    !>  Store linear solver info. These will be use for the final report
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/30/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine store_solver_info(self,ifunc,istep,niter,time,err)
        class(chidg_adjoint_t),     intent(inout)   :: self
        integer(ik),                intent(in)      :: ifunc
        integer(ik),                intent(in)      :: istep
        integer(ik),                intent(in)      :: niter
        real(rk),                   intent(in)      :: time
        real(rk),                   intent(in)      :: err
        
        self%solver_iter(ifunc,istep) = niter
        self%solver_time(ifunc,istep) = time
        self%solver_err( ifunc,istep) = err

    end subroutine store_solver_info
    !************************************************************************************






    !>  Store linear solver info. These will be use for the final report
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/30/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine report(self,funcs_name)
        class(chidg_adjoint_t),     intent(inout)   :: self
        type(svector_t),            intent(in)      :: funcs_name


        integer(ik)          :: istep, ifunc, matrix_iter
        real(rk)             :: functional_time, overall_time, matrix_time, lin_sol_err 
        
        
        ! Accumulate total time for adjoint computation of all the functionals
        overall_time = 0._rk


        ! Adjoint linear solver header
        call write_line(' ')
        call write_line('---------------------------------  Adjoint Linear Solver Report  ----------------------------------')

        
        do ifunc = 1,size(self%v,1)
           
            ! write functional name
            call write_line(' ')
            call write_line('Functional:                ', string_to_upper(funcs_name%data_(ifunc)%get()), columns=.True., column_width=20)
            call write_line(' ')
            
            ! Reset functional computational time
            functional_time = 0._rk

            ! Print per-iteration report
            call write_line('istep', 'Linear residuals', 'Matrix time', 'Matrix iterations', columns=.True., column_width=20)


            ! Loop through stored data and print for each linear iteration the data for each istep 
            do istep = 1,size(self%v,2)

                matrix_iter = self%solver_iter(ifunc,istep)
                matrix_time = self%solver_time(ifunc,istep)
                lin_sol_err = self%solver_err( ifunc,istep)
                
                call write_line(istep, lin_sol_err, matrix_time, matrix_iter, delimiter=', ', columns=.True., column_width=20)
                
                ! Accumulate functional time, for each unstedy step
                functional_time = functional_time + matrix_time
            
            end do
        
            ! Accumulate over all time
            overall_time = overall_time + functional_time 
        
            !call write_line(' ')
            !call write_line('Total computational time:  ', functional_time, columns=.True., column_width=20)
            !call write_line(' ')

        end do !ifunc

        call write_line(' ')
        call write_line(' ')
        call write_line('Overall adjoint run time: ', overall_time, columns=.True., column_width=30)
        call write_line('-------------------------------------------------------------------------------------------------')

    end subroutine report
    !************************************************************************************







    !>  Release allocated data 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   9/13/2017
    !!
    !------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_adjoint_t),     intent(inout)   :: self

        if (allocated(self%v))      deallocate(self%v)
        if (allocated(self%v_in))   deallocate(self%v_in)
        if (allocated(self%q_time)) deallocate(self%q_time)
        if (allocated(self%Jq))     deallocate(self%Jq)
        if (allocated(self%Rd))     deallocate(self%Rd)
        if (allocated(self%vRd))    deallocate(self%vRd)

        if (allocated(self%solver_iter)) deallocate(self%solver_iter)
        if (allocated(self%solver_time)) deallocate(self%solver_time)
        if (allocated(self%solver_err))  deallocate(self%solver_err)

    end subroutine release
    !************************************************************************************




end module type_chidg_adjoint 
