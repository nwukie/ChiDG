module type_chidg_adjointbc
#include<messenger.h>
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_chidg_vector,      only: chidg_vector_t
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: ZERO
    use type_mesh,              only: mesh_t
    use type_nvector,           only: nvector_t
    use type_node,              only: node_t
    use type_storage_flags,     only: storage_flags_t
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM, GLOBAL_MASTER
    use mpi_f08,                only: MPI_BCast, MPI_ISend, MPI_Recv, MPI_INTEGER4,MPI_REAL8,   &
                                      MPI_ISend, MPI_AllReduce, MPI_statuses_ignore, MPI_Request, &
                                      MPI_Status, MPI_Status_size, MPI_Waitall, MPI_status_ignore
    use type_mpi_request_vector,only: mpi_request_vector_t
    implicit none


    !>  Storage for sensitivities of funtionals wrt boundary condition's property
    !!
    !!
    !!  dL/da = v*dR/da + w*dA/da
    !!
    !!      v = primary adjoint variables
    !!      w = auxiliary adjoint variables
    !!      R = residuals of primal problem
    !!      A = residuals of auxiliary problem
    !!      a = boundary condition property (ie 'Total Pressure','Value','AoA',ecc)
    !!
    !!  If 'a' is a BC of the primal problem dA/da will be zero. Whereas, if 'a' is a BC of 
    !!  the auxiliary problem dR/da will be zero.
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    type,   public  :: chidg_adjointbc_t

        ! dR/dX matrix
        type(chidg_vector_t)                    :: Ra
        
        ! Resultant vector from [v]*[dR/da], v being adjoint variables realtive to each 
        ! functional at a particular istep
        real(rk)            ,   allocatable     :: vRa(:)                   ! [ifunc]

        ! Boundary condition sensitivities for auxiliary problem dJ/dBC = wAa 
        real(rk)            ,   allocatable     :: wAa(:)                   ! [ifunc]

        ! Sum of functional total derivatives contribution at each istep
        real(rk)            ,   allocatable     :: Ja_unsteady(:)           ! [ifunc]




        ! Initialization completed
        logical         ::  adjointbc_initialized = .false.


    contains
        
        ! Initialization
        procedure       :: init
        procedure       :: init_containers

        ! Clear 
        procedure       :: vRa_clear
        procedure       :: release

    end type chidg_adjointbc_t
    !************************************************************************************


contains


    !>  Allocate storage for adjointbc derivatives
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !!  param[in]   nstep       Number of steps of the solution (nstep = 1 for steady and HB)
    !!  param[in]   nfunc       Number of functionals
    !!
    !------------------------------------------------------------------------------------
    subroutine init(self,nfunc,nstep,sflags)
        class(chidg_adjointbc_t), intent(inout)   :: self
        integer(ik),              intent(in)      :: nstep
        integer(ik),              intent(in)      :: nfunc
        type(storage_flags_t),    intent(in)      :: sflags

        integer(ik)     :: ierr

        !
        ! Initialize vRx array of chidg_vectors
        !
        if (sflags%vRa) then
            if (allocated(self%vRa)) deallocate (self%vRa)
            allocate( self%vRa(nfunc), stat=ierr )
            if (ierr/=0) call AllocationError
        end if

        !
        ! Initialize Jx array of chidg_vectors
        !
        if (sflags%wAa) then
            if (allocated(self%wAa)) deallocate (self%wAa)
            allocate( self%wAa(nfunc), stat=ierr )
            if (ierr/=0) call AllocationError
        end if

        !
        ! Initialize Junsteady array of chidg_vectors
        !
        if (sflags%Ja_unsteady) then
            if (allocated(self%Ja_unsteady)) deallocate (self%Ja_unsteady)
            allocate( self%Ja_unsteady(nfunc), stat=ierr )
            if (ierr/=0) call AllocationError
        end if

        self%adjointbc_initialized = .true.

    end subroutine init
    !************************************************************************************






    !>  Allocate storage necessary for adjointbc chidg_vectors
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine init_containers(self,mesh,ntime,sflags)
        class(chidg_adjointbc_t),   intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(inout)   :: ntime  
        type(storage_flags_t),      intent(in)      :: sflags
             
        integer(ik)                 :: ifunc,istep,nfuncs


        !
        ! Initialize vector
        !
        if (sflags%Ra)  call self%Ra%init(mesh,ntime)

        !
        ! Set Ja = 0
        !
        if (sflags%wAa) then
            do ifunc = 1,size(self%wAa)
                self%wAa(ifunc) = ZERO
            end do
        end if

        !
        ! Set Ja_unsteady = 0
        !
        if (sflags%Ja_unsteady) then
            do ifunc = 1,size(self%Ja_unsteady)
                self%Ja_unsteady(ifunc) = ZERO
            end do
        end if

    end subroutine init_containers
    !************************************************************************************
    






    !>  Set all vRa to ZERO 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   5/28/2019
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine vRa_clear(self)
        class(chidg_adjointbc_t),     intent(inout)   :: self

        integer(ik)     :: ifunc

        do ifunc = 1,size(self%vRa)
            self%vRa(ifunc) = ZERO
        end do

    end subroutine vRa_clear
    !************************************************************************************    
    
    
    
    
    
    
    
    !>  Release allocated data 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   11/26/2018
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine release(self)
        class(chidg_adjointbc_t),     intent(inout)   :: self
        
        if (allocated(self%vRa))                deallocate (self%vRa)
        if (allocated(self%wAa))                deallocate (self%wAa)
        if (allocated(self%Ja_unsteady))        deallocate (self%Ja_unsteady)

        call self%Ra%release()

    end subroutine release
    !************************************************************************************




end module type_chidg_adjointbc
