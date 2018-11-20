!>
!!
!!
!!  [ A  B ] [x] = [a]
!!  [ C  D ] [y] = [b]
!!
!!  Schur complement
!!
!!  [A - B(D^-1)C] x = a - B(D^-1) b
!!
!!  alpha   = [A - B(D^-1)C]    ! [block]
!!  beta(:) = B(D^-1)           ! [block, block, block, block, ...]
!!
!!
!!  @author Nathan A. Wukie(AFRL)
!!  @date   11/14/2018
!!
!!
!------------------------------------------------------------------------
module precon_schur_element
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, ONE, ZERO
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK
    use mod_io,                 only: verbosity

    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_chidg_vector,      only: chidg_vector_t
    use type_densematrix,       only: densematrix_t
    use type_element,           only: element_t
    use type_mesh,              only: mesh_t
    use type_ivector,           only: ivector_t
    implicit none



    !>  Schur complement decomposition
    !!
    !!  For local element coupling:
    !!  
    !!      .___.
    !!      | D |
    !!  .___.___.___.
    !!  | D | A | D |
    !!  .---.---.---.
    !!      | D |
    !!      .---.
    !!
    !!  A = Coupling of center element with self
    !!  D = Coupling of neighbor elements with selves
    !!
    !!  B = Coupling of center with neighbors
    !!  C = Coupling of neighbors with center
    !!
    !!  Partitioned system:
    !!  [ A  B ] [x] = [a]
    !!  [ C  D ] [y] = [b]
    !!
    !!  Schur complement:
    !!  [A - B(D^-1)C] x = a - B(D^-1) b
    !!
    !!  Storage:
    !!  alpha   = [A - B(D^-1)C]    ! [block]
    !!  beta(:) = B(D^-1)           ! [block, block, block, block, ...]
    !!
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   11/13/2018
    !!
    !--------------------------------------------------------------------
    type, public :: element_local_preconditioner_t

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l

        type(densematrix_t)              :: alpha    
        type(densematrix_t), allocatable :: beta(:) 

        character, allocatable   :: beta_source(:)

    contains
        procedure   :: init   => init_element
        procedure   :: update => update_element
!        procedure   :: apply  => apply_element
    end type element_local_preconditioner_t
    !********************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   11/13/2018
    !!
    !----------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_schur_element_t

        type(chidg_matrix_t) :: D     ! inverse of block diagonal, (ndom,maxelems,ntime)
         
        type(element_local_preconditioner_t), allocatable :: elems(:)

    contains

        procedure   :: init   => init_preconditioner
        procedure   :: update => update_preconditioner
        procedure   :: apply  => apply_preconditioner

    end type precon_schur_element_t
    !**********************************************************************************


contains

    
    !>  Initialize schur_element preconditioner.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   11/13/2018
    !!
    !------------------------------------------------------------------------
    subroutine init_preconditioner(self,data)
        class(precon_schur_element_t),  intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data

        integer(ik) :: nelements, ierr, idom, ielem, iprecon

        ! Initialize storage for inverted diagonal
        call self%D%init(data%mesh,mtype='Diagonal')

        ! Accumulate total number of elements. 
        nelements = data%mesh%nelements()

        ! Allocate element-local preconditioner storage for each
        allocate(self%elems(nelements),stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Call initialization for each allocated element preconditioner
        iprecon = 1
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelements()
                call self%elems(iprecon)%init(data,data%mesh%domain(idom)%elems(ielem)%idomain_g,  &
                                                   data%mesh%domain(idom)%elems(ielem)%idomain_l,  &
                                                   data%mesh%domain(idom)%elems(ielem)%ielement_g, &
                                                   data%mesh%domain(idom)%elems(ielem)%ielement_l)
                iprecon = iprecon + 1
            end do
        end do

        ! Store spectral matrix and indicate initialization
        self%initialized = .true.

    end subroutine init_preconditioner
    !*****************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   6/24/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine init_element(self,data,idomain_g,idomain_l,ielement_g,ielement_l)
        class(element_local_preconditioner_t),  intent(inout)   :: self
        type(chidg_data_t),                     intent(in)      :: data
        integer(ik),                            intent(in)      :: idomain_g
        integer(ik),                            intent(in)      :: idomain_l
        integer(ik),                            intent(in)      :: ielement_g
        integer(ik),                            intent(in)      :: ielement_l

        integer(ik) :: itime, nlocal, noverset, ncoupled, ierr, idiag, ibeta, ilocal, ioverset

        ! Initialize element indices
        self%idomain_g  = idomain_g
        self%idomain_l  = idomain_l
        self%ielement_g = ielement_g
        self%ielement_l = ielement_l

        itime = 1

        ! Compute number of directly-coupled elements
        ! Subtract 1 from locally coupled elements that is the diagonal block
        nlocal   = data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%size()
        noverset = data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%size()
        ncoupled = (nlocal - 1) + noverset


        ! Allocate beta storage for each coupled element
        allocate(self%beta(ncoupled), self%beta_source(ncoupled), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Initialize alpha, beta(:)
        idiag = data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%get_diagonal()
        call self%alpha%init(data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%nterms,         &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%nfields,        &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%dparent_g(),    &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%dparent_l(),    &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%eparent_g(),    &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%eparent_l(),    &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%parent_proc(),  &
                             data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(idiag)%tparent())
                
        

        ! Initialize beta(:)
        ibeta = 1
        do ilocal = 1,nlocal
            if (ilocal /= idiag) then
                call self%beta(ibeta)%init(data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%nterms,        &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%nfields,       &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%dparent_g(),   &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%dparent_l(),   &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%eparent_g(),   &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%eparent_l(),   &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%parent_proc(), &
                                           data%sdata%lhs%dom(idomain_l)%lblks(ielement_l,itime)%data_(ilocal)%tparent())

                self%beta_source(ibeta) = 'L'
                ibeta = ibeta + 1
            end if
        end do !ineighbor


        do ioverset = 1,noverset
            call self%beta(ibeta)%init(data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%nterms,        &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%nfields,       &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%dparent_g(),   &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%dparent_l(),   &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%eparent_g(),   &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%eparent_l(),   &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%parent_proc(), &
                                       data%sdata%lhs%dom(idomain_l)%chi_blks(ielement_l,itime)%data_(ioverset)%tparent())

            self%beta_source(ibeta) = 'O'
            ibeta = ibeta + 1
        end do !ioverset



    end subroutine init_element
    !*****************************************************************************

    






    !>  Update the schur_element preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine update_preconditioner(self,A,b)
        class(precon_schur_element_t),  intent(inout)   :: self
        type(chidg_matrix_t),       intent(in)      :: A
        type(chidg_vector_t),       intent(in)      :: b

        integer(ik) :: idom, ielem, itime, idiag, iprecon

        ! Invert the block diagonal from A and store
        do idom = 1,size(A%dom)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do itime = 1,size(A%dom(idom)%lblks,2)
                    idiag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                    self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = A%dom(idom)%lblks(ielem,itime)%dmat(idiag)
                    self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = inv(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat)
                end do
            end do
        end do
        
        ! Compute beta = B D^-1
        do iprecon = 1,size(self%elems)
            call self%elems(iprecon)%update(A,b,self%D)
        end do !iP

    end subroutine update_preconditioner
    !*****************************************************************************





    !>  Update element
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   11/14/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine update_element(self,A,b,D)
        class(element_local_preconditioner_t),  intent(inout)   :: self
        class(chidg_matrix_t),                  intent(in)      :: A
        class(chidg_vector_t),                  intent(in)      :: b
        class(chidg_matrix_t),                  intent(in)      :: D

        integer(ik) :: ibeta, ddomain_g, ddomain_l, delement_g, delement_l, itime, B_ind, C_ind, idiag

        itime = 1


        ! Update beta
        do ibeta = 1,size(self%beta)
            ddomain_g  = self%beta(ibeta)%dparent_g() 
            ddomain_l  = self%beta(ibeta)%dparent_l() 
            delement_g = self%beta(ibeta)%eparent_g()
            delement_l = self%beta(ibeta)%eparent_l()

            ! B D^-1
            if (self%beta_source(ibeta) == 'L') then
                B_ind = A%dom(self%idomain_l)%lblks(self%ielement_l,itime)%loc(ddomain_g,delement_g,itime)
                self%beta(ibeta)%mat = matmul(A%dom(self%idomain_l)%lblks(self%ielement_l,itime)%data_(B_ind)%mat,D%dom(ddomain_l)%lblks(delement_l,itime)%data_(1)%mat)
            else if (self%beta_source(ibeta) == 'O') then
                B_ind = A%dom(self%idomain_l)%chi_blks(self%ielement_l,itime)%loc(ddomain_g,delement_g,itime)
                self%beta(ibeta)%mat = matmul(A%dom(self%idomain_l)%chi_blks(self%ielement_l,itime)%data_(B_ind)%mat,D%dom(ddomain_l)%lblks(delement_l,itime)%data_(1)%mat)
            end if

        end do !ibeta 


        ! Update alpha (after beta, since we use beta here)
        idiag = A%dom(self%idomain_l)%lblks(self%ielement_l,itime)%get_diagonal()
        self%alpha%mat = A%dom(self%idomain_l)%lblks(self%ielement_l,itime)%data_(idiag)%mat
        do ibeta = 1,size(self%beta)

            ddomain_g  = self%beta(ibeta)%dparent_g() 
            ddomain_l  = self%beta(ibeta)%dparent_l() 
            delement_g = self%beta(ibeta)%eparent_g()
            delement_l = self%beta(ibeta)%eparent_l()


            ! Find C entry to compute: B D^-1 C = beta C, which is subtracted from A and stored as alpha
            if (self%beta_source(ibeta) == 'L') then
                C_ind = A%dom(ddomain_l)%lblks(delement_l,itime)%loc(self%idomain_g,self%ielement_g,itime)
                self%alpha%mat = self%alpha%mat - matmul(self%beta(ibeta)%mat,A%dom(ddomain_l)%lblks(delement_l,itime)%data_(C_ind)%mat)
            else if (self%beta_source(ibeta) == 'O') then
                C_ind = A%dom(ddomain_l)%chi_blks(delement_l,itime)%loc(self%idomain_g,self%ielement_g,itime)
                self%alpha%mat = self%alpha%mat - matmul(self%beta(ibeta)%mat,A%dom(ddomain_l)%chi_blks(delement_l,itime)%data_(C_ind)%mat)
            end if

        end do

        ! Invert alpha
        self%alpha%mat = inv(self%alpha%mat)

    end subroutine update_element
    !*****************************************************************************









    !>  Apply the schur_element preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/25/2018
    !!
    !-----------------------------------------------------------------------------
    function apply_preconditioner(self,A,v,z_old) result(z)
        class(precon_schur_element_t),  intent(inout)   :: self
        type(chidg_matrix_t),           intent(in)      :: A
        type(chidg_vector_t),           intent(in)      :: v
        type(chidg_vector_t),           intent(in), optional :: z_old

        type(chidg_vector_t)    :: z

        integer(ik) :: itime, iprecon, ibeta, ddomain_l, delement_l, idiag, icoupling, couple_domain, couple_element, C_ind

        real(rk),   allocatable :: vec(:), b(:)


        itime = 1

        ! Allocate/clear z
        z = v


        do iprecon = 1,size(self%elems)

            ! Assemble right-hand side: a
            vec = v%dom(self%elems(iprecon)%idomain_l)%vecs(self%elems(iprecon)%ielement_l)%vec


            
            ! Assemble right-hand side: a - BD^-1 b
            do ibeta = 1,size(self%elems(iprecon)%beta)

                ! Get access
                ddomain_l  = self%elems(iprecon)%beta(ibeta)%dparent_l()
                delement_l = self%elems(iprecon)%beta(ibeta)%eparent_l()

                ! Store b
                b = v%dom(ddomain_l)%vecs(delement_l)%vec

!                ! Subtract explicit coupling: b - E y
!                idiag = A%dom(ddomain_l)%lblks(delement_l,itime)%get_diagonal()
!                C_ind = A%dom(ddomain_l)%lblks(delement_l,itime)%loc(self%elems(iprecon)%idomain_g,self%elems(iprecon)%ielement_g,itime)
!                do icoupling = 1,A%dom(ddomain_l)%lblks(delement_l,itime)%size()
!                    if (icoupling /= idiag .and. icoupling /= C_ind) then
!
!                        couple_domain  = A%dom(ddomain_l)%lblks(delement_l,itime)%data_(icoupling)%dparent_l()
!                        couple_element = A%dom(ddomain_l)%lblks(delement_l,itime)%data_(icoupling)%eparent_l()
!
!                        b = b - matmul(A%dom(ddomain_l)%lblks(delement_l,itime)%data_(icoupling)%mat,v%dom(couple_domain)%vecs(couple_element)%vec)
!
!                    end if
!                end do !icoupling

                ! a - BD^-1 b
                vec = vec - matmul(self%elems(iprecon)%beta(ibeta)%mat,b)

            end do !ibeta

            ! Apply preconditioner: z = alpha^-1 (a - BD^-1 b)
            z%dom(self%elems(iprecon)%idomain_l)%vecs(self%elems(iprecon)%ielement_l)%vec = matmul(self%elems(iprecon)%alpha%mat,vec)

        end do !iprecon



    end function apply_preconditioner
    !*****************************************************************************






end module precon_schur_element
