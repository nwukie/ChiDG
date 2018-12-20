module type_rbf_interpolation_matrix
#include <messenger.h>
    use mod_kinds,                  only: ik, rk
    use mod_inv,                    only: inv
    implicit none


    type, public :: rbf_interpolation_matrix_t
        
        ! We assume the form 
        ! |A 0||cb| = |sb|
        ! |B C||ce|   |se|
        !
        ! where A and B are base and C is lower triangular.
        !
        ! We then can solve the system via
        !
        ! 1. Solve the dense system, A*cb = sb to obtain
        !       cb = Ainv*sb
        ! 2. Solve the lower triangular system, C*ce = se-B*cb, via forward substitution to obtain
        !       ce = Cinv*(se-B*cb)
        ! 


        real(rk), allocatable                   :: A(:,:)
        real(rk), allocatable                   :: Ainv(:,:)
        real(rk), allocatable                   :: B(:,:)
        real(rk), allocatable                   :: C(:,:)

        logical                                 :: is_initialized = .false.
        logical                                 :: is_filled = .false.

    contains

        procedure           :: init 
        procedure           :: solve

    end type rbf_interpolation_matrix_t

contains

    subroutine init(self, nnodes_base, nnodes_explicit) 
        class(rbf_interpolation_matrix_t),      intent(inout)           :: self
        integer(ik),                            intent(in)              :: nnodes_base
        integer(ik),                            intent(in), optional    :: nnodes_explicit

        allocate(self%A(nnodes_base, nnodes_base))
        if (present(nnodes_explicit)) then
            if (nnodes_explicit>0) then
                allocate(self%B(nnodes_explicit, nnodes_base))
                allocate(self%C(nnodes_explicit, nnodes_explicit))
            end if
        end if

    end subroutine



    function solve(self, source_val) result(rbf_coeff)
        class(rbf_interpolation_matrix_t),      intent(inout)           :: self
        real(rk),                               intent(in)              :: source_val(:)

        real(rk), allocatable :: rbf_coeff(:)
        real(rk), allocatable :: rbf_coeff_base(:), rbf_coeff_explicit(:), temp(:), rhs(:)
        real(rk), allocatable :: source_val_base(:)
        real(rk), allocatable :: source_val_explicit(:) 

        integer(ik) :: inode_sol, inode_term, nnodes_base, nnodes_explicit

        if (.not. self%is_filled) call chidg_signal(FATAL, 'rbf_interpolation_matrix%solve : matrix is not filled')

        nnodes_base = size(self%A,1)
        ! Invert the base block
        if (.not. allocated(self%Ainv)) then

            self%Ainv = inv(self%A)

        end if

        ! Solve for the basely coupled RBF coefficients
        allocate(rbf_coeff_base(nnodes_base), source_val_base(nnodes_base))
        source_val_base = source_val(1:nnodes_base)
        rbf_coeff_base = matmul(self%Ainv, source_val_base)
        !rbf_coeff_base = solve_dense(self%A, source_val_base)
        




        ! Now, check if we have explicit RBF nodes, and solve if so
        if ((allocated(self%C))) then
            nnodes_explicit = size(self%C,1)
            allocate(rbf_coeff_explicit(nnodes_explicit), source_val_explicit(nnodes_explicit))
            
            temp = matmul(self%B, rbf_coeff_base)

            source_val_explicit = source_val(nnodes_base+1:nnodes_base+nnodes_explicit)
            rhs = source_val_explicit-temp


            rbf_coeff_explicit(1) = rhs(1)/self%C(1,1)
            do inode_sol = 2, size(source_val_explicit,1)

                rbf_coeff_explicit(inode_sol) = rhs(inode_sol)
                do inode_term = 1, (inode_sol-1)
                    rbf_coeff_explicit(inode_sol) = rbf_coeff_explicit(inode_sol) -&
                    self%C(inode_sol,inode_term)*rbf_coeff_explicit(inode_term)

                end do

                rbf_coeff_explicit(inode_sol) = rbf_coeff_explicit(inode_sol)/self%C(inode_sol, inode_sol)

            end do

            
            allocate(rbf_coeff(nnodes_base+nnodes_explicit))
            rbf_coeff(1:nnodes_base) = rbf_coeff_base
            rbf_coeff(nnodes_base+1:(nnodes_base+nnodes_explicit)) = rbf_coeff_explicit
        else
            rbf_coeff = rbf_coeff_base

        end if

    end function solve


end module type_rbf_interpolation_matrix
