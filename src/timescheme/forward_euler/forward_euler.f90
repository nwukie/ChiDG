module forward_euler
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO
    use atype_time_scheme,      only: time_scheme_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_domain,            only: domain_t
    use type_dict,              only: dict_t
    !use type_expansion
    use type_blockvector

    use mod_spatial,            only: update_space

    use mod_tecio,              only: write_tecio_variables
    implicit none
    private



    !>  Solution advancement via the forward-euler method
    !!
    !! Given the system of partial differential equations consisting of the time-derivative
    !! of the solution vector and a spatial residual as
    !!
    !! \f$ \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !! the time derivative is discretized by a one-sided finite-difference approximation as
    !!
    !! \f$ \frac{Q^{n+1} - Q^{n}}{\Delta t} + R(Q^n) = 0 \f$
    !!
    !! The solution at the next time level is then computed as
    !!
    !! \f$ Q^{n+1} = Q^{n} - \Delta t R(Q^n) \f$
    !!
    !! or
    !!
    !! \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !! where \f$ \Delta Q \f$ is defined as
    !!
    !! \f$ \Delta Q = -\Delta t R(Q) \f$
    !!
    !!
    !! This routine computes \f$ \Delta Q \f$ and updates the solution as \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------
    type, extends(time_scheme_t), public :: forward_euler_t


    contains
        procedure   :: solve

        final :: destructor
    end type forward_euler_t
    !-----------------------------------------------------------

contains



    !> Solve for update 'dq'
    !!
    !! \f$ \delta Q = - \delta t R(Q) \f$
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine solve(self,domain,matrixsolver,preconditioner)
        class(forward_euler_t),             intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t),   optional,  intent(inout)   :: matrixsolver
        class(preconditioner_t), optional,  intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, ieqn
        real(rk), allocatable   :: vals(:)


        wcount = 1
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'
            do itime = 1,self%nsteps
                print*, "Step: ", itime


                !call q%dump()

                ! Update Spatial Residual and Linearization (rhs, lin)
                call update_space(domain)


                ! Multiply RHS by mass matrix 
                do ielem = 1,domain%mesh%nelem
                    do ieqn = 1,domain%eqnset%neqns
                        vals = matmul(domain%mesh%elems(ielem)%invmass, rhs%lvecs(ielem)%getvar(ieqn))
                        call rhs%lvecs(ielem)%setvar(ieqn,vals)
                    end do
                end do


                !call rhs%dump()
                !read(*,*)

                ! Compute update vector
                dq = (-dt) * rhs

                ! Advance solution with update vector
                q  = q + dq


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(domain,trim(filename),itime+1)
                    wcount = 0
                end if



                ! Clear residual and linearization storage
                call rhs%clear()
                call lin%clear()

                wcount = wcount + 1
            end do

        end associate

    end subroutine solve







    
    subroutine destructor(self)
        type(forward_euler_t),      intent(in) :: self

    end subroutine




end module forward_euler










