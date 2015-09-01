module solver_forward_euler
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO
    use atype_solver,       only: solver_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_domain,        only: domain_t
    use type_dict,          only: dict_t
    use type_expansion
    use type_blockvector

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables
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
    type, extends(solver_t), public :: forward_euler_s

        real(rk)        :: dt = 0.001_rk    !< Time-step increment
        integer(ik)     :: nsteps = 10000     !< Number of time steps to compute
        integer(ik)     :: nwrite = 100      !< Write data every 'nwrite' steps

    contains
        procedure   :: init
        procedure   :: solve

        final :: destructor
    end type forward_euler_s
    !-----------------------------------------------------------

contains


    !> Solver initialization
    !!  - set solver member data
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine  init(self,domain,options)
        class(forward_euler_s),     intent(inout)   :: self
        type(domain_t),             intent(inout)   :: domain
        type(dict_t), optional,     intent(inout)   :: options

        ! If the options type is passed, use it to set the following data.
        ! Else, the default values will be used.
        if (present(options)) then
            call options%get('dt',self%dt)
            call options%get('nsteps',self%nsteps)
            call options%get('nwrite',self%nwrite)
        end if

    end subroutine init








    !> Solve for update 'dq'
    !!
    !! \f$ \delta Q = - \delta t R(Q) \f$
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine solve(self,domain,matrixsolver)
        class(forward_euler_s),             intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t), optional,    intent(inout)   :: matrixsolver

        character(100)  :: filename
        integer(ik)     :: itime, nsteps, ielem, wcount, iblk


        wcount = 1
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'
            do itime = 1,self%nsteps
                print*, "Step: ", itime


                ! Update Spatial Residual and Linearization (rhs, lin)
                call update_space(domain)


                ! Multiply RHS by mass matrix 
                do ielem = 1,domain%mesh%nelem
                    rhs%lvecs(ielem)%vec = matmul(domain%mesh%elems(ielem)%invmass, rhs%lvecs(ielem)%vec)
                end do


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
        type(forward_euler_s),      intent(in) :: self

    end subroutine




end module solver_forward_euler










