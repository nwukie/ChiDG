module solver_backward_euler
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    use atype_solver,   only: solver_t
    use type_domain,    only: domain_t
    use type_dict,      only: dict_t
    use type_expansion

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------
    type, extends(solver_t), public :: backward_euler_s

        real(rk)        :: dt = 0.001_rk                    !< Time-step increment
        integer(ik)     :: nsteps = 10000                   !< Number of time steps to compute
        integer(ik)     :: nwrite = 100                     !< Write data every 'nwrite' steps

    contains
        procedure   :: init
        procedure   :: solve

        final :: destructor
    end type backward_euler_s
    !-----------------------------------------------------------

contains


    !> Solver initialization
    !!  - set solver member data
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine  init(self,domain,options)
        class(backward_euler_s),     intent(inout)   :: self
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
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine solve(self,domain,matrixsolver)
        class(backward_euler_s),    intent(inout)   :: self
        type(domain_t),             intent(inout)   :: domain
        class(matrixsolver_t),      intent(inout)   :: matrixsolver

        character(100)  :: filename
        integer(ik)     :: itime, nsteps, ielem, wcount, iblk


        wcount = 1
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'
            do itime = 1,self%nsteps
                print*, "Step: ", itime



                do 



                ! Update Spatial Residual and Linearization (rhs, lin)
                call update_space(domain)










                ! Multiply RHS by mass matrix for explicit time-integration
                !do ielem = 1,domain%mesh%nelem
                !    rhs(ielem)%vec = matmul(domain%mesh%elems(ielem)%invmass, rhs(ielem)%vec)
                !end do




                ! We need to solve the matrix system Ax=b for the update vector dq
                call matrixsolver%solver(lin,dq,rhs)





                ! Compute update vector
                dq = dt * rhs

                ! Advance solution with update vector
                q  = q + dq


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(domain,trim(filename),itime+1)
                    wcount = 0
                end if
                wcount = wcount + 1



                ! Clear spatial residual
                do ielem = 1,domain%mesh%nelem
                    rhs(ielem)%vec = ZERO

                    !do iblk = 1,7
                    !    lin%lblks(ielem,iblk)%mat = ZERO
                    !end do
                end do

                call lin%clear()    ! Clear block-sparse matrix containing linearization of the spatial scheme


            end do

        end associate

    end subroutine solve







    
    subroutine destructor(self)
        type(backward_euler_s),      intent(in) :: self

    end subroutine




end module solver_backward_euler










