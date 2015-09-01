module solver_backward_euler
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, DIAG
    use atype_solver,   only: solver_t
    use type_domain,    only: domain_t
    use type_dict,      only: dict_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_blockvector

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

        real(rk)        :: dt = 0.1_rk                    !< Time-step increment
        integer(ik)     :: nsteps = 10000                   !< Number of time steps to compute
        integer(ik)     :: nwrite = 1                     !< Write data every 'nwrite' steps

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
        class(backward_euler_s),            intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t), optional,    intent(inout)   :: matrixsolver

        character(100)      :: filename
        integer(ik)         :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner
        real(rk)            :: resid
        type(blockvector_t) :: b, qn, qold, qnew, dtau
      


        wcount = 1
        ninner = 10
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'

            !
            ! TIME STEP LOOP
            !
            do itime = 1,self%nsteps
                print*, "Step: ", itime


                ! Store qn, since it will be operated on in the inner loop
                qn = q


                !
                ! NONLINEAR CONVERGENCE INNER LOOP
                !
                resid  = ONE    ! Force inner loop entry
                ninner = 1      ! Initialize inner loop counter
                do while ( resid > 1.0e-8_rk )
                    print*, "   ninner: ", ninner


                    ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                    qold = q




                    ! Update Spatial Residual and Linearization (rhs, lin)
                    call update_space(domain)





                    ! Add mass/dt to block diagonal in dR/dQ
                    do ielem = 1,domain%mesh%nelem
                        iblk = DIAG
                       
                        if (allocated(lin%lblks(ielem,iblk)%mat)) then

                            ! Add mass matrix divided by dt to the block diagonal
                            lin%lblks(ielem,iblk)%mat  =  lin%lblks(ielem,iblk)%mat  +  domain%mesh%elems(ielem)%mass/self%dt

                        end if

                    end do



                    ! Divide pseudo-time derivative by dt and multiply by mass matrix
                    dtau = (qold - qn)/self%dt
                    do ielem = 1,domain%mesh%nelem
                        dtau%lvecs(ielem)%vec = matmul(domain%mesh%elems(ielem)%mass,dtau%lvecs(ielem)%vec)
                    end do



                    ! Assign rhs to b, which should allocate storage
                    !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                    b = (-ONE)*dtau - rhs




                    ! We need to solve the matrix system Ax=b for the update vector x (dq)
                    call matrixsolver%solve(lin,dq,b)



                    ! Advance solution with update vector
                    qnew = qold + dq


                    ! Compute residual of nonlinear iteration
                    resid = dq%norm()


                    ! Clear working storage
                    call rhs%clear()
                    call dq%clear()
                    call lin%clear()



                    ! Store updated solution vector (qnew) to working solution vector (q)
                    q = qnew
                    ninner = ninner + 1

                    print*, "   DQ - Norm: ", resid
                end do ! ninner



                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(domain,trim(filename),itime+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do  ! itime

        end associate

    end subroutine solve







    
    subroutine destructor(self)
        type(backward_euler_s),      intent(in) :: self

    end subroutine




end module solver_backward_euler










