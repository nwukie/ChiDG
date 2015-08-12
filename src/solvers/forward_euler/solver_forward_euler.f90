module solver_forward_euler
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    use atype_solver,   only: solver_t
    use type_domain,    only: domain_t
    use type_expansion

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables
    implicit none
    private


    type, extends(solver_t), public :: forward_euler_s

        real(rk)        :: dt

    contains
        procedure   :: init
        procedure   :: solve

        final :: destructor
    end type forward_euler_s

contains


    !> Solver initialization
    subroutine  init(self,domain)
        class(forward_euler_s),     intent(inout)   :: self
        type(domain_t),             intent(inout)   :: domain

        !> Call any other specialized initialization requirements
        self%dt = 0.005_rk
    end subroutine



    !> Solve for update 'dq'
    subroutine solve(self,domain)
        class(forward_euler_s), intent(inout)   :: self
        type(domain_t),         intent(inout)   :: domain

        character(100)  :: filename
        integer(ik)     :: itime, ntime, ielem, wcount, iblk


        ntime = 20000
        wcount = 1
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'
            do itime = 1,ntime
                print*, "Step: ", itime


                !> Update Spatial Residual and Linearization (rhs, lin)
                call update_space(domain)


                !> Multiply RHS by mass matrix for explicit time-integration
                do ielem = 1,domain%mesh%nelem
                    rhs(ielem)%vec = matmul(domain%mesh%elems(ielem)%invmass, rhs(ielem)%vec)
                end do


                !> Compute update vector
                dq = dt * rhs

                !> Advance solution with update vector
                q  = q + dq


                if (wcount == 50) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(domain,trim(filename),itime+1)
                    wcount = 0
                end if



                !> Clear spatial residual
                do ielem = 1,domain%mesh%nelem
                    rhs(ielem)%vec = ZERO

                    do iblk = 1,7
                        lin%lblks(ielem,iblk)%mat = ZERO
                    end do
                end do

                wcount = wcount + 1
            end do

        end associate

    end subroutine



    
    subroutine destructor(self)
        type(forward_euler_s),      intent(in) :: self
    end subroutine

end module solver_forward_euler
