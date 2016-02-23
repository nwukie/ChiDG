module mod_chidg_interpolate
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use type_chidg,         only: chidg_t
    use mod_grid_tools_two, only: compute_element_donor
    use mod_grid_operators, only: solution_point
    use mod_io,             only: nterms_s
    implicit none












contains



    !>  Interpolate solution from one grid to another grid.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/21/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine chidg_interpolate(sourcefile, targetfile)
        character(*),   intent(in)  :: sourcefile
        character(*),   intent(in)  :: targetfile



        type(chidg_t)           :: chidg_source
        type(chidg_t)           :: chidg_target

        real(rk)                :: xi, eta, zeta
        real(rk),   allocatable :: vals(:), val_modes(:)
        integer(ik)             :: idom, ielem, ivar, inode, idom_d, ielem_d, ierr
        type(point_t)           :: node, point_comp

        nterms_s = 3*3*3


        !
        ! Initialize ChiDG environment. Actually should only need to be called once.
        !
        call chidg_source%init('env')


        !
        ! Read grid data from files.
        !
        print*, 'Reading grids: ', trim(sourcefile), trim(targetfile)
        call chidg_source%read_grid(trim(sourcefile))
        call chidg_target%read_grid(trim(targetfile))

        call chidg_source%data%init_sdata()
        call chidg_target%data%init_sdata()


        !
        ! Read solution from source
        !
        print*, 'Reading solution: ', trim(sourcefile)
        call chidg_source%read_solution(trim(sourcefile))











       !
       ! Loop through elements in mesh and call function projection
       !
       do idom = 1,chidg_target%data%ndomains()

           do ielem = 1,chidg_target%data%mesh(idom)%nelem
               print*, 'Domain ', idom, 'Element ', ielem


               do ivar = 1,chidg_target%data%eqnset(idom)%item%neqns


                   !
                   ! Interpolate solution from source to target at integration points for projection.
                   !
                   if ( allocated(vals) ) deallocate(vals)
                   allocate(vals(size(chidg_target%data%mesh(idom)%elems(ielem)%quad_pts)), stat=ierr )
                   if (ierr /= 0) call AllocationError

                   do inode = 1,size(chidg_target%data%mesh(idom)%elems(ielem)%quad_pts)

                       node = chidg_target%data%mesh(idom)%elems(ielem)%quad_pts(inode)

                       !
                       ! Find donor domain/element in source chidg instance.
                       !
                       call compute_element_donor(chidg_source%data%mesh, node, idom_d, ielem_d, point_comp)

                       
                       !
                       ! Get solution at node from source chidg instance
                       !
                       xi   = point_comp%c1_
                       eta  = point_comp%c2_
                       zeta = point_comp%c3_
                       vals(inode) = solution_point(chidg_source%data%sdata%q%dom(idom_d)%lvecs(ielem_d),ivar,xi,eta,zeta)


                   end do !inode

                    
                    !
                    ! Multiply by quadratre weights
                    !
                    vals = vals * chidg_target%data%mesh(idom)%elems(ielem)%gq%vol%weights

                   val_modes = matmul(transpose(chidg_target%data%mesh(idom)%elems(ielem)%gq%vol%val), vals) / chidg_target%data%mesh(idom)%elems(ielem)%gq%vol%dmass



                   !
                   ! Store the projected modes to the solution expansion
                   !
                   call chidg_target%data%sdata%q%dom(idom)%lvecs(ielem)%setvar(ivar,val_modes)


               end do ! ivar

           end do ! ielem

       end do ! idomain






        !
        ! Write interpolated/projected solution to file
        !
        call chidg_target%write_solution(targetfile)

    
    end subroutine chidg_interpolate
    !*********************************************************************************************






end module mod_chidg_interpolate
