module mod_chidg_post
    use mod_chidg_post_hdf2tec,         only: chidg_post_hdf2tec
    use mod_chidg_post_hdf2vtk,         only: chidg_post_hdf2vtk
    use mod_chidg_post_hdf2matplotlib,  only: chidg_post_hdf2matplotlib

    implicit none




contains


    !>  Interface for file conversion.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_post(grid_file,solution_file)
        character(*),   intent(in)  :: grid_file
        character(*),   intent(in)  :: solution_file


        call chidg_post_hdf2tec(grid_file,solution_file)


    end subroutine chidg_post
    !**************************************************************************



    
    !>  Interface for file conversion.
    !!
    !!  @author Mayank Sharma
    !!  @date   11/06/2016
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_post_vtk(grid_file,solution_file)
        character(*),   intent(in)  :: grid_file
        character(*),   intent(in)  :: solution_file


        call chidg_post_hdf2vtk(grid_file,solution_file)


    end subroutine chidg_post_vtk
    !**************************************************************************





    !>  Interface for file conversion.
    !!
    !!  @author Mayank Sharma
    !!  @date   3/24/2017
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_post_matplotlib(grid_file,solution_file)
        character(*),   intent(in)  :: grid_file
        character(*),   intent(in)  :: solution_file


        call chidg_post_hdf2matplotlib(grid_file,solution_file)


    end subroutine chidg_post_matplotlib
    !**************************************************************************



















end module mod_chidg_post
