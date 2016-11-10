module mod_chidg_post
    use mod_chidg_post_hdf2tec, only: chidg_post_hdf2tec
    use mod_chidg_post_hdf2vtk, only: chidg_post_hdf2vtk

    implicit none




contains


    !>  Interface for file conversion.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_post(filename)
        character(*),   intent(in)  :: filename


        call chidg_post_hdf2tec(filename)


    end subroutine chidg_post
    !**************************************************************************



    
    !>  Interface for file conversion.
    !!
    !!  @author Mayank Sharma
    !!  @date   11/06/2016
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_post_vtk(filename)
        character(*),   intent(in)  :: filename


        call chidg_post_hdf2vtk(filename)


    end subroutine chidg_post_vtk
    !**************************************************************************



















end module mod_chidg_post
