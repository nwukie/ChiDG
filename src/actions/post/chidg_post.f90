module mod_chidg_post
    use mod_chidg_post_hdf2tec, only: chidg_post_hdf2tec
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





end module mod_chidg_post
