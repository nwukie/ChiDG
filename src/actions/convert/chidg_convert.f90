module mod_chidg_convert
    use mod_chidg_convert_p3d_hdf5,     only: chidg_convert_p3d_hdf5
    implicit none




contains


    !>  Interface for file conversion.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------
    subroutine chidg_convert(filename)
        character(*),   intent(in)  :: filename


        call chidg_convert_p3d_hdf5(filename)

!        call chidg_convert_gmsh_hdf5(filename)



    end subroutine chidg_convert
    !**************************************************************************





end module mod_chidg_convert
