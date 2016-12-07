module type_chidg_file
    use mod_kinds,  only: rk, ik
    use mod_hdfio,  only: open_file_hdf, close_file_hdf
    use hdf5
    use h5lt
    implicit none





    !>  A class for handling ChiDG-formatted files.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, public :: chidg_file_t

        ! HDF identifier for the file
        integer(HID_T)  :: file_id

        ! HDF identifier for a domain that is open
        integer(HID_T)  :: domain_id

    contains

        procedure   :: open
        procedure   :: close

    end type chidg_file_t
    !*****************************************************************************





contains





    !>  Open a ChiDG file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine open(self,filename)
        class(chidg_file_t),    intent(inout)   :: self
        character(*),           intent(in)      :: filename


        self%file_id = open_file_hdf(filename)


    end subroutine open
    !****************************************************************************************





    !>  Close the ChiDG file
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine close(self)
        class(chidg_file_t),    intent(inout)   :: self


        call close_file_hdf(self%file_id)


    end subroutine close
    !****************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie 
    !!  @date   10/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine convert(self,selector)
        class(chidg_file_t),    intent(inout)   :: self
        character(*),           intent(in)      :: selector




    end subroutine convert
    !****************************************************************************************


end module type_chidg_file
