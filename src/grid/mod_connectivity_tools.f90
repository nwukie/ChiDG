module mod_connectivity_tools
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    implicit none




contains


    !>  Return the element index from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function connectivity_get_ielem(connectivity) result(ielem)
        integer(ik),    intent(in)  :: connectivity(:)

        integer(ik) :: ielem

        ielem = connectivity(1)

    end function connectivity_get_ielem
    !********************************************************************************









    !>  Return the element mapping from a connectivity
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    function connectivity_get_mapping(connectivity) result(mapping)
        integer(ik),    intent(in)  :: connectivity(:)

        integer(ik) :: mapping

        mapping = connectivity(2)

    end function connectivity_get_mapping
    !********************************************************************************








    !>  Return element node index
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/23/2016
    !!
    !!
    !--------------------------------------------------------------------------------
    function connectivity_get_node(connectivity,ipt) result(inode)
        integer(ik),    intent(in)  :: connectivity(:)
        integer(ik),    intent(in)  :: ipt

        integer(ik) :: inode

        inode = connectivity(2+ipt)

    end function connectivity_get_node
    !********************************************************************************



end module mod_connectivity_tools
