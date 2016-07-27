module mod_metis
    implicit none





interface
    !subroutine METIS_PartMeshNodal( ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgt, opts, objval, epart, npart) bind(c)
    function METIS_PartMeshNodal( ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgt, opts, objval, epart, npart) bind(c) result(test)
      use iso_c_binding,    only: c_int, c_ptr
      implicit none

      integer(c_int),                  intent(in)  :: ne, nn
      integer(c_int), dimension(*),    intent(in)  :: eptr, eind
      type(c_ptr),                          value   :: vwgt, vsize   
      integer(c_int),                  intent(in)  :: nparts
      type(c_ptr),                          value   :: tpwgt
      type(c_ptr),                          value   :: opts
      integer(c_int),                  intent(out) :: objval
      integer(c_int), dimension(*),    intent(out) :: epart
      integer(c_int), dimension(*),    intent(out) :: npart
      integer(c_int)                               :: test

    !end subroutine METIS_PartMeshNodal  
    end function METIS_PartMeshNodal  
end interface








end module mod_metis
