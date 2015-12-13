module mod_tecio_interface
    use iso_c_binding
    use mod_kinds,      only: rk,ik,TEC
    use mod_io,         only: OUTPUT_RES
    use type_mesh,      only: mesh_t


    implicit none

#include "tecio.f90"

contains

    ! This opens a new tecplot binary file for writing
    ! and initializes the type, and number of data fields to accept
    subroutine init_tecio_file(title,variables,filename,filetype)
        character(*)    :: title
        character(*)    :: variables
        character(*)    :: filename
        integer(TEC)   :: filetype

        integer(4)     :: tecstat
        character           :: NULLCHAR = char(0)
        integer(TEC)   :: fileformat = 0           ! 0 = tecplot binary .plt   ::   1 = tecplot subzone loadable .szplt
!        integer(kind=TEC)   :: filetype   = 0           ! 0 = full  :: 1 = grid  :: 2 = solution
        integer(TEC)   :: debug      = 0           ! 0 = debug off
        integer(TEC)   :: isdouble   = 1           ! 0 = single precision  :: 1 = double precision

        tecstat = TECINI142(trim(title)//NULLCHAR,      &
                            trim(variables)//NULLCHAR,  &
                            trim(filename)//NULLCHAR,      &
                            '.'//NULLCHAR, &
                            fileformat, &
                            filetype,   &
                            debug,      &
                            isdouble)
        ! Test file initialization
        if (tecstat /= 0) stop "Error: Initializing TECINI142"
    end subroutine


    ! This begins a new zone in the current opened file. Must be called
    ! after init_tecplot_file, because it needs an open binary file.
    ! If multiple files are open, you can switch between them with the
    ! TECFIL142 call, as long as the they can be identified by integer values
    subroutine init_tecio_zone(zonetitle,mesh,writetype,timeindex)
            type(mesh_t) , intent(in)        :: mesh
            integer(ik)  , intent(in)        :: writetype  ! tells us if we are writing a mesh or solution file
            integer(ik)  , intent(in)        :: timeindex

            character(*)   :: zonetitle
            integer(TEC)   :: zonetype                 = 0     ! 0 - ordered
            integer(TEC)   :: imax
            integer(TEC)   :: jmax
            integer(TEC)   :: kmax
            integer(TEC)   :: icellmax                 = 0     ! not used
            integer(TEC)   :: jcellmax                 = 0     ! not used
            integer(TEC)   :: kcellmax                 = 0     ! not used
            real(rk)       :: solutiontime             = 0._rk
            integer(TEC)   :: strandid                 = 0
            integer(TEC)   :: parentzone               = 0
            integer(TEC)   :: isblock                  = 1
            integer(TEC)   :: nfconns                  = 0
            integer(TEC)   :: fnmode                   = 0
            integer(TEC)   :: totalnumfacenodes        = 1
            integer(TEC)   :: totalnumbndryfaces       = 1
            integer(TEC)   :: totalnumbndryconnections = 1
            integer(TEC)   :: passivevars(3)           = 0        ! null = all varaibles active
            integer(TEC)   :: vallocation(3)           = 1        ! null = all variables node-centered
            integer(TEC)   :: sharvarfrom(3)           = 0        ! null = no data shared between zones
            integer(TEC)   :: sharconnfrom             = 0

            integer(4)     :: tecstat
            integer(TEC), pointer :: NullPtr(:) => null()    ! Null pointer array


            solutiontime = real(timeindex,rk)
            strandid = timeindex
            ! Compute the imax, jmax, kmax for the ordered zone. This comes from
            ! multiplying the number of elements in each direction by the output
            ! resolution, output_res, that the continuous modal polynomials are
            ! sampled at
            if (writetype == 0) then        !Write mesh
                imax = (OUTPUT_RES)*mesh%nelem_xi   + 1
                jmax = (OUTPUT_RES)*mesh%nelem_eta  + 1
                kmax = (OUTPUT_RES)*mesh%nelem_zeta + 1
            elseif (writetype == 1) then    ! Write solution
                imax = (OUTPUT_RES+1)*mesh%nelem_xi
                jmax = (OUTPUT_RES+1)*mesh%nelem_eta
                kmax = (OUTPUT_RES+1)*mesh%nelem_zeta
            else
                stop "Error: invalid writetype for init_tecio_zone"
            end if


            tecstat = TECZNE142(trim(zonetitle)//char(0),   &
                                zonetype,                   &
                                imax,                       &
                                jmax,                       &
                                kmax,                       &
                                icellmax,                   &
                                jcellmax,                   &
                                kcellmax,                   &
                                solutiontime,               &
                                strandid,                   &
                                parentzone,                 &
                                isblock,                    &
                                nfconns,                    &
                                fnmode,                     &
                                totalnumfacenodes,          &
                                totalnumbndryfaces,         &
                                totalnumbndryconnections,   &
                                NullPtr,                    &
                                NullPtr,                    &
                                NullPtr,                    &
                                sharconnfrom)
            if(tecstat /= 0) stop "Error in TECZNE initialization"

    end subroutine


    subroutine finalize_tecio()
        integer(kind=TEC) :: tecstat

        tecstat = TECEND142()
        if (tecstat /= 0) stop "Error: TECIO finalization error"

    end subroutine



end module mod_tecio_interface
