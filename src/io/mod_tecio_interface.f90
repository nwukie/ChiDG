module mod_tecio_interface
#include<messenger.h>
    use iso_c_binding
    use mod_kinds,      only: rk,ik,rdouble,TEC
    use mod_constants,  only: OUTPUT_RES
    use type_domain,    only: domain_t


    implicit none

#include "tecio.f90"

contains



    !>  This opens a new tecplot binary file for writing
    !!  and initializes the type, and number of data fields to accept
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  title       Title of the dataset.
    !!  @param[in]  variables   Comma-separated list of variables that will be written.
    !!  @param[in]  filename    Name of the file to be written.
    !!  @param[in]  filetype    Indicating grid or solution file.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_tecio_file(title,variables,filename,filetype)
        character(*)    :: title
        character(*)    :: variables
        character(*)    :: filename
        integer(TEC)    :: filetype

        integer(4)      :: tecstat
        character       :: NULLCHAR = char(0)
        integer(TEC)    :: fileformat = 0       ! 0 = .plt         1 = subzone loadable .szplt
        integer(TEC)    :: isdouble   = 1       ! 0 = single prec  1 = double prec
        integer(TEC)    :: debug      = 0       ! 0 = debug off

        tecstat = TECINI142(trim(title)//NULLCHAR,      &
                            trim(variables)//NULLCHAR,  &
                            trim(filename)//NULLCHAR,   &
                            '.'//NULLCHAR,              &
                            fileformat,                 &
                            filetype,                   &
                            debug,                      &
                            isdouble)

        if (tecstat /= 0) call chidg_signal(FATAL,"init_tecio_file: Error in TecIO file initialization.")

    end subroutine init_tecio_file
    !*****************************************************************************************










    !>  This begins a new zone in the current opened file. Must be called
    !!  after init_tecplot_file, because it needs an open binary file.
    !!  If multiple files are open, you can switch between them with the
    !!  TECFIL142 call, as long as the they can be identified by integer values
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]  zonetitle   Name for the zone being initialized.
    !!  @param[in]  mesh        mesh_t containing the mesh description to be initialized.
    !!  @param[in]  timeindex   Integer index of time strand.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init_tecio_zone(zonetitle,domain,timeindex)
            character(*),   intent(in)  :: zonetitle
            type(domain_t), intent(in)  :: domain
            integer(ik),    intent(in)  :: timeindex

            integer(TEC)   :: zonetype                  = 5    ! 0 - ordered, 5 - FEBRICK
            integer(TEC)   :: numpts
            integer(TEC)   :: numelements
            integer(TEC)   :: numfaces                  = 0    ! not used
            integer(TEC)   :: icellmax                  = 0    ! not used
            integer(TEC)   :: jcellmax                  = 0    ! not used
            integer(TEC)   :: kcellmax                  = 0    ! not used
            real(rk)       :: solutiontime              = 0._rk
            integer(TEC)   :: strandid                  = 0
            integer(TEC)   :: parentzone                = 0
            integer(TEC)   :: isblock                   = 1
            integer(TEC)   :: nfconns                   = 0
            integer(TEC)   :: fnmode                    = 0
            integer(TEC)   :: totalnumfacenodes         = 1
            integer(TEC)   :: totalnumbndryfaces        = 1
            integer(TEC)   :: totalnumbndryconnections  = 1
            integer(TEC)   :: passivevars(3)            = 0    ! null = all vars active
            integer(TEC)   :: vallocation(3)            = 1    ! null = all vars node-centered
            integer(TEC)   :: sharvarfrom(3)            = 0    ! null = zones share no data
            integer(TEC)   :: sharconnfrom              = 0

            integer(4)              :: tecstat
            integer(TEC),   pointer :: NullPtr(:) => null()    ! Null pointer array


            solutiontime = real(timeindex,rk)
            strandid = timeindex

            numpts      = (OUTPUT_RES+1)*(OUTPUT_RES+1)*(OUTPUT_RES+1) * domain%nelem
            numelements = (OUTPUT_RES*OUTPUT_RES*OUTPUT_RES) * domain%nelem



            tecstat = TECZNE142(trim(zonetitle)//char(0),   &
                                zonetype,                   &
                                numpts,                     &
                                numelements,                &
                                numfaces,                   &
                                icellmax,                   &
                                jcellmax,                   &
                                kcellmax,                   &
                                real(solutiontime,rdouble), &
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

            if(tecstat /= 0) call chidg_signal(FATAL,"init_tecio_zone: Error in TecIO zone initialization.")

    end subroutine init_tecio_zone
    !****************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine finalize_tecio()
        integer(kind=TEC) :: tecstat

        tecstat = TECEND142()
        if (tecstat /= 0) call chidg_signal(FATAL,"finalize_tecio: Error in TecIO file end.")

    end subroutine finalize_tecio
    !*****************************************************************************************



end module mod_tecio_interface
