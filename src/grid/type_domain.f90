module type_domain
    use mod_types,          only: rk,ik, TEC
    use mod_constants,      only: HALF, TWO, ONE, DIAG, NFACES, &
                                  XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use mod_tecio,          only: init_tecio_file, init_tecio_zone, finalize_tecio
    use mod_io,             only: output_res

    use type_mesh,          only: mesh_t
    use type_point,         only: point_t
    use atype_equationset,  only: equationset_t
    use type_bcset,         only: bcset_t
    implicit none

    include '../../tecio/tecio.f90'
    private

    !====================================================
    type, public :: domain_t
        type(mesh_t)                    :: mesh
        type(bcset_t)                   :: bcs
        class(equationset_t), pointer   :: eqnset


    contains
        procedure       :: init
        procedure       :: write
        procedure       :: set_qref
        procedure       :: set_rhsref
        procedure       :: filter_solution
        procedure       :: compute_rhs
        procedure       :: compute_lhs
        procedure       :: zero_working_vectors

        final           :: destructor

    end type domain_t
    !=====================================================
contains
    
    subroutine init(self,eqnset,nterms_sol,nterms_mesh,points)
        class(domain_t),        intent(inout)       :: self
        class(equationset_t),   intent(in), target  :: eqnset
        type(point_t),          intent(in)          :: points(:,:,:)
        integer(kind=ik),       intent(in)          :: nterms_sol
        integer(kind=ik),       intent(in)          :: nterms_mesh

        ! Set domain equation set
        self%eqnset => eqnset

        ! Initialize mesh
        call self%mesh%init(eqnset%neqns,nterms_sol,nterms_mesh,points)

        ! Initialize boundary conditions
        call self%bcs%init(self%eqnset,self%mesh)

    end subroutine



    subroutine write(self,writetype,filename,timeindex)
        class(domain_t),  intent(in) :: self
        character(len=*), intent(in) :: writetype
        character(len=*), intent(in) :: filename
        integer(kind=ik), intent(in) :: timeindex

        integer(kind=ik)        :: uind

        integer(kind=ik)        :: nelem_xi, nelem_eta, nelem_zeta
        integer(kind=ik)        :: ielem_xi, ielem_eta, ielem_zeta
        integer(kind=ik)        :: npts_xi,  npts_eta,  npts_zeta
        integer(kind=ik)        :: ipt_xi,   ipt_eta,   ipt_zeta
        integer(kind=ik)        :: xilim,    etalim,    zetalim
        integer(kind=ik)        :: npts
        integer(kind=4)         :: tecstat

        real(kind=rk)           :: xval(1),yval(1),zval(1),val(1)
        real(kind=TEC)          :: xeq(1),yeq(1),zeq(1),valeq(1)
        equivalence                (xeq(1),   xval(1))
        equivalence                (yeq(1),   yval(1))
        equivalence                (zeq(1),   zval(1))
        equivalence                (valeq(1), val(1))
        real(kind=rk)           :: xi,eta,zeta

        ! Store element indices for current block
        nelem_xi   = self%mesh%nelem_xi
        nelem_eta  = self%mesh%nelem_eta
        nelem_zeta = self%mesh%nelem_zeta

        ! using (output_res+1) so that the skip number used in tecplot to
        ! correctly display the element surfaces is the same as the number
        ! specified in the input file
        npts = output_res+1


        !=====================================================
        !
        !   Write mesh routine
        !
        !=====================================================
        if (writetype == 'mesh') then
            ! Initialize TECIO binary for mesh file
            call init_tecio_file('meshfile','X Y Z',filename,1)
            call init_tecio_zone('meshzone',self%mesh,0,1)



            ! Write x-coordinate
            do ielem_zeta = 1,nelem_zeta
                if (ielem_zeta == nelem_zeta) then
                    zetalim = npts
                else
                    zetalim = npts-1
                end if
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        if (ielem_eta == nelem_eta) then
                            etalim = npts
                        else
                            etalim = npts-1
                        end if
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                if (ielem_xi == nelem_xi) then
                                    xilim = npts
                                else
                                    xilim = npts-1
                                end if

                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get x-mesh point
                                    xval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%x(xi,eta,zeta)
                                    tecstat = TECDAT142(1,xeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

            ! Write y-coordinate
            do ielem_zeta = 1,nelem_zeta
                if (ielem_zeta == nelem_zeta) then
                    zetalim = npts
                else
                    zetalim = npts-1
                end if
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        if (ielem_eta == nelem_eta) then
                            etalim = npts
                        else
                            etalim = npts-1
                        end if

                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                if (ielem_xi == nelem_xi) then
                                    xilim = npts
                                else
                                    xilim = npts-1
                                end if

                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get y-mesh point
                                    yval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%y(xi,eta,zeta)
                                    tecstat = TECDAT142(1,yeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do


            ! Write z-coordinate
            do ielem_zeta = 1,nelem_zeta
                if (ielem_zeta == nelem_zeta) then
                    zetalim = npts
                else
                    zetalim = npts-1
                end if

                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        if (ielem_eta == nelem_eta) then
                            etalim = npts
                        else
                            etalim = npts-1
                        end if

                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                if (ielem_xi == nelem_xi) then
                                    xilim = npts
                                else
                                    xilim = npts-1
                                end if

                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    zval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%z(xi,eta,zeta)
                                    tecstat = TECDAT142(1,zeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

            call finalize_tecio()


        !=====================================================
        !
        !   Write solution routine
        !
        !=====================================================
        elseif (writetype == 'soln') then
            ! Initialize TECIO binary for solution file
            call init_tecio_file('solnfile', 'X Y Z RHO RHOU RHOV RHOW RHOE',filename,0)
!            call init_tecio_file('solnfile', 'X Y Z RHO',filename,0)

            call init_tecio_zone('solnzone',self%mesh,1,timeindex)

            xilim   = npts
            etalim  = npts
            zetalim = npts

!            uind = self%eqnset%get_var("u")

           ! Write x-coordinate
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get x-mesh point
                                    xval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%x(xi,eta,zeta)
                                    tecstat = TECDAT142(1,xeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

            ! Write y-coordinate
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get y-mesh point
                                    yval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%y(xi,eta,zeta)
                                    tecstat = TECDAT142(1,yeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do


            ! Write z-coordinate
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    zval = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%z(xi,eta,zeta)
                                    tecstat = TECDAT142(1,zeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do


            ! Write rho-variable
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    val = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%solution_point(1,xi,eta,zeta)
!                                    print*, val
!                                    print*, valeq
!                                    val = 1.0_rk
!                                    val = val + 1.e-16
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do


            ! Write rhou-variable
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    val = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%solution_point(2,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do



            ! Write rhov-variable
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    val = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%solution_point(3,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

            ! Write rhow-variable
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    val = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%solution_point(4,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do

            ! Write rhoE-variable
            do ielem_zeta = 1,nelem_zeta
                do ipt_zeta = 1,zetalim
                    zeta = (((real(ipt_zeta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                    do ielem_eta = 1,nelem_eta
                        do ipt_eta = 1,etalim
                            eta = (((real(ipt_eta,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                            do ielem_xi = 1,nelem_xi
                                do ipt_xi = 1,xilim
                                    xi = (((real(ipt_xi,rk)-ONE)/(real(npts,rk)-ONE)) - HALF)*TWO

                                    ! Get z-mesh point
                                    val = self%mesh%elements(ielem_xi,ielem_eta,ielem_zeta)%solution_point(5,xi,eta,zeta)
                                    tecstat = TECDAT142(1,valeq,1)

                                end do
                            end do

                        end do
                    end do

                end do
            end do


            call finalize_tecio()

        else
            stop "Error: current available write types are - 'mesh'"
        end if





    end subroutine write



    !===============================================
    !
    !   Set reference solution
    !
    !===============================================
    subroutine set_qref(self)
        class(domain_t),    intent(inout)   :: self

        integer(kind=ik)    :: ixi,ieta,izeta,ieq

        do izeta = 1,self%mesh%nelem_zeta
            do ieta = 1,self%mesh%nelem_eta
                do ixi = 1,self%mesh%nelem_xi

                    associate (elem => self%mesh%elements(ixi,ieta,izeta))

                    ! COPY REF SOLUTION TO WORKING ARRAYS
                    elem%q_ref%vals = elem%q%vals

                    end associate

                end do
            end do
        end do

    end subroutine



    !===============================================
    !
    !   Set reference RHS
    !
    !===============================================
    subroutine set_rhsref(self)
        class(domain_t),    intent(inout)   :: self

        integer(kind=ik)    :: ixi,ieta,izeta,ieq

        do izeta = 1,self%mesh%nelem_zeta
            do ieta = 1,self%mesh%nelem_eta
                do ixi = 1,self%mesh%nelem_xi

                    associate (elem => self%mesh%elements(ixi,ieta,izeta))

                    ! COPY REF SOLUTION TO WORKING ARRAYS
                    elem%rhs_ref%vals = elem%rhs%vals

                    end associate

                end do
            end do
        end do

    end subroutine




    !===============================================
    !
    !   Compute RHS spatial integrals
    !
    !===============================================
    subroutine compute_rhs(self,mask)
        class(domain_t),    intent(inout)   :: self
        type(mask_t),       intent(in)      :: mask

        ! COMPUTE BOUNDARY SOLUTION AND FLUX
        call self%bcs%compute_bc_flux()


        ! COMPUTE LIFT MODES
!        call self%eqnset%compute_lift_all(self%mesh)

        ! COMPUTE BOUNDARY AND VOLUME FLUXES
        call self%eqnset%compute_boundary_flux_all(self%mesh,mask)
        call self%eqnset%compute_volume_flux(self%mesh,mask)
        call self%eqnset%compute_volume_source(self%mesh,mask)


    end subroutine





    !==================================================
    !
    !   Reset all temp RHS and lift data to zero
    !
    !==================================================
    subroutine zero_working_vectors(self)
        class(domain_t),    intent(inout)   :: self


        integer(kind=ik)    :: ixi,ieta,izeta,iface


        ! Perturb a given dof by checking it out of q_ref, perturb, and set to q_m
        do izeta = 1,self%mesh%nelem_zeta
            do ieta = 1,self%mesh%nelem_eta
                do ixi = 1,self%mesh%nelem_xi

                    associate (elem    => self%mesh%elements(ixi,ieta,izeta))

                        elem%rhs%vals = 0._rk
                        elem%lift_g%vals = 0._rk

                        do iface = 1,NFACES
                            elem%lift_l(iface)%vals = 0._rk
                        end do

                    end associate

                end do
            end do
        end do


    end subroutine



    !===============================================
    !
    !   Compute LHS
    !
    !===============================================
    subroutine compute_lhs(self)
        use mod_inv,    only: inv
        class(domain_t), intent(inout)  :: self

        type(mask_t)        :: mask
        integer(kind=ik)    :: i,j,ixi,ieta,izeta
        integer(kind=ik)    :: idof,ieqn,offset
        real(kind=rk)       :: spert,lpert,pert
        logical             :: isBoundary
        real(kind=rk),SAVE  :: pseudotime = 1.e-8_rk

        spert = 3.e-8_rk
        lpert = 3.e-8_rk

        ! Set mask to include all elements
        call mask%set_all(self%mesh)

        ! COMPUTE RHS_REF
        call self%bcs%compute_bc_solution()
        call self%compute_rhs(mask)
        call self%set_rhsref()
        call self%set_qref()
        call self%zero_working_vectors()


        !========================================================================
        !
        !                       Compute Linearization Blocks
        !
        !           - Loop through elements and compute linearization via
        !           - finite differencing of the residual vector
        !
        !========================================================================

        ! Perturb a given dof by checking it out of q_ref, perturb, and set to q
        do izeta = 1,self%mesh%nelem_zeta
            do ieta = 1,self%mesh%nelem_eta
                do ixi = 1,self%mesh%nelem_xi

                    ! Check if boundary element
                    if (ixi   == 1 .or. ixi   == self%mesh%nelem_xi    .or. &
                        ieta  == 1 .or. ieta  == self%mesh%nelem_eta   .or. &
                        izeta == 1 .or. izeta == self%mesh%nelem_zeta) then
                        isBoundary = .true.
                    else
                        isBoundary = .false.
                    end if

                    ! Set mask to operate only on the current element
                    call mask%set_elem(ixi,ieta,izeta)



                    ! Associate current element and working variables
                    element: associate (elem    => self%mesh%elements(ixi,ieta,izeta),         &
                                        rhs_ref => self%mesh%elements(ixi,ieta,izeta)%rhs_ref, &
                                        rhs     => self%mesh%elements(ixi,ieta,izeta)%rhs)

                    !==================================================================
                    !
                    !                        Block DIAGONAL
                    !
                    !==================================================================
                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(DIAG), &
                                          elem_p  => self%mesh%elements(ixi,ieta,izeta))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol

                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs
                            rhs_ref%vals = rhs%vals

                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)






                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn

                    block%vals = inv(block%vals)


                    end associate direction

                    !==================================================================
                    !
                    !                        Block XI_MIN
                    !
                    !==================================================================
                    ! If xi_max boundary then neighbor does not exist
                    if (ixi /= 1) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(XI_MIN), &
                                          elem_p  => self%mesh%elements(ixi-1,ieta,izeta))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol

                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)

                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)




                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn



                    end associate direction

                    end if

                    !==================================================================
                    !
                    !                        Block XI_MAX
                    !
                    !==================================================================
                    ! If xi_max boundary then neighbor does not exist
                    if (ixi /= self%mesh%nelem_xi) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(XI_MAX), &
                                          elem_p  => self%mesh%elements(ixi+1,ieta,izeta))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol

                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)






                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn



                    end associate direction


                    end if

                    !==================================================================
                    !
                    !                        Block ETA_MIN
                    !
                    !==================================================================
                    ! If xi_max boundary then neighbor does not exist
                    if (ieta /= 1) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(ETA_MIN), &
                                          elem_p  => self%mesh%elements(ixi,ieta-1,izeta))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol


                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)






                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn



                    end associate direction

                    end if
                    !==================================================================
                    !
                    !                        Block ETA_MAX
                    !
                    !==================================================================
                    ! If ieta boundary then neighbor does not exist
                    if (ieta /= self%mesh%nelem_eta) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(ETA_MAX), &
                                          elem_p  => self%mesh%elements(ixi,ieta+1,izeta))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol


                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)


                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn



                    end associate direction
                    end if

                    !==================================================================
                    !
                    !                        Block ZETA_MIN
                    !
                    !==================================================================
                    ! If xi_max boundary then neighbor does not exist
                    if (izeta /= 1) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(ZETA_MIN), &
                                          elem_p  => self%mesh%elements(ixi,ieta,izeta-1))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol


                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)






                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn




                    end associate direction
                    end if
                    !==================================================================
                    !
                    !                        Block ZETA_MAX
                    !
                    !==================================================================
                    ! If xi_max boundary then neighbor does not exist
                    if (izeta /= self%mesh%nelem_zeta) then


                    ! Associate correct block and element to be perturbed
                    direction: associate (block   => self%mesh%elements(ixi,ieta,izeta)%linrow%blocks(ZETA_MAX), &
                                          elem_p  => self%mesh%elements(ixi,ieta,izeta+1))

                    pert = spert
                    do ieqn = 1,elem%neqns
                        if (ieqn == 5) then
                            pert = lpert
                        end if
                        ! loop through dof's for each equation
                        do idof = 1,elem%nterms_sol


                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) + pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Store rhs in rhs_ref
                            rhs_ref%vals = rhs%vals
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
                            ! Reset rhs vector
                            call self%zero_working_vectors()



                            ! Set perturbed variable
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn) - pert

                            ! Update boundary conditions
                            if (isBoundary) then
                                call self%bcs%compute_bc_solution()
                            end if

                            ! Compute RHS
                            call self%compute_rhs(mask)
                            ! Replace perturbed solution with reference value
                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)






                            ! Compute linearization vector
                            offset = elem%nterms_sol*(ieqn-1)
                            block%vals(:,idof + offset) = -reshape(((rhs%vals - rhs_ref%vals)/(2._rk*pert)), shape(block%vals(:,ieqn)))

!                            ! Replace perturbed solution with reference value
!                            elem_p%q%vals(idof,ieqn) = elem_p%q_ref%vals(idof,ieqn)
!
!
                            ! Reset rhs vector
                            call self%zero_working_vectors()
                        end do ! idof
                    end do ! ieqn



                    end associate direction
                    end if



!                    print*, "DIAG"
!                    print*, elem%linrow%blocks(DIAG)%vals
!
!                    print*, "XI_MIN"
!                    print*, elem%linrow%blocks(XI_MIN)%vals
!
!                    print*, "XI_MAX"
!                    print*, elem%linrow%blocks(XI_MAX)%vals
!
!                    print*, "ETA_MIN"
!                    print*, elem%linrow%blocks(ETA_MIN)%vals
!
!                    print*, "ETA_MAX"
!                    print*, elem%linrow%blocks(ETA_MAX)%vals
!
!                    print*, "ZETA_MIN"
!                    print*, elem%linrow%blocks(ZETA_MIN)%vals
!
!                    print*, "ZETA_MAX"
!                    print*, elem%linrow%blocks(ZETA_MAX)%vals
!
!                    read(*,*)
                    !=========================================================

                    end associate element

                end do
            end do
        end do





    end subroutine




    !===============================================
    !
    !   Set temporary working arrays
    !
    !===============================================
    subroutine filter_solution(self)
        use mod_filter, only: filter
        class(domain_t),    intent(inout)   :: self

        integer(kind=ik)    :: ixi,ieta,izeta,ieq

        do izeta = 1,self%mesh%nelem_zeta
            do ieta = 1,self%mesh%nelem_eta
                do ixi = 1,self%mesh%nelem_xi

                    associate (elem => self%mesh%elements(ixi,ieta,izeta))
!                        do ieq = 1,self%eqnset%neqns    ! Proabably could just do entire matrix operations, don't need loop

                            ! COPY REF SOLUTION TO WORKING ARRAYS
                            call filter(elem%q%vals)
                            call filter(elem%q_ref%vals)
!                            elem%q%vals(:,ieq) = elem%q_ref%vals(:,ieq)

!                        end do
                    end associate

                end do
            end do
        end do

    end subroutine




    !===============================================
    !
    !   Destructor
    !
    !===============================================
    subroutine destructor(self)
        type(domain_t), intent(in) :: self
    end subroutine

end module type_domain
