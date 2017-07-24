module type_oscillator_model
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO, ONE, TWO, PI
    use mod_rigid_body_motion, only: rigid_body_motion_disp_old, rigid_body_motion_disp_new, &
                                    rigid_body_motion_vel, rigid_body_t0, rigid_body_t1
    use mod_chidg_mpi,      only: IRANK, GLOBAL_MASTER
    implicit none

    type    :: oscillator_model_t
        ! Linearly damped oscillator ODE model
        ! mass*x'' + damping_coeff*x' + stiffness_coeff*x = external_force
        real(rk)    :: mass = ONE
        real(rk)    :: damping_coeff(3) = ONE
        real(rk)    :: stiffness_coeff(3) = ONE              
        real(rk)    :: external_forces(3) = ZERO

        real(rk)    :: undamped_angular_frequency(3)
        real(rk)    :: undamped_natural_frequency(3)
        real(rk)    :: damping_factor(3)
        real(rk)    :: minimum_stable_timestep

        character(:), allocatable :: damping_type

        real(rk)    :: eq_pos(3)                        ! Equilibrium position
        real(rk)    :: pos(3)                         ! Displaced position

        real(rk)    :: disp(2,3)                        ! Displacement
        real(rk)    :: vel(2,3)                         ! Velocity


        real(rk), allocatable, dimension(:,:)   :: history_pos(:,:), history_vel(:,:), history_force(:,:)

    contains

        procedure :: init

        procedure :: set_external_forces
        procedure :: update_disp
        procedure :: update_vel
        procedure :: update_oscillator_step
        procedure :: update_oscillator_subcycle_step

    end type oscillator_model_t

contains

    subroutine init(self, mass, damping_coeff, stiffness_coeff)
        class(oscillator_model_t), intent(inout)    :: self
        real(rk),intent(in),optional         :: mass
        real(rk),intent(in),optional         :: damping_coeff(3)
        real(rk),intent(in),optional         :: stiffness_coeff(3)

        real(rk) :: tol

        tol = 1.0e-14

        self%disp = ZERO
        self%pos = ZERO
        self%vel = ZERO
        self%external_forces = ZERO

        if (present(mass)) then
            self%mass = mass 
        end if

        if (present(damping_coeff)) then
            self%damping_coeff = damping_coeff 
        end if

        if (present(stiffness_coeff)) then
            self%stiffness_coeff = stiffness_coeff 
        end if

        self%undamped_angular_frequency = sqrt(self%stiffness_coeff/self%mass)
        self%undamped_natural_frequency = self%undamped_angular_frequency/(TWO*PI)
        self%damping_factor = self%damping_coeff/(TWO*sqrt(self%mass*self%stiffness_coeff))

        !
        ! Compute the minimum stable timestep size for the lepfrog algorithm.
        ! dt < 2/ang_freq
        !

        self%minimum_stable_timestep = TWO/(maxval(self%undamped_angular_frequency))

        if (maxval(self%damping_factor) > ONE + tol) then
            self%damping_type = 'overdamped'
        else if (maxval(self%damping_factor) < ONE-tol) then
            self%damping_type = 'underdamped'
        else
            self%damping_type = 'critically damped'
        end if
        
        print *, 'Oscillaing cylinder damping type:'
        print *, self%damping_type

        print *, 'Oscillating cylinder minimum stable time step size:'
        print *, self%minimum_stable_timestep
    end subroutine init

    subroutine set_external_forces(self, external_forces)
        class(oscillator_model_t)    :: self
        real(rk)                :: external_forces(3)

        self%external_forces = external_forces

    end subroutine set_external_forces

   
    subroutine update_disp(self,dt_struct)
        class(oscillator_model_t) :: self
        real(rk)                :: dt_struct

        

        self%disp(2,:) = self%disp(1,:) + dt_struct*self%vel(1,:)
        self%disp(1,:) = self%disp(2,:)

    end subroutine update_disp

    subroutine update_vel(self,dt_struct)
        class(oscillator_model_t) :: self
        real(rk)                :: dt_struct

        

        real(rk) :: gam(3), mass
        !
        ! Special version of the Leapfrog algorithm for linear damping
        !

        gam = self%damping_coeff 
        mass = self%mass 

        self%vel(2,:) = ((ONE - gam*dt_struct/(TWO*mass))*self%vel(1,:) + &
            dt_struct*(self%external_forces-self%stiffness_coeff*self%disp(1,:))/mass)/ &
            (ONE + gam*dt_struct/(TWO*mass))

        self%vel(1,:) = self%vel(2,:)

    end subroutine update_vel

    
    subroutine update_oscillator_subcycle_step(self, dt_struct, external_forces)
        class(oscillator_model_t) :: self
        real(rk)                :: dt_struct
        real(rk)                :: external_forces(3)

        call self%set_external_forces(external_forces)
        call self%update_vel(dt_struct)
        call self%update_disp(dt_struct)

    end subroutine update_oscillator_subcycle_step


    subroutine update_oscillator_step(self, dt_fluid, t0_in, external_forces)
        class(oscillator_model_t) :: self
        real(rk)                :: dt_fluid
        real(rk)                :: t0_in
        real(rk)                :: external_forces(3)

        real(rk)                :: dt_struct

        integer(ik)             :: nsteps, istep, max_steps


        integer(ik) :: myunit
        logical :: exists
        !
        ! Write to files
        !
        if (IRANK == GLOBAL_MASTER) then
        inquire(file="viv_output.txt", exist=exists)
        if (exists) then
            open(newunit=myunit, file="viv_output.txt", status="old", position="append",action="write")
        else
            open(newunit=myunit, file="viv_output.txt", status="new",action="write")
        end if
        write(myunit,*) t0_in, rigid_body_motion_disp_new(1), rigid_body_motion_disp_new(2),external_forces(1),external_forces(2)
        close(myunit)
        end if

        !
        ! Perform update
        !

        rigid_body_t0 = t0_in
        rigid_body_t1 = t0_in+dt_fluid
        rigid_body_motion_disp_old = rigid_body_motion_disp_new
        ! Check stability of the initial time step and decrease it until it becomes stable
        dt_struct = dt_fluid
        nsteps = 1
        max_steps = 1000

        do while ((dt_struct> 0.1_rk*self%minimum_stable_timestep) .and. (nsteps < max_steps))
            dt_struct = dt_struct/TWO
            nsteps = nsteps*2
        end do
        

        do istep = 1, nsteps

            call self%update_oscillator_subcycle_step(dt_struct, external_forces)

        end do

        !Force one DOF motion
        self%disp(1,1) = ZERO
        self%disp(1,3) = ZERO
        self%vel(1,1) = ZERO
        self%vel(1,3) = ZERO

        rigid_body_motion_disp_new = self%disp(1,:)
        rigid_body_motion_vel = self%vel(1,:)

        
    end subroutine update_oscillator_step
end module type_oscillator_model
