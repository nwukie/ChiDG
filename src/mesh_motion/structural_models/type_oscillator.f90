module type_oscillator_model
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO, ONE, TWO, PI
    use mod_rigid_body_motion, only: rigid_body_motion_disp_old, rigid_body_motion_disp_new, rigid_body_motion_vel
    implicit none

    type    :: oscillator_model_t
        integer(ik) :: oscillator_ID                    !Identifies the current oscillator, allowing for multiple oscillators
        integer(ik) :: structural_time_integrator_ID    !Identifies the structureal time integrator to use

        integer(ik) :: ndof                             !Number of degrees of freedom - 1, 2, or 3

        ! Linearly damped oscillator ODE model
        ! mass*x'' + damping_coeff*x' + stiffness_coeff*x = external_force
        real(rk)    :: mass = ONE
        real(rk)    :: damping_coeff = ONE
        real(rk)    :: stiffness_coeff = ONE              
        real(rk)    :: external_forces(3)               

        real(rk)    :: undamped_angular_frequency
        real(rk)    :: damping_factor
        real(rk)    :: minimum_stable_timestep

        character(:), allocatable :: damping_type

        real(rk)    :: eq_pos(3)                        ! Equilibrium position
        real(rk)    :: pos(3)                         ! Displaced position

        real(rk)    :: disp(2,3)                        ! Displacement
        real(rk)    :: vel(2,3)                         ! Velocity


        real(rk), allocatable, dimension(:,:)   :: history_pos(:,:), history_vel(:,:), history_force(:,:)

    contains
        procedure :: set_external_forces

        procedure :: get_position
        procedure :: get_displacement
        procedure :: get_velocity

        procedure :: update_disp
        procedure :: update_vel
        procedure :: update_oscillator_subcyle_step
        procedure :: update_oscillator_step
        procedure :: store_motion

    end type oscillator_model_t

    subroutine init(self, mass, damping_coeff, stiffness_coeff)
        type(oscillator_model_t), intent(inout)    :: self
        real(rk),intent(in),optional         :: mass
        real(rk),intent(in),optional         :: damping_coeff
        real(rk),intent(in),optional         :: stiffness_coeff

        real(rk) :: tol

        tol = 1.0e-14

        self%disp = ZERO
        self%pos = ZERO
        self%vel = ZERO
        self%force = ZERO

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
        self%undamped_natural_freqency = self%undamped_angular_frequency/(TWO*PI)
        self%damping_ratio = self%damping_coeff/(TWO*sqrt(self%mass*self%stiffness_coeff)

        !
        ! Compute the minimum stable timestep size for the lepfrog algorithm.
        ! dt < 2/ang_freq
        !

        self%minimum_stable_timestep = TWO/(self%undamped_angular_frequency)<F6>

        if (self%damping_ratio > ONE + tol) then
            self%damping_type = 'overdamped'
        else if (self%damping_ratio < ONE-tol) then
            self%damping_type = 'underdamped'
        else
            self%damping_type = 'critically damped'
        end if
    end subroutine init

    subroutine set_external_forces(self, external_forces)
        type(oscillator_model_t)    :: self
        real(rk)                :: external_forces(3)

        self%external_forces = external_forces

    end subroutine set_external_forces

    
    function get_position(self) result(pos)
        type(oscillator_model_t)    :: self

        real(rk)                :: pos(3)

        self%pos = self%eq_pos + self%disp(1,:) 

        pos = self%pos

    end function get_position

    function get_displacement(self) result(disp)
        type(oscillator_model_t)    :: self

        real(rk)                :: disp(3)


        disp = self%disp(1,:)

    end function get_displacement

    function get_velocity(self) result(velocity)
        type(oscillator_model_t)    :: self

        real(rk)                :: velocity(3)


        velocity = self%vel(1,:)

    end function get_position


    subroutine update_disp(self,dt_struct)
        type(oscillator_model_t) :: self
        real(rk)                :: dt_struct

        
        real(rk)                :: dt_struct
        integer(ik)             :: nsteps

        self%disp(2,:) = self%disp(1,:) + dt_struct*self%vel(1,:)
        self%disp(1,:) = self%disp(2,:)

    end subroutine update_disp

    subroutine update_vel(self,dt_struct)
        type(oscillator_model_t) :: self
        real(rk)                :: dt_struct

        
        real(rk)                :: dt_struct
        integer(ik)             :: nsteps

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
        type(oscillator_model_t) :: self
        real(rk)                :: dt_struct
        real(rk)                :: external_forces(3)

        call self%set_external_forces(external_forces)
        call self%update_vel(dt_struct)
        call self%update_disp(dt_struct)

    end subroutine


    subroutine update_oscillator_step(self, dt_fluid, external_forces, t0_in)
        type(oscillator_model_t) :: self
        real(rk)                :: dt_fluid
        real(rk)                :: external_forces(3)

        real(rk)                :: dt_struct

        integer(ik)             :: nsteps


        t0 = t0_in
        t1 = t0_in+dt_fluid
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

        rigid_body_motion_disp_new = self%disp(1,:)
        rigid_body_motion_vel = self%vel(1,:)

    end subroutine update_oscillator_step
end module type_oscillator_model
