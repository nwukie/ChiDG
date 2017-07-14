module mod_rigid_body_motion
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO
    implicit none

    real(rk) :: rigid_body_t0 = ZERO
    real(rk) :: rigid_body_t1 = ZERO
    real(rk) :: rigid_body_motion_disp_old(3) = ZERO
    real(rk) :: rigid_body_motion_disp_new(3) = ZERO
    real(rk) :: rigid_body_motion_vel(3) = ZERO

end module mod_rigid_body_motion

