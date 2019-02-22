!=============================================================================!
!=============================================================================!
!                        Running Parameters for HW6 Code                      !
!=============================================================================!
!=============================================================================!
!------------------------------------------------------------------------------
! Numerical Method choices are made by setting variable (meth):
!    Possible selections are:
!    1 Euler's method
!    2 Leapfrog method
!------------------------------------------------------------------------------
module hw6_params
  implicit none
  !Parameters
  integer, parameter :: meth=1           !Timestep Algorithm choice
  integer, parameter :: nn=100000        !Number of timesteps
  ! Initial Conditions
  real, dimension(1:3) :: x0=(/ 0., 1., 0. /)    !Initial position
  real, dimension(1:3) :: v0=(/ 1., 0., 0. /)    !Initial velocity
  real, dimension(1:2) :: tspan                  !Initial and final times
  

  public :: meth,nn
  public :: x0,v0
  public :: tspan

!contains
!------------------------------------------------------------------------------
end module hw6_params
