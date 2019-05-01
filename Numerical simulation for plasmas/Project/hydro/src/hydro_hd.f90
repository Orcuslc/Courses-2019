!=============================================================================!
!=============================================================================!
!                           HYDRO HYDRODYNAMIC EQUATIONS
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_hd
  implicit none
  private

  public :: dqdt_hd

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Function to calculate the time derivative of fluid variables at one point
!NOTE:  irho = 1             !Density
!       iux  = 2             !x-velocity
!       iuy  = 3             !y-velocity
!       ip   = 4             !Pressure
  function dqdt_hd(q,x,y,gamma)
    use hydro_var, only: imax
    implicit none
    real, dimension(1:imax) :: dqdt_hd     !Time derivatives
    real, dimension(1:,-1:,-1:) :: q   !Fluid Variables
    real, dimension(-1:) :: x,y        !Positions 
    real :: gamma                               !Adiabatic Index
    !Intermediate calculations
    real :: divu             !Div u = du_x/dx + du_y/dy

    !Calculate Divergence of u
    divu=(q(2,1,0)-q(2,-1,0))/(x(1)-x(-1)) + (q(3,0,1)-q(3,0,-1))/(y(1)-y(-1))

    !d rho/dt
    dqdt_hd(1)=-q(2,0,0)* (q(1,1,0)-q(1,-1,0)) / (x(1)-x(-1))   &
               -q(3,0,0)* (q(1,0,1)-q(1,0,-1)) / (y(1)-y(-1))  &
               -q(1,0,0)* divu
    !d ux/dt
    dqdt_hd(2)=-q(2,0,0)* (q(2,1,0)-q(2,-1,0)) / (x(1)-x(-1))   &
               -q(3,0,0)* (q(2,0,1)-q(2,0,-1)) / (y(1)-y(-1))  &
               -(q(4,1,0)-q(4,-1,0)) / (x(1)-x(-1))/ q(1,0,0)
    !d uy/dt
    dqdt_hd(3)=-q(2,0,0)* (q(3,1,0)-q(3,-1,0)) / (x(1)-x(-1))   &
               -q(3,0,0)* (q(3,0,1)-q(3,0,-1)) / (y(1)-y(-1))  &
               -(q(4,0,1)-q(4,0,-1)) / (y(1)-y(-1)) / q(1,0,0)
    !d p/dt
    dqdt_hd(4)=-q(2,0,0)* (q(4,1,0)-q(4,-1,0)) / (x(1)-x(-1))   &
               -q(3,0,0)* (q(4,0,1)-q(4,0,-1)) / (y(1)-y(-1))  &
               -gamma*q(4,0,0)* divu

  end function dqdt_hd
!------------------------------------------------------------------------------


end module hydro_hd
