!=============================================================================!
!=============================================================================!
!!*                                 HW6                                     *!!
!!                        Single Particle Motion                             !!
!!                    Euler and Leapfrog Timstepping                         !!
!!                                                                           !!
!!                PHYS:5905 Numerical Simulation of Plasmas                  !!
!!                         Professor Gregory Howes                           !!
!!                              20 FEB 2019                                  !!
!!*                                                                         *!!
!=============================================================================!
!=============================================================================!
! NOTES:
! SPM: Driver for Single Particle Motion in given E and B fields
! Dimensionless equations:
!    dx'/dt' =v'
!    dv'/dt'= E' + v' x B'
! Normalizations:
!    x'=x/r_L
!    v'=v'/vperp
!    t'= Omega t
!    B'=B/B_0
!    E'=E/(vperp B_0)
!    Relation:   r_L=vperp/Omega
! Numerical Method choices are made by setting variable (meth):
!    Possible selections are:
!    1 Euler's method
!    2 Leapfrog method
! Input x0, v0, B0, E0,
! Compute x and v as a function of time

!------------------------------------------------------------------------------
program hw6
  !Import Module Information
  use hw6_params, only: meth, nn
  use hw6_params, only: x0,v0,tspan
  use hw6_funcs, only: lorentz,magfield,elecfield
  implicit none

  !Declare variables
  real :: pi            
  real, dimension(1:3,0:nn) :: x,v     !Position and velocity arrays
  real, dimension(0:nn) :: t           !Time array
  real :: dt                           !Timestep
  integer :: i
  real, dimension(1:3) :: dxdt,dvdt    !Time derivatives
  real, dimension(1:3) :: B            !Vector magnetic field
  real, dimension(1:3) :: E            !Vector magnetic field
  !Analytical solution
  real, dimension(1:3,0:nn) :: xt      !Theory Position array
  real :: ve                           !E x B velocity
  real :: Bmag                         !Magnitude of B, |B|
  real :: rl                           !Larmor radius at t=0
  real :: err                          !Error at end of integration
  
  !Calculate pi
  pi=4.0*atan(1.0)

  !Set initial and final times for the simulation
  tspan(1)=0.
  tspan(2)=10.* 2.*pi
  !Calculate timestep
  dt=(tspan(2)-tspan(1))/real(nn)
  !Set time variable
  do i=0,nn
     t(i)=real(i)*dt   
  enddo

  !Initialize x and v with initial conditions x0 and v0
  x(1:3,0)=x0(1:3)
  v(1:3,0)=v0(1:3)

  select case(meth)
  case(1) !Euler's method (first-order)
     write(*,'(a)')'Euler method'
     ! Loop over timesteps
     do i=1,nn
        !Calculate time rate of change of x and v
        call lorentz(x(:,i-1),v(:,i-1),t(i-1),dxdt,dvdt)
        !Update x and v using dxdt and dvdt
        x(1:3,i)= x(1:3,i-1) + dt*dxdt(1:3)
        v(1:3,i)= v(1:3,i-1) + dt*dvdt(1:3)
     enddo !End Euler timestep loop

  case(2) !Leapfrog method (second-order)
    write(*,'(a)')'Leapfrog Method'
      !This is HW#6b Prob 2
      call lorentz(x(:, 0), v(:, 0), t(0), dxdt, dvdt)
      x(1:3, 1) = x(1:3, 0) + dt*dxdt(1:3)
      v(1:3, 1) = v(1:3, 0) + dt*dvdt(1:3)
      do i=1,nn
        call lorentz(x(:, i), v(:, i), t(i), dxdt, dvdt)
         x(1:3, i+1) = x(1:3, i-1) + 2*dt*dxdt(1:3)
         v(1:3, i+1) = v(1:3, i-1) + 2*dt*dvdt(1:3)
      enddo
       
  case default
     write(*,*)'ERR: Timestepping algorithm parameter meth unrecognized'
     stop
  end select

  !Compute Analytic Solution: For case with x0=y0, v0=vx0, E=Ey, B=Bz
  call analytical(t,xt)

  ! Compute error as norm of numerical vs. theoretical line
  ! NOTE: function norm2 gives magnitude (2-norm) of argument
  err=norm2(x(1:3,nn)-xt(1:3,nn))

  write(*,'(a,i10,a,es12.4)')"Nstep: ",nn,"        Error in solution: ",err
  
  !Output data to file
  open(27,file='hw6_2.dat',status='replace')
  do i=0,nn
     write(27,'(f12.4,9es12.4)')t(i),x(1:3,i),v(1:3,i),xt(1:3,i)
  enddo
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  subroutine analytical(t1,xt1)
    !NOTE: This routine solves for the total motion (Larmor plus E x B drift)
    !       ASSUMING case with x0=y0, v0=vx0, and constant fields E=Ey, B=Bz
    use hw6_params, only: nn
    use hw6_params, only: x0
    use hw6_funcs, only: magfield,elecfield
    implicit none
    !Passed
    real, intent(in), dimension(0:nn) :: t1          !Time
    real, intent(out), dimension(1:3,0:nn) :: xt1  !Theory Position array
    !Local
    real :: ve                           !E x B velocity
    real :: Bmag                         !Magnitude of B, |B|
    real :: rl                           !Larmor radius at t=0
    real, dimension(1:3) :: B   !Vector magnetic field
    real, dimension(1:3) :: E   !Vector magnetic field
    integer :: j
   
    !Get electric and magnetic fields at x0 and t=0
    E=elecfield(x0,t1(0))
    B=magfield(x0,t1(0))
    Bmag=sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))

    !Subtract ve=ExB velocity from initial perp velocity
    !NOTE: In dimensionless units,  ve' = E' B'
    ve=E(2)*B(3)/Bmag**2.
    !NOTE: In dimensionless units,  r_l/(vperp/Omega) = (v0'-ve')/B'
    rl=(v0(1)-ve)/(B(3));
 
    !NOTE: In dimensionless units,  t' = Omega t
    do j=0,nn
       xt1(1,j)=rl*sin(t1(j)) + (E(2)*B(3)-E(3)*B(2))/Bmag**2.*t1(j) + x0(1)
       xt1(2,j)=rl*cos(t1(j)) + (E(3)*B(1)-E(1)*B(3))/Bmag**2.*t1(j) + (x0(2)-rl)
    enddo
    
  end subroutine analytical
!------------------------------------------------------------------------------
end program hw6
