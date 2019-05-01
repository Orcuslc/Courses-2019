!=============================================================================!
!=============================================================================!
!                           HYDRO TIME ADVANCE
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_dt
  use hydro_var, only: imax
  implicit none
  private

  !AB3 Constants
  real :: aa,bb,cc                              !AB3 Constants

  public :: init_ab3_constants,adv_ab3
  public :: init_twosteps
  public :: adv_euler

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Compute  Third-Order Adam's Bashforth constants
  subroutine init_ab3_constants
    implicit none
    
    !Initialize AB3 constants
    aa=23./12.
    bb=-4./3.
    cc=5./12.
    
  end subroutine init_ab3_constants
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! First-Order Euler Time Stepping
!NOTE:  For dq/dt=f, this is given by
!    (q_(n+1) - q_(n))/dt = f_(n)
  function adv_euler(q,dt,f0)
    implicit none
    real, dimension(1:imax) :: adv_euler !Next timestep value
    real, dimension(:) :: q          !Fluid Variables
    real, dimension(:) :: f0         !Time Derivatives (n)
    real :: dt                                !Timestep

    adv_euler(:)=q(:) + dt * f0

  end function adv_euler
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Third-Order Adam's Bashforth Time Stepping
!NOTE:  For dq/dt=f, this is given by
!    (q_(n+1) - q_(n))/dt = aa*f_(n) + bb*f_(n-1) + cc*f_(n-2)
!    aa=23/12, bb=-4/3, cc=5/12
  function adv_ab3(q,dt,f0,f1,f2)
    implicit none
    real, dimension(1:imax) :: adv_ab3              !Next timestep value
    real, dimension(:) :: q           !Fluid Variables
    real, dimension(:) :: f0,f1,f2   !Time Derivatives (n,n-1,n-2)
    real :: dt                                !Timestep

    adv_ab3(:)=q(:) + dt * (aa*f0 + bb*f1 + cc*f2)

  end function adv_ab3
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Initialize the Timestepping scheme by determining dq/dt at n-2 and n-1 steps
! NOTE: This is done using two first-order Euler timesteps using dt/16.,
!       then stepping up to dt and 2.*dt using AB3
  subroutine init_twosteps
    use hydro_var, only: qq,dqdt,dqdt1,dqdt2,imax,tt,rx,ry
    use hydro_param, only: dt,ny
    use hydro_grid, only: lnx0,lnx1
    use hydro_bcs, only: pass_bcs
    use hydro_ic, only: gamma
    use hydro_hd, only: dqdt_hd
    use hydro_mpi_var, only: proc0
    implicit none
    integer :: it1                !Intermediate step number
    real :: dt1                   !Small intermediate timesteps to initialize
    integer :: ix,iy      !Counting indices for grid
    real, dimension(:,:,:), pointer :: tdqdt1   !d/dt (Fluid Variables) (n-1)
    real, dimension(:,:,:), pointer :: tdqdt2   !d/dt (Fluid Variables) (n-2)

    !Allocate temporary timestep information to save 0.*dt and 1.*dt steps
    allocate(tdqdt1(1:imax,lnx0:lnx1,0:ny-1)); tdqdt1(:,:,:)=0.
    allocate(tdqdt2(1:imax,lnx0:lnx1,0:ny-1)); tdqdt2(:,:,:)=0.
    
    !Set smaller timstep to begin 
    dt1=dt/16.

!    if (proc0) write(*,'(a)')'First timestep' 
    !Get dqdt at t=0 and save in tdqdt2----------------------------------------
    !NOTE: BC's need not be passed because initialization has taken care of this
!    if (proc0) write(*,'(a)')'Calc dqdt' 
    do ix=lnx0,lnx1 !Loop over all local ix's
       do iy=0,ny-1
          dqdt(1:imax,ix,iy)=dqdt_hd(qq(1:imax,ix-1:ix+1,iy-1:iy+1), &
               rx(ix-1:ix+1),ry(iy-1:iy+1),gamma)
       enddo
    enddo
    tdqdt2(:,:,:)=dqdt(:,:,:) !Save Rates at step n=0 (dt)
    
    !Take first-order Euler timestep to t=dt1----------------------------------
!    if (proc0) write(*,'(a)')'Euler advance' 
    do ix=lnx0,lnx1 !Loop over all local ix's
       do iy=0,ny-1
          qq(1:imax,ix,iy)=adv_euler(qq(1:imax,ix,iy),dt1,dqdt(1:imax,ix,iy))
       enddo
    enddo
    tt=tt+dt1
    dqdt1(:,:,:)=dqdt(:,:,:)

    !Take first-order Euler timestep to t=2*dt1----------------------------------
!    if (proc0) write(*,'(a)')'Second timestep' 
!    if (proc0) write(*,'(a)')'Pass BCS' 
    call pass_bcs
!    if (proc0) write(*,'(a)')'Advance' 
    !Get dqdt at t=dt1
    do ix=lnx0,lnx1 !Loop over all local ix's
       do iy=0,ny-1
          dqdt(1:imax,ix,iy)=dqdt_hd(qq(1:imax,ix-1:ix+1,iy-1:iy+1), &
               rx(ix-1:ix+1),ry(iy-1:iy+1),gamma)
       enddo
    enddo
    do ix=lnx0,lnx1 !Loop over all local ix's
       do iy=0,ny-1
          qq(1:imax,ix,iy)=adv_euler(qq(1:imax,ix,iy),dt1,dqdt(1:imax,ix,iy))
       enddo
    enddo
    tt=tt+dt1
!    dqdt2(:,:,:)=dqdt1(:,:,:)
!    dqdt1(:,:,:)=dqdt(:,:,:)
    
    !Loop over small initializing time steps to get to t=2*dt------------------
    do it1=3,32
!       if (proc0) write(*,'(a,i2,a)')'Initialization Iteration ',it1,'/16'
       !Pass BCs
       call pass_bcs
       !Move old timesteps back one
       dqdt2=dqdt1
       dqdt1=dqdt
       !Get dqdt
       do ix=lnx0,lnx1 !Loop over all local ix's
          do iy=0,ny-1
             dqdt(1:imax,ix,iy)=dqdt_hd(qq(1:imax,ix-1:ix+1,iy-1:iy+1), &
                  rx(ix-1:ix+1),ry(iy-1:iy+1),gamma)
          enddo
       enddo
       !Timestep forward using AB3
       do ix=lnx0,lnx1 !Loop over all local ix's
          do iy=0,ny-1
             qq(1:imax,ix,iy)=adv_ab3(qq(1:imax,ix,iy),dt,dqdt(1:imax,ix,iy), &
                  dqdt1(1:imax,ix,iy),dqdt2(1:imax,ix,iy))
          enddo
       enddo
       !Advance time
       tt=tt+dt1
       
       !Save value at t=dt=16.*dt1 for AB3 initialization
       if (it1 .eq. 16) then
          if (proc0) write(*,'(a,i2,a,es12.4)')'Initialization Timestep ',it1/16,': t= ',tt
          tdqdt1(:,:,:)=dqdt(:,:,:) !Save Rates at step n=1 (dt)
       endif

       if (it1 .eq. 32 .and. proc0) write(*,'(a,i2,a,es12.4)')'Initialization Timestep ',it1/16,': t= ',tt

    enddo
!    if (proc0) write(*,'(a)')'Done with initial steps' 

    !Put old dqdt values into correct places for AB3 timesteps using timestep dt
    dqdt1=tdqdt2
    dqdt=tdqdt1
    
!    if (proc0) write(*,'(a)')'dqdt copied' 
    
    !Deallocate temporary variables
    deallocate(tdqdt1,tdqdt2)
    
  end subroutine init_twosteps

!------------------------------------------------------------------------------
end module hydro_dt
