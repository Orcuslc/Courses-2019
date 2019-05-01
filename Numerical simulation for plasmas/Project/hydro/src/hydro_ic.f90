!=============================================================================!
!=============================================================================!
!                           HYDRO Initial Condition Set Up
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_ic
  private
  include 'mpif.h'
  
  integer :: lnx0,lnx1       !Limits of x index on local processor
  integer, dimension(:), pointer :: iproc_ix   !Value of iproc with the data for row ix
  integer, dimension(:), pointer :: tmp        !Temp array for all reduce
  !Sound Wave Initialization
  real :: rho0            !Equilibrium density (kg/m^3)
  real :: p0              !Equilibrium pressure (N/m^2 = kg/(s^2 m))
  real :: ux0             !Equilibrium velocity in x (m/s)
  real :: uy0             !Equilibrium velocity in y (m/s)
  real :: gamma           !Adiabatic Index
  integer :: nkx0         !Integral Sound Wavenumber in x (Number of wavelengths per lx)
  integer :: nky0         !Integral Sound Wavenumber in y (Number of wavelengths per ly)
  real :: p1              !Pressure perturbation (N/m^2)
  real :: cs0             !Equilibrium Sound Speed


  public :: setup_ic
  public :: rho0,p0,ux0,uy0,gamma,nkx0,nky0,p1,cs0
  

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Set up initial conditions
  subroutine setup_ic
    use hydro_param, only: init,pi
    use hydro_mpi_var, only: ierror,proc0
    implicit none
    
    !Initialize pi
    pi = 4.0*ATAN(1.0)

    !Choose initial conditions
    select case(init)
    case(0) !Linear Sound Wave with amplitude p1 and wavenumber kx and ky
       call init_sound_wave
    case default
       if (proc0) then
          write(*,'(a,i6)')'ERROR: Unknown initial condition: init= ',init
          call mpi_finalize (ierror)
          stop
       endif
    end select

  end subroutine setup_ic
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Set up initial conditions
  subroutine init_sound_wave
    use hydro_param, only: runname,pi,lx,ly,cfl,dt
    use hydro_var, only: rx,ry,qq,irho,iux,iuy,ip,tt
    use hydro_grid, only: lnx0,lnx1
    use hydro_mpi_var, only: ierror,proc0
    implicit none
    integer :: in_unit        !Unit file number
    character(50) :: f1       !File name
    integer :: ix             !Counter
    real :: kx,ky             !Wavenumber x and y components
    real :: k                 !Wavenumber
    real :: dr                !Minimum grid spacing
    real :: dtmax             !Max dt for CFL=1.0

    namelist /sound/ rho0,p0,ux0,uy0,gamma,nkx0,nky0,p1
        
    !Read parameters for the Sound Wave initialization
    if (proc0) then        
       !Read parameters on proc0 only-------------------------------------------
       !Open the input file
       in_unit=27
       f1=trim(runname)//".in" !Construct input file name
       open (unit=in_unit,file=f1,status="old",action="read")
       
       !Set default parameter values
       rho0 = 5.0         !Equilibrium density (kg/m^3)
       p0   = 3.0         !Equilibrium pressure (N/m^2 = kg/(s^2 m))
       ux0  = 0.0         !Equilibrium velocity in x (m/s)
       uy0  = 0.0         !Equilibrium velocity in y (m/s)
       gamma = 1.666667   !Adiabatic Index
       nkx0 = 1           !Integral Sound Wavenumber in x (Number of wavelengths per lx)
       nky0 = 0           !Integral Sound Wavenumber in y (Number of wavelengths per ly)
       p1  = 0.003        !Pressure perturbation (N/m^2)
       
       !Read parameters specified in input file namelist &system
       read (unit=in_unit, nml=sound)
       close(in_unit)
       !Summarize Plasma parameters to output
       if (.true.) then
          write(*,'(a)')'====INITIAL CONDITIONS========================================'
          write(*,'(a,g12.4)')   'rho0=         ',rho0
          write(*,'(a,g12.4)')   'p0=           ',p0
          write(*,'(a,g12.4)')   'ux0=          ',ux0
          write(*,'(a,g12.4)')   'uy0=          ',uy0
          write(*,'(a,g12.4)')   'gamma=        ',gamma
          write(*,'(a,i8)')      'nkx0=         ',nkx0
          write(*,'(a,i8)')      'nky0=         ',nky0
          write(*,'(a,g12.4)')   'p1=           ',p1
          write(*,'(a)')'============================================================='
       endif
    endif

    !Broadcast parameters to all procs--------------------------------------
    call mpi_bcast(rho0, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(p0, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(ux0, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(uy0, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(gamma, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(nkx0, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(nky0, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(p1, 1, MPI_REAL8, 0, mpi_comm_world, ierror)

    !Calculate equilibirum sound speed
    cs0=sqrt(gamma*p0/rho0)

    !Calculate inital conditions for linear sound wave eigenfunction------------
    tt=0.
    kx=real(nkx0)*2.*pi/lx
    ky=real(nky0)*2.*pi/ly
    k=sqrt(kx*kx+ky*ky)
    !NOTE: This calculation on each processor includes the boundary values
    do ix=lnx0-1,lnx1+1
       qq(ip,ix,:)=p0+p1*cos(kx*rx(ix)+ky*ry(:))
       qq(irho,ix,:)=rho0+(qq(ip,ix,:)-p0)/(cs0*cs0)
       qq(iux,ix,:)=ux0+cs0*kx/k*(qq(ip,ix,:)-p0)/(gamma*p0)
       qq(iuy,ix,:)=uy0+cs0*ky/k*(qq(ip,ix,:)-p0)/(gamma*p0)
    enddo
    
    !Calculate CFL timestep stability condition---------------------------------
    dr=min(rx(1)-rx(0),ry(1)-ry(0))
    cfl=(dt*cs0)/dr
    dtmax=dr/cs0
    if (proc0) then
       write(*,'(a,es12.4)')'Hydrodynamic variables initialized at t=',tt
       write(*,'(a)')' '
       write(*,'(a,f10.6,a,f10.6)')'CFL Timestep stability ratio= ',cfl,'   Maximum Timestep dt_max=',dtmax
    endif
    if (cfl .gt. 1.) then
       if (proc0) then
          write(*,'(a)')'ERROR:  Timestep too large for  CFL stability'
          write(*,'(a)')   'Finishing HYDRO==================================' 
       endif
       call mpi_finalize (ierror)
       stop
    endif
       

  end subroutine init_sound_wave
!------------------------------------------------------------------------------

end module hydro_ic
