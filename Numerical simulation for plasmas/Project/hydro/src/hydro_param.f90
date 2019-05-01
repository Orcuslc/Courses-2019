!=============================================================================!
!=============================================================================!
!                         Cascade Parameter Module
!=============================================================================!
!=============================================================================!
! Reads in parameters from the input namelist
!------------------------------------------------------------------------------
!                             Copyright,  2007
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_param
  implicit none
  private

  !Parameters
  character(500) :: runname   !Root of input file name
  integer :: nx,ny            !Number of gridpoints
  real :: lx,ly               !Simulation box side lengths
  real :: dt                  !Timestep
  integer :: nsteps           !Total number of timesteps to run
  integer :: nsave            !Save output at each nsave steps
  logical :: save0            !T=Save initial conditions
  integer :: init             !Initialization, 0 = sound wave initial conditions 
  real :: pi                  !PI=3.14159...
  real :: cfl                 !CFL Timestep Stability criterion

  public :: init_param
  public :: runname
  public :: nx,ny,lx,ly
  public :: dt,nsteps,nsave,save0,cfl
  public :: init
  public :: pi

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2007
!------------------------------------------------------------------------------
! Read plasma namelist from input file to get parameters
  subroutine init_param
    use hydro_mpi_var, only: ierror,proc0
    implicit none
    include 'mpif.h'
    integer :: in_unit        !Unit file number
    character(50) :: f1                     !File name
    
    namelist /system/ nx,ny,lx,ly,dt,nsteps,nsave,save0,init
    
    if (proc0) then        
       !Read parameters on proc0 only-------------------------------------------
       !Get the root of the input file and return as runname
       call get_runname(runname)
       
       !Open the input file
       in_unit=27
       f1=trim(runname)//".in" !Construct input file name
       open (unit=in_unit,file=f1,status="old",action="read")
       
       !Set default parameter values
       nx = 256          !Number of gridpoints in x
       ny = 256          !Number of gridpoints in y
       lx = 1.0          !Length of simulation domain in x
       ly = 1.0          !Length of simulation domain in y
       dt = 0.01         !Timestep (s)
       nsteps = 1000     !Total number of timesteps to run
       nsave = 100       !Save output at each nsave steps
       save0=.false.     !F= do not save initial conditions
       init = 0          ! 0 = sound wave initial conditions 
       
       !Read parameters specified in input file namelist &system
       read (unit=in_unit, nml=system)
       close(in_unit)
       !Summarize Plasma parameters to output
       if (.true.) then
          write(*,'(a)')'====SYSTEM PARAMETERS========================================'
          write(*,'(a,i8)')      'nx=           ',nx
          write(*,'(a,i8)')      'ny=           ',ny
          write(*,'(a,es12.4)')  'lx=           ',lx
          write(*,'(a,es12.4)')  'ly=           ',ly
          write(*,'(a,es12.4)')  'dt=           ',dt
          write(*,'(a,i8)')      'nsteps=       ',nsteps
          write(*,'(a,i8)')      'nsave=        ',nsave
          write(*,'(a,l8)')      'save0=        ',save0
          write(*,'(a,i8)')      'init=         ',init
          write(*,'(a)')'============================================================='
       endif
    endif

    !Broadcast parameters to all procs--------------------------------------
    call mpi_bcast (runname, len(runname), MPI_CHARACTER, 0, mpi_comm_world, ierror)
    call mpi_bcast(nx, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(ny, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(lx, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(ly, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(dt, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
    call mpi_bcast(nsteps, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(nsave, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)
    call mpi_bcast(save0, 1, MPI_LOGICAL, 0, mpi_comm_world, ierror)
    call mpi_bcast(init, 1, MPI_INTEGER, 0, mpi_comm_world, ierror)

  end subroutine init_param
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Get runname for output files from input argument
  subroutine get_runname(runname)
    implicit none
    integer       :: l
    character(500) :: arg
    character(500), intent(out) :: runname

    !Get the first argument of the program execution command
    call getarg(1,arg)

    !Check if this is the input file and trim .in extension to get runname
    l = len_trim (arg)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(1:l-3)
    end if
  end subroutine get_runname
!------------------------------------------------------------------------------
end module hydro_param
  
