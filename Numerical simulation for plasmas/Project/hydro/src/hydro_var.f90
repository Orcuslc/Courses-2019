!=============================================================================!
!=============================================================================!
!                           HYDRO HYDRODYNAMIC VARIABLES
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_var
  implicit none
  private

  !Fluid variables
  real, dimension(:,:,:), pointer :: qq      !Fluid Variables
  real, dimension(:,:,:), pointer :: dqdt    !d/dt (Fluid Variables) (n)
  real, dimension(:,:,:), pointer :: dqdt1   !d/dt (Fluid Variables) (n-1)
  real, dimension(:,:,:), pointer :: dqdt2   !d/dt (Fluid Variables) (n-2)
  !Spatial dimensions
  real, dimension(:), pointer :: rx          !x-position
  real, dimension(:), pointer :: ry          !y-position
  !Time
  real :: tt                                 !Time
  
  !Index of variables
  integer, parameter :: irho = 1             !Density
  integer, parameter :: iux  = 2             !x-velocity
  integer, parameter :: iuy  = 3             !y-velocity
  integer, parameter :: ip   = 4             !Pressure
  integer, parameter :: imax = 4             !Number of variables
  
  public :: qq,dqdt,dqdt1,dqdt2
  public :: rx,ry,tt
  public :: irho,iux,iuy,ip,imax

 end module hydro_var
