!=============================================================================!
!=============================================================================!
!                HYDRO Boundary Condition Passing Among Processors
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  NOTE: Boundary Conditions for the whole simulation box are periodic
!------------------------------------------------------------------------------
module hydro_bcs
  implicit none
  private

  public :: pass_bcs

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Pass required information at boundaries of subdomains
  subroutine pass_bcs
    use hydro_mpi_var, only: nproc,iproc,proc0,ierror
    use hydro_var, only: qq
    use hydro_param, only: nx,ny
    use hydro_grid, only: lnx0,lnx1,iproc_ix
   implicit none
   include 'mpif.h'                  !Include MPI library variables
   integer :: sproc, rproc           !Sending and receiving processors
   integer :: sinx,rinx              !ix index to pass first
   integer, dimension (MPI_STATUS_SIZE) :: status

   !Pass periodic information locally in the y-dimension
   qq(:,lnx0:lnx1,-1)=qq(:,lnx0:lnx1,ny-1)
   qq(:,lnx0:lnx1,ny)=qq(:,lnx0:lnx1,0)
   
   !Pass the boundary information to the neighboring processor
   !NOTE: First pass lower boundary, then upper boundary
   !NOTE: mod's below handle wrap around from iproc=(nproc-1) to iproc=0
   !Lower boundary information-------------------------------------------------
   rinx=mod(lnx0-1 + nx, nx) !inx needed by this processor
   sinx=lnx1                 !inx needed by next higher processor 
   !Get processor numbers that have rinx and that need sinx
   rproc=iproc_ix(rinx)
   sproc=mod(iproc+1,nproc)

   !Pass boundary information
   !NOTE:To avoid deadlock, the even iproc always passes first
   if (mod(iproc,2) .eq. 0) then

      call mpi_send(qq(:,sinx,:), size(qq(:,sinx,:)), MPI_REAL8, sproc, 0, mpi_comm_world, ierror)
      call mpi_recv(qq(:,lnx0-1,:), size(qq(:,lnx0-1,:)), MPI_REAL8, rproc, 0, mpi_comm_world, status, ierror)
   else
      call mpi_recv(qq(:,lnx0-1,:), size(qq(:,lnx0-1,:)), MPI_REAL8, rproc, 0, mpi_comm_world, status, ierror)
      call mpi_send(qq(:,sinx,:), size(qq(:,sinx,:)), MPI_REAL8, sproc, 0, mpi_comm_world, ierror)
   endif

   !Upper boundary information-------------------------------------------------
   rinx=mod(lnx1+1, nx)      !inx needed by this processor
   sinx=lnx0                 !inx needed by next higher processor 
   !Get processor numbers that have rinx and that need sinx
   rproc=iproc_ix(rinx)
   sproc=mod(iproc-1+nproc,nproc)
   !Pass boundary information
   !NOTE:To avoid deadlock, the even iproc always passes first
   if (mod(iproc,2) .eq. 0) then
      call mpi_send(qq(:,sinx,:), size(qq(:,sinx,:)), MPI_REAL8, sproc, 0, mpi_comm_world, ierror)
      call mpi_recv(qq(:,lnx1+1,:), size(qq(:,lnx1+1,:)), MPI_REAL8, rproc, 0, mpi_comm_world, status, ierror)
   else
      call mpi_recv(qq(:,lnx1+1,:), size(qq(:,lnx1+1,:)), MPI_REAL8, rproc, 0, mpi_comm_world, status, ierror)
      call mpi_send(qq(:,sinx,:), size(qq(:,sinx,:)), MPI_REAL8, sproc, 0, mpi_comm_world, ierror)
   endif

  end subroutine pass_bcs
!------------------------------------------------------------------------------


end module hydro_bcs
