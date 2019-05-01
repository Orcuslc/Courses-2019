!=============================================================================!
!=============================================================================!
!                           HYDRO Computational Grid Set Up
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_grid
  private
  include 'mpif.h'
  
  integer :: lnx0,lnx1       !Limits of x index on local processor
  integer, dimension(:), pointer :: iproc_ix   !Value of iproc with the data for row ix
  integer, dimension(:), pointer :: tmp        !Temp array for all reduce


  public :: setup_grid,finalize_grid,iproc_ix
  public :: lnx0,lnx1
  

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Set up computational grid
! NOTE: The lower limit of position in x and y for the box is ASSUMED to be 0
  subroutine setup_grid
    use hydro_var, only: qq,dqdt,dqdt1,dqdt2,imax,rx,ry
    use hydro_param, only: nx,ny,lx,ly
    use hydro_mpi_var, only: nproc,iproc,ierror,proc0
    implicit none
    integer :: i            !Counter
    integer :: dnx          !Increment in nx
    
    
    !Allocate position variables for grid------------------------------------------
    !NOTE: These are calculated beyond the grid to allow for IC initialization
    allocate(rx(-1:nx)); rx(:)=0.
    allocate(ry(-1:ny)); ry(:)=0.

    !$omp parallel do
    do i=-1,nx
       rx(i)=real(i)/real(nx)*lx
    enddo
    !$omp end parallel do
    !$omp parallel do
    do i=-1,ny
       ry(i)=real(i)/real(ny)*ly
    enddo
    !$omp end parallel do


    !SIMPLE Domain Decomposition in x direction-------------------------------------
    if (mod(nx,nproc) .eq. 0) then
       dnx=nx/nproc
    else
       dnx=int(real(nx)/real(nproc)+1.)
    endif
    !Calculate minimum local x index lnx0 and maximum lnx1
    lnx0=(iproc)*dnx 
    lnx1=min((iproc+1)*dnx-1,nx-1)

    !Construct the variable for determining which processor 
    !    has the data for row ix
    allocate(iproc_ix(0:nx-1)); iproc_ix(:)=0
    !Set value for the local information
    !NOTE: This handles the case that no points are on a processor (lnx0>lnx1)
    !      and will not cause array out of bounds
    !$omp parallel do
    do i=lnx0,lnx1
       iproc_ix(i)=iproc
    enddo
    !$omp end parallel do

    !Gather all of the information on proc0
    allocate(tmp(0:nx-1)); tmp=iproc_ix
    call mpi_allreduce(tmp, iproc_ix, size(iproc_ix), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierror)
    deallocate(tmp)
    
    !Allocate variables to contain the local information ---------------------------
    !NOTE: Boundary values are included for fluid variable information qq in x and y
    allocate(qq(1:imax,lnx0-1:lnx1+1,-1:ny)); qq(:,:,:)=0.
    allocate(dqdt(1:imax,lnx0:lnx1,0:ny-1)); dqdt(:,:,:)=0.
    allocate(dqdt1(1:imax,lnx0:lnx1,0:ny-1)); dqdt1(:,:,:)=0.
    allocate(dqdt2(1:imax,lnx0:lnx1,0:ny-1)); dqdt2(:,:,:)=0.

    
    !DEBUG 
    if (.false.) then
       do i=0,nproc-1
          if (i .eq. iproc ) write(*,'(3(a,i4))') 'iproc= ',iproc,'   lnx0= ',lnx0,'   lnx1= ', lnx1
!         call mpi_barrier(mpi_comm_world)
       enddo
       call mpi_barrier(mpi_comm_world)
!       if (proc0) then
!          do i=0,nx-1
!             write(*,'(2(a,i4))') 'ix= ',i,'   iproc= ',iproc_ix(i)
!          enddo
!       endif
    endif


  end subroutine setup_grid
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Deallocate variables from grid 
  subroutine finalize_grid
    use hydro_var, only: qq,dqdt,dqdt1,dqdt2,imax,rx,ry
    implicit none
    
    !Deallocate variables
    deallocate(rx,ry)
    deallocate(iproc_ix)
    deallocate(qq)
    deallocate(dqdt,dqdt1,dqdt2)

  end subroutine finalize_grid
!------------------------------------------------------------------------------

end module hydro_grid
