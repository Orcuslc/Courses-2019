!=============================================================================!
!=============================================================================!
!                           HYDRO Diagnostic Routines
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module hydro_diag
  private
  
  public :: save_output
  

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Save output
  subroutine save_output(it)
    use hydro_var, only: qq,rx,ry,tt,imax
    use hydro_mpi_var, only: nproc,iproc,ierror,proc0
    use hydro_param, only: runname,nsteps,nsave,nx,ny
    use hydro_grid, only: lnx0,lnx1,iproc_ix
    implicit none
    include 'mpif.h'
    !Passed
    integer, intent(in) :: it  !Timestep counter
    !Local
    integer :: out_unit        !Unit file number
    character(50) :: f1        !File name
    integer :: ix,iy           !Counting indices for grid
    integer :: ip              !iproc counter
    real, dimension(:,:,:), allocatable :: tmp, tmp2 !Structure for data
    integer ::  rproc          !Receiving processors
    integer :: outmeth         !Output Method: 0=reduce
                               !               1=pass to proc0
    real, dimension(:,:,:), allocatable :: tqq !Structure for data
    integer :: rnx0,rnx1     !Limits of x index on remote processor
    integer, dimension (MPI_STATUS_SIZE) :: status

    !Output whole data set to one file
    !First, create the data file using proc0
    out_unit=21
    write(f1,'(a5,a,a5,i4.4)')'data/',trim(runname),'.out.',int(it/nsave)
    if (proc0) write(*,'(a,a,a,es12.4)')'Output saved to file: ',trim(f1),' at t=',tt
!    if (proc0) then
!       open (unit=out_unit,file=f1,status="replace",action="write")
!       close(out_unit)
!    endif
    
    !Choose method for output
    outmeth=0
    select case(outmeth)
    case(0) !Reduce (requires more memory)
       !Loop over each processor to collect data on proc0 for output to a single file
       !NOTE: Gather cannot be used because processors may have different amounts of data
       allocate(tmp(1:imax,0:nx-1,0:ny-1)); tmp=0.
       !Copy data from local processor into tmp 
       tmp(:,lnx0:lnx1,0:ny-1)=qq(:,lnx0:lnx1,0:ny-1)
       !Reduce data from all processors
       rproc=0

       if (proc0) then
          allocate(tmp2(1:imax,0:nx-1,0:ny-1)); tmp2=tmp
          call mpi_reduce(tmp2, tmp, size(tmp2), MPI_REAL8, MPI_SUM, rproc, mpi_comm_world, ierror)
          deallocate(tmp2)
       else
          call mpi_reduce(tmp, tmp, size(tmp), MPI_REAL8, MPI_SUM, rproc, mpi_comm_world, ierror)
       endif
       
       !Output data to file
       if (proc0) then
          open (unit=out_unit,file=f1,status="replace",action="write")
          do ix=0,nx-1
             do iy=0,ny-1
                !             write(*,'(es12.4,2i8,2es12.4,4x,4es14.6,i6)')tt,ix,iy,rx(ix),ry(iy),tmp(:,ix,iy),iproc
                write(out_unit,'(es12.4,2i8,2es12.4,4x,4es14.6,i6)')tt,ix,iy,rx(ix),ry(iy),tmp(:,ix,iy),iproc
             enddo
          enddo
          close(out_unit)
       endif
       
       deallocate(tmp)
       
    case(1) !Pass to proc0 and output (Better for large problems
       !ERROR: THIS SEEMS NOT TO WORK YET!!!!!!
       !Output proc0 information first
       if (proc0) then
          open (unit=out_unit,file=f1,status="old",action="write",position="append")
          !Write data to file
          do ix=lnx0,lnx1
             do iy=0,ny-1
                write(*,'(es12.4,2i8,2es12.4,4x,4es12.4,i6)')tt,ix,iy,rx(ix),ry(iy),qq(:,ix,iy),iproc
                write(out_unit,'(es12.4,2i8,2es12.4,4x,4es12.4,i6)')tt,ix,iy,rx(ix),ry(iy),qq(:,ix,iy),iproc
             enddo
          enddo
       endif
       !Loop over iproc>0, pass data to proc0, and write to file
       do ip=1,nproc-1
          if (proc0) then
             !Get limits on iproc=ip on proc0, rnx0,rnx1
             rnx0=nx+1
             rnx1=-1
             do ix=0,nx-1
                if (iproc_ix(ix) .eq. ip) then
                   if (ix .gt. rnx1) rnx1=ix
                   if (ix .lt. rnx0) rnx0=ix
                endif
             enddo
             !Allocate temporary variable on proc0 to hold passed info
             allocate(tqq(1:imax,rnx0:rnx1,0:ny-1)); tqq=0.
             !Get data from iproc
             write(*,'(a,i4,2(a,i6))')'iproc:',iproc,' rnx0=',rnx0,' rnx1=',rnx1
             call mpi_recv(tqq(:,:,:), size(tqq), MPI_REAL8, ip, 0, mpi_comm_world, status, ierror)
             !Write data to file
             do ix=rnx0,rnx1
                do iy=0,ny-1
                   write(*,'(es12.4,2i8,2es12.4,4x,4es12.4,i6)')tt,ix,iy,rx(ix),ry(iy),tqq(:,ix,iy),ip
                   write(out_unit,'(es12.4,2i8,2es12.4,4x,4es12.4,i6)')tt,ix,iy,rx(ix),ry(iy),tqq(:,ix,iy),ip
                enddo
             enddo
             !Deallocate
             deallocate(tqq)
             
          elseif(ip .eq. iproc) then
             !Pass data to proc0
             write(*,'(a,i4,2(a,i6))')'iproc:',iproc,' lnx0=',lnx0,' lnx1=',lnx1
             call mpi_send(qq(:,lnx0:lnx1,0:ny-1), size(qq(:,lnx0:lnx1,0:ny-1)), MPI_REAL8, 0, &
                  0, mpi_comm_world, ierror)
          endif
          call mpi_barrier(mpi_comm_world)
       enddo
          
       close(out_unit)
    end select
    
  end subroutine save_output
!------------------------------------------------------------------------------
end module hydro_diag
