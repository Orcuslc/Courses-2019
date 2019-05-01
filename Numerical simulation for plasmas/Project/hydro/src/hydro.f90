!------------------------------------------------------------------------------
!                             Copyright,  2009
!                                Greg Howes
!------------------------------------------------------------------------------
!          HYDRO: SAMPLE PARALLEL CODE FOR HYDRODYNAMIC SIMULATION
!------------------------------------------------------------------------------
! NOTES:
program hydro
  use hydro_param, only: init_param,runname,nsteps,nsave,save0,dt,ny
  use hydro_mpi_var, only: nproc,iproc,proc0,ierror
  use hydro_grid, only: setup_grid,finalize_grid,lnx0,lnx1
  use hydro_ic, only: setup_ic,gamma,cs0
  use hydro_dt, only: init_ab3_constants,adv_ab3,init_twosteps
  use hydro_var, only: qq,dqdt,dqdt1,dqdt2,rx,ry,imax,tt
  use hydro_hd, only: dqdt_hd
  use hydro_bcs, only: pass_bcs
  use hydro_diag, only: save_output
  implicit none
  include 'mpif.h'      !Include MPI library variables
  integer :: it         !Step number
  integer :: ix,iy      !Counting indices for grid
  logical :: finish=.false.     !Flag to quit program
  !Timing variables
  real, dimension(1:2) :: t_init=0., t_adv=0., t_comm=0., t_loc=0., t_save=0., t_tot=0.
  character(50) :: f1        !File name

  !Initialize MPI message passing---------------------------------------------
  call mpi_init (ierror)
  call mpi_comm_size (mpi_comm_world, nproc, ierror)
  call mpi_comm_rank (mpi_comm_world, iproc, ierror)
  !Set logical proc0=true if iproc=0
  proc0= (iproc == 0)

  call timing(t_tot)
  call timing(t_init)

  if (proc0) write(*,'(a)')'Starting HYDRO===================================' 
  
  !Check to be sure nproc is even, otherwise shutdown
  if (mod(nproc,2) .ne. 0) then
     if (proc0) then
        write(*,'(a,i6)')'ERROR: Number of processors must be even: nproc= ',nproc
        write(*,'(a)')   'Finishing HYDRO==================================' 
     endif
     call mpi_finalize (ierror)
     stop
  endif

  !Read parameters------------------------------------------------------------
  call init_param
  if (proc0) write(*,'(a,a)') 'Input File name: ',runname

  !TEST
  if (.false.) then
     write(*,'(a,i6)') 'Hello, I am processor number ',iproc
  endif

  !Set up grid----------------------------------------------------------------
  call setup_grid

  !Set up initial conditions--------------------------------------------------
  call setup_ic
  !Output initial conditions
  if (save0) call save_output(0)

  !Initialize Timestepping for AB3 time advance (step=0 and step =1)----------
  if (proc0) write(*,'(a)')'>>> BEGIN: Initialization of timesteps-----------' 
  call init_ab3_constants
  call init_twosteps
  if (proc0) write(*,'(a)')'<<< COMPLETE: Initialization of timesteps--------' 
  
  call timing(t_init)

  if (.true.) then
  !MAIN PROGRAM TIMESTEP LOOP=================================================

  do it=3,nsteps
     call timing(t_adv)
     call timing(t_comm)
     !Pass Boundary Condition information among processors--------------------
     call pass_bcs
     call timing(t_comm)

     call timing(t_loc)
     !Calculate Time derivatives of fluid variables---------------------------
     !Move old timesteps back one
     dqdt2=dqdt1
     dqdt1=dqdt

     !$omp parallel do
     !Loop over all local points to calculate current derivatives dqdt
     do ix=lnx0,lnx1 !Loop over all local ix's
        do iy=0,ny-1
           dqdt(1:imax,ix,iy)=dqdt_hd(qq(1:imax,ix-1:ix+1,iy-1:iy+1), &
                rx(ix-1:ix+1),ry(iy-1:iy+1),gamma)
        enddo
     enddo
     !$omp end parallel do

     !$omp parallel do
     !Time step fluid variables forward using AB3-----------------------------
     do ix=lnx0,lnx1 !Loop over all local ix's
        do iy=0,ny-1
           qq(1:imax,ix,iy)=adv_ab3(qq(1:imax,ix,iy),dt,dqdt(1:imax,ix,iy), &
                dqdt1(1:imax,ix,iy),dqdt2(1:imax,ix,iy))
        enddo
     enddo
     !$omp end parallel do
     
     !Advance time
     tt=tt+dt
     if (proc0) write(*,'(a,i6,a,es12.4)')'Timestep ',it,': t=',tt
     call timing(t_loc)
     
     !Check if we need to save------------------------------------------------
     if (mod(it,nsave) .eq. 0) then
        call timing(t_adv)
        call timing(t_save)
        call save_output(it)
        call timing(t_save)
        call timing(t_adv)
     endif

     !Check if we have signaled to exit the program---------------------------
     call checkstop(finish)

     call timing(t_adv)
     if (finish .or. it .eq. nsteps) then
!        call final_diagnostics
        exit !Jump out of main timestep loop
     endif
     
  enddo !END MAIN PROGRAM LOOP================================================
  endif
  
!  call mpi_barrier(mpi_comm_world)

  if (proc0) write(*,'(a)')'Finishing HYDRO==================================' 
  call finalize_grid

  call timing(t_tot)

  !Output timing statistics
  if (proc0) then
     write(f1,'(a,a)')trim(runname),'.time'
     open (unit=37,file=f1,status="replace",action="write")
     write(37,'(a,a,a,i6,a)')'Output from run    ',trim(runname),'.in      using  ',nproc,' processing cores'
     write(37,'(a,es14.6,a,f10.4,a,f10.3,a)') &
          'Initialization     ',t_init(1), ' sec     (',t_init(1)/60., ' min)     ',t_init(1)/t_tot(1)*100.,' %'
     write(37,'(a,es14.6,a,f10.4,a,f10.3,a)') &
          'Save               ',t_save(1), ' sec     (',t_save(1)/60., ' min)     ',t_save(1)/t_tot(1)*100.,' %'
     write(37,'(a,es14.6,a,f10.4,a,f10.3,a)') &
          'Advance            ',t_adv(1),  ' sec     (',t_adv(1)/60.,  ' min)     ',t_adv(1)/t_tot(1)*100.,' %'
     write(37,'(a,es14.6,a,f10.4,a,f10.3,a)') &
          '    Local          ',t_loc(1),  ' sec     (',t_loc(1)/60.,  ' min)     ',t_loc(1)/t_tot(1)*100.,' %'
     write(37,'(a,es14.6,a,f10.4,a,f10.3,a)') &
          '    Communication  ',t_comm(1), ' sec     (',t_comm(1)/60., ' min)     ',t_comm(1)/t_tot(1)*100.,' %'
     write(37,'(a)')' '
     write(37,'(a,es14.6,a,f10.4,a)') &
          'Total              ',t_tot(1),  ' sec     (',t_tot(1)/60.,  ' min)'
      
    close(37)
  endif

  !Finalize MPI message passing
  call mpi_finalize (ierror)

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Check if the file runname.stop has been created. If so, finish1=T.
  subroutine checkstop(finish1)
    implicit none
    logical, intent (in out) :: finish1
    character (len=300) :: filename
    
    logical :: exit_local
    
    ! If .stop file has appeared, set finish1 flag
    filename=trim(runname)//".stop"
    
    if (proc0) then
       inquire(file=filename,exist=exit_local)
       finish1 = finish1 .or. exit_local
    end if

    call mpi_bcast (finish1, 1, MPI_LOGICAL, 0, mpi_comm_world, ierror)

  end subroutine checkstop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                           Greg Howes, 2009
!------------------------------------------------------------------------------
! Compile timing statistics for run
! At first call, targ(2) is set to current cpu timer in seconds
! At second call, targ(1) is incremented with time since last call, tnew-targ(2)
!    and targ(2) is reset to zero.
  subroutine timing(targ)
    implicit none
    !Passed
    real, dimension(1:2), intent(inout) :: targ ! tsum and told
    !Local
    real :: tnew

    !Get CPU time in seconds
    tnew=mpi_wtime()
!    call cpu_time(tnew)

    !If first call, set targ(2) to tnew
    !Otherwise, add to targ(1) the time difference between this call and last call
    !    and set targ(2)=0    
    if (targ(2) .eq. 0.) then
       targ(2) = tnew
    else
       targ(1)=targ(1)+tnew-targ(2)
       targ(2)=0.
    end if

  end subroutine timing
  !------------------------------------------------------------------------------
end program hydro
