!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
program cola
! Description:
!  A 3D unstructured finite volume solver for Computational Fluid Dynamics.   
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use fieldManipulation
  use sparse_matrix
  use temperature
  use concentration
  use mhd
  use utils, only: timestamp

  implicit none

  integer :: k, ios
  integer :: iter
  integer :: itimes, itimee
  real(dp):: source
  real :: start, finish, t_start, t_end
  character( len = 9) :: timechar
!                                                                       
!******************************************************************************
!

!  Check command line arguments
  ! narg=command_argument_count()
  call get_command_argument(1,input_file)
  call get_command_argument(2,restart_file)
  call get_command_argument(3,monitor_file)

  ! Open simulation log file or Print.
  if(len_trim(monitor_file).NE.0)then
    open(unit=6,file=monitor_file,iostat=ios)
    if(ios.EQ.0)then
        rewind 6
    endif
  endif

  ! Print logo in the header of monitor file
    write(6,'(a)')' '
    write(6,'(a)')'  Notice:'
    WRITE(6,'(2x,a)') '------------------------------------------------------------------------------'
    write(6,'(2xa)')'You are running code Cola. Cola is an opensource computational fluid dynamic '
    write(6,'(2xa)')'software distributed under GPL license. It is also a modified version of free'
    write(6,'(2xa)')'Cappuccino, that can be got from https://github.com/nikola-m/freeCappuccino. ' 
    write(6,'(2xa)')'Cola is designed simplely and itself is an experiment tool for new CFD techs. ' 
    write(6,'(2xa)')'The initial idea of codes comes from caffe3d, a code in the book of Ferziger.'
    call timestamp
  ! Record start time.
  call cpu_time(t_start)

!-----------------------------------------------------------
!  Initialization, mesh definition, sparse matrix allocation
!-----------------------------------------------------------
  call read_input_file

  call read_mesh

  call create_CSR_matrix

  call create_fields

  call init

  !
  !===============================================
  ! Time loop: 
  !===============================================
  
  itimes = itime+1
  itimee = itime+numstep

  time_loop: do itime=itimes,itimee


    ! Update variables - shift in time: 
    call time_shift
   
    ! Set inlet boundary conditions
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    call CourantNo

    iteration_loop: do iter=1,maxit 

      write(6,'(2(2x,a,i0))') 'Timestep: ',itime,',Iteration: ',iter

      call cpu_time(start)

      ! Calculate velocities.
      if(lcal(iu)) call calcuvw  

      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(lcal(ip)) then

        if(SIMPLE)   call calcp_simple
        if(PISO)     call calcp_piso

      endif 

      ! Turbulence
      if(lturb)    call modify_viscosity

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field

      ! if(lcal(ivart)) call calculate_temperature_variance_field
      
      if(lcal(icon))  call calculate_concentration_field

      if(lcal(iep))  call calculate_electric_potential

     
      ! Log scaled residuals
      write(6,'(2x,a)') 'Scaled residuals:'
      write(6,'(2x,11(a,4x))')     (chvarSolver(k), k=1,nphi)
      write(6,'(a,11(1PE9.3,2x))') '> ',(resor(k), k=1,nphi)

      call cpu_time(finish)
      write(timechar,'(f9.5)') finish-start
      write(6,'(3a)') '  ExecutionTime = ',adjustl(timechar),' s'
      write(6,*)

      !---------------------------------------------------------------
      !  Simulation management - residuals, loop control, output, etc.
      !---------------------------------------------------------------

      ! Check residual and stop program if residuals diverge
      source = max(resor(iu),resor(iv),resor(iw)) 

      if( source.gt.slarge ) then
          write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
          stop
      endif

      ! If residuals fall to level below tolerance level - simulation is finished.
      if( .not.ltransient .and. source.lt.sormax ) then
        call write_restart_files
        call writefiles
        exit time_loop
      end if

      if(ltransient) then 

        ! Has converged within timestep or has reached maximum no. of SIMPLE iterations per timetstep:
        if( source.lt.sormax .or. iter.ge.maxit ) then 

          ! Correct driving force for a constant mass flow rate simulation:
          if(const_mflux) then
            call constant_mass_flow_forcing
          endif

          ! Write values at monitoring points and recalculate time-average values for statistics:
          call writehistory
          call calc_statistics 

          ! Write field values after nzapis iterations or at the end of time-dependent simulation:
          if( mod(itime,nzapis).eq.0  .or. itime.eq.numstep ) then
            call write_restart_files
            call writefiles
          endif

          cycle time_loop
       
        endif

      end if 

    end do iteration_loop

    ! Write field values after nzapis iterations or at the end of false-time-stepping simulation:
    if(.not.ltransient .and. ( mod(itime,nzapis).eq.0 .or. (itime-itimes+1).eq.numstep ) ) then
      call write_restart_files
      call writefiles
    endif

    if(ltransient) call flush(6)

  end do time_loop

  ! Stop Record. give cpu time used.
  call cpu_time(t_end)
  write(6,'(a,I0,1x,a)')'>>MISSION COMPLETE, ', INT(t_end - t_start), ' SECONDS OF CPU TIME CONSUMED.'
  call timestamp
end program