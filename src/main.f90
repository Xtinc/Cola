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
  use utils, only: timestamp,CLS_CMD_Progress,CMD_PROGRESS_ABSOLUTE,WriteHead,ConverTime

  implicit none

  integer :: k, ios
  integer :: total_cycles = 0
  integer :: local_cycles = 0
  integer :: iter
  integer :: itimes, itimee
  real(dp):: source, res_ratio
  real :: start, finish, t_start, t_end
  logical :: diverge = .false.
  logical :: last_iter_diverge = .false.
  logical :: fulloop = .true.
  character( len = 9) :: timechar
  type( CLS_CMD_Progress ) ::Progress
!                                                                       
!******************************************************************************
!
! Print logo in the header of monitor file
  call WriteHead
  call timestamp  
  call read_input_file

  open(unit=6,file=monitor_file,iostat=ios)

  ! Record start time.
  call cpu_time(t_start)

!-----------------------------------------------------------
!  Initialization, mesh definition, sparse matrix allocation
!-----------------------------------------------------------

  call read_mesh

  call create_CSR_matrix

  call create_fields

  call init

  !
  !===============================================
  ! Time loop: 
  !===============================================
  
  itimes = itime+1
  itimee = numstep
  !we won't check res in first 5% cycle
  
  WRITE(*,*)''
  WRITE(*,'(2x,a)')'Looping:'
  WRITE(*,'(2x,a)') '------------------------------------------------------------------------------'

  time_loop: do itime=itimes,itimee
    
    total_cycles=total_cycles+local_cycles  
    
    !ProgressBar
    call Progress%Set( N = numstep , L = 60 )        
    call Progress%Put(itime,CMD_PROGRESS_ABSOLUTE)
    ! Update variables - shift in time: 
    call time_shift
   
    ! Set inlet boundary conditions
    if(itime.eq.itimes) call bcin

    ! Courant number report:
    call CourantNo
    
    iteration_loop: do iter=1,maxit 
        
      local_cycles=iter
  
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

      call RecordRes(total_cycles+iter)
      
      if(EnableMenu) call exec_menu(total_cycles+iter)

      call cpu_time(finish)
      write(timechar,'(f9.5)') finish-start
      write(6,'(3a)') '  ExecutionTime = ',adjustl(timechar),' s'
      write(6,*)

      !---------------------------------------------------------------
      !  Simulation management - residuals, loop control, output, etc.
      !---------------------------------------------------------------

      res_ratio = source
      ! Check residual and stop program if residuals diverge
      source = max(resor(iu),resor(iv),resor(iw),resor(ip)) 
      !res_ratio=source/lastres
      res_ratio = source/res_ratio

      if(source.gt.slarge) then
          if(last_iter_diverge.and.(res_ratio.gt.10)) then
              ! Residual is growing and derivative is more than 10(should be user specified).
              write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
              diverge=.true. 
              fulloop=.false.
              exit time_loop
          else
              last_iter_diverge = .true.
          end if
      endif
      

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
          exit iteration_loop   
        endif
      else
          ! If residuals fall to level below tolerance level - simulation is finished.
          if(source.lt.sormax)then
              call write_restart_files
              call writefiles
              fulloop=.false.
              exit time_loop
          end if
          
      end if 

    end do iteration_loop

    if(ltransient) then
        call flush(6)
    else
        ! Write field values after nzapis iterations or at the end of false-time-stepping simulation:
        if(mod(itime,nzapis).eq.0 .or. (itime-itimes+1).eq.numstep)then
            call write_restart_files
            call writefiles
        end if
    end if    
    
  end do time_loop
  

  ! Stop Record. give cpu time used.
  call cpu_time(t_end)
  
  if(.not.fulloop)then
      write(*,*)''
      if(diverge)then
          write(*,"(2x,a,a,a)") "Program Terminated -  Iterations Diverge, ", ConverTime(INT(t_end - t_start)), 'CPU Time Consumed.'
      end if
      
  else
      write(*,'(2x,a,a,a)')'Mission Complete, ', ConverTime(INT(t_end - t_start)), 'CPU Time Consumed.'
  end if
  
  WRITE(*,'(2x,a)') '------------------------------------------------------------------------------'

  call timestamp

end program
    
        
subroutine RecordRes(iter)
    use title_mod
    use parameters,only:resor,nphi
    use variables,only:plt
    
    implicit none
    integer,intent(in)::iter
    integer::k,m
    
    m=0
    ! Log scaled residuals
    write(6,'(2x,a)') 'Scaled residuals:'
    write(6,'(2x,11(a,4x))')     (chvarSolver(k), k=1,nphi)
    write(6,'(a,11(1PE9.2,2x))') '> ',(resor(k), k=1,nphi)
    
    ! get Residual for plot
    if(iter.gt.plt%NominalSize())then
        call plt%Resize
    end if
    if(mod(iter-1,2**plt%Level).eq.0)then
        m=(iter-1)/2**plt%Level
        do k=1,plt%col()
            plt%arraytmp(m+1,k)=resor(k)
        end do
    end if
    
end subroutine