!***********************************************************************
!
subroutine read_input_file
!
!***********************************************************************
!
! Open & Read and Process Input File
!
!***********************************************************************
  use types
  use parameters
  use gradients, only: lstsq, lstsq_qr, lstsq_dm, gauss, limiter
  use interpolation
  use title_mod
  use utils, only: formatSTR,rc1,rc2,lc

  implicit none

  integer :: i,imon
  character(len=2) :: trpn
  character(len=25) :: convective_scheme
  character(len=:), allocatable :: varnames
  character(len=8) :: len8str
!
!***********************************************************************
!

  OPEN(UNIT=5,FILE=input_file,status='old')
! on *inx platform, open function will creat a new file with given name 
! whether state specified. Maybe it has more relations with gcc. But now
! status quo. 
  REWIND 5

  READ(5,'(a)') TITLE
  READ(5,'(a)') mesh_format 
  READ(5,*) LREAD,LWRITE,LTEST
  READ(5,*) (LCAL(I),I=1,NPHI)
  READ(5,*) monCell,pRefCell,MPoints
  READ(5,*) SLARGE,SORMAX
  READ(5,*) DENSIT,VISCOS
  READ(5,*) PRANL,TREF,BETA
  READ(5,*) LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ
  READ(5,*) roughWall,EROUGH,ZZERO
  READ(5,*) FACNAP,FACFLX
  READ(5,*) LTRANSIENT,BDF,BDF2,BDF3,CN
  READ(5,*) LEVM,LASM,LLES
  READ(5,*) LSGDH,LGGDH,LAFM
  READ(5,*) TurbModel
  READ(5,*) convective_scheme
  READ(5,*) limiter
  READ(5,*) (GDS(I),I=1,NPHI)
  READ(5,*) (URF(I),I=1,NPHI)
  READ(5,*) (SOR(I),I=1,NPHI)
  READ(5,*) (NSW(I),I=1,NPHI)
  READ(5,*) NUMSTEP,TIMESTEP,NZAPIS,MAXIT
  READ(5,*) lstsq, lstsq_qr, lstsq_dm, gauss
  READ(5,*) NPCOR, NIGRAD
  READ(5,*) SIMPLE,PISO,ncorr
  READ(5,*) const_mflux,magUBar
  READ(5,*) CoNumFix, CoNumFixValue
!.END: READ INPUT FILE.............................................!
  CLOSE (5)
!.Check the input.  
  convective_scheme = judge_convection_scheme(convective_scheme)
  limiter = judge_convection_limiter(limiter)
!.Create an input file reading log:
  WRITE(*,'(2x,a)') formatSTR('Input log:',.false.,lc)
  WRITE(*,'(2x,a)')     '------------------------------------------------------------------------------'
  WRITE(*,'(2x,a,a1,a)')'#Basic________________________________________________________________________'
  WRITE(*,'(2x,a,a1,a)')formatSTR('Title',.false.,lc),'|',formatSTR(Title,.true.,20)
  WRITE(*,'(2x,a,a1,a)')formatSTR('Mesh Format',.false.,lc),'|',formatSTR(mesh_format,.true.,20)
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Read RESTART File',.false.,lc),'|',LREAD
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Write RESTART File',.false.,lc),'|',LWRITE
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Print Linear Solver Residual',.false.,lc),'|',LTEST
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Monitor Cell',.false.,lc),'|',monCell
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Pressure Reference Cell',.false.,lc),'|',pRefCell
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Monitor Cell Number',.false.,lc),'|',MPoints
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Density           ',.false.,lc),'|',DENSIT
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Dynamic Viscosity  ',.false.,lc),'|',VISCOS
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Prandtl Number        ',.false.,lc),'|',PRANL
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Reference Temperature',.false.,lc),'|',TREF
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Thermal Expansion Coefficient',.false.,lc),'|',BETA
  WRITE(*,'(2x,a,a1,f20.7)')formatSTR('Gravity X Compoment',.false.,lc),'|',GRAVX
  WRITE(*,'(2x,a,a1,f20.7)')formatSTR('Gravity Y Compoment',.false.,lc),'|',GRAVY
  WRITE(*,'(2x,a,a1,f20.7)')formatSTR('Gravity Z Compoment',.false.,lc),'|',GRAVZ
  WRITE(*,'(2x,a,a1,a)')'#Solver__________________________________________________','|','____________________'
  WRITE(*,'(2x,a,a1,a)')formatSTR('Convective Scheme:',.false.,lc),'|',formatSTR(trim(convective_scheme),.true.,rc2)
  varnames=''
  if(lstsq)then
      varnames = 'LSQR'
  elseif(lstsq_qr)then
      varnames = 'LSQR_QR'
  elseif(lstsq_dm)then
      varnames = 'Weighted LSQR'
  else
      varnames = 'Green-Gauss'
  endif
  WRITE(*,'(2x,a,a1,a)')formatSTR('Gradient Reconstruction Method',.false.,lc),'|',formatSTR(varnames,.true.,rc2)
  WRITE(*,'(2x,a,a1,a)')formatSTR('Gradient Limiter',.false.,lc),'|',formatSTR(limiter,.true.,rc2)
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Iterations Of Gradient Calculation',.false.,lc),'|',NIGRAD
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Iterations Of Pressure Corrections',.false.,lc),'|',NPCOR
  if(simple)then
      varnames = 'SIMPLE'
  elseif(piso)then
      varnames = 'PISO'
  else
      simple = .true.
      varnames = 'SIMPLE'
  endif
  WRITE(*,'(2x,a,a1,a)')formatSTR('Pressure Velocity Coupled Algrothrim',.false.,lc),'|',formatSTR(varnames,.true.,rc2)
  if(piso)then
      WRITE(*,'(2x,a,a1,i20)')formatSTR('Iterations Of PISO Corrections',.false.,lc),'|',ncorr
  endif
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Convergence Residual In A Timestep',.false.,lc),'|',SORMAX
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Divergence Residual In A Timestep',.false.,lc),'|',SLARGE
  if(LTRANSIENT)then
      if(BDF)then 
          varnames = 'BDF in 1st order'
      elseif(BDF2)then
          varnames = 'BDF in 2nd order'
      elseif(BDF3)then
          varnames = 'BDF in 3rd order'
      elseif(CN)then
          varnames = 'Crank-Nicolson'
      else
          BDF = .true.
          varnames = 'BDF in 1st order'
      endif
  else
      varnames = 'None'
  endif
  WRITE(*,'(2x,a,a1,a)')formatSTR('Temporal Scheme',.false.,lc),'|',formatSTR(varnames,.true.,rc2)  
  if(LTRANSIENT)then
      WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Timestep Size',.false.,lc),'|',TIMESTEP
  endif
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Total Number Of Cycles',.false.,lc),'|',NUMSTEP
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Output Intervals (Cycles Number)',.false.,lc),'|',NZAPIS
  WRITE(*,'(2x,a,a1,i20)')formatSTR('Maximum Cycles In A Timestep',.false.,lc),'|',MAXIT
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Constant Flow Rate Hypothesis',.false.,lc),'|',const_mflux
  WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Magnitude Of The Bulk Velocity',.false.,lc),'|',magUBar
  if(CoNumFix)then
      WRITE(*,'(2x,a,a1,es20.7)')formatSTR('Fixed Value For Courant Number',.false.,lc),'|',CoNumFixValue
  endif
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Bouyancy Effects In Momentum & Turbulence Eqns:',.false.,lc),'|',LBUOY
  WRITE(*,'(2x,a,a1,L20)')formatSTR('Boussinesq Hypothesis:',.false.,lc),'|',BOUSSINESQ
  WRITE(*,'(2x,a,a1,a)')'#Turbulence______________________________________________','|','____________________'
  ! Turbulent flow computation condition:
  ! Current Turb Models are mainly eddy visc models.
  ! In fact, I think input checking should be more strict. like LSGDH,LGGDH,LAFM and some others have
  ! relations. but since a new input form will be given by nicloa mirkov, just stay the status quo.
  lturb = levm.or.lasm.or.lles
  if(lturb)then
      if(lasm)then
          WRITE(*,'(2x,a)')formatSTR('Algebraic Turbulence Stress Model',.true.,rc1)
      endif
      if(levm.or.lles)then
          WRITE(*,'(2x,a)')formatSTR(trim(judge_turbulence_model(TurbModel)),.true.,rc1)
      endif
      varnames=''
      if(LSGDH)then
          varnames = formatSTR('Simple Gradient Hypothesis',.true.,rc2)
      elseif(LGGDH)then
          varnames = formatSTR('Generalized Gradient Hypothesis',.true.,rc2)
      elseif(LAFM)then
          varnames = formatSTR('Algerbaic Flux Model',.true.,rc2)
      else
          varnames = formatSTR('None',.true.,rc2)
      endif
      WRITE(*,'(2x,a,a1,a)')formatSTR('Turbulent Heat Flux Calculation Method',.false.,lc),'|',varnames
      WRITE(*,'(2x,a,a1,f20.7)') formatSTR('Underelaxation Factor For Reynolds Stresses',.false.,lc),'|',FACNAP
      WRITE(*,'(2x,a,a1,f20.7)') formatSTR('Underelaxation Factor For Turbulent Heat Fluxes',.false.,lc),'|',FACFLX
  else
      WRITE(*,'(2x,a)')formatSTR('None',.true.,rc1)
  endif
  WRITE(*,'(2x,a)')'------------------------------------------------------------------------------'

  WRITE(*,*)' '
  WRITE(*,'(2x,a)')'Solving For:'
  WRITE(*,'(2x,a)')'------------------------------------------------------------------------------'
  varnames = ''
  do i =1,NPHI
      if(LCAL(I))then
          write(len8str,'(a1,a7)')'|',trim(adjustl(chvar(i)))
          varnames = varnames//len8str
      end if
  end do
  WRITE(*,'(2x,a,a)')'Var.Name    ',varnames
  varnames = ''
  do i =1,NPHI
      if(LCAL(I))then
          write(len8str,'(a1,f7.2)')'|',GDS(i)
          varnames = varnames//len8str
      end if
  end do
  WRITE(*,'(2x,a,a)')'DC.Factor   ',varnames
  varnames = ''
  do i =1,NPHI
      if(LCAL(I))then
          write(len8str,'(a1,f7.2)')'|',URF(i)
          varnames = varnames//len8str
      end if
  end do
  WRITE(*,'(2x,a,a)')'Relax.Factor',varnames
  varnames = ''
  do i =1,NPHI
      if(LCAL(I))then
          write(len8str,'(a1,f7.2)')'|',SOR(i)*100
          varnames = varnames//len8str
      end if
  end do
  WRITE(*,'(2x,a,a)')'Max.Res.%   ',varnames
  varnames = ''
  do i =1,NPHI
      if(LCAL(I))then
          write(len8str,'(a1,f7.2)')'|',NSW(i)
          varnames = varnames//len8str
      end if
  end do
  WRITE(*,'(2x,a,a)')'Max.Iters.  ',varnames
  WRITE(*,'(2x,a)')'------------------------------------------------------------------------------'
  !
  ! Convective scheme:
  !
  call set_convective_scheme(convective_scheme)  


  !
  ! Open files for data at monitoring points 
  !
  if( ltransient .and. mpoints>0 ) then
    open(unit=89,file='transient_monitoring_points')
    rewind 89
    do imon=1,mpoints
      write(trpn,'(i2)') imon
      open(91+imon,file="transient_monitor_point_"//trpn, access='append')
      if(.not.lread) rewind(91+imon)
    end do
  end if

end subroutine