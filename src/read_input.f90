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
  use m_config

  implicit none

  integer :: i,imon
  character(len=2) :: trpn
  character(len=25) :: convective_scheme
  character(len=:), allocatable :: varnames
  character(len=8) :: len8str
  type(CFG_t) :: my_cfg
  logical::LCONTROL(3) = .false.
  real(dp)::GRAVITY(3) = 0.0_DP
  character(30)::inputstr = ''
!
!***********************************************************************
!
  
  CALL CFG_ADD(MY_CFG,"TITLE","LAMINAR","TITLE")
  CALL CFG_ADD(MY_CFG,"MESHFORMAT","NATIVEMESH","FOAMMESH/NATIVEMESH")
  CALL CFG_ADD(MY_CFG,"LCONTROL",[.FALSE.,.TRUE.,.FALSE.],"READ RESTART FILE/WRITE RESTART FILE/PRINT LINEAR SOLVER RESIDUAL", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"LCAL",[.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.],"SOLVING FOR U V W P...", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"MONCELL",1,"MONITOR CELL")
  CALL CFG_ADD(MY_CFG,"PREFCELL",1,"PRESSURE REFERENCE CELL")
  CALL CFG_ADD(MY_CFG,"IPREFPROCESSMPOINTS",0,"MONITOR CELL NUMBER")
  CALL CFG_ADD(MY_CFG,"SLARGE",1.0_DP,"DIVERGENCE RESIDUAL")
  CALL CFG_ADD(MY_CFG,"SORMAX",0.01_DP,"CONVERGENCE RESIDUAL")
  CALL CFG_ADD(MY_CFG,"DENSIT",1.0_DP,"DENSITY")
  CALL CFG_ADD(MY_CFG,"VISCOS",0.00001_DP,"DYNAMIC VISCOSITY")
  CALL CFG_ADD(MY_CFG,"PRANL",0.71_DP,"PRANDTL NUMBER")
  CALL CFG_ADD(MY_CFG,"TREF",20.0_DP,"REFERENCE TEMPERATURE")
  CALL CFG_ADD(MY_CFG,"BETA",0.0001_DP,"THERMAL EXPANSION COEFFICIENT")
  CALL CFG_ADD(MY_CFG,"LBUOY",.FALSE.,"BOUYANCY EFFECTS IN MOMENTUM & TURBULENCE EQUATIONS")
  CALL CFG_ADD(MY_CFG,"GRAVITY",[0.0_DP,-9.81_DP,0.0_DP],"GRAVITY")
  CALL CFG_ADD(MY_CFG,"BOUSSINESQ",.TRUE.,"BOUSSINESQ HYPOTHESIS")
  CALL CFG_ADD(MY_CFG,"FACNAP",1.0_DP,"UNDERELAXATION FACTOR FOR REYNOLDS STRESSES")
  CALL CFG_ADD(MY_CFG,"FACFLX",1.0_DP,"UNDERELAXATION FACTOR FOR TURBULENT HEAT FLUXES")
  CALL CFG_ADD(MY_CFG,"LTRANSIENT",.TRUE.,"TRANSIENT")
  CALL CFG_ADD(MY_CFG,"TEMPROAL","BDF","TEMPROAL SCHEME")
  CALL CFG_ADD(MY_CFG,"TURBMODEL1","NONE","LEVM/LASM/LLES")
  CALL CFG_ADD(MY_CFG,"TURBMODEL2",0,"TURBMODEL2")
  CALL CFG_ADD(MY_CFG,"TURBGRAD","NONE","LSGDH/LGGDH/LAFM/NONE")
  CALL CFG_ADD(MY_CFG,"CONVECTIVESCHEME","LINEARUPWIND","LINEARUPWIND/")
  CALL CFG_ADD(MY_CFG,"GRADIENTLIMITER","NONE","VENKATAKRISHNAN/NONE")
  CALL CFG_ADD(MY_CFG,"GDS",[1.0_DP,1.0_DP,1.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP],"GAMMA BLENDING FACTOR [0,1] FOR DEFFERED CORRECTION FOR CONVECTION TERMS", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"URF",[1.0_DP,1.0_DP,1.0_DP,1.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP],"UNDERRELAXATION FACTOR", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"SOR",[0.001_DP,0.001_DP,0.001_DP,0.00001_DP,0.001_DP,0.001_DP,0.001_DP,0.001_DP,0.001_DP,0.001_DP,0.001_DP],"TOLERANCE LEVEL FOR RESIDUAL IN LINEAR SOLVER FOR EACH VARIABLE", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"NSW",[5,5,5,100,10,10,5,5,5,5,5],"NUMBER OF ALLOWED ITERATIONS IN LINEAR SOLVER FOR EACH VARIABLE", dynamic_size=.true.)
  CALL CFG_ADD(MY_CFG,"NUMSTEP",1000,"NUMSTEP")
  CALL CFG_ADD(MY_CFG,"TIMESTEP",0.001_DP,"TIMESTEP")
  CALL CFG_ADD(MY_CFG,"NZAPIS",100,"RECORD RESULT EVERY NZAPIS CYCLES")
  CALL CFG_ADD(MY_CFG,"MAXIT",5,"MAXIMUM NUMBER OF ITERATIONS IN TIMESTEP")
  CALL CFG_ADD(MY_CFG,"GRADIENTRECONSTRUCTION","GAUSS","LSTSQ/GAUSS")
  CALL CFG_ADD(MY_CFG,"NPCOR",1,"NO. OF PRESSURE-CORRECTIONS; NON-ORTHOGONALITY CORRECTORS")
  CALL CFG_ADD(MY_CFG,"NIGRAD",2,"NO. OF ITERS. FOR ITERATIVE CELL-CENTERED GRADIENT CALCULATION")
  CALL CFG_ADD(MY_CFG,"PVCOUPLE","PISO","PVCOUPLE")
  CALL CFG_ADD(MY_CFG,"NCORR",1,"PISO CONTROL PARAMETER: NO. OF PISO CORRECTIONS")
  CALL CFG_ADD(MY_CFG,"CONSTMFLUX",.FALSE.,"CONTROL FOR CONSTANT FLOW RATE")
  CALL CFG_ADD(MY_CFG,"MAGUBAR",0.0_DP,"MAGNITUDE OF THE BULK VELOCITY")
  CALL CFG_ADD(MY_CFG,"CONUMFIX",.TRUE.,"FIXED VALUE FOR COURANT NUMBER")
  CALL CFG_ADD(MY_CFG,"CONUMFIXVALUE",1.0_DP,"COURANT NUMBER")
  CALL CFG_ADD(MY_CFG,"ENABLEMENU",.TRUE.,"ENABLE MENU INVOKED BY KEYBOARD WHEN PROGRAM RUNS")
  CALL CFG_ADD(MY_CFG,"GNUPLOTPOINTS",1250,"MAX SAMPLE POINTS FOR GNUPLOT WHEN PLOTTING RESIDUAL")
  CALL CFG_ADD(MY_CFG,"MONITORFILE","MONITOR","MONITOR FILE NAME")
  CALL CFG_ADD(MY_CFG,"RESTARTFILE","RESTART","RESTART FILE NAME")
  
  call CFG_read_file(my_cfg,"cola_settings.ini")
  call CFG_update_from_arguments(my_cfg)
  
  call CFG_get(my_cfg,"TITLE",TITLE)
  call CFG_get(my_cfg,"MESHFORMAT",mesh_format)
  call CFG_get_dynamic(my_cfg,"LCONTROL",LCONTROL)
  LREAD =LCONTROL(1)
  LWRITE=LCONTROL(2)
  LTEST =LCONTROL(3)
  call CFG_get_dynamic(my_cfg,"LCAL",LCAL)
  call CFG_get(my_cfg,"MONCELL",monCell)
  call CFG_get(my_cfg,"PREFCELL",pRefCell)
  call CFG_get(my_cfg,"IPREFPROCESSMPOINTS",MPoints)
  call CFG_get(my_cfg,"SLARGE",SLARGE)
  call CFG_get(my_cfg,"SORMAX",SORMAX)
  call CFG_get(my_cfg,"DENSIT",DENSIT)
  call CFG_get(my_cfg,"VISCOS",VISCOS)
  call CFG_get(my_cfg,"PRANL",PRANL)
  call CFG_get(my_cfg,"TREF",TREF)
  call CFG_get(my_cfg,"BETA",BETA)
  call CFG_get(my_cfg,"LBUOY",LBUOY)
  call CFG_get_dynamic(my_cfg,"GRAVITY",GRAVITY)
  GRAVX=GRAVITY(1)
  GRAVY=GRAVITY(2)
  GRAVZ=GRAVITY(3)
  call CFG_get(my_cfg,"BOUSSINESQ",BOUSSINESQ)
  call CFG_get(my_cfg,"FACNAP",FACNAP)
  call CFG_get(my_cfg,"FACFLX",FACFLX)
  call CFG_get(my_cfg,"LTRANSIENT",LTRANSIENT)
  BDF=.false.
  BDF2=.false.
  BDF3=.false.
  CN=.false.
  call CFG_get(my_cfg,"TEMPROAL",inputstr)
  select case (trim(adjustl(inputstr)))
  case('BDF')
      BDF=.true.
  case('BDF2')
      BDF2=.true.
  case('BDF3')
      BDF3=.true.
  case('CN')
      CN=.true.
  case default
      BDF=.true.
  end select
  LEVM=.false.
  LASM=.false.
  LLES=.FALSE.
  call CFG_get(my_cfg,"TURBMODEL1",inputstr)
  select case (trim(adjustl(inputstr)))
  case('LEVM')
      LEVM=.TRUE.
  case('LASM')
      LASM=.TRUE.
  case('LLES')
      LLES=.TRUE.
  case default
  end select
  LSGDH=.false.
  LGGDH=.false.
  LAFM=.false.
  call CFG_get(my_cfg,"TURBGRAD",inputstr)
  select case (trim(adjustl(inputstr)))
  case('LSGDH')
      LSGDH=.TRUE.
  case('LGGDH')
      LGGDH=.TRUE.
  case('LAFM')
      LAFM=.TRUE.
  case default
  end select
  call CFG_get(my_cfg,"TURBMODEL2",TurbModel)
  call CFG_get(my_cfg,"CONVECTIVESCHEME",inputstr)
  convective_scheme = judge_convection_scheme(inputstr)
  call CFG_get(my_cfg,"GRADIENTLIMITER",inputstr)
  limiter = judge_convection_limiter(inputstr)
  call CFG_get_dynamic(my_cfg,"GDS",GDS)
  call CFG_get_dynamic(my_cfg,"URF",URF)
  call CFG_get_dynamic(my_cfg,"SOR",SOR)
  call CFG_get_dynamic(my_cfg,"NSW",NSW)
  call CFG_get(my_cfg,"NUMSTEP",NUMSTEP)
  call CFG_get(my_cfg,"TIMESTEP",TIMESTEP)  
  call CFG_get(my_cfg,"NZAPIS",NZAPIS)  
  call CFG_get(my_cfg,"MAXIT",MAXIT)  
  lstsq=.false.
  lstsq_qr=.false.
  lstsq_dm=.false.
  gauss=.false.
  call CFG_get(my_cfg,"GRADIENTRECONSTRUCTION",inputstr)
  select case(trim(adjustl(inputstr)))
  case('lstsq')
      lstsq=.true.
  case('lstsq_qr')
      lstsq_qr=.true.
  case('lstsq_dm')
      lstsq_dm=.true.
  case('gauss')
      gauss=.true.
  case default
      gauss=.true.
  end select
  call CFG_get(my_cfg,"NPCOR",NPCOR)  
  call CFG_get(my_cfg,"NIGRAD",NIGRAD) 
  SIMPLE=.false.
  PISO=.false.
  call CFG_get(my_cfg,"PVCOUPLE",inputstr)
  select case(trim(adjustl(inputstr)))
  case('SIMPLE')
      SIMPLE=.true.
  case('PISO')
      PISO=.true.
  case default
      PISO=.true.
  end select
  call CFG_get(my_cfg,"NCORR",ncorr)
  call CFG_get(my_cfg,"CONSTMFLUX",const_mflux)
  call CFG_get(my_cfg,"MAGUBAR",magUBar)
  call CFG_get(my_cfg,"CONUMFIX",CoNumFix)  
  call CFG_get(my_cfg,"CONUMFIXVALUE",CoNumFixValue)   
  call CFG_get(my_cfg,"ENABLEMENU",EnableMenu)
  call CFG_get(my_cfg,"GNUPLOTPOINTS",GNUPLOTPOINTS)
  call CFG_get(my_cfg,"MONITORFILE",monitor_file)
  call CFG_get(my_cfg,"RESTARTFILE",restart_file)
!.END: READ INPUT FILE.............................................!


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
  
  call CFG_write(my_cfg, "cola_settings.ini")
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