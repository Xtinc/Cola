module types
!
! define precision for floating-point numbers
!
  ! double precision    
  integer, parameter :: dp = kind(1.0d0) 
end module types

module parameters  
  use types

  integer, parameter :: nphi=11 ! Max. no. of variables to solve
  integer, parameter :: iu=1    ! Variable identifiers
  integer, parameter :: iv=2
  integer, parameter :: iw=3
  integer, parameter :: ip=4
  integer, parameter :: ite=5
  integer, parameter :: ied=6
  integer, parameter :: ien=7
  integer, parameter :: ivis=8
  integer, parameter :: ivart=9
  integer, parameter :: icon=10
  integer, parameter :: iep=11
  integer, parameter :: ibmag=12

  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: small = 1e-20
  real(dp), parameter :: great = 1e+20
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: twothirds = 2./3._dp
  real(dp), parameter :: onethird = 1./3._dp
  
  ! Law of the wall parameters
  real(dp), parameter :: CAPPA = 0.41_dp 
  real(dp), parameter :: ELOG = 8.432_dp
  real(dp), parameter :: ctrans = 11.63


  integer :: monCell   ! Monitoring point for values, usually for log file
  integer :: pRefCell  ! Pressure reference cell
  integer :: mpoints   ! No. of monitoring points

  real(dp) :: flomas   ! mass flow at inlet
  real(dp) :: flomom   ! momentum of flow at inlet
  real(dp) :: prt1     ! Used in calcheatflux, not set at the moment TODO!
  real(dp) :: densit   ! Fluid density
  real(dp) :: viscos   ! Molecular dynamic viscosity
  real(dp) :: sormax   ! Residual toleance for SIMPLE
  real(dp) :: slarge   ! Upper limit for residuals before simulation blowup
  real(dp) :: facnap   ! Underelaxation factor for Reynolds stresses
  real(dp) :: facflx   ! Underelaxation factor for turbulent heat fluxes
  real(dp) :: uin,vin,win,tein,edin,tin,vartin,conin ! Inlet values (assumed constant across inlet)

  real(dp) :: pranl     ! (= 0.7_dp for air, 7.0_dp for water, read it from input file.)
  logical :: lbuoy      ! Bouyancy effects - are they included in momentum and turbulence eqns. If yes we calculate heat fluxes.
  logical :: boussinesq ! Is Boussinesq hypothesis evoked yes=1/no=0, read from input file.
  real(dp) :: tref      ! Reference temperature, read from input file
  real(dp) :: beta      ! Thermal expansion coefficient, read from input file
  real(dp) :: phit      ! Parameter used for GGDH and AFM in calcheatflux subroutine, read from input file.
  real(dp) :: sksi      ! Parameter used in calcheatflux, read from input file
  real(dp) :: eta       ! Parameter used in calcheatflux, read from input file
  real(dp) :: gravx,gravy,gravz ! Components of gravity acceleration vector, read from input file
  real(dp) :: facvis            ! Under-relaxation for viscosity
  real(dp) :: erough            ! E parameter for rough wall
  real(dp) :: zzero             ! z0 - wall roughness [m] for aerodynamically rough walls

  !
  ! Mesh format type
  !
  character( len=10 ) :: mesh_format
  character(len = 20) :: linear_solver
  character(len = 20) :: symmetry_linear_solver

  !
  ! Timesteping control
  !
  integer :: numstep   ! Total number of timesteps
  integer :: itime     ! Current timestep value
  integer :: nzapis    ! After how many timesteps we write backup and output files
  integer :: maxit     ! Maximum number of iterations in timestep, also max. number of iterations for SIMPLE iteration
  real(dp) :: timestep ! Timestep size
  real(dp) :: time     ! Total elapsed time
  real(dp) :: btime    ! Coefficient for Backward Euler time stepping algorithm, if btime = 1. => BDF2 Three time implicit (2nd order), if btime=0. Backward Euler (1st order)
  
  logical :: CoNumFix         ! Is Courant no. fixed during time-stepping
  real(dp) :: CoNumFixValue   ! Fixed value for Courant number - set in modinp for now - may be read in input
  real(dp) :: CoNum,meanCoNum ! Courant number.  
  ! character(len=9) :: timechar! A char string to write current timestep

    
  ! Logicals, mostly read from simulation-input file:
  logical :: lturb,lread,lwrite,ltest             ! turbulent simulation, read restart file, write restart file, print residual of the linear solver,.,..      
  logical :: ltransient                           ! LTRANSIENT is TRUE for transient (non-stationary) simulations              
  logical :: levm,lasm,lles,lsgdh,lggdh,lafm      ! eddy-viscosity, algebraic stress model or LES, simple gradient or generalized gradient hypothesis, algerbaic flux model
  logical :: bdf,bdf2,bdf3,cn                     ! control for the time-stepping algorithm
  logical :: simple,piso                          ! control for the velocity-pressure coupling algorithm
  logical :: const_mflux                          ! control for constant flow rate 
  logical :: solveOmega, solveEpsilon, SolveTKE   ! Selfexplanatory, used in 'init'


  integer :: ncorr                   ! PISO control parameter: no. of Piso corrections.
  integer :: icorr                   ! PISO iteration no.: icorr=1..ncorr
  integer :: npcor                   ! No. of pressure-corrections; non-orthogonality correctors
  integer :: ipcorr                  ! Iteration no.: ipcorr=1..npcor
  integer :: nigrad                  ! No. of iters. for iterative cell-centered gradient calculation
  integer, parameter :: nipgrad = 2  ! No. of stages for 'pressure extrapolation at boundary + pressure gradient update' loop

  integer :: TurbModel ! Turbulence model case selector

  logical :: roughWall                                            ! Is aerodynamically rough wall assumed
  real(dp) :: magUbar, gradPcmf                                   ! Magnitude of the bulk velocity,and pressure grad that will drive the constant mass-flux flow (cmf)
  real(dp) :: sumLocalContErr, globalContErr, cumulativeContErr   ! Continuity errors


  ! Those with nphi are related to each field that we calculate U,V,W,P,T,TKE,ED...
  logical,  dimension(nphi) :: lcal      ! Logical do we calculate that particular field
  integer,  dimension(nphi) :: nsw       ! Number of allowed iterations in linear solver for each variable
  real(dp), dimension(nphi) :: sor       ! Tolerance level for residual in linear solver for each variable
  real(dp), dimension(nphi) :: resor     ! Residual
  real(dp), dimension(nphi) :: prtinv    ! Inverse Prandtl numbers, (see diffusion term discretisation)
  real(dp), dimension(nphi) :: urf       ! Underrelaxation factor
  real(dp), dimension(nphi) :: urfr      ! Recipr. value of urf: 1/urf
  real(dp), dimension(nphi) :: urfm      ! 1-urf
  real(dp), dimension(nphi) :: gds       ! Gamma blending factor [0,1] for deffered correction for convection terms: uds + gds*(uhigh-uds), 0-UDS, 1-Higher order diff.scheme

  !menu control hardcoded
  logical :: EnableMenu=.true.
  integer :: Gnuplotpoints = 1250
  
end module parameters

module variables
!%%%%%%%%%%%%%%
  use types

    ! These are cellwise defined variables, that is - the fields
    real(dp), dimension(:), allocatable :: u,v,w              ! Velocity components
    real(dp), dimension(:), allocatable :: flmass             ! Mass fluxes trough inner faces
    real(dp), dimension(:), allocatable :: p,pp               ! Pressure, Press. correction,  
    real(dp), dimension(:), allocatable :: te,ed              ! Turb. kin. energy, Dissipation,
    real(dp), dimension(:), allocatable :: vis                ! Effective viscosity
    real(dp), dimension(:), allocatable :: den                ! Density
    real(dp), dimension(:), allocatable :: gen                ! Turb. kin. energy generation
    real(dp), dimension(:), allocatable :: t                  ! Temperature, given in Degree Celsius.
    real(dp), dimension(:), allocatable :: vart               ! Temperature variance
    ! real(dp), dimension(:), allocatable :: edd              ! Dissipation of temperature variance
    real(dp), dimension(:), allocatable :: utt,vtt,wtt        ! Turbulent heat fluxes
    ! real(dp), dimension(:), allocatable :: ret              ! Turbulent Re number
    real(dp), dimension(:), allocatable :: con                ! Concentration
    real(dp), dimension(:), allocatable :: uu,vv,ww,uv,uw,vw  ! Reynolds stress tensor components
    real(dp), dimension(:,:), allocatable :: bij              ! Reynolds stress anisotropy tensor
    ! real(dp), dimension(:), allocatable :: fmi, fmo           ! Mass fluxes trough boundary faces
    real(dp), dimension(:), allocatable :: visw,ypl,tau       ! Effective visc. for boundary face, the y+ non-dimensional distance from wall, tau - wall shear stress
  
    ! values from n-1 timestep
    real(dp), dimension(:), allocatable :: uo, vo, wo
    real(dp), dimension(:), allocatable :: po
    real(dp), dimension(:), allocatable :: flmasso
    real(dp), dimension(:), allocatable :: teo
    real(dp), dimension(:), allocatable :: edo
    real(dp), dimension(:), allocatable :: to
    real(dp), dimension(:), allocatable :: varto
    real(dp), dimension(:), allocatable :: cono

    ! values from n-2 time step
    real(dp), dimension(:), allocatable :: uoo,voo,woo
    real(dp), dimension(:), allocatable :: poo
    real(dp), dimension(:), allocatable :: flmassoo
    real(dp), dimension(:), allocatable :: teoo
    real(dp), dimension(:), allocatable :: edoo
    real(dp), dimension(:), allocatable :: too
    real(dp), dimension(:), allocatable :: vartoo
    real(dp), dimension(:), allocatable :: conoo 

    ! values from n-3 time step
    real(dp), dimension(:), allocatable :: uooo,vooo,wooo
    real(dp), dimension(:), allocatable :: pooo
    real(dp), dimension(:), allocatable :: flmassooo

    ! Gradients
    real(dp),dimension(:,:), allocatable :: dUdxi,dVdxi,dWdxi
    real(dp),dimension(:,:), allocatable :: dPdxi
    real(dp),dimension(:,:), allocatable :: dTEdxi
    real(dp),dimension(:,:), allocatable :: dEDdxi
    real(dp),dimension(:,:), allocatable :: dTdxi
    real(dp),dimension(:,:), allocatable :: dCondxi
    real(dp),dimension(:,:), allocatable :: dVartdxi

    real(dp), dimension(:), allocatable :: magStrain          ! Strain magnitude
    real(dp), dimension(:), allocatable :: Vorticity          ! Vorticity magnitude
    
    type :: PlotArray
        integer,private::ArraySize=10000
        integer,private::ArrayCol=11
        integer,public ::Level =0                             ! hardcoded now
        real,dimension(:,:),allocatable:: arraytmp           ! for plot
    contains
        procedure::NominalSize
        procedure::setSize
        procedure::ReSize
        procedure::Row
        procedure::Col
        procedure::Mealloc
        procedure::Dealloc
    end type
    type(PlotArray)::plt
    contains
    
    function NominalSize(this)
    class(plotArray)::this
    integer::NominalSize
    NominalSize=this%ArraySize*2**this%Level
    end function
    
    subroutine setSize(this,M,N)
    class(PlotArray)::this
    integer,intent(in)::M,N
    this%ArraySize=min(this%ArraySize,M)
    this%ArrayCol=min(this%ArrayCol,N)
    end subroutine
    
    subroutine ReSize(this)
    class(PlotArray)::this
    integer::i,j
    real,dimension(this%ArraySize)::tmp
    do j=1,this%ArrayCol
        tmp=0.0
        do i=1,floor(this%ArraySize/2.0)
            tmp(i)=this%arraytmp(2*i-1,j)
        end do
        this%arraytmp(:,j)=tmp(:)
    end do
    this%Level=this%Level+1
    end subroutine
    
    function Row(this)
    class(PlotArray)::this
    integer::Row
    Row=this%ArraySize
    end function
    
    function Col(this)
    class(PlotArray)::this
    integer::Col
    Col=this%ArrayCol
    end function
    
    subroutine Mealloc(this)
    class(PlotArray)::this
    allocate(this%arraytmp(this%ArraySize,this%ArrayCol))
    end subroutine
    
    subroutine Dealloc(this)
    class(PlotArray)::this
    deallocate(this%arraytmp)
    end subroutine
    
end module variables


module title_mod
!%%%%%%%%%%%%
!
! Used for output to monitor file.
!
use parameters, only: nphi

  character(len=70) :: title
  character(len=4), dimension(nphi) ::  chvar = &
  (/'  U ', '  V ', '  W ', '  P ',' TE ', ' ED ', '  T ', ' VIS', 'VART', ' CON', 'EPOT' /)
  character(len=7), dimension(nphi) ::  chvarSolver = &
  (/'U      ','V      ','W      ','p      ','k      ','epsilon','Temp   ','Visc   ','VarTemp','Conc   ','Epot   ' /)
  character(len=100):: inlet_file,grid_file,monitor_file,restart_file
    contains
    
  function judge_convection_scheme(scheme)
  character(len=*),intent(in)::scheme
  character(len=25)::judge_convection_scheme
    select case( trim(scheme) ) 

    case ('cds')
      judge_convection_scheme = 'cds'
    case ('cdscorr')
      judge_convection_scheme = 'cdscorr'
    case ('luds')
      judge_convection_scheme = 'luds'
    case('kappa')
      judge_convection_scheme = 'kappa'
    case ('smart')
      judge_convection_scheme = 'smart'
    case ('avl-smart')
      judge_convection_scheme = 'avl-smart'
    case ('muscl')
      judge_convection_scheme = 'muscl'
    case ('umist')
      judge_convection_scheme = 'umist'
    case ('koren')
      judge_convection_scheme = 'koren'
    case ('charm')
      judge_convection_scheme = 'charm'
    case ('ospre')
      judge_convection_scheme = 'ospre'
    case ('minmod')
      judge_convection_scheme = 'minmod'
    case ('central')
      judge_convection_scheme = 'central'
    case ('linearUpwind')
      judge_convection_scheme = 'linearUpwind'
    case ('boundedLinearUpwind')
      judge_convection_scheme = 'boundedLinearUpwind'
    case ('boundedLinearUpwind02')
      judge_convection_scheme = 'boundedLinearUpwind02'
    case ('boundedCentral')
      judge_convection_scheme = 'boundedCentral'
    case default
      judge_convection_scheme = '2nd order upwind'
    end select
    return
  end function
  
  function judge_convection_limiter(scheme)
  character(len=*),intent(in)::scheme
  character(len=20)::judge_convection_limiter
    select case( trim(scheme) ) 
    case ('Barth-Jespersen')
      judge_convection_limiter = 'Barth-Jespersen'
    case ('Venkatakrishnan')
      judge_convection_limiter = 'Venkatakrishnan'
    case ('MDL')
      judge_convection_limiter = 'Multidimensional'
    case default
      judge_convection_limiter = 'None'
    end select
    return
  end function
  
  function judge_turbulence_model(scheme)
    integer,intent(in)::scheme
    character(78)::judge_turbulence_model
    select case(scheme) 
    case (1)
      judge_turbulence_model = 'Standard k-epsilon'
    case (2)
      judge_turbulence_model = 'Renormalization Group k-epsilon'
    case (3)
      judge_turbulence_model = 'Shear Stress Transport k-Omega'
    case (4)
      judge_turbulence_model = 'Standard k-Omega'
    case (5)
      judge_turbulence_model = 'Spalart Allmaras'
    case (6)
      judge_turbulence_model = 'Yoshizawa'
    case (7)
      judge_turbulence_model = '2-Layer Enchanced Wall k-epsilon'
    case (8)
      judge_turbulence_model = 'Realizable k-epsilon'
    case (9)
      judge_turbulence_model = '2-layer Enchanced Wall Realizable k-epsilon'
    case (10)
      judge_turbulence_model = 'Detached Eddy Model with Standard k-Omega'
    case (11)
      judge_turbulence_model = 'Improved Delayed Detached Eddy Model'
    case (12)
      judge_turbulence_model = 'Smagorinsky SGS.'
    case default
      judge_turbulence_model = 'None'
    end select
    return
  end function
  
  function judge_linear_solver(scheme,symmetry)
  character(len=*),intent(in)::scheme
  logical,intent(in)::symmetry
  character(len=20)::judge_linear_solver
    select case( trim(scheme) ) 
    case ('gauss-seidel')
      judge_linear_solver = 'gauss-seidel'
    case('bicgstab_ilu')
      judge_linear_solver = 'bicgstab_ilu'
    case('pmgmres_ilu')
      judge_linear_solver = 'pmgmres_ilu'
    case('dcg')
      if(symmetry)then
          judge_linear_solver = 'dcg'
      else
          judge_linear_solver = 'bicgstab_ilu'
      end if
    case('iccg')
        if(symmetry)then
          judge_linear_solver = 'iccg'
        else
          judge_linear_solver = 'bicgstab_ilu'  
        end if
    case default
      judge_linear_solver = 'bicgstab_ilu'
    end select
    return
  end function
  
end module title_mod


module statistics
! Variables for collecting statistics
! these are ensamble averaged values over time steps
! used in t_rans-urans-hybrid rans/les approach
! these are saved in tecplot_stat file
  use types

    integer :: n_sample,istat,ifreq

    real(dp), dimension(:), allocatable :: u_aver,v_aver,w_aver
    real(dp), dimension(:), allocatable :: uu_aver,vv_aver,ww_aver, uv_aver,uw_aver,vw_aver
    real(dp), dimension(:), allocatable :: te_aver
    real(dp), dimension(:), allocatable :: t_aver
    real(dp), dimension(:), allocatable :: ut_aver,vt_aver,wt_aver
    real(dp), dimension(:), allocatable :: tt_aver
                                             
end module statistics


module hcoef
!%%%%%%%%%%%%%
  use types  
    ! Arrays used in piso algorithm 
    real(dp), dimension(:), allocatable :: h,rU,rV,rW
end module hcoef