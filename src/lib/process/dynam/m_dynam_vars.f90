MODULE M_Dynam_Vars
  USE M_Kinds
  USE M_T_Component,ONLY: T_Component
  USE M_T_KinModel, ONLY: T_KinModel
  USE M_T_Phase,    ONLY: T_Phase
  
  IMPLICIT NONE
  
  PUBLIC
  
  LOGICAL:: CoupledCoores= .FALSE. !-> coores
  !
  ! if(AdjustVolBox) then all calculations done for VolBox= VBox0
  ! else calcualtions done with user's VolBox
  LOGICAL:: AdjustVolBox=  .TRUE.
  !
  !!<new 200912>
  LOGICAL,PARAMETER:: TimeIsSeconds= .FALSE. ! .TRUE.  !
  REAL(dp):: TimeScale
  !!</new 200912>
  !
  LOGICAL,PARAMETER:: OnlyActiv= .FALSE.
  !
  TYPE(T_Component),ALLOCATABLE:: vCpnInj(:)
  TYPE(T_Component),ALLOCATABLE:: vCpnBox(:)
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vQTotInj !1:nCp, flow rates of components by time unit (mol.s-1)

  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vTotF,   &  !1:nCp, total amounts of components (incl. redox balance) in the FLUID at time T
  & vTotInj, &  !1:nCp, nr moles of components in given amount of injected fluid
  & vMolSec, &  !1:nAs, nr moles of aqu'second'species
  & vMolF       !1:nAq, nr moles of aqu'species in the box
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  ! used in DynamCalc, Residual
  & vLnAct, & !Ln(Activity) aqu'species in box
  & vLnGam, & !Ln(Activ'Coeff) aqu'species in box 
  & vLnBuf    !Ln(Activity) buffering species
  
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  ! used in DynamCalc, Residual
  & vDG_As, & !DeltaGibbsFreeEnergy of reaction, for Second'Aqu'Species
  & vDG_Sp, & !DeltaGibbsFreeEnergy of reaction, for all species
  & vDLGam_As
  !
  !-----------------------------------------kinetic phases, allocatables
  CHARACTER,DIMENSION(:),ALLOCATABLE:: & !
  & vStatusK       !1:nMk, saturation state
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vMolK,     &   !1:nMk, nr moles of kin'species in the box
  & vSurfK,    &   !1:nMk, surface of kin'species in the box
  & vKinQsk        !1:nMk, saturation index
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vKinMinim, &   !1:nMk, minimal nr moles of kin'species in the box
  & vMolarVol, &   !1:nMk, Molar volume (=vFas(vKinFas(1:nMk)%iFas)%V)
  & vDG_Kin        !1:nMk,
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vVmAct,    &   !1:nMk, kinetic constant * (activators/inhibitors)
  & vVmQsK         !1:nMk, affinity factor (e.g. QsK-1)
  !!& vVm0         !1:nMk, value of vVm at begin time step (END former time step)
  !                !-> vVm(iMk)= vSurfK(iMk) *vVmAct(iMk) *vVmQsK(iMk)
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vSurfK0,   &   !1:nMk, surface of phase at ref'time
  & vMolK0         !1:nMk, mole number of phase at ref'time
  !
  LOGICAL,ALLOCATABLE:: vLKinActiv(:)
  LOGICAL,ALLOCATABLE:: vLEquActiv(:)
  INTEGER,ALLOCATABLE:: vKinPrm(:)
  !
  REAL(dp),ALLOCATABLE:: tStoikioAqu(:,:) ! stoikio aqu'sp/elements
  REAL(dp),ALLOCATABLE:: tStoikioKin(:,:) ! stoikio kin'sp/elements
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vWeitSp, &  !1:nAq, vSpc(vOrdAq(1:nAq))%WeitKg
  & vWeitCp     !1:nCp, vEle(vOrdCp(1:nCp))%WeitKg
  !
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: & !
  & tNu_Kin,  &   !1:nMd,1:nCp
  & tAlfKin       !1:nCp,1:nMd
  !
  TYPE(T_KinModel),ALLOCATABLE:: vKinMod(:)
  !----------------------------------------/kinetic phases, allocatables
  !
  !---------------------------------------------------------------------
  REAL(dp):: QsK_Iota = Zero !1.D-9 !  1.D-6 !
  !
  CHARACTER(LEN=7):: sModelSurf= "SPHERE" !"CRUNCH" !"SPHERE",
  !
  REAL(dp):: PhiInert !volume fraction occupied by inactive minerals
  !
  !------------------default texture parameters for secondary minerals--
  REAL(dp):: RadiusMinim= 1.0D-9  !radius (meter), 1.0D-9= 1 nm
  REAL(dp):: DensitMinim= 1.0D+18  !number of grains / m^3 fluid volume
  !---------------------------------------------------------------------
  ! vKinMinim(:)= lower limit for mole number in the box
  ! -> to insure that the mineral is always present (potentially)
  ! its mole number cannot decrease below vKinMinim(:)
  !!!! vKinMinim currently not used as constant !!!
  !!!! we use RadiusMinim & DensitMinim to compute VolMinim : 
  !!!! -> cf Dynam_Init_KinFas_Texture
  !
  !=< used in TP_Changed, to check whether TP-dependent param's must be updated
  REAL(dp):: TdgK0,Pbar0
  REAL(dp):: TolTdgK= One
  REAL(dp):: TolPbar= One
  !---/
  !
  INTEGER:: fSavTime= 0 !file for saving results at (nearly) fixed time intervals (dTSav)
  INTEGER:: fSavRate= 0 !file for details on the rate parameters of minerals
  !
  LOGICAL:: SteadyState_Stop= .FALSE. !.TRUE.
  !
  !-------------------------------------- box volume is FIXED or FREE --
  !------- FIXED- default, FREE- for closed system simulations mainly --
  LOGICAL:: UpdateMassFluid
  !---/
  !
  !---------------------------------------------------------------------
  TYPE:: T_DynBox  ! lumping all box related parameters
    REAL(dp)::   & !
    & VBox,      & ! volume of the "box", m^3
    & dX,        & ! length of the "box", m
    & Fout,      & ! flow rate output, m^3.s-1
    & PhiF,      & ! fluid fraction, i.e. porosity
    & RhoF,      & ! fluid density
    & UDarcy       ! Darcy velocity of fluid
    INTEGER:: nCell
    LOGICAL:: VFixed
  ENDTYPE T_DynBox
  !
  TYPE(T_DynBox):: DynBox,DynBoxUser
  !
  !---------------------------------------------------------------------
  REAL(dp)::  &
  & VBox,     & ! volume of the "box", m^3
  & dX,       & ! length of the "box", m
  & Fout        ! flow rate output, m^3.s-1
  !
  ! volume of the "box", m^3, default value
  REAL(dp):: VBox0= 1.D-3
  !
  REAL(dp)::  &  !
  & UDarcy,   &  !Darcy velocity of fluid
  & PhiF,     &  !fluid vol'fraction, i.e. porosity (and the gas ?)
  & PhiF0,    &  !idem at ref'time
  & RhoF,     &  !fluid density, kg.m^-3
  & pH_
  !
  !---------------------------------------------------------------------
  TYPE:: T_DynColumn
    INTEGER::  nCell
    INTEGER::  Method
    REAL(dp):: UDarcy, Disp
    REAL(dp):: dt, dx, duration
    REAL(dp):: time_save
  ENDTYPE T_DynColumn
  !
  TYPE(T_DynColumn):: DynColumn
  !---------------------------------------------------------------------
  TYPE:: T_DynTime
    CHARACTER(LEN=6):: TUnit !"SECOND","MINUTE","HOUR","DAY","YEAR"
    REAL(dp):: TimeFactor
    REAL(dp)::  & !
    & Time,     & ! current time
    & TFinal,   & ! end of simulation
    & dTime,    & ! current time step
    & dTmin,    & ! minimal time step
    & dTMax,    & ! maximal time step
    & dTSav       ! time laps between two records on x_time.tab
  ENDTYPE T_DynTime
  !
  TYPE(T_DynTime):: DynTime
  !
  CHARACTER(LEN=6):: TUnit !"DAY","YEAR","HOUR","SECOND","MINUTE"
  REAL(dp):: TimeFactor != the length, in seconds, of current time unit
  !
  REAL(dp)::  & !
  & Time,     & !current time
  & dTime,    & !current time step
  & TFinal      !end of simulation
  !
  REAL(dp)::  &
  & dTmin,    & !minimal time step
  & dTMax,    & !maximal time step
  & dTSav       !time laps between two records on x_time.tab
  !
  !---------------------------------------------------------------------
  LOGICAL:: Extrapole= .TRUE. ! .FALSE. !
  !
  LOGICAL:: LogForAqu=.TRUE. ! .FALSE. !
  !
  LOGICAL:: LogForMin=.FALSE.
  !if .TRUE. then works on LOG(MoleNr) instead of MoleNr for kin'phases
  !
  LOGICAL:: DirectSub=.FALSE.
  !
  ! implicitation of rate parameters
  ! CAVEAT: as Activity Coeffs are not Implicited, Rate_Act is not Implicited !!!
  LOGICAL:: Implicit_ActivFactor !.TRUE.  !used only in Residual
  LOGICAL:: Implicit_Surface     !.FALSE. !used only in residual
  !
  !------------------------------------------- parameters for SOLVER ---
  !
  !INTEGER :: iMethod
  CHARACTER(LEN=30):: cMethod
  ! default value given in Dynam_Zero_Numeric (=7 -> Newton_Walker)
  ! can be changed at run time (keyword METHOD in DYNAMIC.NUMERIC)
  ! (read by Dynam_ReadNumeric)
  INTEGER :: iCtrlTime
  REAL(dp):: &
  & Time_Decrease= 0.5_dp, & !timestep decrease factor
  & Time_Increase= 2.0_dp    !timestep increase factor
  !& Time_Increase= 1.5_dp    !timestep increase factor
  !
  INTEGER:: nEvalFunc
  !
  INTEGER:: &
  & NewtMaxIts,  &
  & NewtIterMax, & !IF(iter>NewtIterMax) dT0=dT0*Time_Decrease
  & NewtIterMin    !IF(iter<NewtIterMin) dT0=dT0*Time_Increase
  
  REAL(dp):: &
  & NewtTOLF,       & ! convergence criterion on function values
  & NewtTolF_Equil, & ! idem for aqu'species equilibrium
  & NewtTOLMIN,     & ! criterion for spurious convergence
  & NewtTOLX          ! convergence criterion on dx
  !---------------------------------------------------------------------
  !
  LOGICAL:: &
  & DebNewt, &
  & DebJacob, &
  & TestJacob, &
  & TestMax, &
  & bFinDIF
  !
  INTEGER,PUBLIC::   &
  & Dynam_nStep,     &
  & Dynam_nNewtIter, &
  & Dynam_nTotalNewtIter
  !
ENDMODULE M_Dynam_Vars
