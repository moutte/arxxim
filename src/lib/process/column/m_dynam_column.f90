MODULE M_Dynam_Column
  USE M_Kinds
  USE M_Trace,      ONLY: iDebug,Pause_,T_
  USE M_T_Component,ONLY: T_Component
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dynam_Column
  !
  ! PUBLIC:: TUnit, Tfinal, Dtime, Time
  ! PUBLIC:: PhiF, Fout, vQTotInj
  ! PUBLIC:: CoupledCoores
  
CONTAINS

SUBROUTINE Dynam_Column
  !
  USE M_IOTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut
  !
  USE M_Equil_Write, ONLY: Equil_Write_Detail
  !
  USE M_Dynam_Tools   !Dynam_Zero_X,Dynam_Close
  USE M_Dynam_Read    !Dynam_Read
  USE M_Dynam_Tools   !Dynam_Init
  USE M_Dynam_Files   !Dynam_Files_Init,Dynam_Files_Close
  USE M_Dynam_Calc    !Dynam_Calc,DTimeTooSmall,PorosTooSmall,...
  !
  USE M_Dynam_Cell
  !
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_Global_Vars, ONLY: vKinFas,vFas
  !
  USE M_Dynam_Vars,  ONLY: DynTime,DynBox,DynBoxUser
  USE M_Dynam_Vars,  ONLY: vTotF,vTotInj,vMolarVol
  USE M_Dynam_Vars,  ONLY: VBox
  USE M_Dynam_Vars,  ONLY: DynColumn
  !
  USE M_Dynam_Vars, ONLY: & !
  !-----------------------import from M_Dynam_Vars to publish in M_Dynam
  & TUnit, Tfinal, Dtime, & !
  & Time,                 & !
  & PhiF, Fout, vQTotInj, & !
  & CoupledCoores
  !---------------------------------------------------------------------
  LOGICAL:: Ok
  CHARACTER(LEN=80):: Msg
  !
  INTEGER:: nCell
  TYPE(T_Cell_State), ALLOCATABLE:: vCell(:)
  LOGICAL:: Time_To_Small
  !
  INTEGER:: nCp,nAq,nMk
  INTEGER:: ix,i,j,k
  !
  real(dp):: CpuBegin,CpuEnd
  !
  !------------------------------------------------------------rt1d vars
  integer :: T1D_Method
  real(dp):: Retard,Disp,VV
  real(dp):: T1D_dt,dx,dt
  logical :: T1D_error
  !
  real(dp), allocatable:: pConc(:,:)
  real(dp), allocatable:: Conc(:,:)
  real(dp), allocatable:: tMolK_Init(:,:)
  !
  real(dp):: time_rt, Duration !, time_0, time_1,
  real(dp):: time_print, Time_Save, T1D_TSave
  integer :: iter
  integer:: fCp,fMk,fVar,fQsk ! file indexes
  logical:: time_refine = .FALSE. !.TRUE.
  !integer:: nt, iter, Pulse_0, Pulse_1
  !-----------------------------------------------------------/rt1d vars
  
  !-------------------------------------------------- Dynamic calc. init
  CALL Dynam_Zero_Numeric
  CALL Dynam_Zero_Time(DynTime)
  CALL Dynam_Zero_Box(DynBox)
  !
  CALL Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  IF(.NOT. Ok) THEN
    PRINT *,"Error in Dynam_Read !!"
    CALL Pause_
    RETURN !------------------------------------------------------RETURN
  ENDIF
  !
  CALL Dynam_Alloc
  CALL Dynam_Init(TdgK,Pbar,DynTime,DynBox)
  !! print '(20G15.7)',(vTotInj(i),i=1,size(vTotInj))   ;   pause
  !
  CALL Dynam_ReStart_Box(DynBox,DynTime%TimeFactor)
  CALL Dynam_ReStart_Time(DynTime)
  !
  CALL Dynam_Calc_Init
  !
  IF(iDebug>2) CALL Dynam_Files_Init
  !
  CALL Dynam_Save
  !--------------------------------------------------/Dynamic calc. init
  !
  write(*,*) "Dynam_Column"
  write(*,*) "TUnit=", DynTime%TUnit
  write(*,*) "VBox=", VBox
  write(*,*) "pHIf=", PhiF
  write(*,*) "TotF="
  do i=1,size(vTotF(:))
    write(*,'(A15,G12.3)') TRIM(vCpn(i)%NamCp),vTotF(i)
  end do
  write(*,*) "TotInj="
  do i=1,size(vTotInj(:))
    write(*,'(A15,G12.3)') TRIM(vCpn(i)%NamCp),vTotInj(i)
  end do
  if (iDebug>2) call pause_
  !
  CALL GetUnit(fCp)  ;  OPEN(fCp, FILE=TRIM(DirOut)//"_cells_fluid.tab")
  CALL GetUnit(fMk)  ;  OPEN(fMk, FILE=TRIM(DirOut)//"_cells_minerals.tab")
  CALL GetUnit(fVar) ;  OPEN(fVar,FILE=TRIM(DirOut)//"_cells_min_delta.tab")
  CALL GetUnit(fQsk) ;  OPEN(fQsk,FILE=TRIM(DirOut)//"_cells_min_qsk.tab")
  !
  !------------------------------------------------------------rt1d init
  !
  !--default parameters
  DynColumn%Method   = 2
  DynColumn%UDarcy   = 0.5D0
  DynColumn%Disp     = 0.001D0
  DynColumn%dx       = 0.1D0
  DynColumn%dt       = 0.1D0
  DynColumn%Duration = 1.D0
  DynColumn%Time_Save= DynColumn%Duration /10.D0
  DynColumn%nCell    = 20
  !
  !--read parameters from file
  CALL Dynam_ReadColumn(DynColumn,Ok,Msg)
  !
  if(DynColumn%nCell>200) DynColumn%nCell= 200
  !IF(Ok) THEN
    vv        = DynColumn%UDarcy   
    Disp      = DynColumn%Disp     
    dx        = DynColumn%dx       
    T1D_dt    = DynColumn%dt       
    T1D_Method= DynColumn%Method   
    Duration  = DynColumn%Duration 
    T1D_TSave = DynColumn%Time_Save
    nCell     = DynColumn%nCell    
  !ELSE
  !  PRINT *,TRIM(Msg)  ;  PAUSE
  !ENDIF
  !
  ! call rt1d_check_peclet(disp, vv, dx)
  call rt1d_check_courant(vv, dx, T1D_dt)
  !
  !print *,"dt=",dt  ;  pause
  !dt= 0.01D0
  Retard= 1.0D0
  Time_Save=  T1D_TSave !/100.D0
  time_print= Time_Save
  !-----------------------------------------------------------/rt1d init
  !
  !-----------------------------------------------------------init cells
  ALLOCATE(vCell(nCell))
  !
  CALL Dynam_Column_Initialize(nCp,nAq,nMk,vCell)
  !
  ALLOCATE(tMolK_Init(nCell,nMk))
  DO ix= 1,nCell
    tMolK_Init(ix,:)=vCell(ix)%vMolK(:)
  END DO
  !! write(12,'(200G12.3)') (vCell(i)%vMolK(1), i=1,nCell)
  !! write(13,'(200G12.3)') (vCell(i)%vMolK(2), i=1,nCell)
  !--------------------------------------------------------/init cells--
  !
  !-----------------------------------------------------------rt1d run--
  allocate(Conc(nCp,nCell+1))
  allocate(pConc(nCp,nCell+1))
  !
  call conc_initial !!(nCell,Conc)
  !
  time_rt= 0.D0
  !
  call output
  !
  call pause_
  !
  ! call Dynam_Column_Reaction(Time_To_Small)
  !
  call CPU_TIME(CpuBegin)
  !
  do ! time loop
    !
    !iter= iter +1
    !print *,"iter=",iter
    if(time_refine) then
      if(time_rt<Duration/1000.D0) then
        dt= T1D_dt/100.D0
        !Time_Save= T1D_TSave /100.D0
      else if(time_rt<Duration/100.D0) then
        dt= T1D_dt/10.D0
        !Time_Save= T1D_TSave /10.D0
      else
        dt= T1D_dt
        !Time_Save= T1D_TSave
      end if
    else
      dt= T1D_dt
    end if
    !
    !----------------------------------------------------------transport
    Conc(:,nCell+1)= Conc(:,nCell)
    !--save previous time step concentrations
    pConc(:,:)= Conc(:,:)
    !--set the  boundary condition for the time step
    call setboundary ! (iter, Pulse_0, Pulse_1) !!, nCp, Conc)
    !
    call Dynam_Column_Transport( &
    & T1D_Method,  &
    & nCp,nCell,   &
    & vv,dt,dx,    &
    & Disp,Retard, &
    & pConc, Conc)
    !---------------------------------------------------------/transport
    !
    !-----------------------------------------------------------reaction
    call Dynam_Column_Reaction(dt,conc(:,:),vCell(:),Time_To_Small)
    if(Time_To_Small) then
      ! dt= dt /2.D0
      ! cycle
      exit
    end if
    !----------------------------------------------------------/reaction
    !
    !!   if(DTimeTooSmall) then   ;  dt= dt/10.-dp
    !!   else                     ;  exit
    !!   end if
    !!   !
    !! end do
    !
    ! if (modulo(iter,nt/NPRINT)==0) call output
    if (time_rt>=time_print) then
      call output
      time_print= time_rt + Time_Save != time for next print
    end if
    !
    time_rt= time_rt + dt
    print '(A,G15.8)',"Time= ",time_rt
    !
    if (time_rt>Duration) exit
    !
  end do
  !------------------------------------------------------------/rt1d run
  !
  call CPU_TIME(CpuEnd)
  print '(A,G15.6,/)',"CPU_TIME=",(CpuEnd-CpuBegin)
  !
  CALL Dynam_Calc_Close
  CALL Dynam_Save
  !
  !--------------------------------------write speciation of final fluid
  !! CALL Equil_Write_Detail("DYN",TdgK,Pbar,vCpn)
  !
  CALL Dynam_Files_Close
  !
  CALL Dynam_Clean
  !-------------------------------------------------------/Dynamic calc.
  !
  close(fCp)
  close(fMk)
  close(fVar)
  close(fQsk)
  
CONTAINS

  SUBROUTINE conc_initial
    INTEGER:: ix
    !-------------------------------------------------------------------
    !print *,size(conc,2), size(vCell)  ;  pause
    DO ix=1,nCell
      ! conc(:,ix)= vCell(ix)%vTotF(:) /vCell(ix)%PhiF
      conc(:,ix)= vCell(ix)%vTotF(:) / (vCell(ix)%PhiF *VBox *1.D3)
      !-> components' concentr'ns in mole /litre of aq'phase
    END DO
    !
  END SUBROUTINE conc_initial
  
  SUBROUTINE setboundary !(iter, Pulse_0, Pulse_1)
    ! INTEGER,INTENT(IN):: iter, Pulse_0, Pulse_1
    !!!-----------------------------------------------------------------
    !!if (iter <= Pulse_0) then
    !!  Conc(:, 1)= 1.D-3
    !!else if (iter <= Pulse_1) then
    !!  Conc(:, 1)= 1.D0
    !!else
    !!  Conc(:, 1)= 1.D-3
    !!end if
    !!!-----------------------------------------------------------------
    if (vv>0.D0) conc(:,1)= vTotInj(:)
    ! write(12,'(10G12.3)') (Conc(i,1), i=1,SIZE(conc,1))
  END SUBROUTINE setboundary
  
  SUBROUTINE output
    INTEGER:: i,ix
    !
    DO i=1,nCp
      WRITE(fCp,'(A,1X,G15.8,200G15.8)') &
      & trim(vCpn(i)%NamCp), &
      & time_rt, &
      & (vCell(ix)%vTotF(i)/vCell(ix)%vTotF(1), ix=1,nCell)
    END DO
    WRITE(fMk,'(A,1X,G15.8,200G15.4)') &
    & "POROSITY", &
    & time_rt, &
    & (vCell(ix)%PhiF*VBox*1.D3, ix=1,nCell)
    DO i=1,nMk
      WRITE(fMk,'(A,1X,G15.8,200G15.4)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vMolK(i)*vMolarVol(i)*1.D3, ix=1,nCell)
    END DO
    DO i=1,nMk
      WRITE(fVar,'(A,1X,G15.8,200G15.4)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vMolK(i)-tMolK_Init(ix,i), ix=1,nCell)
    END DO
    DO i=1,nMk
      WRITE(fQsk,'(A,1X,G15.8,200G15.8)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vKinQsk(i), ix=1,nCell)
    END DO
    !
  END SUBROUTINE output
  
END SUBROUTINE Dynam_Column

SUBROUTINE Dynam_Column_Transport( &
& T1D_Method,  &
& nCp,nCell,   &
& vv,dt,dx,    &
& Disp,Retard, &
& pConc, Conc)
  !---------------------------------------------------------------------
  use M_Dynam_Transport
  !---------------------------------------------------------------------
  integer, intent(in)   :: T1D_Method
  integer, intent(in)   :: nCp,nCell
  real(dp),intent(in)   :: vv,dt,dx
  real(dp),intent(in)   :: Disp,Retard
  real(dp),intent(in)   :: pConc(:,:)
  real(dp),intent(inout):: Conc(:,:)
  !---------------------------------------------------------------------
  integer:: i
  !---------------------------------------------------------------------
  select case(T1D_Method)
  case(0) ! explicit advection
    do i= 1,nCp
      call T1D_Adv_Explicit(vv,dt,dx, Retard, nCell, pConc(i,:), Conc(i,:))
      if (Disp > 0.D0) call T1D_Dispersion(dt,dx,Disp, Retard, nCell, Conc(i,:))
    end do
  case(1) ! implicit advection
    do i= 1,nCp
      call T1D_Adv_Dis_Implicit(vv,dt,dx,Disp,Retard, nCell, Conc(i,:))
    end do
  case(2) ! TVD & Implicit dispersion
    do i= 1,nCp
      call T1D_Adv_TVD(vv,dt,dx, Retard, nCell, pConc(i,:), Conc(i,:))
      if (Disp > 0.D0) &
      & call T1D_Dispersion(dt,dx,Disp, Retard, nCell, Conc(i,:))
    end do
  end select
  !
END SUBROUTINE Dynam_Column_Transport

SUBROUTINE Dynam_Column_Initialize(nCp,nAq,nMk,vCell)
  USE M_Global_Vars, ONLY: vKinFas,vFas
  USE M_Dynam_Cell,  ONLY: T_Cell_State
  USE M_Dynam_Vars,  ONLY: vTotF,vMolF,vLnAct,vLnGam,vLnBuf !,vTotInj
  USE M_Dynam_Vars,  ONLY: vMolK,vKinQsk,vSurfK,vStatusK,vMolarVol
  USE M_Dynam_Vars,  ONLY: PhiF
  !
  INTEGER,           INTENT(OUT)  :: nCp,nAq,nMk
  TYPE(T_Cell_State),INTENT(INOUT):: vCell(:)
  !
  INTEGER :: ix,i
  !
  nCp= SIZE(vTotF)
  nAq= SIZE(vMolF)
  nMk= SIZE(vMolK)
  !
  DO i= 1,SIZE(vCell)
    CALL vCell(i)%New_(nCp,nAq,nMk)
  END DO
  !
  DO ix= 1,SIZE(vCell)
    CALL vCell(ix)%Set_( &  !
    & vMolK,    & !-> mole nr's of kin'phases in cell
    & vSurfK,   & !-> surfaces of kin'phases in cell
    & vKinQsk,  & !-> sat'states of kin'phases
    & vStatusK, & !
    & vTotF,    & !
    & vMolF,    & ! mole nr's of aqu'species in cell
    & vLnAct,   & !
    & vLnGam,   & !
    & vLnBuf,   & !
    & PhiF      & !-> vol'fraction of fluid
    & )
  END DO
  !
  print *,"Dynam_Column_Initialize"
  do i=1,SIZE(vKinFas)
    print *,"KinFas= ",vFas(vKinFas(I)%iFas)%NamFs,vCell(1)%vMolK(i)
  end do
  if (iDebug>2) call pause_
  !
END SUBROUTINE Dynam_Column_Initialize

SUBROUTINE Dynam_Column_Reaction( &
& dt,        &
& conc,      &
& vCell,     &
& error)
  !USE M_Dynam_Calc    !Dynam_Calc,DTimeTooSmall,PorosTooSmall,...
  USE M_Dynam_Cell,  ONLY: T_Cell_State
  USE M_Dynam_Vars,  ONLY: vTotF,vMolF,vLnAct,vLnGam,vLnBuf !,vTotInj
  USE M_Dynam_Vars,  ONLY: vMolK,vKinQsk,vSurfK,vStatusK,vMolarVol
  USE M_Dynam_Vars,  ONLY: PhiF,VBox,vQTotInj
  USE M_Dynam_Vars,  ONLY: Time,TFinal
  USE M_Dynam_Calc,  ONLY: DTimeTooSmall
  USE M_Dynam_Calc,  ONLY: Dynam_Calc_Main
  !
  REAL(dp),          INTENT(IN)   :: dt
  REAL(dp),          INTENT(INOUT):: conc(:,:)
  TYPE(T_Cell_State),INTENT(INOUT):: vCell(:)
  LOGICAL,           INTENT(OUT)  :: error
  !
  REAL(dp):: tmpTime
  INTEGER :: ix,i
  INTEGER :: nRepeat
  INTEGER,PARAMETER :: maxRepeat = 1000
  !---------------------------------------------------------------------
  error= .FALSE.
  !
  vQTotInj= 0.D0
  !
  DO ix= 1,SIZE(vCell)
    !
    nRepeat= 1
    DO_Repeat: DO
    
      DO i=1,nRepeat
        !
        !---------retrieve the properties of cell ix from last time step
        CALL vCell(ix)%Get_( &  !
        & vMolK,    & !-> mole nr's of kin'phases in cell
        & vSurfK,   & !-> surfaces of kin'phases in cell
        & vKinQsk,  & !-> sat'states of kin'phases
        & vStatusK, & !
        & vTotF,    & !
        & vMolF,    & ! mole nr's of aqu'species in cell
        & vLnAct,   & !
        & vLnGam,   & !
        & vLnBuf,   & !
        & PhiF      & !-> vol'fraction of fluid
        & )
        !--------------------------------------------------------------/
        !
        !------retrieve fluid mole nrs resulting from transport operator
        vTotF(:)= conc(:,ix) *(vCell(ix)%PhiF*VBox*1.D3)
        !--------------------------------------------------------------/
        Time=   0.D0
        TFinal= dt/nRepeat !DynTime%TFinal
        CALL Dynam_Calc_Main
        !
        IF (DTimeTooSmall) THEN
          PRINT *,"Reaction, Cell Nr",iX," nRepeat=",nRepeat
          nRepeat= 5 *nRepeat
          CYCLE DO_Repeat
          IF(nRepeat>maxRepeat) THEN
            error= .TRUE.
            RETURN
          END IF
        ELSE
          EXIT DO_Repeat
        END IF
        !
      END DO
    
    END DO DO_Repeat
    !
    !----------------update the fluid composition for transport operator
    conc(:,ix)= vTotF(:) /(PhiF *VBox *1.D3)
    !------------------------------------------------------------------/
    
    !----------------update the properties of cell ix for next time step
    CALL vCell(ix)%Set_( &  !
    & vMolK,    & !-> mole nr's of kin'phases in cell
    & vSurfK,   & !-> surfaces of kin'phases in cell
    & vKinQsk,  & !-> sat'states of kin'phases
    & vStatusK, & !
    & vTotF,    & !
    & vMolF,    & ! mole nr's of aqu'species in cell
    & vLnAct,   & !
    & vLnGam,   & !
    & vLnBuf,   & !
    & PhiF      & !-> vol'fraction of fluid
    & )
    !------------------------------------------------------------------/
    !
  END DO
  
  RETURN
END SUBROUTINE Dynam_Column_Reaction
  
subroutine rt1d_check_courant( & !
& vv,dx, & !
& dt)
  !
  real(8),intent(in)    :: vv,dx
  real(8),intent(inout) :: dt
  !
  real(8):: Courant
  
  Courant= vv * dt / dx
  print *,"Courant number= ",Courant
  if (Courant>1.D0) then
    dt= dx / vv *0.5D0
    print *,"New dt= ",dt
  end if
  call pause_
  !
end subroutine rt1d_check_courant

subroutine rt1d_check_peclet( & !
& disp,vv, & !
& dx)
  !
  real(8),intent(in)    :: disp,vv
  real(8),intent(inout) :: dx
  !
  real(8):: Pec

  Pec= vv * dx / disp
  print *,"Peclet number= ",Pec
  if (Pec > 2.D0) then
    print *,"Peclet number > 2 !!!"
    dx= disp / vv *1.5D0
    print *,"New dx= ", dx
  end if
  call pause_
  !
end subroutine rt1d_check_peclet

END MODULE M_Dynam_Column

