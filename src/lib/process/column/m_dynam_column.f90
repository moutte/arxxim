module M_Dynam_Column
  use M_Kinds
  use M_Trace,      only: iDebug,Pause_,T_
  use M_T_Component,only: T_Component
  !
  implicit none
  !
  private
  !
  public:: Dynam_Column
  !
  ! public:: TUnit, Tfinal, Dtime, Time
  ! public:: PhiF, Fout, vQTotInj
  ! public:: CoupledCoores
  
contains

subroutine Dynam_Column
  !
  use M_IOTools,only: GetUnit
  use M_Files,  only: DirOut
  !
  use M_Equil_Write, only: Equil_Write_Detail
  !
  use M_Dynam_Tools   !Dynam_Zero_X,Dynam_Close
  use M_Dynam_Read    !Dynam_Read
  use M_Dynam_Tools   !Dynam_Init
  use M_Dynam_Files   !Dynam_Files_Init,Dynam_Files_Close
  use M_Dynam_Calc    !Dynam_Calc,DTimeTooSmall,PorosTooSmall,...
  !
  use M_Dynam_Cell
  !
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_Global_Vars, only: vKinFas,vFas
  !
  use M_Dynam_Vars,  only: DynTime,DynBox,DynBoxUser
  use M_Dynam_Vars,  only: vTotF,vTotInj,vMolarVol
  use M_Dynam_Vars,  only: VBox
  use M_Dynam_Vars,  only: DynColumn
  !
  use M_Dynam_Vars, only: & !
  !-----------------------import from M_Dynam_Vars to publish in M_Dynam
  & TUnit, Tfinal, Dtime, & !
  & Time,                 & !
  & PhiF, Fout, vQTotInj, & !
  & CoupledCoores
  !---------------------------------------------------------------------
  logical:: Ok
  character(len=80):: Msg
  !
  integer:: nCell
  type(T_Cell_State), allocatable:: vCell(:)
  logical:: Time_To_Small
  !
  integer:: nCp,nAq,nMk
  integer:: ix,i,j,k
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
  logical:: time_refine = .false. !.true.
  !integer:: nt, iter, Pulse_0, Pulse_1
  !-----------------------------------------------------------/rt1d vars
  !
  !-------------------------------------------------- Dynamic calc. init
  call Dynam_Zero_Numeric
  call Dynam_Zero_Time(DynTime)
  call Dynam_Zero_Box(DynBox)
  !
  call Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  if(.not. Ok) then
    print *,"Error in Dynam_Read !!"
    call Pause_
    return !------------------------------------------------------return
  end if
  !
  call Dynam_Alloc
  call Dynam_Init(TdgK,Pbar,DynTime,DynBox)
  !! print '(20G15.7)',(vTotInj(i),i=1,size(vTotInj))   ;   pause
  !
  call Dynam_ReStart_Box(DynBox,DynTime%TimeFactor)
  call Dynam_ReStart_Time(DynTime)
  !
  call Dynam_Calc_Init
  !
  if(iDebug>2) call Dynam_Files_Init
  !
  call Dynam_Save
  !--------------------------------------------------/Dynamic calc. init
  !
  write(*,*) "Dynam_Column"
  write(*,*) "TUnit=", DynTime%TUnit
  write(*,*) "VBox=", VBox
  write(*,*) "PhiF=", PhiF
  write(*,*) "TotF="
  do i=1,size(vTotF(:))
    write(*,'(A15,G12.3)') trim(vCpn(i)%NamCp),vTotF(i)
  end do
  write(*,*) "TotInj="
  do i=1,size(vTotInj(:))
    write(*,'(A15,G12.3)') trim(vCpn(i)%NamCp),vTotInj(i)
  end do
  if (iDebug>2) call pause_
  !
  call GetUnit(fCp)  ;  open(fCp, file=trim(DirOut)//"_cells_fluid.tab")
  call GetUnit(fMk)  ;  open(fMk, file=trim(DirOut)//"_cells_minerals.tab")
  call GetUnit(fVar) ;  open(fVar,file=trim(DirOut)//"_cells_min_delta.tab")
  call GetUnit(fQsk) ;  open(fQsk,file=trim(DirOut)//"_cells_min_qsk.tab")
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
  call Dynam_ReadColumn(DynColumn,Ok,Msg)
  !
  if(DynColumn%nCell>200) DynColumn%nCell= 200
  !if(Ok) then
    vv        = DynColumn%UDarcy   
    Disp      = DynColumn%Disp     
    dx        = DynColumn%dx       
    T1D_dt    = DynColumn%dt       
    T1D_Method= DynColumn%Method   
    Duration  = DynColumn%Duration 
    T1D_TSave = DynColumn%Time_Save
    nCell     = DynColumn%nCell    
  !else
  !  print *,trim(Msg)  ;  pause
  !end if
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
  allocate(vCell(nCell))
  !
  call Dynam_Column_Initialize(nCp,nAq,nMk,vCell)
  !
  allocate(tMolK_Init(nCell,nMk))
  do ix= 1,nCell
    tMolK_Init(ix,:)=vCell(ix)%vMolK(:)
  end do
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
    ! if (modulo(iter,nt/Nprint)==0) call output
    if (time_rt>=time_print) then
      call output
      time_print= time_rt + Time_Save != time for next print
    end if
    !
    time_rt= time_rt + dt
    print '(A,G15.8)',"TIME= ",time_rt
    !
    if (time_rt>Duration) exit
    !
  end do
  !------------------------------------------------------------/rt1d run
  !
  call CPU_TIME(CpuEnd)
  print '(A,G15.6,/)',"CPU_TIME=",(CpuEnd-CpuBegin)
  !
  call Dynam_Calc_Close
  call Dynam_Save
  !
  !--------------------------------------write speciation of final fluid
  !! call Equil_Write_Detail("DYN",TdgK,Pbar,vCpn)
  !
  call Dynam_Files_Close
  !
  call Dynam_Clean
  !-------------------------------------------------------/Dynamic calc.
  !
  close(fCp)
  close(fMk)
  close(fVar)
  close(fQsk)
  
contains

  subroutine conc_initial
    integer:: ix
    !-------------------------------------------------------------------
    !print *,size(conc,2), size(vCell)  ;  pause
    do ix=1,nCell
      ! conc(:,ix)= vCell(ix)%vTotF(:) /vCell(ix)%PhiF
      conc(:,ix)= vCell(ix)%vTotF(:) / (vCell(ix)%PhiF *VBox *1.D3)
      !-> components' concentr'ns in mole /litre of aq'phase
    end do
    !
  end subroutine conc_initial
  
  subroutine setboundary !(iter, Pulse_0, Pulse_1)
    ! integer,intent(in):: iter, Pulse_0, Pulse_1
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
    ! write(12,'(10G12.3)') (Conc(i,1), i=1,size(conc,1))
  end subroutine setboundary
  
  subroutine output
    integer:: i,ix
    !
    do i=1,nCp
      write(fCp,'(A,1X,G15.8,200G15.8)') &
      & trim(vCpn(i)%NamCp), &
      & time_rt, &
      & (vCell(ix)%vTotF(i)/vCell(ix)%vTotF(1), ix=1,nCell)
    end do
    write(fMk,'(A,1X,G15.8,200G15.4)') &
    & "POROSITY", &
    & time_rt, &
    & (vCell(ix)%PhiF*VBox*1.D3, ix=1,nCell)
    do i=1,nMk
      write(fMk,'(A,1X,G15.8,200G15.4)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vMolK(i)*vMolarVol(i)*1.D3, ix=1,nCell)
    end do
    do i=1,nMk
      write(fVar,'(A,1X,G15.8,200G15.4)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vMolK(i)-tMolK_Init(ix,i), ix=1,nCell)
    end do
    do i=1,nMk
      write(fQsk,'(A,1X,G15.8,200G15.8)') &
      & trim(vFas(vKinFas(i)%iFas)%NamFs), &
      & time_rt, &
      & (vCell(ix)%vKinQsk(i), ix=1,nCell)
    end do
    !
  end subroutine output
  
end subroutine Dynam_Column

subroutine Dynam_Column_Transport( &
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
end subroutine Dynam_Column_Transport

subroutine Dynam_Column_Initialize(nCp,nAq,nMk,vCell)
  use M_Global_Vars, only: vKinFas,vFas
  use M_Dynam_Cell,  only: T_Cell_State
  use M_Dynam_Vars,  only: vTotF,vMolF,vLnAct,vLnGam,vLnBuf !,vTotInj
  use M_Dynam_Vars,  only: vMolK,vKinQsk,vSurfK,vStatusK,vMolarVol
  use M_Dynam_Vars,  only: PhiF
  !
  integer,           intent(out)  :: nCp,nAq,nMk
  type(T_Cell_State),intent(inout):: vCell(:)
  !
  integer :: ix,i
  !
  nCp= size(vTotF)
  nAq= size(vMolF)
  nMk= size(vMolK)
  !
  do i= 1,size(vCell)
    call vCell(i)%New_(nCp,nAq,nMk)
  end do
  !
  do ix= 1,size(vCell)
    call vCell(ix)%Set_( &  !
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
  end do
  !
  print *,"Dynam_Column_Initialize"
  do i=1,size(vKinFas)
    print *,"KinFas= ",vFas(vKinFas(I)%iFas)%NamFs,vCell(1)%vMolK(i)
  end do
  if (iDebug>2) call pause_
  !
end subroutine Dynam_Column_Initialize

subroutine Dynam_Column_Reaction( &
& dt,        &
& conc,      &
& vCell,     &
& error)
  !use M_Dynam_Calc    !Dynam_Calc,DTimeTooSmall,PorosTooSmall,...
  use M_Dynam_Cell,  only: T_Cell_State
  use M_Dynam_Vars,  only: vTotF,vMolF,vLnAct,vLnGam,vLnBuf !,vTotInj
  use M_Dynam_Vars,  only: vMolK,vKinQsk,vSurfK,vStatusK,vMolarVol
  use M_Dynam_Vars,  only: PhiF,VBox,vQTotInj
  use M_Dynam_Vars,  only: Time,TFinal
  use M_Dynam_Calc,  only: DTimeTooSmall
  use M_Dynam_Calc,  only: Dynam_Calc_Main
  !
  real(dp),          intent(in)   :: dt
  real(dp),          intent(inout):: conc(:,:)
  type(T_Cell_State),intent(inout):: vCell(:)
  logical,           intent(out)  :: error
  !
  real(dp):: tmpTime
  integer :: ix,i
  integer :: nRepeat
  integer,parameter :: maxRepeat = 1000
  !---------------------------------------------------------------------
  error= .false.
  !
  vQTotInj= 0.D0
  !
  do ix= 1,size(vCell)
    !
    nRepeat= 1
    do_Repeat: do
    
      do i=1,nRepeat
        !
        !---------retrieve the properties of cell ix from last time step
        call vCell(ix)%Get_( &  !
        & vMolK,    & !-> mole nr's of kin'phases in cell
        & vSurfK,   & !-> surfaces of kin'phases in cell
        & vKinQsk,  & !-> sat'states of kin'phases
        & vStatusK, & !
        & vTotF,    & !
        & vMolF,    & ! mole nr's of aqu'species in cell
        & vLnAct,   & !
        & vLnGam,   & !
        & vLnBuf,   & !
        & PhiF      ) !-> vol'fraction of fluid
        !--------------------------------------------------------------/
        !
        !------retrieve fluid mole nrs resulting from transport operator
        vTotF(:)= conc(:,ix) *(vCell(ix)%PhiF*VBox*1.D3)
        !--------------------------------------------------------------/
        Time=   0.D0
        TFinal= dt/nRepeat !DynTime%TFinal
        call Dynam_Calc_Main
        !
        if (DTimeTooSmall) then
          print *,"Reaction, Cell Nr",iX," nRepeat=",nRepeat
          nRepeat= 5 *nRepeat
          cycle do_Repeat
          if(nRepeat>maxRepeat) then
            error= .true.
            return
          end if
        else
          exit do_Repeat
        end if
        !
      end do
    
    end do do_Repeat
    !
    !----------------update the fluid composition for transport operator
    conc(:,ix)= vTotF(:) /(PhiF *VBox *1.D3)
    !------------------------------------------------------------------/
    
    !----------------update the properties of cell ix for next time step
    call vCell(ix)%Set_( &  !
    & vMolK,    & !-> mole nr's of kin'phases in cell
    & vSurfK,   & !-> surfaces of kin'phases in cell
    & vKinQsk,  & !-> sat'states of kin'phases
    & vStatusK, & !
    & vTotF,    & !
    & vMolF,    & ! mole nr's of aqu'species in cell
    & vLnAct,   & !
    & vLnGam,   & !
    & vLnBuf,   & !
    & PhiF      ) !-> vol'fraction of fluid
    !------------------------------------------------------------------/
    !
  end do
  
  return
end subroutine Dynam_Column_Reaction
  
subroutine rt1d_check_courant( & !
& vv,dx, & !
& dt)
  !
  real(dp),intent(in)    :: vv,dx
  real(dp),intent(inout) :: dt
  !
  real(dp):: Courant
  
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
  real(dp),intent(in)    :: disp,vv
  real(dp),intent(inout) :: dx
  !
  real(dp):: Pec

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

end module M_Dynam_Column

