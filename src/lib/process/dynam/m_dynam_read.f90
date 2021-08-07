module M_Dynam_Read
  use M_Kinds
  use M_Trace
  implicit none
  !
  private
  public:: Dynam_Read
  public:: Dynam_ReadColumn

contains

subroutine Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  use M_Numeric_Const,only: LN10
  use M_KinModel_Read,only: KinModel_Init
  use M_Dynam_Tools,  only: Dynam_TimeFactor
  !
  use M_Global_Vars,  only: vSpc,vKinModel
  use M_Dynam_Vars,   only: T_DynTime,T_DynBox,TimeFactor,TUnit
  use M_Dynam_Vars,   only: sModelSurf,Implicit_Surface,VBox0
  use M_Dynam_Vars,   only: LogForAqu,DirectSub,bFinDif
  use M_Dynam_Vars,   only: AdjustVolBox
  !
  type(T_DynTime),intent(inout):: DynTime
  type(T_DynBox), intent(inout):: DynBox
  type(T_DynBox), intent(out)  :: DynBoxUser
  logical,        intent(out)  :: Ok
  !
  logical:: ReadOk
  character(len=255):: Msg
  real(dp):: A
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_Read"
  !
  Ok=.true.
  !
  !-----------------------------------------------read time parameters--
  call Dynam_ReadTime(DynTime,ReadOk)
  !-> TUnit,TFinal,dTime,dTMin,dTMax,dTSav,TAgain
  if (.not.ReadOk) then
    call Stop_("Block DYNAMIC.TIME not found")
    Ok = .false.
    return
  end if
  DynTime%TimeFactor= Dynam_TimeFactor(DynTime%TUnit)
  !--------------------------------------------------------------------/
  !
  !------------------------------------------------read box parameters--
  call Dynam_ReadBox(DynBox,ReadOk,Msg)
  !-> VBox,dX,etc.
  if(.not. ReadOk .and. iDebug>2) print '(A)',trim(Msg)
  !
  DynBoxUser= DynBox
  !
  if(AdjustVolBox .and. DynBox%VBox /= VBox0) then
    A= (DynBox%VBox/VBox0)**(1.0D0/3.0D0)
    DynBox%VBox=   VBox0
    DynBox%DX=     DynBox%DX     /A
    DynBox%UDarcy= DynBox%UDarcy /A
  end if
  !
  if(iDebug>2) then
    write(fTrc,'(/,A)')     "<  BOX VOL/DX/UDARCY ADJUSTED==>"
    write(fTrc,'(A,2G15.6)') "  VBox  =", DynBoxUser%VBox,   DynBox%VBox
    write(fTrc,'(A,2G15.6)') "  dX    =", DynBoxUser%dX,     DynBox%dX
    write(fTrc,'(A,2G15.6)') "  UDarcy=", DynBoxUser%UDarcy, DynBox%UDarcy
    write(fTrc,'(A,/)')     "</ BOX VOL/DX/UDARCY ADJUSTED==>"
  end if
  !---------------------------------------------------------------------
  !
  !------------------------------------------ read numeric parameters --
  call Dynam_ReadNumeric(ReadOk)
  !
  !if(.not. LogForAqu) DirectSub= .true.
  !if(.not. LogForAqu) bFinDif=   .true.
  if(DirectSub) bFinDif=   .true.
  !
  if(sModelSurf=="CRUNCH") Implicit_Surface=.false.
  !---------------------------------------------------------------------
  !
  !---------------------------------------------- read kinetic models --
  call KinModel_Init(vSpc)
  !! ReadOk= size(vKinModel)>0
  !! if(.not. ReadOk) call Stop_("NO Kinetic Data Found") !______stop
  if(size(vKinModel)>0 .and. iDebug>2) call Check_KinRate_CalcQsK
  !---------------------------------------------/ read kinetic models --
  !
  !------------------------------------- read box mineral composition --
  call Dynam_ReadRock(ReadOk)
  !! if(.not.ReadOk)  call Stop_("NO Kinetic Minerals Found")
  !---------------------------------------------------------------------
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_Read"
  !
  return
end subroutine Dynam_Read

subroutine Dynam_ReadTime(T,Ok)
!--
!-- scan the input file, read block DYNAMIC.TIME
!-- -> temporal conditions of the run
!-- type:: T_DynTime
!--   character(len-6):: TUnit !-- "DAY","YEAR","HOUR"
!--   real(dp)::  &
!--   & Time,     & !-- current time
!--   & dTime,    & !-- current time step
!--   & dTmin,    & !-- minimal time step
!--   & dTMax,    & !-- maximal time step
!--   & dTSav       !-- time laps between two records on x_time.tab
!-- end type T_DynTime
!--
  use M_Files,only: NamFInn
  use M_IOTools
  use M_Dynam_Vars, only: SteadyState_Stop,TimeIsSeconds
  use M_Dynam_Vars, only: T_DynTime,TUnit,TFinal,dTime,dTMin,dTMax,dTSav,TimeFactor
  use M_Dynam_Tools,only: Dynam_TimeFactor
  !
  type(T_DynTime),intent(inout):: T
  logical,        intent(out)  :: Ok
  !
  character(len=255):: L
  character(len=80) :: W,W1,W2,WW
  logical           :: EoL
  integer           :: F,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_ReadTime"
  !
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))
  Ok=.false.
  !
  Do01: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
    call LinToWrd(L,W,EoL)
    !
    if(W(1:1)=='!')   cycle Do01 !skip comment lines
    if(W=="END" .and. .not. Eol ) then !when end followed by keyword, append
      call LinToWrd(L,WW,EoL); W=trim(W)//trim(WW)
    end if
    !
    select case(trim(W))
    !
    case("ENDINPUT"); exit  Do01
    !
    case("DYNAMIC.TIME","DYNAMIC")
      !
      Ok=.true.
      !
      Do02: do
        !
        read(F,'(A)',iostat=ios) L  ;  if(ios/=0) exit Do01
        call LinToWrd(L,W1,EoL)     ;  call AppendToEnd(L,W1,EoL)
        if(W1(1:1)=='!') cycle Do02
        call LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        select case(trim(W1))
        !
        !----------------------------------------------timing parameters
        case("ENDINPUT")                                    ;  exit Do01
        case("END","ENDDYNAMIC.TIME","ENDDYNAMIC")          ;  exit Do02
        !
        case("STEADY")
          SteadyState_Stop= (trim(W2)=="stop")
        case("TUNIT")
          select case(trim(W2))
          case("SECOND")  ;  T%TUnit="SECOND"
          case("MINUTE")  ;  T%TUnit="MINUTE"
          case("HOUR")    ;  T%TUnit="HOUR"
          case("DAY" )    ;  T%TUnit="DAY"
          case("YEAR")    ;  T%TUnit="YEAR"
          case default
            Ok= .false.
            call Stop_ &
            & (trim(W1)//" <-unknown TUNIT in DYNAMIC.TIME / DYNAMIC")
          end select
        !
        case("TFIN")   ;  call WrdToReal(W2,T%TFinal) !duration of simulation!
        case("TFINAL") ;  call WrdToReal(W2,T%TFinal) !duration of simulation!
        case("DTIME")  ;  call WrdToReal(W2,T%dTime)  !initial time step     !
        case("DTMIN")  ;  call WrdToReal(W2,T%dTMin)  !minimal time step     !
        case("DTMAX")  ;  call WrdToReal(W2,T%dTMax)  !maximal time step     !
        case("DTSAV")  ;  call WrdToReal(W2,T%dTSav)  !edition time step     !
        !! case default
        !!   Ok= .false.
        !!   call Stop_ &
        !!   & (trim(W1)//" <-unknown KeyWord in DYNAMIC.TIME or DYNAMIC")
        !
        !---------------------------------------------/timing parameters
        end select
        !
      end do Do02
    !_case("DYNAMIC")
    end select
  end do Do01
  call closeFILE(F)
  !
  if(T%dTSav==Zero) T%dTSav= T%TFinal/2.E3
  !
  !------------------------ if user given max time step is too small ---
  if(T%dTMax>T%TFinal/100._dp) T%dTMax= T%TFinal/100._dp
  !-------------------- max time step set to 1/100 of total duration ---
  if(T%dTMax==Zero) T%dTMax= T%TFinal/100._dp
  !
  T%TimeFactor= Dynam_TimeFactor(T%TUnit)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_ReadTime"
  !
  return
end subroutine Dynam_ReadTime

subroutine Dynam_ReadRock_(Ok)
!--
!-- read block DYNAMIC.ROCK
!--
  use M_Global_Vars, only: vSpc,vFas,vKinFas,vKinModel
  use M_KinFas_Read, only: T_LnkKin,KinFas_BuildLnk,KinFas_LnkToVec
  use M_Dynam_Vars,  only: sModelSurf
  !
  logical,intent(out):: Ok
  !
  type(T_LnkKin),pointer:: LnkKin
  real(dp):: x
  integer :: I,N
  !
  !------------------------------------------------------read ROCK block
  !
  !if(allocated(vKinFas)) deallocate(vKinFas)
  !allocate(vKinFas(100))
  !
  call KinFas_BuildLnk( &  !
  & vFas,vKinModel,sModelSurf, & !IN
  & N,LnkKin)                    !OUT
  !
  if(N>0) then !mod 12/06/2008 18:15
    Ok=.true.
    if(allocated(vKinFas)) deallocate(vKinFas)
    allocate(vKinFas(N))
    call KinFas_LnkToVec(LnkKin,vFas,vKinFas) !IN,OUT
  else
    Ok=.false.
    return !------------------------------------------------------return
  end if
  !
  if(iDebug>2) then
    !
    print '(/,A)',"Checking input from DYNAMIC.ROCK:"
    !
    if(size(vKinFas)>0) then
      
      print '(5A16)', &
      & "Kinetic_Phase___", &
      & "Thermo_Model____", &
      & "Kinetic_Model___", &
      & "__-___Radius/m__", &
      & "________Surf/g__"
      
      do I=1,size(vKinFas)
        ! if(vKinFas(I)%iFas<1) &
        ! & call Stop_("problem with "//trim(vKinFas(I)%Name))
        if(vKinFas(I)%iKin>0) then
          print '(3(A15,1X),2(G15.3,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & vKinModel(vKinFas(I)%iKin)%Name, &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg
        else
          print '(3(A15,1X),2(G15.6,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & "EQUILIBRIUM", &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg/1.0D3
        end if
      end do
    
    else
      
      print *,"EMPTY BOX ... !!"
    
    end if
    !
    call Pause_
  end if
  !
  !-----------------------------------------------------/read ROCK block
  !
  !------------------------------- normalize PhiM to sum(VolFract)-1 ---
  x= sum(vKinFas(:)%Dat%PhiM) !, mask=vKinFas(:)%Dat%cSat /= "MINIMAL")
  if(x>Zero) then
    do I=1,size(vKinFas)
      if(vKinFas(I)%Dat%cSat /= "M") then           ! M(INIMAL
        vKinFas(I)%Dat%PhiM= vKinFas(I)%Dat%PhiM /x
      else
        vKinFas(I)%Dat%PhiM= Zero
      end if
    end do
  else
    call Stop_("sum(vKinFas(1:nMk)%PhiM <0 ???") !------------------stop
  end if
  !---------------------------------------------------------------------
  !
  return
end subroutine Dynam_ReadRock_

subroutine Dynam_ReadRock(Ok)
!--
!-- read block DYNAMIC.ROCK
!--
  use M_Global_Vars, only: vSpc,vFas,vKinFas,vKinModel
  use M_KinFas_Read, only: T_LnkKin,KinFas_BuildLnk,KinFas_LnkToVec
  use M_Dynam_Vars,  only: sModelSurf
  !
  logical,intent(out):: Ok
  !
  type(T_LnkKin),pointer:: LnkKin
  real(dp):: x
  integer :: I,N
  !
  !------------------------------------------------------read ROCK block
  !
  call KinFas_BuildLnk( &  !
  & vFas,vKinModel,sModelSurf, & !IN
  & N,LnkKin)                    !OUT
  !
  if(N>0) then !mod 12/06/2008 18:15
    Ok=.true.
    if(allocated(vKinFas)) deallocate(vKinFas)
    allocate(vKinFas(N))
    call KinFas_LnkToVec(LnkKin,vFas,vKinFas) !IN,OUT
  else
    Ok=.false.
    return !------------------------------------------------------return
  end if
  !
  if(iDebug>2) then
    !
    print '(/,A)',"Checking input from DYNAMIC.ROCK:"
    !
    if(size(vKinFas)>0) then
      
      print '(5A16)', &
      & "Kinetic_Phase___", &
      & "Thermo_Model____", &
      & "Kinetic_Model___", &
      & "__-___Radius/m__", &
      & "________Surf/g__"
      
      do I=1,size(vKinFas)
        ! if(vKinFas(I)%iFas<1) &
        ! & call Stop_("problem with "//trim(vKinFas(I)%Name))
        if(vKinFas(I)%iKin>0) then
          print '(3(A15,1X),2(G15.3,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & vKinModel(vKinFas(I)%iKin)%Name, &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg
        else
          print '(3(A15,1X),2(G15.6,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & "EQUILIBRIUM", &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg/1.0D3
        end if
      end do
    
    else
      
      print *,"EMPTY BOX ... !!"
    
    end if
    !
    call Pause_
  end if
  !
  !-----------------------------------------------------/read ROCK block
  !
  !------------------------------- normalize PhiM to sum(VolFract)-1 ---
  x= sum(vKinFas(:)%Dat%PhiM) !, mask=vKinFas(:)%Dat%cSat /= "MINIMAL")
  if(x>Zero) then
    do I=1,size(vKinFas)
      if(vKinFas(I)%Dat%cSat /= "M") then           ! M(INIMAL
        vKinFas(I)%Dat%PhiM= vKinFas(I)%Dat%PhiM /x
      else
        vKinFas(I)%Dat%PhiM= Zero
      end if
    end do
  else
    call Stop_("sum(vKinFas(1:nMk)%PhiM <0 ???") !------------------stop
  end if
  !---------------------------------------------------------------------
  !
  return
end subroutine Dynam_ReadRock

subroutine Dynam_ReadBox(Box,Ok,Msg)
!--
!-- read block DYNAMIC.BOX -> spatial conditions
!--
  use M_Files,only: NamFInn
  use M_IOTools
  use M_Dynam_Vars,only: T_DynBox,sModelSurf
  !
  type(T_DynBox),intent(inout):: Box
  logical,       intent(out)  :: Ok
  character(*),  intent(out)  :: Msg
  !
  character(len=255):: L
  character(len=80) :: W,W1,W2
  logical           :: EoL
  integer           :: F,ios
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_ReadBox"
  !
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))
  Ok=.false.
  Msg= "Ok"
  !
  Do01: do
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
    call LinToWrd(L,W,EoL); call AppendToEnd(L,W,EoL)
    if(W(1:1)=='!')   cycle Do01 !skip comment lines
    !
    select case(trim(W))
    !
    case("ENDINPUT"); exit  Do01
    !
    case("DYNAMIC.BOX","DYNAMIC")
      Ok=.true.
      Do02: do
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
        call LinToWrd(L,W1,EoL); call AppendToEnd(L,W1,EoL)
        if(W1(1:1)=='!') cycle Do02
        call LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        select case(trim(W1))
        !
        case("ENDINPUT"); exit Do01
        case("END","ENDDYNAMIC.BOX","ENDDYNAMIC"); exit Do02
        !
        case("MODEL","SURFACE")
          select case(trim(W2))
          case("SPHERE") ; sModelSurf= "SPHERE" ! (= the default model)
          case("CRUNCH") ; sModelSurf= "CRUNCH"
          case default
            Ok= .false.
            Msg="In DYNAMIC / DYNAMIC.BOX, "//trim(W2)//" is Unknown Keyword"
          end select
          !
        case("VOLUME")
          select case(trim(W2))
            case("FREE");  Box%VFixed= .false. !UpdateMassFluid=.false.
            case("FIXED"); Box%VFixed= .true.  !UpdateMassFluid=.true.
            case default
              Ok= .false.
              Msg= "In DYNAMIC / DYNAMIC.BOX, VOLUME is either FREE or FIXED"
          end select
        !---------------------------------------------------------------------!
        !--------------------------------------------------spatial parameters.!
        case("DX");       call WrdToReal(W2,Box%dX)     !length of box........!
        case("VOLBOX");   call WrdToReal(W2,Box%VBox)   !volume of box........!
        case("UDARCY");   call WrdToReal(W2,Box%UDarcy) !flux rate,length/time!
        case("POROSITY"); call WrdToReal(W2,Box%PhiF)   !initial porosity.....!
        !---------------------------------------------------------------------!
        !
        case default
          Ok= .false.
          Msg= "In DYNAMIC / DYNAMIC.BOX, "//trim(W1)//" is unknown KeyWord !!"
        !
        end select
      end do Do02
    !_case("DYNAMIC")
    end select
  end do Do01
  call closeFILE(F)
  !
  if(.not. Ok) Msg= "Block DYNAMIC.BOX not found, using default values !!!"
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_ReadBox"
  !
end subroutine Dynam_ReadBox

subroutine Dynam_ReadColumn(Column,Ok,Msg)
!--
!-- read block DYNAMIC.1D
!--
  use M_Files,only: NamFInn
  use M_IOTools
  use M_Dynam_Vars,only: T_DynColumn
  !
  type(T_DynColumn),intent(inout):: Column
  logical,          intent(out)  :: Ok
  character(*),     intent(out)  :: Msg
  !
  character(len=255):: L
  character(len=80) :: W,W1,W2
  logical           :: EoL
  integer           :: F,ios
  
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_ReadColumn"
  !
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))
  !
  Ok=.false.
  Msg= "Ok"
  !
  Do01: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
    call LinToWrd(L,W,EoL); call AppendToEnd(L,W,EoL)
    if(W(1:1)=='!') cycle Do01 !skip comment lines
    !
    select case(trim(W))
    !
    case("ENDINPUT"); exit  Do01
    !
    case("DYNAMIC.COLUMN")
      Ok=.true.
      Do02: do
        !
        read(F,'(A)',iostat=ios) L
        if(ios/=0) exit Do01
        !
        call LinToWrd(L,W1,EoL)
        call AppendToEnd(L,W1,EoL)
        !
        if(W1(1:1)=='!') cycle Do02
        !
        call LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        select case(trim(W1))
        !
        case("ENDINPUT"); exit Do01
        case("END","ENDDYNAMIC.COLUMN"); exit Do02
        !
        !-------------------------------------------spatial parameters--
        case("DX")      ;   call WrdToReal(W2,Column%dX)
        case("DT")      ;   call WrdToReal(W2,Column%dT)
        case("DISP")    ;   call WrdToReal(W2,Column%Disp)
        case("UDARCY")  ;   call WrdToReal(W2,Column%UDarcy)
        case("TFIN")    ;   call WrdToReal(W2,Column%Duration)
        case("Tsave")   ;   call WrdToReal(W2,Column%Time_Save)
        case("NCELL")   ;   call WrdToInt (W2,Column%nCell)
        case("METHOD")  ;   call WrdToInt (W2,Column%Method)
        !---------------------------------------------------------------
        !
        case default
          Ok= .false.
          Msg= "In DYNAMIC.COLUMN, "//trim(W1)//" is unknown KeyWord !!"
        !
        end select
      end do Do02
    !
    end select
    !
  end do Do01
  !
  call closeFILE(F)
  
  if(.not. Ok) Msg= "Block DYNAMIC.COLUMN not found !!!"
  
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_ReadColumn"
  
  return
end subroutine Dynam_ReadColumn

subroutine Dynam_ReadNumeric(Ok)
!--
!-- read block DYNAMIC.NUMERIC -> numerical options
!--
  use M_Files,only: NamFInn
  use M_IOTools
  use M_Dynam_Vars,only: &
  & iCtrlTime,cMethod,Implicit_Surface, &
  & LogForMin,LogForAqu,DirectSub, &
  & bFinDif,DebNewt,TestJacob,TestMax, &
  & NewtMaxIts,NewtIterMax,NewtIterMin, &
  & NewtTOLF,NewtTOLX
  !
  logical,intent(out)::Ok
  !
  character(len=255):: L
  character(len=255) :: W0,W1,W2
  logical           :: EoL
  integer           :: F,ios
  !
  call GetUnit(F)
  call OpenFile (F,file=trim(NamFInn))
  Ok=.false.
  
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_ReadNumeric"
    
  Do01: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
    call LinToWrd(L,W0,EoL); call AppendToEnd(L,W0,EoL)
    if(W0(1:1)=='!')   cycle Do01 !skip comment lines
    !
    select case(trim(W0))
    !
    case("ENDINPUT")  ;  exit Do01
    !
    !------------------------------------- param's for numeric method --
    case("DYNAMIC.NUMERIC","DYN.NUMERIC")
      !
      Ok=.true.
      !
      DoReadNum: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit Do01
        !read(F,'(A)') L
        call LinToWrd(L,W1,EoL)
        if(W1(1:1)=='!') cycle DoReadNum
        call AppendToEnd(L,W1,EoL)
        !
        ! print *,'Dynam_ReadNumeric,W1=',trim(W1)   ;  pause
        
        select case(trim(W1))
        case("ENDINPUT")
          exit Do01
        case("END","ENDDYNAMIC.NUMERIC","ENDDYN.NUMERIC")
          exit DoReadNum
        end select
        !
        call LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        ! print *,'Dynam_ReadNumeric,W2=',trim(W2)  ;  pause
        
        select case(trim(W1))
        !!case("STEADY.stop");
        case("SURFACE")
          select case(trim(W2))
            case("IMPLICIT")  ; Implicit_Surface= .true.
            case("EXPLICIT")  ; Implicit_Surface= .false.
            case default
              call Stop_(trim(W2)//"-> unknown keyword (implicit,EXPLICIT)")
          end select
        !
        case("MAXITER")  ; call WrdToInt(W2,NewtMaxIts)
        case("ITERMAX")  ; call WrdToInt(W2,NewtIterMax)
        case("ITERMIN")  ; call WrdToInt(W2,NewtIterMin)
        !
        !if(iter>NewtIterMax) dT0=dT0/2.0D0
        !if(iter>NewtIterMax) dT0=dT0/2.0D0
        !if(iter<NewtIterMin) dT0=dT0*2.0D0
        !
        case("CONTROL")  ; call WrdToInt(W2,iCtrlTime)   !
        !
        case("METHOD")  !!then
          select case(trim(W2))
          case("NEWTON")       ;   cMethod= trim(W2)  ! iMethod= 1
          case("NEWTLNSRCH")   ;   cMethod= trim(W2)  ! iMethod= 2
          case("NEWTONPRESS")  ;   cMethod= trim(W2)  ! iMethod= 2
          case("BROYDEN")      ;   cMethod= trim(W2)  ! iMethod= 3
          case("NEWTONCHESS")  ;   cMethod= trim(W2)  ! iMethod= 4
          case("TENSOLVE_1")   ;   cMethod= trim(W2)  ! iMethod= 5
          case("NEWTONKELLEY") ;   cMethod= trim(W2)  ! iMethod= 6
          case("NEWTONWALKER") ;   cMethod= trim(W2)  ! iMethod= 7 !-> default value
          case default         ;   call Warning(W2)
          end select
        !_
        case("NEWTTOLF")
          !convergence criterion on function values
          call WrdToReal(W2,NewtTolF)
        !case("NEWTTOLMIN")
        !  !criterion for spurious convergence
        !  call WrdToReal(W2,NewtTolF)
        !case("NEWTTOLX")
        !  !convergence criterion on dx
        !  call WrdToReal(W2,NewtTolX)
        !
        !------------------------------------------------- obsolete --
        case("NOTLOG")  ;       LogForMin=.false.
        case("LOGMIN")  ;       LogForMin=.true.
        case("FINDif")  ;       bFinDif=  .true.
        !
        case("TESTMAX")  ;      TestMax=  .true.
        !------------------------------------------------/ obsolete --
        !
        case("LOGFORAQU")
          select case(trim(W2))
            case("TRUE", "T", "YES")  ;  LogForAqu= .true.
            case("FALSE", "F", "NO")  ;  LogForAqu= .false.
          end select

        case("LOGFORMIN")
          select case(trim(W2))
          case("TRUE", "T", "YES")  ;  LogForMin= .true.
          case("FALSE", "F", "NO")  ;  LogForMin= .false.
          end select

        case("DIRECT")
          select case(trim(W2))
          case("TRUE", "T", "YES")  ;  DirectSub= .true.
          case("FALSE", "F", "NO")  ;  DirectSub= .false.
          end select

        case("JACOBIAN")
          select case(trim(W2))
          case("NUMERIC")   ;  bFinDif= .true.
          case("ANALYTIC")  ;  bFinDif= .false.
          end select

        case default
          call Warning(W1)

        end select
        
        ! print *,"W1=",trim(W1)  ;  pause

      end do DoReadNum
    !_case("DYNAMIC.NUMERIC")
    !---/
    end select
    !
  end do Do01
  !
  call closeFILE(F)
  !
  DebNewt=   (iDebug>2)
  TestJacob= (iDebug>2)
  !
  if(.not. &
  &       (NewtMaxIts>NewtIterMax &
  & .and. NewtIterMax>NewtIterMin)) then
    if(idebug>1) write(fTrc,'(A)') &
    & "Must Have MAXITER>ITERMAX .and. ITERMAX>ITERMIN"
    call Stop_("Must Have MAXITER>ITERMAX .and. ITERMAX>ITERMIN")
  end if
  !
  if(iCtrlTime<1 .or. iCtrlTime>2) iCtrlTime=1
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_ReadNumeric"
  !
  return
end subroutine Dynam_ReadNumeric

subroutine Warning(W)
  character(len=*),intent(in):: W
  if(iDebug>2) &
  & write(fTrc,'(A)') "!!WARNING!! "//trim(W)//" <-unknown KeyWord !!"
  if(iDebug>2) then
    print '(A)',"!!WARNING!! "//trim(W)//" <-unknown KeyWord !!"
    call Pause_
  end if
  return
end subroutine Warning

subroutine Check_KinRate_CalcQsK
!--
!-- check kinetic models
!--
  use M_IoTools,   only: GetUnit, Closefile, OpenFile
  use M_Files,     only: DirLog,Files_Index_Write
  use M_Global_Vars,only: vKinModel, nAq
  use M_KinRate,   only: KinRate_CalcQsKFactor
  !
  real(dp):: Aff,Qsk,VmQsk
  integer :: I,J,F
  real(dp),dimension(1:nAq):: vDX,vX
  !
  call GetUnit(F)
  call OpenFile(F,file=trim(DirLog)//"kinrate_vmqsk.log")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirLog)//"kinrate_vmqsk.log",&
  & "DYNAMIC/LOG: check KinModels")
  !
  write(F,'(2(A,A1))',advance="NO") "Ord",T_,"Aff",T_
  do I=1,size(vKinModel)
    write(F,'(A,A1)',advance="NO") vKinModel(I)%Name,T_
  end do
  write(F,*)
  !
  Aff= Zero
  vX=  Zero
  J=1
  do
    write(F,'(I3,A1,G15.6,A1)',advance="NO") J,T_,Aff,T_
    QsK=exp(-Aff)
    do I=1,size(vKinModel)
      !
      call KinRate_CalcQsKFactor(&
      & "D",          & !IN was "DISSOLU"
      & vKinModel(I), & !IN: Kinetic model
      & QsK,          & !IN
      & vX,           & !IN !dQsKdLnXi
      & VmQsK,        & !OUT
      & vDX)            !OUT dVmQdLnX_M
      !
      write(F,'(G15.6,A1)',advance="NO") VmQsK,T_
    end do
    write(F,*)
    !
    Aff= Aff +0.1_dp; J=J+1
    if(Aff>6.0_dp) exit
    !
  end do
  !
  call closeFILE(F)
  !
end subroutine Check_KinRate_CalcQsK

end module M_Dynam_Read

