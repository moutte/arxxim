module M_TPCond_Read
  use M_Kinds
  use M_Trace

  implicit none

  private
  !
  public:: TPpath_Read
  public:: TPgrid_Read
  public:: TPgrid_Build
  !
contains

subroutine TPpath_Read(TdgK0,Pbar0)
!--
!-- scan the input file, reads block TP.TABLE, if present
!-- -> allocate & build vTPpath
!--
  use M_IoTools !,only:GetUnit,dimV
  use M_Files,     only: NamFInn
  !
  use M_Dtb_Const, only: T_CK,Tref,Pref
  use M_Dtb_Calc,  only: Dtb_TP_Check
  use M_T_Tpcond,  only: T_TPCond
  use M_Fluid_Calc,only: Eos_H2O_Psat
  !
  use M_Dtb_Vars,  only: DtbFormat,DtbLogK_vTPCond,Psat_Auto
  use M_Path_Vars, only: vTPpath,DimPath
  !
  real(dp),intent(in),optional:: TdgK0,Pbar0
  != default values when no TP.TABLE or TP.PATH is found
  !
  logical :: EoL,Ok
  integer :: i,mDum,ios,f,N !,m1,m2
  real(dp):: TdgK,Pbar
  character(len=512):: L,W,W1
  character(len=80) :: Msg
  !
  real(dp):: vX(dimV)
  type(T_TPCond):: vCond(dimV)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< TPpath_Read"
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  Ok=.false.
  !
  vCond(:)%Name= "x"
  vCond(:)%Pbar= Pref
  !
  DoFile: do
    !
    read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit  DoFile
    !
    case("TP.TABLE","TPTABLE","TP.PATH")

      if(iDebug>0 .and. &
      &   DtbFormat=="LOGKTBL" .and. &
      &  (trim(W)=="TP.TABLE" .or. trim(W)=="TPTABLE")) then
        call Warning_("with logK table, use TP.PATH, to distinguish from TP.TABLE")
      end if

      Ok=  .true.
      N= dimV
      DoTPtable: do

        read(f,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,Eol)
        if(W(1:1)=='!') cycle DoTPtable
        !write(fTrc,'(A)') L
        call AppendToEnd(L,W,EoL)

        select case(W)

        !! case("ENDINPUT"); exit DoFile
        case("END","ENDTPTABLE","ENDTP.TABLE","ENDTP.PATH")
          exit DoTPtable
        !---------------!

        case("TDGC", "TDGK", "PBAR","PMPA", "TDEGC","TDEGK")
          !
          !! call ReadRValsV(L,mDum,vX)
          !
          i=0
          do
            call LinToWrd(L,W1,EoL)
            i=i+1
            if(i>DimV) exit
            if(trim(W1)=="PSAT") then
              if(vCond(i)%TdgC>Zero) then
                call Eos_H2O_Psat(vCond(i)%TdgC+T_CK,vX(i))
                !print *,"TdgC= ", vCond(i)%TdgC
                !print *,"psat= ", vX(i)  ;  call pause_
              else
                call Stop_("error in TP.TABLE: invalid TdgC for Psat computation")
              end if
            else
              call WrdToReal(W1,vX(i))
            end if
            if(EoL) exit
          end do
          !
          mDum= i
          N=MIN(N,mDum)
          if(iDebug>0) &
          & write(fTrc,'(A,A1,A,2I3)') trim(W),T_," DIM= ", mDum, N
          !
          select case(trim(W))
          !
          case("TDGC");  vCond(:)%TdgC= vX(:)
          case("TDGK");  vCond(:)%TdgC= vX(:)-T_CK
          !
          case("PBAR");  vCond(:)%Pbar= vX(:)       !pressure in bar
          case("PMPA");  vCond(:)%Pbar= vX(:)*1.D-1 !pressure in MegaPascal
          !
          case("TDEGC"); vCond(:)%TdgC=vX(:)
            call Warning_("Depreciated Keyword TDEGC")
          case("TDEGK"); vCond(:)%TdgC=vX(:)-T_CK
            call Warning_("Depreciated Keyword TDEGK")
          !
          end select
        !---------------!

        case("NAME")
          i=0
          DoName: do
            call LinToWrd(L,W,EoL)
            i=i+1
            if(I>DimV) exit DoName
            vCond(i)%Name=trim(W)
            if(EoL) exit DoName
          end do DoName
        !---------------!

        case default
          call Stop_(trim(W)//"<<unknown Keyword...")

        end select !end if
        !end if
      end do DoTPtable
    end select
    !
  end do DoFile
  close(f)
  !
  if(Ok) then

    do I=1,N
      TdgK= vCond(I)%TdgC +T_CK
      Pbar= vCond(I)%Pbar
      !
      call Dtb_TP_Check(DtbFormat,DtbLogK_vTPCond,Psat_Auto,TdgK,Pbar,Ok,Msg)
      !
      if(.not. Ok) call Stop_(trim(Msg))
      !
      vCond(I)%TdgC= TdgK -T_CK
      vCond(I)%Pbar= Pbar
      !
      if(iDebug>0) &
      & write(fTrc,'(2(G12.3,1X))') vCond(I)%TdgC, vCond(I)%Pbar !!i,vCond(i)%Name
    end do

  else

    N=1
    if(present(TdgK0)) then ;  vCond(1)%TdgC= TdgK0 -T_CK
    else                    ;  vCond(1)%TdgC= Tref  -T_CK
    end if
    if(present(Pbar0)) then ;  vCond(1)%Pbar= Pbar0
    else                    ;  vCond(1)%Pbar= Pref
    end if

  end if
  !
  if(allocated(vTPpath)) deallocate(vTPpath)
  allocate(vTPpath(1:N))
  vTPpath(1:N)= vCond(1:N)
  !
  DimPath= N
  !
  !!if(.not.Ok) call Stop_("Block TP.PATH not Found ...!!!")
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ TPpath_Read"
  !
end subroutine TPpath_Read

subroutine TPgrid_Read( &
& Ok, &
& T_Min,T_Max,T_ratio,T_delta, &
& P_Min,P_Max,P_ratio,P_delta)
  !
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_Dtb_Const,only: T_CK
  use M_Files,    only: NamFInn
  !
  logical, intent(out):: Ok
  real(dp),intent(out):: T_Min,T_Max,T_ratio,T_delta,P_Min,P_Max,P_ratio,P_delta
  !
  character(len=512):: L,W,W1
  logical :: EoL
  integer :: F,ios
  real(dp):: rBegin,rFinal,rRatio,rDelta
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Dtb_Read_TPgrid"
  !
  Ok= .false.
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  !
  DoFile: do

    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines

    call AppendToEnd(L,W,EoL)
    select case(W)
    !
    case("ENDINPUT"); exit DoFile
    !
    case("TP.GRID","TPGRID") !!!!!!canevas for reading one "block"
      !... I=0
      Ok= .true.
      !
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        select case(W)
        !
        case("ENDINPUT"); exit DoFile
        !
        case("END","ENDTP.GRID","ENDTPGRID"); exit DoBlock
        !
        case default
          call Pause_("< WARNING!!! "//trim(W)//"=UNKNOWN keyword in TPgrid")
        !
        case("TDGC","TDGK","PBAR")
          !
          rBegin=Zero; rFinal=Zero; rRatio=One; rDelta=Zero
          !
          DoLine: do
            if(EoL) exit DoLine
            call LinToWrd(L,W1,EoL)
            select case(W1)
            case("INITIAL") ; call LinToWrd(L,W1,EoL); call WrdToReal(W1,rBegin)
            case("FINAL")   ; call LinToWrd(L,W1,EoL); call WrdToReal(W1,rFinal)
            case("RATIO")   ; call LinToWrd(L,W1,EoL); call WrdToReal(W1,rRatio)
            case("DELTA")   ; call LinToWrd(L,W1,EoL); call WrdToReal(W1,rDelta)
            case default    ; call Stop_(W1//"= unknown keyword in TPgrid !!!") !stop
            end select
          end do DoLine
          !
          select case(W)
          case("TDGC")
            T_Min=   rBegin ;      T_Max=   rFinal
            T_ratio= rRatio ;      T_delta= rDelta
          case("TDGK")
            T_Min=   rBegin-T_CK;  T_Max=   rFinal-T_CK
            T_ratio= rRatio ;      T_delta= rDelta
          case("PBAR")
            P_Min=   rBegin ;      P_Max=   rFinal
            P_ratio= rRatio ;      P_delta= rDelta
          end select !case(W)
        !
        end select !case(W)
        !
      end do DoBlock
    !
    end select
  end do DoFile
  close(F)

  if(iDebug>0) write(fTrc,'(A,/)') "</ Dtb_Read_TPgrid"

  return
end subroutine TPgrid_Read

subroutine TPgrid_Build(Ok)
  use M_T_TPcond,   only: T_TPCond
  use M_Dtb_Const,  only: T_CK
  use M_Fluid_Calc, only: Eos_H2O_Psat
  !
  use M_Path_Vars,  only: vTPpath
  !
  logical,intent(out):: Ok
  !
  real(dp):: T_Min, T_Max, T_delta, T_ratio
  real(dp):: P_Min, P_Max, P_delta, P_ratio
  real(dp):: TdgC,TdgK,Pbar,Psat
  !
  integer:: N
  !
  call TPgrid_Read( Ok, &
  & T_Min,T_Max,T_ratio,T_delta, &
  & P_Min,P_Max,P_ratio,P_delta)
  !
  ! when a specific TP.GRID block is not found,
  ! use a simple TP.TABLE
  if(.not. Ok) then
    != if TP.GRID block not found
    !N= size(vTPcond)
    !if(allocated(vTPgrid)) deallocate(vTPgrid)
    !allocate(vTPgrid(1:N))
    !vTPgrid= vTPcond
    return
  end if
  !
  if(iDebug==4) then
    print *,"TPgrid_Build"
    print '(A,4G15.6)',"T_Min,T_Max,T_ratio,T_delta",T_Min,T_Max,T_ratio,T_delta
    print '(A,4G15.6)',"P_Min,P_Max,P_ratio,P_delta",P_Min,P_Max,P_ratio,P_delta
    call Pause_
  end if

  !-------------------------------- first, count number of T,P points --
  N= 0
  TdgC= T_Min
  do
    Pbar= P_Min
    do
      TdgK= TdgC +T_CK
      Psat= Zero
      if(TdgK<=647.25D0) call Eos_H2O_Psat(TdgK,Psat)
      if(Pbar>Psat) N= N+1
      if(P_ratio>One) then ; Pbar= Pbar * P_ratio
      else                 ; Pbar= Pbar + P_delta
      end if
      if(Pbar>P_Max) exit
    end do
    TdgC= TdgC + T_delta
    if(TdgC>T_max) exit
  end do
  !
  if(allocated(vTPpath)) deallocate(vTPpath)
  allocate(vTPpath(1:N))

  !----------------------------------------------- then, fill vTPpath --
  TdgC= T_Min
  N= 0
  do
    Pbar= P_Min
    do
      TdgK= TdgC +T_CK
      !print '(F7.2,1X,F7.2)',TdgC,Pbar
      Psat= Zero
      if (TdgK<=647.25D0) call Eos_H2O_Psat(TdgK,Psat)
      if(Pbar>Psat) then
        N= N+1
        vTPpath(N)%TdgC= TdgC
        vTPpath(N)%Pbar= Pbar
      end if
      if(P_ratio>One) then ; Pbar= Pbar * P_ratio
      else                 ; Pbar= Pbar + P_delta
      end if
      if(Pbar>P_Max) exit
    end do
    TdgC= TdgC + T_delta
    if(TdgC>T_max) exit
  end do

  return
end subroutine TPgrid_Build

end module M_TPcond_Read
