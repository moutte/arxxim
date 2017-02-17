module M_Dtb_Read
!--
!-- reading databases
!--
  use M_Kinds
  use M_IoTools
  use M_Trace,        only: iDebug,T_,fTrc,fHtm,Stop_,Pause_, Warning_
  !
  use M_T_DtbAquHkf,  only: T_DtbAquHkf
  use M_T_DtbMinHkf,  only: T_DtbMinHkf
  use M_T_DtbMinThr,  only: T_DtbMinThr
  use M_T_DtbLogKTbl, only: T_DtbLogKTbl
  !
  implicit none
  !
  private
  !
  public:: Dtb_Read_Species
  public:: Dtb_Read_Format
  !
contains

subroutine Dtb_Read_Species(vEle)
!--
!-- read species block
!-- -> if species logk exit
!--    else build vDtbAquHkf, vDtbMinThr, etc.
!--
  use M_T_Element,only: T_Element
  use M_Files,    only: NamFInn
  !
  use M_Dtb_Read_Vars
  use M_Dtb_Read_DtbMinHkf
  use M_Dtb_Read_DtbAquHkf
  use M_Dtb_Read_DtbMinThr
  use M_Dtb_Read_DtbLogKTbl
  use M_Dtb_Read_DtbLogKAnl
  !
  use M_Dtb_Vars, only: Dtb_Vars_Zero, DtbFormat
  !
  type(T_Element),dimension(:),intent(in):: vEle
  !
  character(len=255):: L,W
  logical           :: EoL
  integer:: F,ios,N
  !
  if(iDebug>0) write(fTrc,"(/,A)") "< Dtb_Read_Species"
  !
  call Dtb_Vars_Zero
  !
  call Dtb_Read_Format
  !
  if(DtbFormat=="LOGKTBL") then
    call DtbLogK_TPCond_Read(N)
    if(N==0) call Stop_("TP table not found for LogKTbl database")
    !! print *,"DtbLogK_Dim=",N ; pause_
  end if
  !
  call GetUnit(F)  ; open(F,file=trim(NamFInn))
  !
  !-- dimensions and lined lists initialized outside loop
  !-- -> should allow reading several successive SPECIES blocks
  !
  nAquHkf=  0 ; nullify(LisAquHkf)
  nMinHkf=  0 ; nullify(LisMinHkf)
  nMinThr=  0 ; nullify(LisMinThr)
  nLogKTbl= 0 ; nullify(LisLogKTbl)
  nLogKAnl= 0 ; nullify(LisLogKAnl)
  !
  DoFile: do 
    
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    
    select case(W)
      
      case("SPECIES")
        
        if(EoL) then
          !LogK TP.table is considered the default database format
          DtbFormat= "LOGKTBL"
        else
          !
          call LinToWrd(L,W,EoL)
          !
          select case(trim(W))
          case("LOGK","LOGK.TABLE")  ;  W= "LOGKTBL"
          case("LOGK.ANALYTIC")      ;  W= "LOGKANL"
          end select
          !
          DtbFormat= trim(W)
          !
        end if
        
        select case(DtbFormat) !each of them will end on "ENDSPECIES"
        !
        case("MIN.THR") ; call DtbMinThr_Read(f,vEle,nMinThr)      !-> LisMinThr
        case("GAS.THR") ; call DtbMinThr_Read(f,vEle,nMinThr)      !-> LisMinThr
        case("HSV.THR") ; call DtbMinThr_Read_Line(f,vEle,nMinThr) !-> LisMinThr
        !---------------------------------------------------------------
        case("AQU.HKF") ; call DtbAquHKF_Read(f,vEle,nAquHkf)      !-> LisAquHkf
        !---------------------------------------------------------------
        case("MIN.HKF") ; call DtbMinHKF_Read_Old(f,vEle,nMinHkf) !-> LisMinHkf
        case("GAS.HKF") ; call DtbMinHKF_Read_Old(f,vEle,nMinHkf) !-> LisMinHkf
        case("HSV.HKF") ; call DtbMinHKF_Read(f,vEle,nMinHkf)     !-> LisMinHkf
        !---------------------------------------------------------------
        case("LOGKTBL") ; call DtbLogKTbl_Read(f,vEle,nLogKTbl) !-> LisLogKTbl
        case("LOGKANL") ; call DtbLogKAnl_Read(f,vEle,nLogKAnl) !-> LisLogKTbl
        !---------------------------------------------------------------
        case default
          call Stop_(trim(DtbFormat)//" <-unknown format !!??")
        !
        end select
        
    end select
  
  end do DoFile
  close(F)
  !
  if(fLin/=0) close(fLin)
  !
  !! if(nAquHkf==0) call Stop_("Found NO Aqu.Species ...") !________stop
  !
  if(nAquHkf>0)  call DtbAquHKF_Build  !-> vDtbAquHkf
  if(nMinHkf>0)  call DtbMinHKF_Build  !-> vDtbMinHkf
  if(nMinThr>0)  call DtbMinThr_Build  !-> vDtbMinThr
  if(nLogKTbl>0) call DtbLogKTbl_Build !-> vDtbLogKTbl
  if(nLogKAnl>0) call DtbLogKAnl_Build !-> vDtbLogKAnl
  !
  !if(iDebug>0) then
  !  if(nAquHkf>0) call DtbAquHKF_ShoStoik(vEle)
  !  if(nMinHkf>0) call DtbMinHKF_ShoStoik(vEle)
  !  if(nMinThr>0) call DtbMinTHR_ShoStoik(vEle)
  !end if
  !
  if(iDebug>2) call Dtb_Read_Check
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Dtb_Read_Species"
  
  return
end subroutine Dtb_Read_Species

subroutine Dtb_Read_Format
!--
!-- scan the input file for all SPECIES block
!-- if at least one SPECIES block is of log format,
!-- then DtbFormat is logK,
!-- and any SPECIES block of different format is discarded
!--

  use M_Files,only: NamFInn
  use M_Dtb_Vars,only: DtbFormat
  
  character(len=255):: L,W
  logical           :: EoL
  integer:: ios, f
  
  call GetUnit(F); open(F,file=trim(NamFInn))
  !
  DoFile: do 
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    
    select case(W)
    !
    case("ENDINPUT")
      exit DoFile
    !
    case("SPECIES")
      if(EoL) then
        DtbFormat= "LOGKTBL"
        ! <-LogK TP.table is considered the default database format
      else
        call LinToWrd(L,W,EoL)
        DtbFormat= trim(W)
        if(trim(DtbFormat)=="LOGK") DtbFormat= "LOGKTBL"
      end if
      if(DtbFormat=="LOGKTBL") then
        close(F)
        return
      end if
    !
    end select
  end do DoFile
  !
  close(F)
  !
  return
end subroutine Dtb_Read_Format

subroutine Dtb_Read_Check
  use M_Dtb_Vars
  use M_Dtb_Write
  !
  if(size(vDtbAquHkf)>0) call DtbAquHkf_Write !check results DtbAquHKF_Read
  if(size(vDtbMinHkf)>0) call DtbMinHkf_Write !check results DtbMinHkf_Read
  if(size(vDtbMinThr)>0) call DtbMinThr_Write !check results DtbMinThr_Read
  !
end subroutine Dtb_Read_Check

end module M_Dtb_Read

!comments from deCapitani files
!
!THE NAMES OF THE GASES CAN NOT BE CHANGED
!AS THEY ARE useD TO RECOGNIZE WHICH GAS ROUTINE SHOULD BE useD.
!ST=  code for standard STATE PROPERTIES - G, H, S, V
!   UNITS: J, J, J/DEG, J/bar
!   NOTE: G=  H - T*S  where S is third law entropy,
!   not entropy of formation; i.e. it"s an apparent G of formation
!CP calculated with eqn:
!   Cp=  K0 + K1/SQRT(T) + K2/T/T + K3/T/T/T + K4/T + K5*T + K6*T*T
!C1=  CODE FOR THE FIRST LINE OF CP TERMS - K0, K1, K2, K3
!C2=  CODE FOR THE 2ND LINE OF CP TERMS - K4, K5, K6
!C3=  CODE FOR MAIER-KELLY CP EQUATION - K0, K5, K2
!   UNITS: J AND DEGREES
!CP(disorder) calculated by:
!   Cp(dis)=  D0 + D1/SQRT(T) + D2/T/T + D3/T + D4*T + D5*T*T
!   To get H, S (disorder), this eqn integrated by TMIN and TMAX
!D1=  CODE FOR THE FIRST LINE OF DISORDER TERMS - D0, D1, D2, D3
!   UNITS: J  AND DEGREES
!D2=  CODE FOR THE SECOND LINE OF DISORDER TERMS - D4, D5, TMIN, TMAX
!   UNITS: J AND DEGREES
!LAMBDA HEAT CAPACITY EQUATION IS:
!   CP(LAMBDA)=  T * (L1 + L2 * T)**2
!T1 (T3,T5 FOR MORE TRANSITIONS)=  CODE FOR THE FIRST LINE OF TRANSITION TERMS
!   T TRANS(K), REFERENCE T FOR INTEGRATION, LAMBDA CP TERM L1, LAMBDA CP TERM L2, HEAT OF TRANSITION
!   UNITS: J AND DEGREES 
!   (SEE BERMAN AND BROWN, 1985, CMP, 89, 168-183)
!T2 (T4,T6 FOR MORE TRANSITIONS)=  CODE FOR THE SECOND LINE OF TRANS TERMS 
!   DT/DP SLOPE OF TRANSITION, parameter FOR HIGHER ORDER FIT OF TRANSITION T AS function OF P, 
!   DV/DT FOR LOW T POLYMORPH, DV/DP FOR LOW T POLYMORPH
!   UNITS: JOULES, BARS, DEGREES
!VOLUME EQUATION:
!   V(P,T)/V(1,298)=  1 + V1(T-298) + V2(T-298)**2 + V3(P-1) + V4(P-1)**2
!V1=  CODE FOR VOLUME parameterS - V1, V2,V3,V4
!   NOTE: V1,V2,V3 terms need to be divided by 10**5, V4 by 10**8
