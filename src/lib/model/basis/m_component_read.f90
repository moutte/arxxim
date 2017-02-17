module M_Component_Read
  use M_Kinds
  use M_T_Component,only: T_Component
  use M_Trace
  !
  implicit none

  private

  public:: Components_Read
  !
contains

subroutine Components_Read( & !
!--
!-- scan the SYSTEM block -> build vCpn, read (T,P) conditions --
!--
& sKeyWord,  & !IN
& vEle,      & !IN
& vSpc,      & !IN
& vMixFas,   & !IN
& TdgK,Pbar, & !INOUT
& N,         & !OUT
& SysType,   & !OUT
& Ok,        & !OUT
& Msg,       & !OUT
& vCpn       & !INOUT
)
  !---------------------------------------------------------------------
  use M_IoTools
  use M_Files,      only: NamFInn
  use M_Dtb_Const,  only: T_CK
  use M_Numeric_Const,only: Ln10
  use M_T_Element,  only: T_Element,Element_Index
  use M_T_Species,  only: T_Species,Species_Index,Species_Rename
  use M_T_MixPhase, only: T_MixPhase,MixPhase_Index
  use M_T_Component,only: Component_Zero,CpnMolMinim
  !
  use M_Dtb_Vars,only: DtbLogK_vTPCond
  !
  use M_System_Vars,only: BufferIsExtern
  !---------------------------------------------------------------------
  character(len=*), intent(in):: sKeyWord
  type(T_Element),  intent(in):: vEle(:)
  type(T_Species),  intent(in):: vSpc(:)
  type(T_MixPhase), intent(in):: vMixFas(:)
  !
  real(dp),         intent(inout):: TdgK,Pbar
  integer,          intent(out)  :: N
  character(len=7), intent(out)  :: SysType
  logical,          intent(out)  :: Ok
  character(*),     intent(out)  :: Msg
  type(T_Component),intent(out)  :: vCpn(:)
  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------development
  logical,parameter:: MoleIsTrue= .false.  !.true. !
  != MoleIsTrue => when Statut is MOLE or GRAM
  != species stoichiometry is considered for
  != computation of element mole numbers
  !---------------------------------------------------------/development
  type(T_Component) :: Cpn
  character(len=255):: L,sListElem
  character(len=80) :: W0,W,V1,V2,Wnum,WType
  character(len=7)  :: Statut
  logical :: EoL
  logical :: bForceH2O =.false.
  integer :: f,ios,iTP
  integer :: I,iEl,iSp,iW,iH_,iO_,nEl
  real(dp):: X1,Eps
  real(dp),allocatable:: vX(:)
  integer, allocatable:: tStoikCp(:,:)
  !---------------------------------------------------------------------
  
  if(iDebug>0) write(fTrc,'(/,A)') "< Read_Input_System"
  !
  Ok= .true.
  Msg= "Ok"
  !
  call GetUnit(f)
  call OpenFile(f,file=trim(NamFInn))
  !
  iTP= 0
  SysType="Z"
  sListElem="" !used to prevent reading same element more than once
  !
  N=   0
  iH_= 0
  iO_= 0
  nEl= size(vEle)
  !
  allocate(tStoikCp(size(vCpn),0:nEl+1))
  tStoikCp(:,:)= 0 ! stoikiometry
  tStoikCp(:,0)= 1 ! divisor
  !
  DoFile: do
    !
    read(f,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W0,EoL)
    if(W0(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W0,EoL)
    !
    if(W0=="ENDINPUT") exit DoFile
    !!select case(W0)
    !!  case("SYSTEM","SYSTEM.AQUEOUS"); SysType="AQUEOUS"
    !!  case("SYSTEM.MOMAS");            SysType="MOMAS"
    !!end select
    if(trim(W0)==trim(sKeyWord)) then
      !
      !-- possibility to have several SYSTEM blocks, each with a Name,
      !-- -> to define several systems, make them react with another, etc.
      if ( trim(sKeyWord)=="SYSTEM"         .or. &
      &    trim(sKeyWord)=="SYSTEM.AQUEOUS" .or. &
      &    trim(sKeyWord)=="SYSTEM.BOX"     .or. &
      &    trim(sKeyWord)=="SYSTEM.MIX"     .or. &
      &    trim(sKeyWord)=="SYSTEM.INJECT")  &
      & SysType="AQUEOUS"
      !
      call Component_Zero(Cpn)
      !
      Statut= "INERT" !default value
      !
      if(SysType=="AQUEOUS") then
        !
        if(Element_Index("O__",vEle)==0) then
          Ok= .false.
          Msg= "O__ not found among elements !!!"
          return
        end if
        !
        if(Element_Index("H__",vEle)==0) then
          Ok= .false.
          Msg= "H__ not found among elements !!!"
          return
        end if
        !
        if(Species_Index("H2O",vSpc)==0) then
          Ok= .false.
          Msg= "H2O not found among species  !!!"
          return
        end if
        !
        iW= Species_Index("H2O",vSpc)
        !
        if(bForceH2O) then
          !---------------make water (solvent) as Cpn Nr1
          !---------------------------ele-"O__",spc-'H2O"
          iO_= 1
          Cpn%iEle=     Element_Index("O__",vEle)
          Cpn%NamCp=    vEle(Cpn%iEle)%NamEl
          Cpn%iSpc=     Species_Index("H2O",vSpc)
          Cpn%Mole=     One /vSpc(iW)%WeitKg
          Cpn%LnAct=    Zero
          Cpn%Statut=   "INERT"
          Cpn%namSol=   "Z"
          !
          sListElem=    "O__"
          N=N+1
          !
          vCpn(N)= Cpn
          !vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
          tStoikCp(N,Cpn%iEle)= 1
          !
        end if
        !
      end if
      !
      DoBlock: do
        !
        read(f,'(A)',iostat=ios) L
        if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL) !-> read 1st word -> W
        if(W(1:1)=='!') cycle DoBlock !-> skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        if(    W=="ENDINPUT" &
        & .or. W=="END"      &
        & .or. W=="END"//trim(sKeyWord)) exit DoFile
        !
        !------------------------------------read Temperature / Pressure
        !if(len_trim(W)>3) then
        if(IsKeyword("_TP.POINT_TPPOINT_TDGK_TDGC_PBAR_PMPA_TDEGK_TDEGC_",W)) &
        & then
          !
          call LinToWrd(L,Wnum,EoL)
          !-> !!!§§§ add something to check input !!
          select case(trim(W))
          !
          case("TP.POINT","TPPOINT")
            ! call Warning_(trim(W)//" is Obsolete, better give TdgC, Pbar !!!")
            !
            call WrdToInt(Wnum,iTP)
            !
            if(iTP>0) then
              !! call Warning_ &
              !! & ("keyword TP.POINT soon Obsolete !! better give TdgC, Pbar in SYSTEM")
              !
              if(iTP>size(DtbLogK_vTPCond)) then
                Ok= .false.
                Msg= "TP.POINT outside range"
                return
              end if
              !
              TdgK= DtbLogK_vTPCond(iTP)%TdgC +T_CK
              Pbar= DtbLogK_vTPCond(iTP)%Pbar
              !
            end if
          !
          case("TDGK","TDGC","TDEGK","TDEGC","PBAR","PMPA")
            call WrdToReal(Wnum,X1)
            !
            select case(trim(W))
            !
            case("TDGK");  TdgK= X1
            case("TDGC");  TdgK= X1 + T_CK
            ! print *,"debuggg Components_Read"
            ! print *,TdgK
            !
            case("PBAR");  Pbar= X1
            case("PMPA");  Pbar= X1 * 1.0D-1
            !
            case("TDEGK")
              TdgK= X1
              call Warning_("Depreciated Keyword TDEGK")
            case("TDEGC")
              TdgK= X1 + T_CK
              call Warning_("Depreciated Keyword TDEGC")
            !
            end select
          !
          case default
            Ok= .false.
            Msg= "Unknown (or Obsolete) Keyword: "//trim(W)
            return
          !
          end select
          !
          cycle DoBlock
          !
        end if
        !-----------------------------------/read Temperature / Pressure
        !
        !----------------------------------------------read element name
        if(.not. IsElement(vEle,W)) then
          Ok= .false.
          Msg=           "Element "//trim(W)//" Not in Element base !!!"
          return !------------------------------------------------return
        end if
        !
        call Str_Append(W,3)
        !
        !--position of current element in vEle
        iEl= Element_Index(trim(W),vEle)
        !
        if(index(sListElem,trim(W))>0) then
          Ok= .false.
          Msg=              "Component "//trim(W)//" Already listed !!!"
          return !------------------------------------------------return
        end if
        !
        N=N+1
        if(N>size(vCpn)) then
          Ok= .false.
          Msg=                                     "Too many components"
          return !------------------------------------------------return
        end if
        !
        Cpn%NamCp= trim(W)
        Cpn%iEle= iEl
        !Cpn%NamCp= trim(vEle(Cpn%iEle)%NamEl)
        !Cpn%vStoikCp(iEl)= 1
        tStoikCp(N,iEl)= 1
        !
        if(W=="H__") iH_= N
        if(W=="O__") iO_= N
        !
        sListElem=trim(sListElem)//trim(W)
        !-> prevent reading same component in subsequent loops
        !
        !---------------------------------------------/read element name
        !
        !------------------------------------------read component status
        call LinToWrd(L,WType,EoL)
        !
        !if(.not. IsKeyword &
        !& ("_INERT_MOLE_GRAM_MOBILE_ACTIVITY_PK_BUFFER_BALANCE_",WType)) then
        !  Ok= .false.
        !  Msg=         trim(WType)//" -> unrecognized species status !!"
        !  return !------------------------------------------------return
        !end if
        !
        select case(trim(WType))
        case("INERT","MOLE","GRAM","PPM")        ;  Statut="INERT"
        case("MOBILE","ACTIVITY","PK")           ;  Statut="MOBILE"
        case("BUFFER")                           ;  Statut="BUFFER"
        case("BALANCE")                          ;  Statut="BALANCE"
        !! case("STOIKIO")                       ;  Statut="STOIKIO"
        case default
          Ok= .false.
          Msg=      trim(WType)//" -> unrecognized species status !!"
          return !---------------------------------------------return
        end select
        !
        !-----------------------------------------------------for Oxygen
        if(N==iO_ .and. Statut/="INERT") then
          Ok= .false.
          Msg= "Oxygen should be INERT !!!"
          return
        end if
        !----------------------------------------------------/for Oxygen
        !
        !-- assign component mobility status
        Cpn%Statut= trim(Statut)
        !-----------------------------------------/read component status
        !
        !----------------------------------------------read species name
        if(EoL) then
          Ok= .false.
          Msg=                   "Species lacking for "//trim(Cpn%NamCp)
          return
        end if
        call LinToWrd(L,V1,EoL)!= species name
        !
        call Species_Rename(V1)
        !
        iSp=Species_Index(V1,vSpc)
        !
        if(iSp==0) then
          Ok= .false.
          Msg=              trim(V1)//" <-Species NOT in database ???!!"
          return
        else
          if(Statut=="INERT"  .and. vSpc(iSp)%Typ/="AQU") then
            Msg=  trim(V1)//" <- an INERT SPECIES should be AQUEOUS !!!"
            return
          end if
          !<BUFFER MODIFIED
          if(BufferIsExtern   .and. &
          &  Statut=="BUFFER" .and. &
          &  vSpc(iSp)%Typ=="AQU") then
            Msg= trim(V1)//" <-a BUFFER SPECIES should NOT be AQUEOUS !!!"
            return
          end if
          !</BUFFER MODIFIED
        end if
        !
        !------------------------------------------------- for Oxygen --
        if(N==iO_ .and. trim(V1)/="H2O") then
          Ok= .false.
          Msg=                   "for Oxygen, species should be H2O !!!"
          return
        end if
        !------------------------------------------------/ for Oxygen --
        !
        Cpn%iSpc=iSp != index of "associated" species
        !
        if(Cpn%Statut=="INERT") then
          if(vSpc(iSp)%vStoikio(iEl)==0) then
            Ok= .false.
            Msg= &
            & "Element "//trim(vEle(iEl)%NamEl)// &
            & " is not contained in species "//trim(vSpc(iSp)%NamSp)
            return
          end if
        end if
        !
        !<new, 200911>
        !-----------------write species'stoichio to component'stoichio--
        select case(trim(WType))
        case("MOLE","GRAM","PPM")
          tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
        end select
        !if(trim(WType)=="MOLE" .or. trim(WType)=="GRAM") then
        !  !Cpn%vStoikCp(0:nEl+1)= vSpc(iSp)%vStoikio(0:nEl+1)
        !  tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
        !end if
        !</new>
        !
        !---------------------------------------------/read species name
        !
        !------------------------- read numeric data (mole, pk, etc.) --
        X1=Zero
        Cpn%namSol= "Z" !trim(vSpc(iSp)%Typ) !default value
        !
        if(Statut/="BALANCE") then
          !
          if(EoL) then
          !----------------------if no numeric data, take default values
            !
            !------------ for mobile pure species, default is saturation
            if(Statut=="MOBILE" .or. Statut=="BUFFER") X1=Zero
            !
            if(Statut=="INERT") X1= CpnMolMinim
            !
          else
            !
            call LinToWrd(L,V2,EoL)
            !
            if(trim(V2)/="SOLUTION" .and. trim(V2)/="MIXTURE") then
              ! read real value -> total amount, colog10(activity), etc
              call WrdToReal(V2,X1)
            else
              !
              if(trim(V2)=="SOLUTION") &
              & call Warning_(trim(V2)//" soon Obsolete, better use MIXTURE !!!")
              !
              if(Cpn%Statut=="INERT") then
                Ok= .false.
                Msg= "SOLUTION/MIXTURE is NOT VALID KEYWORD for INERT components !!!"
                return
              end if
              !-------------------------read name of "controlling phase"
              call LinToWrd(L,V2,EoL)
              if(MixPhase_Index(vMixFas,V2)==0) then
                Ok= .false.
                Msg=          trim(V2)//" <-this phase is unknown !!!!!"
                return
              end if
              !
              Cpn%namSol=trim(V2)
              !!print *, trim(V2)  ;  pause
              !
            end if
            !
          end if
          !
        end if !!if(Statut/="BALANCE")
        !---------------------------------------------/read numeric data
        !
        !-------------------------------------------process numeric data
        select case(trim(Statut))
        !
        case("INERT") !---------------------------------------case INERT
          !
          Cpn%Mole=X1  != tot.amount of element
          !
          Cpn%Factor= 1.D0
          !---------------------------------------------unit conversions
          select case(trim(WType))
          case("GRAM")  ;  Cpn%Factor= vSpc(iSp)%WeitKg *1.0D3
          case("PPM")   ;  Cpn%Factor= vSpc(iSp)%WeitKg *1.0D6
          end select
          !if (trim(WType)=="GRAM") then
          !  Cpn%Factor= vSpc(iSp)%WeitKg *1.0D3
          !else
            !---------------------------------------------------OBSOLETE
            !-------------------------------maintained for compatibility
            if(.not. EoL) then
              !------------------------- read unit: MG/KG__, UG/KG__,...
              call LinToWrd(L,V2,EoL)
              select case(trim(V2))
              case("MG/KG","PPM")
                Cpn%Factor= vSpc(iSp)%WeitKg *1.0D6
              case("UG/KG","PPB")
                Cpn%Factor= vSpc(iSp)%WeitKg *1.0D9
              case default
              !------------------------------------unknown keyword, stop
                Ok= .false.
                Msg= trim(V2)//" = unknown unit..."
                return
              end select
            end if
            !--------------------------------------------------/OBSOLETE
          !end if
          !--------------------------------------------/unit conversions
          Cpn%Mole= Cpn%Mole /Cpn%Factor
          !
        case("MOBILE","BUFFER") !----------------------------case MOBILE
          !
          Cpn%LnAct= -X1 *Ln10 !from colog10 to ln
          !
          !---------------------------------------------unit conversions
          select case(trim(WType))
          !
          case("PK")
            Cpn%LnAct= -X1 *Ln10
          !
          case("ACTIVITY")
            Eps=EPSILON(X1)
            if(ABS(X1)<Eps) then
              Ok= .false.
              Msg= "Activity cannot be Zero !!!"
              return
            end if
            Cpn%LnAct= log(X1)
          !
          end select
          !--------------------------------------------/unit conversions
          !
          !<OBSOLETE>
          !!<keyword BUFFER is now placed differently>
          !!if(.not. EoL) then
          !!  call LinToWrd(L,V2,EoL)
          !!  if(trim(V2)=="BUFFER") Statut="BUFFER"
          !!  !!! BUFFER MODIFIED
          !!  if(BufferIsExtern .and. Statut=="BUFFER" .and. vSpc(iSp)%Typ=="AQU") &
          !!  & call Stop_(trim(V1)//" <-a BUFFER SPECIES should NOT be AQUEOUS !!!!!")
          !!  !!! BUFFER MODIFIED end
          !!end if
          !</OBSOLETE>
          !
        end select
        !--------------------------------------/ process numeric data --
        !
        if (iDebug>2) write(fTrc,'(3(A,A1),2(G15.6,A1))') &
        & Cpn%Statut,          T_, &
        & vEle(Cpn%iEle)%NamEl,T_, &
        & vSpc(Cpn%iSpc)%NamSp,T_, &
        & Cpn%Mole,            T_, &
        & Cpn%LnAct,           T_
        !
        vCpn(N)= Cpn
        !
      end do DoBlock
    end if !
  end do DoFile
  call closeFILE(f)
  !
  if (iDebug>2) write(fTrc,'(A)') &
  & "========================================================================="
  !
  if(N==0) return
  !
  !! if(Check_TP) then
  !!   if(DtbFormat=="LOGK") then
  !!     iTP_= TPcond_IndexT(TdgK,vTPCond)
  !!     if(iTP_==0) then
  !!       print '(A,G15.6)',"TdgC=", TdgK -T_CK
  !!       call Stop_("Temperature not found in TP-series")
  !!     else
  !!       TdgK= vTPCond(iTP_)%TdgC +T_CK
  !!       Pbar= vTPCond(iTP_)%Pbar
  !!     end if
  !!   else
  !!     iTP_= TPcond_IndexTP(TdgK,Pbar,vTPCond)
  !!     if(iDebug==4) call Pause_("TPcond_IndexTP= "//trim(FIntToStr3(iTP_)))
  !!     !
  !!     if(TdgK <Zero) TdgK= Zero
  !!     !
  !!     ! prevent pressure to be below vapor saturation for pure H2O
  !!     if (TdgK<=647.25D0) then; call fluid_h2o_psat(TdgK,PSatBar)
  !!     else                    ; PSatBar= 220.55D0
  !!     end if
  !!     if(Pbar<PSatBar) Pbar=PSatBar
  !!   end if
  !!   if(iDebug==4) print '(/,A,2G15.6,/)',"TdgC,Pbar=",TdgK-T_CK,Pbar
  !! end if
  !
  if(SysType=="AQUEOUS") then
    !---------------------------add O__/H2O if not found in species list
    if(iO_==0) then
      !
      N= N+1
      if(N>size(vCpn)) then
        Ok= .false.
        Msg= "Too many components"
        return
      end if
      !
      iO_= N
      !
      Cpn%iEle=   Element_Index("O__",vEle)
      Cpn%NamCp=  vEle(Cpn%iEle)%NamEl
      Cpn%iSpc=   Species_Index("H2O",vSpc)
      Cpn%Mole=   One /vSpc(iW)%WeitKg
      Cpn%LnAct=  Zero
      Cpn%Statut= "INERT"
      Cpn%namSol= "Z"
      !
      vCpn(N)= Cpn
      !Cpn%vStoikCp(0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      tStoikCp(N,Cpn%iEle)= 1
      !
      sListElem= trim(sListElem)//"O__"
      !
    end if
    !------------------------------------------------------/ add O__/H2O
    !
    !----------------------------add H__/H+ if not found in species list
    if(iH_==0) then !
      N=   N+1
      if(N>size(vCpn)) then
        Ok= .false.
        Msg= "Too many components"
        return
      end if
      !
      iH_= N
      !
      Cpn%iEle=   Element_Index("H__",vEle)
      Cpn%NamCp=  vEle(Cpn%iEle)%NamEl
      Cpn%iSpc=   Species_Index("H+", vSpc)
      !!Cpn%Mole=   2.0D0 /vSpc(iW)%WeitKg
      Cpn%LnAct=  Zero
      Cpn%Statut= "INERT"
      Cpn%namSol= "Z"
      !
      vCpn(N)= Cpn
      tStoikCp(N,Cpn%iEle)= 1
      !Cpn%vStoikCp(0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !
      sListElem=trim(sListElem)//"H__"
      !
    end if
    !-------------------------------------------------------/ add H__/H+
  end if !if(SysType=="AQUEOUS")
  !
  !-------------------------------------compute mole numbers of elements
  if(MoleIsTrue) then
    allocate(vX(N))
    vX(:)= Zero
    !
    do iEl=1,N
      do I=1,N
        if(vCpn(I)%Statut=="INERT") &
        & vX(iEl)= vX(iEl) +vCpn(I)%Mole *tStoikCp(I,vCpn(iEl)%iEle)
      end do
    end do
    !
    !--- trace
    if(iDebug==4) then
      print *,"== Components_Read =="
      do iEl=1,N
        print *,vCpn(iEl)%NamCp,vCpn(iEl)%iEle
        do I=1,N
          print *, &
          & "        ",      &
          & vCpn(I)%NamCp,    &
          & trim(vSpc(vCpn(I)%iSpc)%NamSp), &
          & tStoikCp(I,vCpn(iEl)%iEle)
          !!& vCpn(I)%vStoikCp(vCpn(iEl)%iEle)
        end do
        print '(A,G15.6)',"=======================================tot=",vX(iEl)
      end do
      call Pause_
    end if !--</ trace
    !
    vCpn(1:N)%Mole= vX(1:N)
    !
    deallocate(vX)
  end if
  !-----------------------------------/ compute mole numbers of elements
  !
  deallocate(tStoikCp)
  !
  if(iH_/=0) vCpn(iH_)%Mole= 2.0D0*vCpn(iO_)%Mole
  !
  !--------------------------------------------place Oxygen as species 1
  if(iO_/=0) then
    if(iO_/=1) then
      Cpn=       vCpn(1)
      vCpn(1)=   vCpn(iO_)
      vCpn(iO_)= Cpn
    end if
  end if
  !------------------------------------------/ place Oxygen as species 1
  !
  if (iDebug>0) then
    do I=1,N
      Cpn= vCpn(I)
      write(fTrc,'(3(A,A1),2(G15.6,A1))') &
      & Cpn%Statut,          T_, &
      & vEle(Cpn%iEle)%NamEl,T_, &
      & vSpc(Cpn%iSpc)%NamSp,T_, &
      & Cpn%Mole,            T_, &
      & Cpn%LnAct,           T_
    end do
  end if
  !
  if(iDebug>0) write(fTrc,'(A10,A)') "LISTELEM= ",trim(sListElem)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ ReadInput_System"
  !
end subroutine Components_Read

logical function IsSpecies(vSpc,W)
  use M_T_Species,only: T_Species,Species_Index,Species_Rename
  type(T_Species), intent(in):: vSpc(:)
  character(len=*),intent(in):: W
  !
  IsSpecies= (Species_Index(trim(W),vSpc)/=0)
  !
end function IsSpecies

logical function IsElement(vEle,W)
  use M_T_Element,only: T_Element,Element_Index
  use M_IoTools,  only: Str_Append
  type(T_Element), intent(in):: vEle(:)
  character(len=*),intent(in):: W
  !
  character(len=3):: W2
  !
  W2= W(1:3)
  call Str_Append(W2,3)
  IsElement= (Element_Index(trim(W2),vEle)>0)
  !
end function IsElement

logical function IsKeyword(sList,W)
  character(len=*),intent(in):: sList
  character(len=*),intent(in):: W
  !
  IsKeyword= (index(trim(sList),"_"//trim(W)//"_")>0)
  !
end function IsKeyword

end module M_Component_Read
