module M_Path_Read
!--
!-- read parameters for path calculations --
!--
!
  use M_Kinds
  use M_Trace
  use M_T_MixPhase,only: T_MixPhase
  !
  implicit none
  !
  private
  !
  public:: Path_ReadParam
  public:: Path_ReadParam_new
  public:: Path_ReadMode
  ! 
contains

subroutine Path_ReadMode( &
& NamFInn, &
& PathMode,Ok,Msg)
!--
!-- read path mode --
!--
  use M_IOTools
  !
  character(len=*),intent(in) :: NamFInn
  character(len=*),intent(out):: PathMode
  logical,         intent(out):: Ok
  character(len=*),intent(out):: Msg
  !
  character(len=255):: L
  character(len=80) :: W
  logical:: EoL
  integer:: f, ios
  !
  Ok=.false.
  Msg= ""
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L  ; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!')   cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    if(W=="PATH") then
      !
      if(EoL) then
        Ok=  .false.
        Msg= "Error in PATH: Path type undefined"
        return
      end if
      ! 
      !------------------------------------------- read the path mode --
      call LinToWrd(L,W,EoL)
      !
      Ok=.true.
      !
      select case(trim(W))
      
      case("MIX")         ; PathMode= "MIX"
      case("ADD")         ; PathMode= "ADD"
      case("ADDALL")      ; PathMode= "ADA"
      case("CHANGE")      ; PathMode= "CHG"
      case("GRID")        ; PathMode= "GRD"
      case("LOGK")        ; PathMode= "LGK"
      
      case default
        Ok=  .false.
        Msg= "Error in PATH: Path mode should be CHANGE | ADD | MIX | LOGK"
        return
      
      end select
      !------------------------------------------/ read the path mode --
      !
    end if !Cod=="PATH"
    !
  end do DoFile
  !
  close(f)
  !
  if(.not. Ok) Msg= "PATH block not found !!"
  !
  return
end subroutine Path_ReadMode

subroutine Path_ReadParam( &
& NamFInn,   &  !IN
& PathMode,  &  !IN
& vCpn,      &  !IN
& TdgK,Pbar, &  !IN
& Ok,Msg)       !OUT
!--
!-- read the path parameters --
!--
  use M_IOTools
  use M_Dtb_Const,   only: T_CK
  use M_System_Tools,only: System_TP_Update
  use M_T_Element,   only: Formula_Read,Element_Index
  use M_T_Species,   only: Species_Index
  use M_T_Component, only: T_Component,Component_Index
  use M_Basis,       only: Basis_Change
  use M_T_MixPhase,  only: MixPhase_Index
  !
  use M_Global_Vars, only: vEle,vSpc,vFas,vMixFas,vMixModel
  use M_Basis_Vars,  only: tAlfFs
  use M_Path_Vars,   only: vLPath,vTPpath,tPathData,DimPath
  use M_Path_Vars,   only: vPhasBegin,vPhasFinal
  use M_Path_Vars,   only: iLogK,vPathLogK
  !
  character(len=*), intent(in)   :: NamFInn
  character(len=*), intent(in)   :: PathMode
  type(T_Component),intent(in)   :: vCpn(:)
  real(dp),         intent(in)   :: TdgK,Pbar
  logical,          intent(out)  :: Ok
  character(len=*), intent(out)  :: Msg
  !
  integer,parameter:: DimMax= 255
  !
  character(len=1024):: L
  character(len=80) :: W,W1
  logical :: EoL, fOk, OkMix, PathAdd, BuildTable, TestMixture
  integer :: ios, ieOx, Z, Zxs, nDiv, iMix
  integer :: I, J, K, f, iCp, nCp, nStep
  integer :: fCheck
  real(dp):: rBegin,rFinal,rRatio,rDelta,AddedTot,X
  integer,dimension(1:size(vCpn))::vStoik !for PathAdd
  !
  real(dp):: vTmp(DimMax)
  real(dp),allocatable:: tTmp(:,:)
  integer, allocatable:: vNstep(:)
  integer, allocatable:: vStAdd(:)
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Path_ReadParam"
  !
  !if(iDebug==4) print '(A)',"Path_Read"
  !
  Ok=.false.
  Msg= ""
  TestMixture= .false.
  !
  nCp=size(vCpn)
  !
  allocate(tTmp(1:nCp+2,1:DimMax))
  ! 1:nCp = components
  ! nCp+1 = TdgC
  ! nCp+2 = Pbar
  tTmp(:,:)= Zero
  tTmp(nCp+1,:)= TdgK -T_CK
  tTmp(nCp+2,:)= Pbar
  !
  allocate(vNstep(nCp+2))  ;  vNstep(:)= 0
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  if(iDebug>2) then
    do I=1,size(vCpn)
      print '(I3,1X,3A)', &
      & I,vCpn(I)%NamCp," -> STATUT=",vCpn(I)%Statut
    end do
    call Pause_
  end if
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile !---------------------------------end of file
    !
    call LinToWrd(L,W,EoL)
    !
    if(W(1:1)=='!') cycle DoFile !--------------------skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    BuildTable= .false.
    !
    if(W=="PATH") then 
      !
      Ok=.true.
      !-----------------------------------------------allocate variables
      !---------------------------------------according to the Path mode
      select case(PathMode) 
      case("ADD","ADA")
        allocate(vStAdd(0:nCp))  ; vStAdd=0
        !
      case("CHG")
        allocate(vLPath(1:nCp+2))  ;  vLPath=.false.
        !
        allocate(vPhasBegin(1:nCp))  ;  vPhasBegin=0
        allocate(vPhasFinal(1:nCp))  ;  vPhasFinal=0
        !
      end select
      !---------------------------------------------/ allocate variables
      !
      !---------------------------------------- read the path parameters
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!')   cycle DoBlock !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        if(W=="ENDINPUT") exit  DoFile
        if(W=="END" .or. W=="ENDPATH") exit  DoFile
        !-> read only the first PATH block found !!
        !
        select case(trim(PathMode))
        !-------------------------------------------- case "LGK" (-LOGK)
        case("LGK")
          
          iLogK= Species_Index(trim(W),vSpc)
          if(iLogK==0) then
            Ok=  .false.
            Msg=      "In PATH LOGK, "//trim(W)//"= unknown species !!!"
            return !----------------------------------------------------
          end if
          
          !--------------------------------------------------read values
          rBegin=Zero; rFinal=Zero; rDelta=Zero
          do
            if(EoL) exit
            
            call LinToWrd(L,W1,EoL)
            
            select case(W1)
            case("INITIAL")
              call LinToWrd(L,W1,EoL); call WrdToReal(W1,rBegin)
            case("FINAL")
              call LinToWrd(L,W1,EoL); call WrdToReal(W1,rFinal)
            case("DELTA")
              call LinToWrd(L,W1,EoL); call WrdToReal(W1,rDelta)
            case default
              Ok=  .false.
              Msg= "In PATH LOGK, "//trim(W1)//"= invalid keyword !!!"
              return !--------------------------------------------return
            end select
          end do
          
          if(abs(rDelta)<Epsilon(rDelta)) then
            Ok=  .false.
            Msg=                       "In PATH LOGK, Delta not defined"
            return !----------------------------------------------return
          end if
          
          ! if(abs(rFinal-rBegin) < 0.1) rFinal= rBegin +One
          
          rDelta= SIGN(rDelta,rFinal-rBegin)
          !-- function SIGN(a,b) returns abs(a)*sign(b) --
          !
          !------------------------------------------------/ read values
          !
          !--------------------------------------------- build vPathLogK
          I= 1
          vTmp(I)= rBegin
          do
            if(I>DimMax) exit
            if(abs(vTmp(I)-rBegin)>abs(rFinal-rBegin)) exit
            I=I+1
            vTmp(I)= vTmp(I-1) + rDelta
            print '(A,I3,G15.6)',"logK=",I,vTmp(i)
          end do
          ! call pause_
          !
          DimPath= I
          allocate(vPathLogK(1:I))
          vPathLogK(1:I)= vTmp(1:I)
          !--------------------------------------------/ build vPathLogK
          !
        !-------------------------------------------/ case "LGK" (-LOGK)
        !
        !---------------------------------------------------- case "ADD"
        case("ADD","ADA")
          !
          !------------------------- read stoichio of added material --
          call Formula_Read(W,vEle,Z,nDiv,fOk,vStoik)
          !
          if(Z/=0) then
            Ok=  .false.
            Msg=     "In PATH ADD, "//trim(W)//"= should be neutral !!!"
            return !----------------------------------------------return
          end if
          !
          !---------------------compute stoikio coeff of redox component
          ieOx= Element_Index("OX_",vEle)
          Zxs= dot_product(vStoik(:),vEle(:)%Z)
          !-- -> oxidation state
          if(Zxs/=Z) then
            if(ieOx/=0) then
              vStoik(ieOx)= Z - Zxs
              fOk= .true.
            else
              fOk= .false.
            end if
          end if
          !-------------------/ compute stoikio coeff of redox component
          !
          if(.not.fOk) then
            Ok=  .false.
            Msg= "In PATH ADD, "//trim(W)//"= problem in stoikiometry ?"
            return !----------------------------------------------return
          end if
          !
          vStAdd(1:nCp)= vStoik(vCpn(1:nCp)%iEle)
          vStAdd(0)=     nDiv
          !--------------- vStAdd(1:nCp)/vStAdd(0) is the stoikio'vector
          !
          !----------------------------/ read stoichio of added material
          !
          rBegin=  1.D-6
          rFinal=  1.0D0
          rRatio=  1.5D0
          rDelta=  Zero
          !
          nStep= 1
          AddedTot= Zero
          tTmp(1:nCp,1)= vCpn(:)%Mole
          !
          !------------------------------------------------- read values
          do
            if(EoL) exit
            !
            call LinToWrd(L,W,EoL)
            select case(W)
              
              case("INITIAL","FINAL","RATIO","DELTA")
                BuildTable= .true.
                select case(W)
                  case("INITIAL"); call LinToWrd(L,W,EoL); call WrdToReal(W,rBegin)
                  case("FINAL");   call LinToWrd(L,W,EoL); call WrdToReal(W,rFinal)
                  case("RATIO");   call LinToWrd(L,W,EoL); call WrdToReal(W,rRatio)
                  case("DELTA");   call LinToWrd(L,W,EoL); call WrdToReal(W,rDelta)
                end select
              
              case default
                if(BuildTable) then
                  Ok=  .false.
                  Msg=       trim(W)//"= unknown keyword in PATH CHANGE"
                  return !----------------------------------------return
                else
                  do
                    nStep= nStep+1
                    call WrdToReal(W,X)
                    !
                    do iCp=1,size(vCpn)
                      tTmp(iCp,nStep)= tTmp(iCp,1) &
                      + vStAdd(iCp) /real(vStAdd(0)) *X
                    end do
                    !
                    if(EoL) exit
                    if(nStep==DimMax) exit
                    !
                    call LinToWrd(L,W,EoL)
                  end do
                  !vNstep(I)= nStep
                end if
                
            end select
          end do
          !------------------------------------------------/ read values
          !
          !--------------------------------------------- build tPathData
          PathAdd= (rDelta>Zero) !.false.
          !
          if(iDebug>2) then
            print '(A)',"conditions for PATH ADD"
            print '(A)',"stoikiometry vStAdd" ; call Pause_
            do I=1,nCp
              if(vStAdd(I)/=0) print '(I3,G9.2)',I,vStAdd(I)  
            end do
            !! print '(A,3G15.6,I3)',"rAddBegin,rAddFinal,rAddRatio,rAddDelta", &
            !! & rAddBegin,rAddFinal,rAddRatio,rAddDelta
            call Pause_
          end if
          !
          if(BuildTable) then
            !
            if(iDebug>2) then
              call GetUnit(fCheck)
              open(fCheck,file="debug_pathadd_table.restab")
            end if
            !
            do
              nStep= nStep+1
              !
              do iCp=1,size(vCpn)
                !
                if(vStAdd(iCp)/=0) then
                  !
                  if(PathAdd) then
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp)/real(vStAdd(0)) *rDelta !*(nStep-1)
                  else
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp)/real(vStAdd(0)) *rBegin *rRatio**(nStep-1)
                    !AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
                  end if
                  !
                else
                  tTmp(iCp,nStep)= tTmp(iCp,nStep-1)
                end if
                !
                if(iDebug>2) &
                & write(fCheck,'(G15.6,1X)',advance="no") tTmp(iCp,nStep)
                !
              end do
              !
              if(iDebug>2) write(fCheck,*)
              if(iDebug>2) flush(fCheck)
              !
              if(PathAdd) then
                AddedTot= AddedTot +rDelta
              else
                AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
              end if
              !
              if(AddedTot>=rFinal) exit
              !
              if(nStep>=DimMax) exit
              !
            end do
          end if
          !
          DimPath= nStep
          !
          allocate(tPathData(nCp,DimPath))
          tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
          !
          allocate(vTPpath(DimPath))
          vTPpath(1:DimPath)%TdgC= tTmp(nCp+1,1:DimPath)
          vTPpath(1:DimPath)%Pbar= tTmp(nCp+2,1:DimPath)
          !--------------------------------------------/ build tPathData
          !
        !_ADD
        !---------------------------------------------------/ case "ADD"
        !!! case("MIX")
        !!!   !
        !!!   call Str_Append(W,3)
        !!!   I=Element_Index(W,vEle)
        !!!   if(I<1) then
        !!!     call Stop_(trim(W)//"<-NOT in the current system !!!!")
        !!!   else 
        !!!     call LinToWrd(L,W,EoL)
        !!!     call WrdToReal(W,X1)
        !!!     vSys2(I)=X1
        !!!   end if
        !!! !_MIX
        !
        !---------------------------------------------------- case "CHG"
        case("CHG")
          !
          select case(trim(W))
          !
          case("TDGC")  ;  I= nCp+1
          !
          case("PBAR")  ;  I= nCp+2
          !
          case default
            call Str_Append(W,3)
            !-- from component name (max 3 chars), find index in vCpn
            I= Component_Index(W,vCpn)
            if(I<=0) then
              Ok=  .false.
              Msg= "PATH CHANGE: Component "//trim(W)//" Not in the system .."
              return !--------------------------------------------------
            end if
          !
          end select
          !
          vLPath(I)=.true.
          !
          if(I<=nCp) then
            J=    vCpn(I)%iSpc !-> find related species
            iMix= vCpn(I)%iMix !-> component is an end-member of mixture phase iMix
          else
            J= 0
            iMix= 0
          end if
          !
          !--------------- prim' mobile species is aqueous or pure phase
          if(iMix==0) then
            !
            rBegin= Zero; rFinal= Zero
            rRatio= One;  rDelta= Zero
            !----------------------------------------------- read values
            do
              !
              if(EoL) exit
              call LinToWrd(L,W,EoL)
              !
              select case(W)
              !
              case("INITIAL","FINAL","RATIO","DELTA")
                BuildTable= .true.
                select case(W)
                case("INITIAL")
                  call LinToWrd(L,W,EoL)  ;  call WrdToReal(W,rBegin)
                case("FINAL")
                  call LinToWrd(L,W,EoL)  ;  call WrdToReal(W,rFinal)
                case("RATIO")
                  call LinToWrd(L,W,EoL)  ;  call WrdToReal(W,rRatio)
                case("DELTA")
                  call LinToWrd(L,W,EoL)  ;  call WrdToReal(W,rDelta)
                end select
              !
              case default
                if(BuildTable) then
                  Ok=  .false.
                  Msg=     trim(W)//"= unknown keyword in PATH CHANGE"
                  return !----------------------------------------------
                else
                  nStep= 0
                  do
                    nStep= nStep+1
                    call WrdToReal(W,tTmp(I,nStep))
                    if(EoL) exit
                    if(nStep==DimMax) exit
                    call LinToWrd(L,W,EoL)
                  end do
                  vNstep(I)= nStep
                end if
                  
              end select
              !
            end do
            !----------------------------------------------/ read values
            PathAdd= (rDelta>Zero)
            !
            !------------------------------------------------ BuildTable
            if(BuildTable) then
              if(I<=nCp) then
                !---------------------------- case of chemical change --
                if(vCpn(I)%Statut=="INERT" .and.(.not. PathAdd)) then
                  !->  default case (geometric)
                  if(rRatio < 1.10D0 ) rRatio= 1.10D0
                  if(rRatio > 3.0D0  ) rRatio= 3.00D0
                  if(rFinal < rBegin ) rRatio= One/rRatio
                else !MOBILE or BUFFER or DELTA
                  if(rDelta==Zero)  rDelta= abs(rFinal-rBegin) /20.0D0
                  !-> default value, in case not in input
                  if(rFinal>rBegin) rDelta= abs(rDelta)
                  if(rFinal<rBegin) rDelta=-abs(rDelta)
                end if
                !---/
              else
                !--------------------------------- case of (T,P) changes
                if(.not. PathAdd) then
                  Ok=  .false.
                  Msg=           "PATH CHANGE: for (T,P) only DELTA ..."
                  return !----------------------------------------------
                end if
                if(rDelta==Zero)  rDelta= abs(rFinal-rBegin) /20.0D0
                !-> default value, in case not in input
                if(rFinal>rBegin) rDelta= abs(rDelta)
                if(rFinal<rBegin) rDelta=-abs(rDelta)
                !---/
              end if
              !
              if(idebug>1) write(fTrc,'(I3,A20,4G12.3)') &
              & I," =Begin,Final,Step= ",rBegin,rFinal,rDelta,rRatio
              !
              nStep= 0
              do
                !
                nStep= nStep+1
                if(PathAdd) then
                  tTmp(I,nStep)= rBegin + (nStep-1) *rDelta
                else
                  tTmp(I,nStep)= rBegin *rRatio**(nStep-1)
                end if
                !
                if(rFinal>rBegin) then
                  if(tTmp(I,nStep) > rFinal) exit
                else
                  if(tTmp(I,nStep) < rFinal) exit
                end if
                !
                if(nStep==DimMax) exit
                !
              end do
              vNstep(I)= nStep-1
              !
            end if
            !-----------------------------------------------/ BuildTable
            !
            ! print *,"unit change"
            if(I<=nCp) then
              if(vCpn(I)%Statut=="INERT") &
              & tTmp(I,1:vNstep(I))= tTmp(I,1:vNstep(I)) /vCpn(I)%Factor
            end if
            !
            !----------------------------- add amounts in SYSTEM.ROCK --
            if(I<=nCp) then
              if(vCpn(I)%Statut=="INERT") then
                do J=1,vNstep(I)
                  !! print *,"SIZE(tTmp)=",size(tTmp,1),size(tTmp,2)
                  !! print *,"SIZE(vFas)=",size(vFas)  ;  pause
                  do K=1,size(vFas)
                    tTmp(I,J)= tTmp(I,J) + &
                    & tAlfFs(I,K) &
                    & *vFas(K)%MolFs
                  end do
                end do
              end if
            end if
            !--------------------------------/add amounts in SYSTEM.ROCK
            !
            if(idebug>1) then
              write(fTrc,'(I3,A1)',advance="no") I,t_
              do J= 1,vNstep(I)
                write(fTrc,'(G12.3,A1)',advance="no") tTmp(I,J),t_
              end do
              write(fTrc,*)
            end if
            !
          else
          !------ iMix/-0 -> species activity controlled by non-aqueous phase --
            TestMixture= .true.
            !----------------------------------------------- read values
            do
              !
              if(EoL) exit
              call LinToWrd(L,W,EoL)
              !
              select case(W)
                !
              case("INITIAL","FINAL")
                !
                call LinToWrd(L,W1,EoL) !-> name of initial phase
                K= MixPhase_Index(vMixFas,W1)
                !
                if(K==0) then
                  Ok=  .false.
                  Msg= trim(W1)//" = UNKNOWN PHASE in PATH CHANGE!!!"
                  return
                end if
                !
                if(vMixFas(K)%iModel /= vMixFas(iMix)%iModel)  then
                  Ok=  .false.
                  Msg=                                  "Phase "//trim(W1)// &
                  &   " should be of same model as "//trim(vMixFas(iMix)%Name)
                  return !----------------------------------------return
                end if
                !
                select case(W)
                case("INITIAL")
                  vPhasBegin(I)=K
                  
                  if(idebug>1) write(fTrc,'(4A)') &
                  & "PhasBegin=",vMixFas(K)%Name, &
                  & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
                  
                case("FINAL")
                  vPhasFinal(I)=K
                  
                  if(idebug>1) write(fTrc,'(4A)') &
                  & "PhasFinal=",vMixFas(K)%Name, &
                  & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
                  
                end select
                !
                !! case("RATIO")
                !!   Ok=  .false.
                !!   Msg=                  "only DELTA IN THIS case !!!"
                !!   return !----------------------------------------------
                !!   !call LinToWrd(L,W,EoL); call WrdToReal(W,rRatio)
                !! !
                !! case("DELTA")
                !!   call LinToWrd(L,W,EoL); call WrdToReal(W,rDelta)
                !!   !currently not used ...!!! (always on 100 steps)
                !
                case default
                  Ok=  .false.
                  Msg=           trim(W)//"= unknown keyword in PATH CHANGE !!!"
                  return !----------------------------------------return
                !  
              end select
              !
            end do
            !
            if(iDebug>1) then
              print '(4A)', &
              & "PhasBegin=",vMixFas(vPhasBegin(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
              print '(4A)', &
              & "PhasFinal=",vMixFas(vPhasFinal(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
            end if
            !----------------------------------------------/ read values
            !
            if(vPhasBegin(I)*vPhasFinal(I)==0) then
              Ok= .false.
              Msg=    "Problem in reading PATH conditions on solid solutions ??"
              return !--------------------------------------------return
            end if
            !
          end if
          !-----/ iMix/-0 -> species activity controlled by non-aqueous phase --
        ! end case CHANGE
        !---------------------------------------------------/ case "CHG"
        !
        end select
      end do DoBlock
      !---------------------------------------/ read the path parameters
      !! if(idebug>1) call Pause_
    end if !Cod=="PATH"
    !
  end do DoFile
  close(f)
  !
  !--------------------------------------------------- !! DEVELOPMENT !!
  if(TestMixture) then
    DimPath= 99
    allocate(vTPpath(DimPath))
    vTPpath(:)%TdgC= TdgK -T_CK
    vTPpath(:)%Pbar= Pbar
  end if
  !---------------------------------------------------/!! DEVELOPMENT !!
  !
  if(count(vNstep>0)>0) then
    !
    DimPath= MINVAL(vNstep(:),mask=vNstep(:)>0)
    allocate(tPathData(nCp,DimPath))
    tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
    !! print *,"minval(vNstep,mask=vNstep>0)", minval(vNstep,mask=vNstep>0)
    !! pause
    !
    allocate(vTPpath(DimPath))
    vTPpath(:)%TdgC= tTmp(nCp+1,1:DimPath)
    vTPpath(:)%Pbar= tTmp(nCp+2,1:DimPath)
    !
  end if
  !
  deallocate(tTmp)
  deallocate(vNstep)
  !
  if(allocated(vStAdd))  deallocate(vStAdd)
  !
  if(iDebug>1) call Pause_
  if(idebug>1) write(fTrc,'(A,/)') "</ Path_ReadParam"
  if(idebug>1) flush(fTrc)
  
  ! do i= 1, size(tPathData,2)
  !   do j=1,nCp
  !     write(11,'(G11.3,1X)',advance="no") tPathData(j,i)
  !   end do
  !   write(11,*)
  ! end do

end subroutine Path_ReadParam

subroutine Path_ReadParam_new( &
& NamFInn,   &  !IN
& PathMode,  &  !IN
& vCpn,      &  !IN
& TdgK,Pbar, &  !IN
& Ok,Msg)       !OUT
!--
!-- read the path parameters --
!-- -> initialize variables of M_Path_Vars
!--   vLPath,vTPpath
!--   vPhasBegin,vPhasFinal
!--   iLogK,vPathLogK
!--   tPathData,DimPath
!--
  use M_Dtb_Const,   only: T_CK
  use M_IOTools
  use M_System_Tools,only: System_TP_Update
  use M_T_Element,   only: Formula_Read,Element_Index
  use M_T_Species,   only: Species_Index
  use M_T_Component, only: T_Component,Component_Index
  use M_Basis,       only: Basis_Change
  use M_T_MixPhase,  only: MixPhase_Index
  !
  !--- vars
  use M_Global_Vars, only: vEle,vSpc,vFas,vMixFas,vMixModel
  !
  use M_Path_Vars,   only: vLPath,vTPpath
  use M_Path_Vars,   only: vPhasBegin,vPhasFinal
  use M_Path_Vars,   only: tPathData,DimPath
  use M_Path_Vars,   only: iLogK,vPathLogK
  !---/vars

  character(len=*), intent(in)   :: NamFInn
  character(len=*), intent(in)   :: PathMode
  type(T_Component),intent(in)   :: vCpn(:)
  real(dp),         intent(in)   :: TdgK,Pbar
  logical,          intent(out)  :: Ok
  character(len=*), intent(out)  :: Msg
  !
  integer,parameter:: DimMax= 255
  !
  character(len=1024):: L
  character(len=80) :: W,W1
  !
  logical :: EoL, fOk, PathAdd, BuildTable, TestMixture
  integer :: ios, ieOx, Z, Zxs, nDiv, iMix
  integer :: I, J, K, f, iCp, nCp, nStep
  integer :: fCheck
  real(dp):: rBegin,rFinal,rRatio,rDelta,AddedTot,X
  !
  integer,dimension(1:size(vCpn))::vStoik !for PathAdd
  ! real(dp):: TdgKAdd,PbarAdd
  !
  ! type(T_Component),allocatable:: vCpnAdd(:)
  !
  real(dp):: vTmp(1:DimMax)
  real(dp):: tTmp(1:size(vCpn)+2,1:DimMax)
  
  integer, allocatable:: vNstep(:)
  real(dp),allocatable:: vStAdd(:)
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Path_Read"
  !
  !if(iDebug==4) print '(A)',"Path_Read"
  !
  Ok=.false.
  Msg= ""
  TestMixture= .false.
  !
  nCp=size(vCpn)
  !
  ! 1:nCp = components
  ! nCp+1 = TdgC
  ! nCp+2 = Pbar
  tTmp(:,:)= Zero
  tTmp(nCp+1,:)= TdgK -T_CK
  tTmp(nCp+2,:)= Pbar
  !
  allocate(vNstep(nCp+2))  ;  vNstep(:)= 0
  !
  call GetUnit(f)
  open(f,file=trim(NamFInn))
  !
  if(iDebug>2) then
    do I=1,size(vCpn)
      !print '(I3,1X,3A)',I,vEle(vCpn(I)%iEle)%NamEl," -> STATUT=",vCpn(I)%Statut
      print '(I3,1X,3A)', &
      & I,vCpn(I)%NamCp," -> STATUT=",vCpn(I)%Statut
    end do
    call Pause_
  end if
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile !---------------------------------end of file
    !
    call LinToWrd(L,W,EoL)
    !
    if(W(1:1)=='!') cycle DoFile !--------------------skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    BuildTable= .false.
    !
    if(W=="PATH") then
      !
      Ok=.true.
      !
      !------------------- allocate variables according to the Path mode
      select case(PathMode)
      
      case("ADD","ADA")
        allocate(vStAdd(1:nCp))  ; vStAdd= 0.D0
        !
      case("CHG")
        allocate(vLPath(1:nCp+2))  ;  vLPath=.false.
        !
        allocate(vPhasBegin(1:nCp))  ;  vPhasBegin=0
        allocate(vPhasFinal(1:nCp))  ;  vPhasFinal=0
        !
      end select
      !-------------------------------------------------------/ allocate
      !
      if(iDebug>2) then
        print *,'Path_ReadParam'
        do i=1,size(vCpn)
          print *,i,vCpn(i)%Mole
        end do
        call pause_
      end if

      !! do I=1,size(vCpn)
        !! if(vCpn(I)%Statut=="INERT") tTmp(I,:)= vCpn(:)%Mole
        !! !if() tTmp(I,:)= vCpn(:)%Mole
      !! end do
      !
      !---------------------------------------- read the path parameters
      DoBlock: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!')   cycle DoBlock !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        if(W=="ENDINPUT") exit  DoFile
        if(W=="END" .or. W=="ENDPATH") exit  DoFile
        !-> read only the first PATH block found !!
        !
        select case(trim(PathMode))
        !-------------------------------------------- case "LGK" (-LOGK)
        case("LGK")
          !
          iLogK= Species_Index(trim(W),vSpc)
          if(iLogK==0) then
            Ok=  .false.
            Msg=      "In PATH LOGK, "//trim(W)//"= unknown species !!!"
            if(idebug>1) write(fTrc,'(A)') Msg
            return !----------------------------------------------return
          end if
          !------------------------------------------------- read values
          rBegin=Zero; rFinal=Zero; rDelta=Zero
          do
            if(EoL) exit
            !
            call LinToWrd(L,W1,EoL)
            select case(W1)
              !
              case("INITIAL")
                call LinToWrd(L,W1,EoL); call WrdToReal(W1,rBegin)
              case("FINAL")
                call LinToWrd(L,W1,EoL); call WrdToReal(W1,rFinal)
              case("DELTA")
                call LinToWrd(L,W1,EoL); call WrdToReal(W1,rDelta)
              case default
                Ok=  .false.
                Msg= "In PATH LOGK, "//trim(W1)//"= invalid keyword !!!"
                if(idebug>1) write(fTrc,'(A)') Msg
                return !------------------------------------------return
              !
            end select
            !
          end do
          !
          if(abs(rDelta)<Epsilon(rDelta)) then
            Ok=  .false.
            Msg=                       "In PATH LOGK, Delta not defined"
            if(idebug>1) write(fTrc,'(A)') Msg
            return !----------------------------------------------return
          end if
          if(abs(rFinal-rBegin) < 0.1) rFinal= rBegin +One
          rDelta= SIGN(rDelta,rFinal-rBegin)
          !-- function SIGN(a,b) returns abs(a)*sign(b) --
          !
          !------------------------------------------------/ read values
          !
          !--------------------------------------------- build vPathLogK
          I= 1
          vTmp(I)= rBegin
          do
            if(I>DimMax) exit
            if(abs(vTmp(I)-rBegin)>abs(rFinal-rBegin)) exit
            I=I+1
            vTmp(I)= vTmp(I-1) + rDelta
            !!print '(I3,G15.6)',I,vTmp(i)
          end do
          !!pause
          !
          DimPath= I
          allocate(vPathLogK(1:I))
          vPathLogK(1:I)= vTmp(1:I)
          !--------------------------------------------/ build vPathLogK
          !
        !-------------------------------------------/ case "LGK" (-LOGK)
        !
        !---------------------------------------------------- case "ADD"
        case("ADD","ADA")
          !
          ! if(trim(W)=="MIX") then
          !
          !   allocate(vCpnAdd(1:size(vCpn)))
          !   vCpnAdd= vCpn
          !
          !   call Path_Calc_Fluid( &
          !   & "SYSTEM.MIX",vCpnAdd, &
          !   & TdgKAdd,PbarAdd,Ok)
          !
          !   if(.not. Ok) then
          !     !! call Warning_("SYSTEM.MIX NOT FOUND")
          !     deallocate(vCpnAdd)
          !     Msg= "SYSTEM.MIX NOT FOUND"
          !     return
          !   end if
          !
          !   vStAdd(1:nCp)= vCpnAdd(1:nCp)%Mole
          !
          !   deallocate(vCpnAdd)
          !
          ! else

            !--------------------------- read stoichio of added material
            call Formula_Read(W,vEle,Z,nDiv,fOk,vStoik)
            !
            if(Z/=0) then
              Ok=  .false.
              Msg=     "In PATH ADD, "//trim(W)//"= should be neutral !!!"
              if(idebug>1) write(fTrc,'(A)') Msg
              return !--------------------------------------------return
            end if
            !
            !------------------ compute stoikio coeff of redox component
            ieOx= Element_Index("OX_",vEle)
            Zxs= dot_product(vStoik(:),vEle(:)%Z)
            !-- -> oxidation state
            if(Zxs/=Z) then
              if(ieOx/=0) then
                vStoik(ieOx)= Z - Zxs
                fOk= .true.
              else
                fOk= .false.
              end if
            end if
            !-----------------/ compute stoikio coeff of redox component
            !
            if(.not.fOk) then
              Ok=  .false.
              Msg= "In PATH ADD, "//trim(W)//"= problem in stoikiometry ?"
              if(idebug>1) write(fTrc,'(A)') Msg
              return !--------------------------------------------return
            end if
            !
            vStAdd(1:nCp)= vStoik(vCpn(1:nCp)%iEle) /real(nDiv)
            !------------- vStAdd(1:nCp)/vStAdd(0) is the stoikio'vector
          ! end if
          !----------------------------/ read stoichio of added material
          !
          rBegin=  1.D-6
          rFinal=  1.0D0
          rRatio=  1.5D0
          rDelta=  Zero
          !
          nStep= 1
          AddedTot= Zero
          tTmp(1:nCp,1)= vCpn(:)%Mole
          !
          !------------------------------------------------- read values
          do
            if(EoL) exit
            !
            call LinToWrd(L,W,EoL)
            select case(W)

            case("INITIAL","FINAL","RATIO","DELTA")
              BuildTable= .true.
              select case(W)
                case("INITIAL"); call LinToWrd(L,W,EoL); call WrdToReal(W,rBegin)
                case("FINAL");   call LinToWrd(L,W,EoL); call WrdToReal(W,rFinal)
                case("RATIO");   call LinToWrd(L,W,EoL); call WrdToReal(W,rRatio)
                case("DELTA");   call LinToWrd(L,W,EoL); call WrdToReal(W,rDelta)
              end select

            case default
              if(BuildTable) then
                Ok=  .false.
                Msg=          trim(W)//"= unknown keyword in PATH ADD"
                if(idebug>1) write(fTrc,'(A)') Msg
                return !------------------------------------------return
              else
                do
                  nStep= nStep+1
                  call WrdToReal(W,X)
                  !
                  do iCp=1,size(vCpn)
                    tTmp(iCp,nStep)= tTmp(iCp,1) + vStAdd(iCp) *X
                  end do
                  !
                  if(EoL) exit
                  if(nStep==DimMax) exit
                  call LinToWrd(L,W,EoL)
                end do
                !vNstep(I)= nStep
              end if

            end select
          end do
          !------------------------------------------------/ read values
          !
          !--------------------------------------------- build tPathData
          PathAdd= (rDelta>Zero) !.false.
          !
          !------------------------------------------------------- debug
          if(iDebug>2) then
            print '(A)',"conditions for PATH ADD"
            print '(A)',"stoikiometry vStAdd" ; call Pause_
            do I=1,nCp
              if(vStAdd(I)/=0.D0) print '(I3,G15.6)',I,vStAdd(I)
            end do
            !! print '(A,3G15.6,I3)',"rAddBegin,rAddFinal,rAddRatio,rAddDelta", &
            !! & rAddBegin,rAddFinal,rAddRatio,rAddDelta
            call Pause_
          end if
          !-------------------------------------------------------/debug
          !
          if(BuildTable) then
            !
            if(iDebug>2) then
              call GetUnit(fCheck)
              open(fCheck,file="debug_pathadd_table.restab")
            end if
            !
            do
              nStep= nStep+1
              !
              do iCp=1,size(vCpn)
                !
                if(vStAdd(iCp)/=0.D0) then
                  !
                  if(PathAdd) then
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp) *rDelta !*(nStep-1)
                  else
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp) *rBegin *rRatio**(nStep-1)
                    !AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
                  end if
                  !
                else
                  tTmp(iCp,nStep)= tTmp(iCp,nStep-1)
                end if
                !
                if(iDebug>2) &
                & write(fCheck,'(G15.6,1X)',advance="no") tTmp(iCp,nStep)
                !
              end do
              !
              if(iDebug>2) write(fCheck,*)
              if(iDebug>2) flush(fCheck)
              !
              if(PathAdd) then ;  AddedTot= AddedTot +rDelta
              else             ;  AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
              end if
              !
              if(AddedTot>=rFinal) exit
              !
              if(nStep>=DimMax) exit
              !
            end do
            !
          end if !if(BuildTable)
          !
          DimPath= nStep
          !
          if(allocated(tPathData)) deallocate(tPathData)
          allocate(tPathData(nCp,DimPath))
          tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
          !
          if(allocated(vTPpath)) deallocate(vTPpath)
          allocate(vTPpath(DimPath))
          vTPpath(:)%TdgC= tTmp(nCp+1,:)
          vTPpath(:)%Pbar= tTmp(nCp+2,:)
          !--------------------------------------------/ build tPathData
          !
        !_ADD
        !---------------------------------------------------/ case "ADD"
        !!! case("MIX")
        !!!   !
        !!!   call Str_Append(W,3)
        !!!   I=Element_Index(W,vEle)
        !!!   if(I<1) then
        !!!     call Stop_(trim(W)//"<-NOT in the current system !!!!")
        !!!   else
        !!!     call LinToWrd(L,W,EoL)
        !!!     call WrdToReal(W,X1)
        !!!     vSys2(I)=X1
        !!!   end if
        !!! !_MIX
        !
        !---------------------------------------------------- case "CHG"
        case("CHG")
          !
          select case(trim(W))
          !
          case("TDGC")  ;  I= nCp+1
          !
          case("PBAR")  ;  I= nCp+2
          !
          case default
            call Str_Append(W,3)
            !-- from component name (max 3 chars), find index in vCpn
            I= Component_Index(W,vCpn)
            if(I<=0) then
              Ok=  .false.
              Msg= "PATH CHANGE: Component "//trim(W)//" Not in the system .."
              if(idebug>1) write(fTrc,'(A)') Msg
              return !--------------------------------------------return
            end if
          !
          end select
          !
          vLPath(I)=.true.
          !
          if(I<=nCp) then
            J=    vCpn(I)%iSpc !-> find related species
            iMix= vCpn(I)%iMix !-> component is an end-member of mixture phase iMix
          else
            J= 0
            iMix= 0
          end if
          !
          ! if(iDebug>2) then
          !   if(I<=nCp) print *,vCpn(I)%NamCp
          !   print *,'Path_ReadParam, iMix=',iMix
          ! end if
          !--------------- prim' mobile species is aqueous or pure phase
          if(iMix==0) then
          
            rBegin= Zero; rFinal= Zero
            rRatio= One;  rDelta= Zero
            
            !----------------------------------------------- read values
            do
              if(EoL) exit
              call LinToWrd(L,W,EoL)
              select case(W)

                case("INITIAL","FINAL","RATIO","DELTA")
                  BuildTable= .true.
                  select case(W)
                    case("INITIAL"); call LinToWrd(L,W,EoL); call WrdToReal(W,rBegin)
                    case("FINAL")  ; call LinToWrd(L,W,EoL); call WrdToReal(W,rFinal)
                    case("RATIO")  ; call LinToWrd(L,W,EoL); call WrdToReal(W,rRatio)
                    case("DELTA")  ; call LinToWrd(L,W,EoL); call WrdToReal(W,rDelta)
                  end select

                case default
                  if(BuildTable) then
                    Ok=  .false.
                    Msg=     trim(W)//"= unknown keyword in PATH CHANGE"
                    if(idebug>1) write(fTrc,'(A)') Msg
                    return !--==========================================
                  else
                    nStep= 0
                    do
                      nStep= nStep+1
                      call WrdToReal(W,tTmp(I,nStep))
                      if(EoL) exit
                      if(nStep==DimMax) exit
                      call LinToWrd(L,W,EoL)
                    end do
                    vNstep(I)= nStep
                  end if

              end select
            end do
            !----------------------------------------------/ read values
            !
            PathAdd= (rDelta>Zero)
            !
            !------------------------------------------------ BuildTable
            if(BuildTable) then
              if(I<=nCp) then
                !------------------------------- case of chemical change
                if(vCpn(I)%Statut=="INERT" .and.(.not. PathAdd)) then
                  !->  default case (geometric)
                  if(rRatio < 1.10D0 ) rRatio= 1.10D0
                  if(rRatio > 3.0D0  ) rRatio= 3.00D0
                  if(rFinal < rBegin ) rRatio= One/rRatio
                else !MOBILE or BUFFER or DELTA
                  if(rDelta==Zero)  rDelta= abs(rFinal-rBegin) /20.0D0
                  !-> default value, in case not in input
                  if(rFinal>rBegin) rDelta= abs(rDelta)
                  if(rFinal<rBegin) rDelta=-abs(rDelta)
                end if
                !------------------------------------------------------/
              else
                !--------------------------------- case of (T,P) changes
                if(.not. PathAdd) then
                  Ok=  .false.
                  Msg=           "PATH CHANGE: for (T,P) only DELTA ..."
                  return !--------------------------------------return--
                end if
                if(rDelta==Zero)  rDelta= abs(rFinal-rBegin) /20.0D0
                !-> default value, in case not in input
                if(rFinal>rBegin) rDelta= abs(rDelta)
                if(rFinal<rBegin) rDelta=-abs(rDelta)
                !------------------------------------------------------/
              end if
              !
              if(idebug>1) write(fTrc,'(I3,A20,4G12.3)') &
              & I," =Begin,Final,Step= ",rBegin,rFinal,rDelta,rRatio
              !
              nStep= 0
              do
                !
                nStep= nStep+1
                if(PathAdd) then
                  tTmp(I,nStep)= rBegin + (nStep-1) *rDelta
                else
                  tTmp(I,nStep)= rBegin *rRatio**(nStep-1)
                end if
                !
                if(rFinal>rBegin) then
                  if(tTmp(I,nStep) > rFinal) exit
                else
                  if(tTmp(I,nStep) < rFinal) exit
                end if
                !
                if(nStep==DimMax) exit
                !
              end do
              vNstep(I)= nStep-1
            end if
            !-----------------------------------------------/ BuildTable
            !
            if(I<=nCp) then
              !if(vCpn(I)%Statut=="INERT") then
                do J=1,vNstep(I)
                  tTmp(I,J)= tTmp(I,J) /vCpn(I)%Factor
                  !if (iDebug>2) then
                  !  print *,"unit change"
                  !  print *,tTmp(I,J)
                  !  call pause_
                  !end if
                end do
              !end if
            end if
            !
            !----------------------------- add amounts in SYSTEM.ROCK --
            ! if(I<=nCp .and. vCpn(I)%Statut=="INERT") then
            !   do J=1,vNstep(I)
            !     print *,"SIZE(tTmp)=",size(tTmp,1),size(tTmp,2)
            !     print *,"SIZE(vFas)=",size(vFas)  ;  pause
            !     do K=1,size(vFas)
            !       tTmp(I,J)= tTmp(I,J) + &
            !       & tAlfFs(I,K) &
            !       & *vFas(K)%Mole
            !     end do
            !   end do
            ! end if
            !-----------------------------/add amounts in SYSTEM.ROCK --
            !
            if(idebug>1) then
              write(fTrc,'(I3,A1)',advance="no") I,t_
              do J= 1,vNstep(I)
                write(fTrc,'(G12.3,A1)',advance="no") tTmp(I,J),t_
              end do
              write(fTrc,*)
            end if
            !
          else
          !- iMix/-0 -> species activity controlled by non-aqueous phase
            TestMixture= .true.
            !----------------------------------------------- read values
            do
              !
              if(EoL) exit
              call LinToWrd(L,W,EoL)
              !
              select case(W)
                !
                case("INITIAL","FINAL")
                  !
                  call LinToWrd(L,W1,EoL) !-> name of initial phase
                  K= MixPhase_Index(vMixFas,W1)
                  !
                  if(K==0) then
                    Ok=  .false.
                    Msg= trim(W1)//" = UNKNOWN PHASE in PATH CHANGE!!!"
                    return
                  end if
                  !
                  if(vMixFas(K)%iModel /= vMixFas(iMix)%iModel)  then
                    Ok=  .false.
                    Msg=                                  "Phase "//trim(W1)// &
                    &   " should be of same model as "//trim(vMixFas(iMix)%Name)
                    return !--------------------------------------return
                  end if
                  !
                select case(W)
                  case("INITIAL")
                    vPhasBegin(I)=K

                    if(idebug>1) write(fTrc,'(4A)') &
                    & "PhasBegin=",vMixFas(K)%Name, &
                    & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name

                  case("FINAL")
                    vPhasFinal(I)=K

                    if(idebug>1) write(fTrc,'(4A)') &
                    & "PhasFinal=",vMixFas(K)%Name, &
                    & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name

                end select
                !
                ! case("RATIO")
                !   Ok=  .false.
                !   Msg=                             "only DELTA IN THIS case !!!"
                !   return !----------------------------------------------
                !   !call LinToWrd(L,W,EoL); call WrdToReal(W,rRatio)
                ! !
                ! case("DELTA")
                !   call LinToWrd(L,W,EoL); call WrdToReal(W,rDelta)
                !   !currently not used ...!!! (always on 100 steps)
                ! !
                ! case default
                !   Ok=  .false.
                !   Msg=           trim(W)//"= unknown keyword in PATH CHANGE !!!"
                !   return !----------------------------------------------
                !
              end select
            end do
            !
            if(iDebug>1) then
              print '(4A)', &
              & "PhasBegin=",vMixFas(vPhasBegin(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
              print '(4A)', &
              & "PhasFinal=",vMixFas(vPhasFinal(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
            end if
            !------------------------___-------------------/ read values
            !
            if(vPhasBegin(I)*vPhasFinal(I)==0) then
              Ok= .false.
              Msg= "Problem in reading PATH conditions on mixtures ??"
              return !--------------------------------------------return
            end if
            !
          end if
          !/ iMix/-0 -> species activity controlled by non-aqueous phase
        !_CHANGE
        !---------------------------------------------------/ case "CHG"
        !
        end select
      end do DoBlock
      !---------------------------------------/ read the path parameters
      !! if(idebug>1) call Pause_
    end if !Cod=="PATH"
    !
  end do DoFile
  close(f)
  !
  !--------------------------------------------------- !! DEVELOPMENT !!
  if(TestMixture) then
    DimPath= 99
    allocate(vTPpath(DimPath))
    vTPpath(:)%TdgC= TdgK -T_CK
    vTPpath(:)%Pbar= Pbar
  end if
  !---------------------------------------------------/!! DEVELOPMENT !!
  !
  if(count(vNstep>0)>0) then
    !
    DimPath= MINVAL(vNstep,mask=vNstep>0)
    allocate(tPathData(nCp,DimPath))
    tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
    !! print *,"minval(vNstep,mask=vNstep>0)", minval(vNstep,mask=vNstep>0)
    !! pause
    !
    allocate(vTPpath(DimPath))
    vTPpath(:)%TdgC= tTmp(nCp+1,1:DimPath)
    vTPpath(:)%Pbar= tTmp(nCp+2,1:DimPath)
    !
  end if
  !
  deallocate(vNstep)
  !
  if(allocated(vStAdd))  deallocate(vStAdd)
  !
  if(iDebug>1) call Pause_
  if(idebug>1) write(fTrc,'(A,/)') "</ Path_Read"
  
end subroutine Path_ReadParam_new

end module M_Path_Read
