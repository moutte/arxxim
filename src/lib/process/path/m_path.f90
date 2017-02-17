module M_Path
!--
!-- implements various path calculations
!--
  use M_Kinds
  use M_Trace,only: Stop_,fTrc,T_,iDebug,Pause_,Warning_
  use M_T_Component,only: T_Component

  implicit none
  !
  private
  !
  public:: Path_Execute
  !
  integer,parameter:: SafeStep= 255
  !
  type(T_Component),allocatable:: vCpnMix(:)
  real(dp),         allocatable:: vSys1(:), vSys2(:)
  !
  real(dp):: TdgKMix, PbarMix
  real(dp):: pH_,pO2
  !
contains

subroutine Path_Execute(EquMod5)
  !---------------------------------------------------------------------
  use M_Files,       only: NamFInn
  use M_System_Tools,only: System_TP_Update
  use M_TPcond_Read, only: TPpath_Read, TPgrid_Build
  use M_Path_Read,   only: Path_ReadMode, Path_ReadParam
  use M_Basis,       only: Basis_Change
  use M_Equil,       only: Equil_Calc
  use M_Equil_Read,  only: Equil_Read_PhaseAdd
  !
  use M_Global_Vars, only: vFas
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_Path_Vars,   only: Path_Vars_Clean,vTPpath,vFasFound
  !---------------------------------------------------------------------
  character(len=5), intent(in):: EquMod5
  !
  character(len=3) :: PathMod3
  character(len=3) :: EquMod3
  character(len=80):: Msg
  logical:: Ok,OkGrid,OkMix,Singular
  logical:: TPpath
  integer:: iErr,I
  !---------------------------------------------------------------------
  
  allocate(vFasFound(1:size(vFas)))
  vFasFound= .false.
  !
  TPpath= (EquMod5(4:5)=="TP")
  EquMod3= EquMod5(1:3)
  !
  !------------------------------------------------------- initialize --
  if(TPpath) then
    !
    call TPgrid_Build(OkGrid)
    if(.not. OkGrid) call TPpath_Read(TdgK,Pbar)
    !
    PathMod3= "TP_"
    !
    if(EquMod3(1:2)=="EQ") then
      if(count(vCpn(:)%Statut=="MOBILE")>0 .or. &
      &  count(vCpn(:)%Statut=="BUFFER")>0) then
        call Equil_Calc("SPC")
        call Basis_Change("DYN",vCpn)
      end if
    end if
    !
  else
    !
    call Path_ReadMode( & 
    & NamFInn, &
    & PathMod3,Ok,Msg)
    !
    if(.not.Ok) then
      if(iDebug>0) call Warning_(trim(Msg))
      return !--=================================< i/o error, return ==
    end if
    !
    !--------------------------- compute speciation of current system --
    call Equil_Calc("SPC")
    ! call Equil_Calc(EquMod3)
    !
    select case(PathMod3)
    !
    case("ADD")
      !
      !--- change basis, MOBILE > INERT
      call Basis_Change("EQU",vCpn)
      call System_TP_Update(TdgK,Pbar)
      !
    case("ADA")
      !
      !--- change basis, MOBILE > INERT
      call Basis_Change("EQU",vCpn)
      call System_TP_Update(TdgK,Pbar)
      !
      if(EquMod3(1:2)=="EQ") call Equil_Read_PhaseAdd(vFas)
      !
    case("MIX")
      !
      !--- save composition of system 1 to vSys1
      call System_TP_Update(TdgK,Pbar)
      !
      call Basis_Change("EQU",vCpn)
      !
      allocate(vSys1(1:size(vCpn)))
      vSys1(:)= vCpn(:)%Mole
      !
      !--- compute composition of system 2 and save to vSys2
      allocate(vCpnMix(1:size(vCpn)))
      !
      call Path_Calc_FluidMix(TdgKMix,PbarMix,OkMix)
      !
      allocate(vSys2(1:size(vCpn)))
      vSys2(:)= vCpnMix(:)%Mole
      !
    end select
    !---------------------------/compute speciation of current system --
    !
    if(EquMod3(1:2)=="EQ") call Equil_Read_PhaseAdd(vFas)
    !
    call Path_ReadParam( &
    & NamFInn,   &  !IN
    & PathMod3,  &  !IN
    & vCpn,      &  !IN
    & TdgK,Pbar, &  !IN
    & Ok,Msg)       !OUT
    !
    if(PathMod3=="PHP") then
      pH_= 1.0D0
      pO2= 50.0D0
    end if
    !
    if(.not.Ok) then
      if(iDebug>0) call Warning_(trim(Msg))
      return !--==================================< i/o error, return ==
    end if
    !
  end if
  !------------------------------------------------------/ initialize --
  !
  !----------------------------------------------------- compute path --
  call Path_Single(PathMod3,EquMod3)
  !----------------------------------------------------/ compute path --
  !
  !------------------------------------------------------------ clean --
  !
  if(PathMod3=="MIX") deallocate(vSys1,vSys2,vCpnMix)
  !
  call Path_Vars_Clean
  !-----------------------------------------------------------/ clean --
  !
  return
end subroutine Path_Execute

subroutine Path_Double
  
end subroutine Path_Double

subroutine Path_Single(PathMod3,EquMod3)
  !---------------------------------------------------------------------
  use M_Dtb_Const,   only: T_CK
  use M_System_Tools,only: System_TP_Update
  use M_Basis
  use M_Equil_Tools, only: Equil_Zero,Equil_Restart,Equil_Save
  use M_Equil_Tools, only: Equil_Trace_Init,Equil_Trace_Close
  use M_Equil_Tools, only: Equil_Errors,Equil_Errors_Warning
  use M_Equil_Specia,only: Equil_Specia
  use M_Equil_1,     only: Equil_Eq1
  use M_Equil_2,     only: Equil_Eq2
  use M_Equil_Write
  use M_Path_Write
  !
  use M_Global_Vars, only: vFas
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_Equil_Vars,  only: Equil_Vars_Clean,vYesList
  use M_Path_Vars,   only: vTPpath,vFasFound,DimPath,TotalMixStep,vLPath
  !---------------------------------------------------------------------
  character(len=3),intent(in):: PathMod3
  character(len=3),intent(in):: EquMod3
  !---------------------------------------------------------------------
  real(dp):: CpuBegin,CpuEnd
  integer :: I,nCp,nStep,iErr
  logical :: PathExit,PathRedox
  !---------------------------------------------------------------------
  
  call System_TP_Update(TdgK,Pbar)
  !
  call Basis_Change("SPC",vCpn)
  !
  call Equil_Zero(EquMod3)
  call Equil_Trace_Init
  !
  !----------------------- initial speciation in case PATH is ADD or MIX
  if(PathMod3=="ADD" .or. &
  &  PathMod3=="ADA" .or. &
  &  PathMod3=="MIX") then
    !
    call Equil_Restart
    !
    select case(EquMod3)
    
    case("SPC")
      call Equil_Specia(iErr)
    
    case("EQ1")
      call Equil_Eq1(.true.,iErr)
    
    case("EQ2")
      call Equil_Eq2(.true.,iErr)
      if(iErr<0) call Equil_Eq1(.true.,iErr)
      
    case("EQM")
      call Equil_Eq2(.false.,iErr)
      
    end select
    !
    if(iErr<0) call Equil_Errors(iErr)
    !
    call Equil_Save
    !
  end if
  !----------------------/ initial speciation in case PATH is ADD or MIX
  !
  if (trim(PathMod3)=="MIX") then
    call Basis_Change("EQU",vCpn)
    vSys1(:)= vCpn(:)%Mole
    !
    DimPath= TotalMixStep
    !
    allocate(vTPpath(DimPath))
    do I=1,TotalMixStep
      vTPpath(I)%TdgC= TdgK + (TdgKMix - TdgK)*I/TotalMixStep -T_CK
      vTPpath(I)%Pbar= Pbar
    end do
  end if
  !
  call Path_Write_FasEnTete
  !
  if(iDebug>2) call CPU_TIME(CpuBegin)
  !
  !----------------------------------------------------------- path loop
  nStep=1
  do
    !
    if(nStep>SafeStep) then
      if(iDebug>1) print '(A)',"Reached max' number of steps !!"
      exit !to prevent unterminated loops !!!
    end if
    !
    !if(nStep>size(vTPpath)) exit !---------------------- exit path loop
    !
    if(PathMod3 /= "LGK") then
      if( ABS(vTPpath(nStep)%TdgC +T_CK -TdgK) > 1.D-2 .or. &
      &   ABS(vTPpath(nStep)%Pbar       -Pbar) > 1.D-2) then
        !
        TdgK= vTPpath(nStep)%TdgC +T_CK
        Pbar= vTPpath(nStep)%Pbar
        if(iDebug>1) print &
        & '(A,I3,2G15.6)',"TP changed, at Step ",nStep,TdgK-T_CK,Pbar
        !
        call System_TP_Update(TdgK,Pbar)
        !
      end if
    end if
    !
    !-------------------------------------------- conditions of new step
    select case(PathMod3)
    !
    case("TP_")
    !--- path along the TP table --
      !
      if(nStep>size(vTPpath)) exit !----------------------exit path loop
      !
      !! !--- conditions of new step --
      !! TdgK= vTPpath(nStep)%TdgC +T_CK
      !! Pbar= vTPpath(nStep)%Pbar
      !! call System_TP_Update(TdgK,Pbar)
      !! !---/
      !
    case("PHP")
      call Path_NewStep_pHpE(PathExit)
      if(PathExit) pH_= pH_ +1.0D0
      if(pH_>13.0D0) exit
      !
    case("ADD","ADA","MIX","CHG","EVP","LGK")
    !--- chemical path defined by the PATH block --
      !
      !--- conditions of new step --
      call Path_NewStep(PathMod3,nStep,PathExit)
      !---/
      !
    end select
    !-------------------------------------------/ conditions of new step
    !
    !-------------------------------------------------- compute new step
    call Equil_Restart
    !
    select case(EquMod3)
    !
    case("SPC")
      call Equil_Specia(iErr)
    case("EQM")
      call Equil_Eq2(.false.,iErr)
      call FasFound_Record(vFas,vFasFound)
    case("EQ1")
      call Equil_Eq1(.true.,iErr)
      call FasFound_Record(vFas,vFasFound)
    case("EQ2")
      call Equil_Eq2(.true.,iErr)
      if(iErr<0) call Equil_Eq1(.true.,iErr)
      call FasFound_Record(vFas,vFasFound)
    !~ case("EQ3")
      !~ call Equil_Eq3(iErr)
      !~ call FasFound_Record(vFas,vFasFound)
    !
    end select
    !
    if(iErr<0) then
      call Equil_Errors_Warning(iErr)
      iErr= 0
      if(PathExit) exit !---------------------------------exit path loop
      nStep= nStep+1
      cycle
    end if
    !
    !~ if(iDebug==4) call Basis_Change_Wrk(vCpn)
    !
    call Equil_Save
    !-------------------------------------------------/ compute new step
    !
    select case(EquMod3)
    case("SPC"); call Path_Write_Line(WrCount=nStep,WrCod="SPC")
    case("EQM"); call Path_Write_Line(WrCount=nStep,WrCod="EQM",vYes=vYesList)  
    case("EQ1"); call Path_Write_Line(WrCount=nStep,WrCod="EQ1",vYes=vYesList)  
    case("EQ2"); call Path_Write_Line(WrCount=nStep,WrCod="EQ2",vYes=vYesList)   
    end select
    !
    if(iDebug==4) call Path_Write_Distrib(nStep)
    !
    call Path_Write_FasAff(nStep)
    !if(OkKinetics) call KinRateActiv_Test(nStep)
    !
    if(PathExit) exit !-----------------------------------exit path loop
    !
    nStep= nStep+1
    !
  end do
  !----------------------------------------------------------/ path loop
  if(iDebug>2) call CPU_TIME(CpuEnd)
  if(iDebug>2) print '(A,G15.6)',"CPU=", CpuEnd -CpuBegin
  !
  if(iDebug>0 .and. EquMod3(1:2)=="EQ") call FasFound_Sho(vFas,vFasFound)
  !
  call Equil_Trace_Close
  call Equil_Vars_Clean
  call Path_Files_Close
  !
  return
end subroutine Path_Single

subroutine Path_NewStep_pHpE(PathExit)
!
! 2.H2O=  O2 + 4.H+ + 4.e-
! 2.G(H2O)+lnA_H2O=  G(O2)+lnA_O2  +4.lnA_H+ +4.lnA_e-
! 4.(pH + pE)=  [G(O2) +lnA_O2 -2.G(H2O) -2.lnA_H2O]/ln10
! log10(Act(O2))=4*(pH_+pE_)+ 2[G(H2O) + lnA_H2O] - [G(O2) + lnA_O2]
!
! if redox species is H2aq:
! 2H+ + 2E-=  H2 -> pH + pE=  -( vSpc(iOx)%G0 + vLnAct(iOx)) /2
!
!  pE= & 
!  & (  vSpc(isO2)%G0rt + vSpc(isO2)%Dat%LAct          &
!  &  -(vSpc(isW )%G0rt + vSpc(isW )%Dat%LAct) *Two  ) &
!  & /LN10 /4.0D0 &
!  & - pH
!  
!  pE= vSpc(isO2)%G0rt     -vSpc(isW )%G0rt     *Two &
!  & + vSpc(isO2)%Dat%LAct -vSpc(isW )%Dat%LAct *Two
!  pE= pE /LN10 /4.0D0 - pH
!
  use M_Numeric_Const,only: Ln10
  use M_Global_Vars,  only: vSpc
  use M_Basis_Vars,   only: isH_,isO2
  !
  logical,intent(out):: PathExit
  !
  real(dp):: pHmin,pHmax !,pEmin,pEmax
  real(dp):: pH,pK_O2 !,pE
  !
  vSpc(isH_)%Dat%LAct= -pH_ *LN10
  vSpc(isO2)%Dat%LAct= -pO2 *LN10
  !
  pO2= pO2 + 1.0D0
  PathExit= (pO2>1.0D0)
  !
end subroutine Path_NewStep_pHpE

subroutine FasFound_Sho(vFas,vFasFound)
  use M_T_Phase
  type(T_Phase),intent(in):: vFas(:)
  logical,      intent(in):: vFasFound(:)
  !
  integer:: I
  !
  if(count(vFasFound)>0) then
    if(iDebug>0) then
      print '(/,A,/)',"phases encountered in PATHEQU:"
      do I=1,size(vFas)
        if(vFasFound(I)) print '(A)',vFas(I)%NamFs
      end do
    end if
  end if
  !
end subroutine FasFound_Sho

subroutine FasFound_Record(vFas,vFasFound)
!--
!-- keep record of phases found during whole path --
!--
  use M_T_Phase
  type(T_Phase),intent(in):: vFas(:)
  logical,   intent(inout):: vFasFound(:)
  !
  integer:: iFs
  !
  do iFs=1,size(vFas) 
    if(vFas(iFs)%MolFs>Zero) vFasFound(iFs)=.true.
  end do
  !
  return
end subroutine FasFound_Record

subroutine Path_NewStep(PathMod3,nStep,PathExit)
  use M_Files
  use M_Numeric_Const,only: Ln10
  use M_Dtb_Const,    only: T_CK
  use M_T_Species,    only: Species_Index
  use M_T_MixPhase,   only: T_MixPhase, MixPhase_CalcActivs
  !
  use M_Global_Vars,only: vEle,vSpc,vFas,vMixFas,vMixModel,nAq
  use M_System_Vars,only: vCpn,TdgK,Pbar
  use M_Basis_Vars, only: isW,iH_,MWSv,tAlfFs,iBal
  use M_Path_Vars,  only: &
  & vPhasBegin,vPhasFinal,vFasFound, &
  & vLPath,tPathData,DimPath,TotalMixStep, &
  & iLogK,vPathLogK
  !
  character(len=3),intent(in) :: PathMod3
  integer,         intent(in) :: nStep
  logical,         intent(out):: PathExit
  !
  integer,parameter:: AdjustMode= 1
  !
  integer :: iCp,J,iFs,iP,NPole,nCp
  real(dp):: X
  !
  type(T_MixPhase):: Fas, PhasBegin, PhasFinal
  
  PathExit=.false.
  !
  nCp= size(vCpn)
  !
  if (iDebug>0) write(fTrc,'(/,A,I3)') "< Path_NewStep ==== step= ", nStep
  
  select case(PathMod3)
  
  !---------------------------------------------------------- case "LGK"
  case("LGK")
    if(iDebug>1) print '(A,I3,G15.6)',"nStep",nStep,vPathLogK(nStep)
    !
    vSpc(iLogK)%G0rt= -vPathLogK(nStep) *Ln10
    do J=1,size(vFas)
      if(vFas(J)%iSpc == iLogK) vFas(J)%Grt= vSpc(iLogK)%G0rt
    end do
    PathExit= nStep==size(vPathLogK)
    !print *,"nStep          = ",nStep
    !print *,"SIZE(vPathLogK)= ",size(vPathLogK)
    !pause
  !---------------------------------------------------------/ case "LGK"
  
  !--------------------------------------------- case "EVP" (-EVAPORATE)
  case("EVP")
    X= vCpn(isW)%Mole *0.10D0
    vCpn(isW)%Mole= vCpn(isW)%Mole - X
    vCpn(iH_)%Mole= vCpn(iH_)%Mole - X*2.0D0
    PathExit= nStep==255
    
  !------------------------------------------------ case "CHG" (-CHANGE)
  case("CHG")
  ! case("CHG","EVP")
    !
    if(iDebug>0) print '(A,I3)',"nStep",nStep
    !
    DoComponent: do iCp=1,nCp
      !
      if(vLPath(iCp)) then
      !----------------------------------------- change in component iCp
        J= vCpn(iCp)%iSpc !-> related species
        !
        select case(trim(vCpn(iCp)%Statut))
        
        case("MOBILE","BUFFER")
        !------------------ case of CHANGE on mobile or buffer component
        !-------------------- -> change ACTIVITY -> change vSpc(J)%LnAct
          
          if(vCpn(iCp)%iMix==0) then
          !------------------ mobile species is aqueous or in pure phase
            vCpn(iCp)%LnAct= - tPathData(iCp,nStep) *Ln10
            vSpc(J)%Dat%LAct= vCpn(iCp)%LnAct
            !
            PathExit= (nStep==DimPath)
            !
          else
          !----------------- mobile species in non-aqueous mixture phase
            PhasBegin= vMixFas(vPhasBegin(iCp))
            PhasFinal= vMixFas(vPhasFinal(iCp))
            Fas= PhasBegin
            !
            NPole= vMixModel(Fas%iModel)%NPole
            do iP=1,NPole
              if(Fas%vLPole(iP)) then
                Fas%vXPole(iP)= &
                & (PhasBegin%vXPole(iP)*(100 - nStep)  &
                +  PhasFinal%vXPole(iP)*(nStep)      ) &
                & /100.0D0
              end if
            end do
            !
            call MixPhase_CalcActivs( &
            & TdgK,Pbar,&
            & vMixModel(Fas%iModel),&
            & Fas)
            !
            !----------------------------------------------------- debug
            if(iDebug>2) then 
              write(fTrc,'(A15,A1)',advance="no") "PoleFraction=",T_
              do iP=1,NPole
                if(Fas%vLPole(iP)) &
                & write(fTrc,'(G15.6,A1)',advance="no") Fas%vXPole(iP),T_
              end do
              write(fTrc,*)
              write(fTrc,'(A15,A1)',advance="no") "Activity=",T_
              do iP=1,NPole
                if(Fas%vLPole(iP)) &
                & write(fTrc,'(G15.6,A1)',advance="no") exp(Fas%vLnAct(iP)),T_
              end do
              write(fTrc,*)
            end if
            !-----------------------------------------------------/debug
            !
            vSpc(J)%Dat%LAct= Fas%vLnAct(vCpn(iCp)%iPol) !not necessary ???
            vCpn(iCp)%LnAct=  vSpc(J)%Dat%LAct
            !
            if(iDebug>0) write(fTrc,'(A4,A3,A1, A4,A23,A1, A7,G15.6)') &
            & "Cpn=",vEle(vCpn(iCp)%iEle)%NamEl,T_, &
            & "Spc=",vSpc(J)%NamSp,T_, &
            & "LogAct=", vSpc(J)%Dat%LAct/Ln10
            !
            PathExit=(nStep==99)
            !
          end if
          !--------------------/ mobile species in non-aqueous sol'phase
        !-----------------/ case of CHANGE on mobile or buffer component

        case("INERT") 
        !----------------------------- case of CHANGE on INERT component
        !---- for inert species, change MOLE -> change vCpn(iCp)%Mole
          vCpn(iCp)%Mole= tPathData(iCp,nStep)
          PathExit= (nStep==DimPath)
          !
          !------------------------------ bri-collage !!!§§§->rework !!!
          !-- vCpn will be used for the FLUID system !!!
          !-- it coincides with total system in SPC mode,
          !-- but, in EQU mode, we control here the total system ...
          !-- substract mole nrs from non-fluid phases
          !-- (case of EQn calculation)
          !-------------------------------------------------------------
          do iFs=1,size(vFas)
            if(vFas(iFs)%MolFs>Zero) &
            & vCpn(iCp)%Mole= vCpn(iCp)%Mole - tAlfFs(iCp,iFs) *vFas(iFs)%MolFs
          end do
          !
          do iFs=1,size(vFas)
            if(vFas(iFs)%MolFs>Zero) vFasFound(iFs)=.true.
          end do
          !
        !----------------------------/ case of CHANGE on INERT component
        end select
      end if ! if(vLPath(iCp))
    end do DoComponent ! 
    !
    if(PathMod3=="EVP") then
      do iFs=1,size(vFas)
        if(vFas(iFs)%MolFs>Zero) vFas(iFs)%MolFs= 1.D-12
      end do
    end if
    !
    call NewStep_Adjust_TotO_TotH( & !
    & AdjustMode,        &           !
    & iBal,isW,iH_,vEle, &           !IN
    & vCpn)                          !INOUT
    !
  !----------------------------------------------------------/case "CHG"
  
  !---------------------------------------------------------- case "ADD"
  case("ADD","ADA")
    !
    if(iDebug>0) print '(A,I3)',"nStep",nStep
    !
    vCpn(:)%Mole= tPathData(:,nStep)
    !
    !------------------------------------ bri-collage !!!§§§->rework !!!
    !-- vCpn is the FLUID system !!!
    !-- it coincides with total system in SPC mode,
    !-- but, in EQU mode, we control here the total system ...
    !---substract mole nrs from non-fluid phases (case of EQn calculation)
    do iCp=1,size(vCpn)
      do iFs=1,size(vFas)
        if(vFas(iFs)%MolFs>Zero) &
        & vCpn(iCp)%Mole= vCpn(iCp)%Mole - tAlfFs(iCp,iFs) *vFas(iFs)%MolFs
      end do
    end do
    !------------------------------------------------------------------/
    !
    ! if(iDebug>2) then
    !   write(61,'(A,1X)',advance="no") "Path_NewStep"
    !   do J=1,size(vCpn)
    !     write(61,'(G15.6,1X)',advance="no") vCpn(J)%Mole
    !   end do
    !   write(61,*)
    ! end if
    !
    PathExit= nStep==DimPath
  !----------------------------------------------------------/case "ADD"
  
  !---------------------------------------------------------- case "MIX"
  case("MIX")
    !
    if(nStep==TotalMixStep) then
      PathExit=.true.
    else
      do iCp=1,size(vCpn)
        !
        vCpn(iCp)%Mole= vSys1(iCp) &
        &             + (vSys2(iCp) -vSys1(iCp)) *nStep/real(TotalMixStep)
        !
        if(iDebug>0) write(fTrc,'(A3,1X,3G15.6,1X,A)') &
        & vEle(vCpn(iCp)%iEle)%NamEl, &
        & vSys1(iCp),vSys2(iCp),vCpn(iCp)%Mole,&
        & vCpn(iCp)%Statut
        !
      end do
      !
    end if
  !----------------------------------------------------------/case "MIX"
  
  end select ! case(PathMod3)
  
  !------------------------ keep record of phase found during whole path
  do iFs=1,size(vFas) !
    if(vFas(iFs)%MolFs>Zero) vFasFound(iFs)=.true.
  end do
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ Path_NewStep"
  !
end subroutine Path_NewStep

subroutine NewStep_Adjust_TotO_TotH( & !
& AdjustMode, &                        !
& iBal,ic_W,ic_H,vEle, &               !
& vCpn)
!--
!-- adjust total molar amounts of O and H according 
!-- to molar amounts of other elements
!-- -> values of vTotF
!--
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  !
  integer,          intent(in)   :: AdjustMode
  integer,          intent(in)   :: iBal
  integer,          intent(in)   :: ic_W,ic_H
  type(T_Element),  intent(in)   :: vEle(:)
  type(T_Component),intent(inout):: vCpn(:)
  !
  real(dp):: Zbalance,TotO_
  integer :: I,iBal_
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< NewStep_Adjust_TotO_TotH"
  !
  TotO_= vCpn(ic_W)%Mole
  !
  if(iBal==0) then ;  iBal_= ic_H
  else             ;  iBal_= iBal
  end if

  if(iBal_ /= ic_H .and. vCpn(ic_H)%Statut=="INERT") then
    vCpn(ic_W)%Mole= TotO_
    vCpn(ic_H)%Mole= TotO_*Two
    !return
  end if
  
  Zbalance= Zero
  do I=1,size(vCpn)
    if(vCpn(I)%Statut=="INERT" .and. I/=ic_W .and. I/=ic_H) then
      !print *,vEle(vCpn(I)%iEle)%NamEl, vEle(vCpn(I)%iEle)%Z
      Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
    end if
  end do
  !pause

  select case(AdjustMode)
  
  case(1)
    !--- if excess cation, equiv'Na > equiv'Cl, adjust with OH
    if(Zbalance > Zero) then
      vCpn(ic_W)%Mole=  TotO_     +Zbalance
      vCpn(iBal_)%Mole= TotO_*Two +Zbalance
    !--- if excess anion, Cl>Na, adjust with H
    else
      vCpn(ic_W)%Mole=  TotO_
      vCpn(iBal_)%Mole= (TotO_*Two -Zbalance) /real(vEle(iBal_)%Z)
    end if

  case(2)
    !---------- Oxygen number unchanged, equilibrate using hydrogen only
    vCpn(ic_W)%Mole= TotO_
    vCpn(iBal_)%Mole= TotO_*Two -Zbalance
    
  end select
  
  !! !write(12,'(3G15.6)') Zbalance,vCpn(ic_W)%Mole,vCpn(iBal_)%Mole
  !! 
  !! if(iBal /= ic_H) then
  !!   vCpn(ic_W)%Mole= TotO_
  !!   vCpn(ic_H)%Mole= TotO_*Two
  !!   !return
  !! end if
  !! 
  !! Zbalance= Zero
  !! do I=1,size(vCpn)
  !!   if(vCpn(I)%Statut=="INERT" .and. I/=ic_W .and. I/=ic_H) &
  !!   & Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
  !! end do
  !! !
  !! select case(AdjustMode)
  !!
  !! case(1)
  !!   !--- if excess cation, Na>Cl, adjust with OH
  !!   if(Zbalance > Zero) then
  !!     vCpn(ic_W)%Mole= TotO_     +Zbalance
  !!     vCpn(iBal)%Mole= TotO_*Two +Zbalance
  !!   !--- if excess anion, Cl>Na, adjust with H
  !!   else
  !!     vCpn(ic_W)%Mole=  TotO_
  !!     vCpn(iBal)%Mole= (TotO_*Two -Zbalance) /real(vEle(iBal)%Z)
  !!   end if
  !! 
  !! case(2)
  !!   !--- Oxygen number unchanged, equilibrate using hydrogen only
  !!   vCpn(ic_W)%Mole= TotO_
  !!   vCpn(iBal)%Mole= TotO_*Two -Zbalance
  !! 
  !! end select
  !
  if(iDebug>2) then
    write(fTrc,'(A,G15.6)')  "Zbalance       =", Zbalance
    write(fTrc,'(A,G15.6)')  "vCpn(ic_W)%Mole=", vCpn(ic_W)%Mole
    write(fTrc,'(A,G15.6)')  "vCpn(ic_H)%Mole=", vCpn(ic_H)%Mole
  end if
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ NewStep_Adjust_TotO_TotH"
  !
end subroutine NewStep_Adjust_TotO_TotH
  !
subroutine Path_Calc_FluidMix(TdgKMix,PbarMix,Ok)
!--
!-- read constraints on mix fluid,
!-- then calc' speciation, then make all cpn inert,
!--
  use M_T_Component, only: T_Component
  use M_Basis,       only: Basis_Change
  use M_System_Tools,only: System_Build_Custom,System_TP_Update
  use M_Equil
  !
  use M_Global_Vars, only: vEle,vSpc,vMixFas
  use M_System_Vars, only: vCpn,TdgK,Pbar
  !
  real(dp),intent(out):: TdgKMix,PbarMix
  logical, intent(out):: Ok
  !
  type(T_Component),dimension(:),allocatable:: vC
  integer :: I,J,N
  real(dp):: T0,P0
  !
  allocate(vC(1:size(vCpn)))
  vC= vCpn
  T0= TdgK
  P0= Pbar
  !
  !if(allocated(vCpnMix)) deallocate(vCpnMix)
  vCpnMix= vCpn
  !
  call System_Build_Custom( &
  & "SYSTEM.MIX",vEle,vSpc,vMixFas,vCpn, & !IN
  & TdgK,Pbar, &
  & vCpnMix,Ok)
  !
  if(Ok) then
    vCpn= vCpnMix
    !
    call System_TP_Update(TdgK,Pbar)
    
    call Equil_Calc("SPC")
    !-> calc' speciation for mixing end-member system
    call Basis_Change("EQU",vCpn)
    !-> make all components inert
    vCpnMix= vCpn
  end if
  !
  TdgKMix= TdgK
  PbarMix= Pbar
  !
  vCpn= vC
  TdgK= T0
  Pbar= P0
  !
  !! if(PbarMix/=Pbar) then
  !!   call Warning_("pressure should be the same as master system !!")
  !!   PbarMix= Pbar
  !! end if
  !
  if(Ok) then
    !permute vCpnMix -> its element order must be same as vCpnBox
    !-> use vC for swapping
    vC= vCpnMix
    N= size(vCpn)
    do I=1,N
      do J=1,N; if(vC(J)%iEle==vCpn(I)%iEle) exit; end do
      vCpnMix(I)= vC(J)
    end do
  end if
  deallocate(vC)
  !
  return
end subroutine Path_Calc_FluidMix

subroutine Path_Record
  use M_Global_Vars,only: vFas
  use M_Path_Vars,  only: vFasFound
  !  
  integer:: iFs
  do iFs=1,size(vFas) !keep record of phase found during whole path
    if(vFas(iFs)%MolFs>Zero) vFasFound(iFs)=.true.
  end do
  !
end subroutine Path_Record

end module M_Path

