module M_Equil_Tools
!--
!-- common tools for equilibrium calculations
!--
  use M_Trace,only: iDebug,fTrc,fHtm,T_,Stop_,Pause_,Warning_
  use M_Info_Value
  use M_Kinds
  !
  implicit none
  !
  private
  !
  public:: Equil_Zero      !for a new speciation or equilibrium
  public:: Equil_Restart
  public:: Equil_Save      !store results of Equil_Calc in vCpn, vSpc 
  public:: Equil_Errors
  public:: Equil_Errors_Warning
  !
  public:: Equil_CalcFasAff
  !
  public:: Equil_Trace_Init
  public:: Equil_Trace_Close
  !
  public:: Equil_FasTrace_EnTete
  public:: Equil_FasTrace_Write
  
contains

subroutine Equil_Errors_Warning(iErr)
  integer,intent(in):: iErr
  !
  if (iErr>=0) then
  
    call Info_Value("Error Code", iErr)
    call Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  else
    
    select case(iErr)
    case(-1)  ;  call Warning_("Equil_Error: NewtMaxIts Exceeded in Newton")
    case(-2)  ;  call Warning_("Equil_Error: Singular Matrix in Newton")
    case(-3)  ;  call Warning_("Equil_Error: Roundoff Problem in LineSearch")
    case(-4)  ;  call Warning_("Equil_Error: No Convergence on Gammas")
    case(-5)  ;  call Warning_("Equil_Error: Stationary Point in Newton")
    case(-7)  ;  call Warning_("Equil_Error: MaxIter Exceeded in Outer Loop")
    case default
      call Info_Value("Error Code", iErr)
      call Warning_("Equil_Error: Unknown Code")
    end select
    
  end if
end subroutine Equil_Errors_Warning

subroutine Equil_Errors(iErr)
  integer,intent(in):: iErr
  !
  if (iErr>=0) then
  
    call Info_Value("Error Code", iErr)
    call Warning_("Equil_Errors Called with iErr >= 0 !!")
  
  else
    
    select case(iErr)
    case(-1)  ;  call Stop_("Equil_Error: NewtMaxIts Exceeded in Newton")
    case(-2)  ;  call Stop_("Equil_Error: Singular Matrix in Newton")
    case(-3)  ;  call Stop_("Equil_Error: Roundoff Problem in LineSearch")
    case(-4)  ;  call Stop_("Equil_Error: No Convergence on Gammas")
    case(-5)  ;  call Stop_("Equil_Error: Stationary Point in Newton")
    case(-7)  ;  call Stop_("Equil_Error: MaxIter Exceeded in Outer Loop")
    case default
      call Info_Value("Error Code", iErr)
      call Stop_("Equil_Error: Unknown Code")
    end select
    
  end if
end subroutine Equil_Errors

subroutine Equil_Init_TotO_TotH( & !
& AdjustMode, &                       !
& iBal,ic_W,ic_H,vEle, &              !IN
& vCpn)                               !INOUT
!--
!-- adjust total molar amounts of O and H according 
!-- to molar amounts of other elements
!-- -> values of vTotF
!--
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  !---------------------------------------------------------------------
  integer,          intent(in)   :: AdjustMode
  integer,          intent(in)   :: iBal
  integer,          intent(in)   :: ic_W,ic_H
  type(T_Element),  intent(in)   :: vEle(:)
  type(T_Component),intent(inout):: vCpn(:)
  !---------------------------------------------------------------------
  real(dp):: Zbalance,TotO_
  integer :: I,iBal_
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Init_TotO_TotH"
  !
  TotO_= vCpn(ic_W)%Mole
  !
  if(iBal==0) then ;  iBal_= ic_H
  else             ;  iBal_= iBal
  end if
  !
  if(iBal_ /= ic_H .and. vCpn(ic_H)%Statut=="INERT") then
    vCpn(ic_W)%Mole= TotO_
    vCpn(ic_H)%Mole= TotO_*Two
    !return
  end if
  !
  Zbalance= Zero
  do I=1,size(vCpn)
    if(vCpn(I)%Statut=="INERT" .and. I/=ic_W .and. I/=ic_H) then
      Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
    end if
  end do
  !
  !! select case(AdjustMode)
  !! !
  !! case(1)
  !!   !--- if excess cation, Na>Cl, adjust with OH
  !!   if(Zbalance > Zero) then
  !!     vCpn(ic_W)%Mole= TotO_     +Zbalance
  !!     vCpn(ic_H)%Mole= TotO_*Two +Zbalance
  !!   !--- if excess anion, Cl>Na, adjust with H
  !!   else
  !!     vCpn(ic_W)%Mole= TotO_
  !!     vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !!   end if
  !! !
  !! case(2)
  !!   !--- Oxygen number unchanged, equilibrate using hydrogen only
  !!   vCpn(ic_W)%Mole= TotO_
  !!   vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !! !
  !! end select
  !
  select case(AdjustMode)  
  case(1)
    !----------------------------if excess cation, Na>Cl, adjust with OH
    if(Zbalance > Zero) then
      vCpn(ic_W)%Mole=  TotO_     +Zbalance
      vCpn(iBal_)%Mole= TotO_*Two +Zbalance
    !------------------------------if excess anion, Cl>Na, adjust with H
    else
      vCpn(ic_W)%Mole=  TotO_
      vCpn(iBal_)%Mole= (TotO_*Two -Zbalance) /real(vEle(iBal_)%Z)
    end if
  case(2)
    !-----------Oxygen number unchanged, equilibrate using hydrogen only
    vCpn(ic_W)%Mole= TotO_
    vCpn(iBal_)%Mole= TotO_*Two -Zbalance
  end select
  !
  if(iDebug>2) then !----------------------------------------------trace
    write(fTrc,'(A,G15.6)')  "Zbalance       =", Zbalance
    write(fTrc,'(A,G15.6)')  "vCpn(ic_W)%Mole=", vCpn(ic_W)%Mole
    write(fTrc,'(A,G15.6)')  "vCpn(ic_H)%Mole=", vCpn(ic_H)%Mole
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Init_TotO_TotH"
  !---------------------------------------------------------------/trace
  !  
end subroutine Equil_Init_TotO_TotH

subroutine Equil_Zero(Cod)
!--
!-- for "very initial" speciation
!--
  use M_T_Species,  only: T_SpcData
  use M_Equil_Read ! <-equil_Read_Debug,Equil_Read_Conditions,Equil_Read_PhaseAdd,Equil_Read_YesList
  use M_Equil_Vars, only: Equil_Vars_Alloc
  !
  use M_Global_Vars,only: vSpc,vEle,vFas,nAq,vSpcDat
  use M_Global_Vars,only: SolModel,vMixModel
  use M_System_Vars,only: vCpn,iO_,iH_
  !
  use M_Basis_Vars, only: iBal,initMolNu
  use M_Basis_Vars, only: vLAx,vLMx,vLCi,nAs,tAlfFs
  use M_Equil_Vars, only: DebFormula,DebNewt,DebJacob,vYesList
  use M_Equil_Vars, only: LogForAqu,DirectSub,cMethod,bFinDif
  use M_Equil_Vars, only: vLnAct
  use M_Equil_Vars, only: vTotS
  use M_Equil_Vars, only: nEquFas
  !
  character(len=3),intent(in):: Cod
  !!type(T_Component),dimension(:),intent(inout):: vCpn
  !
  integer,parameter:: AdjustMode= 2 !!1
  !
  logical :: Error,OK
  character(len=80):: ErrorMsg
  integer :: I,iSpc
  integer :: isW
  real(dp):: MWSv
  type(T_SpcData),allocatable:: vDat(:)
  !
  logical:: test_initmolnu= .false.
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Zero"
  !
  isW=  SolModel%iSolvent
  MWSv= SolModel%MolWeitSv
  !------------------------------------------------------ default values
  bFinDif=   .false.
  LogForAqu= .true.
  DirectSub= .false.  !! .true.  !! 
  initMolNu= 1.0D-6
  !-----------------------------------------------------/ default values
  !
  !------------------------------------------------------read CONDITIONS
  call Equil_Read_Debug(DebFormula,DebNewt,DebJacob)
  call Equil_Read_Numeric( &
  & initMolNu,LogForAqu,DirectSub,cMethod,bFinDif,Error,ErrorMsg)
  !
  !if(.not. LogForAqu) DirectSub= .true.
  !-----------------------------------------------------/read CONDITIONS
  !
  call Equil_Init_TotO_TotH( & !
  & AdjustMode,        &          !IN
  & iBal,iO_,iH_,vEle, &          !IN
  & vCpn)                         !INOUT
  !-> calc. mole.nrs of O and H in fluid from mole.nrs of other elements
  !
  vFas(:)%MolFs= Zero
  nEquFas= 0
  !
  if(Cod(1:2)=="EQ") call Equil_Read_PhaseAdd(vFas)
  !
  !!!-> read SYSTEM.ROCK -> modify vFas(:)%Mole
  !!! new deleted 18/10/2007 18:25
  !! do I=1,size(vFas)
  !!   vCpn(:)%Mole= vCpn(:)%Mole + tAlfFs(:,I) *vFas(I)%Mole
  !! end do
  !
  !--------------------------------------------------------- allocations
  call Equil_Vars_Alloc(size(vEle),size(vSpc),nAq,nAs,size(vFas))
  !! print *,"vLnAct(isW)=", vLnAct(isW)  ;  call pause_
  !--------------------------------------------------------/ allocations
  !
  !! if(trim(Cod)/="SPC") call Equil_Read_YesList(vFas,vYesList)
  if(Cod(1:2)=="EQ") call Equil_Read_YesList(size(vFas),vFas,vYesList)
  !
  vSpcDat(:)%LGam= Zero !--------------------------------initial gamma's
  !
  !----------------------------------------initial mole numbers in fluid
  ! initialize vSpcDat(1:nAq)%Mole, used as initial guess for solver
  vSpcDat(1:nAq)%Mole= InitMolNu !initial value for all aqu'species
  vSpcDat(isW)%Mole=   One/MWSv  !except Solvent
  !
  if(test_initmolnu) then
    do I=1,size(vCpn)
      if(I/=iH_) then
        iSpc= vCpn(I)%iSpc
        if(vLCi(iSpc)) vSpcDat(iSpc)%Mole= vCpn(I)%Mole
      end if
    end do
  end if
  !
  !-------------------------------------------------------read from file
  if(iDebug==1) then
    allocate(vDat(size(vSpcDat)))
    call Equil_Read_FromFile(vDat,OK)
    if(OK) then
      do i=1,size(vSpcDat)
        if(.not.(vLAx(i) .or. vLMx(i))) then
          vSpcDat(i)= vDat(i)
        end if
      end do
    end if
    deallocate(vDat)
  end if
  !-----------------------------------------------------------/from file
  !
  != for mobile aqu'species, compute mole numbers estimates
  != from their activities (which are fixed) and the current gamma's
  do I=1,nAq
    if(vLAx(I)) &
    & vSpcDat(I)%Mole= exp(vSpcDat(I)%LAct -vSpcDat(I)%LGam)
    !this gives molality, should multiply by mass of solvent ...
  end do
  !---------------------------------------/initial mole numbers in fluid
  !!
  !call Specia_InitBalance(vTotS) !!NEW201010!!
  if(iBal/=0) vTotS(iBal)= Zero
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Zero"
  !
end subroutine Equil_Zero

subroutine Equil_Restart
!--
!-- restart a speciation, using vCpn%Mole, vSpcDat%LnAct, vSpcDat%LnGam
!--
  use M_IoTools,    only: OutStrVec
  !
  use M_Global_Vars,only: vSpc,nAq,vFas,vSpcDat
  !
  use M_System_Vars,only: vCpn
  !
  use M_Basis_Vars, only: iBal,nCi,nAs,tNuAs,tNuFas,tAlfFs,vOrdPr,vOrdAs
  !
  use M_Equil_Vars, only: vDeltaG_As,vTotS,vLnAct,vLnGam,vMolF
  use M_Equil_Vars, only: vDeltaG_Fs,vFasMole,vAffScale
  !
  integer:: I
  !
  !! DEBUGG !!
  ! write(6,'(A)') "==Equil_Restart=="
  ! do I=1,size(vFas)
  !   write(6,'(A)')  trim(vFas(I)%NamFs)
  ! end do
  ! pause

  !----------------------------------------------------- update deltaG's
  ! update deltaG's of second'aqu'species
  ! (may have changed if TP changed ...)
  do I=1,nAs
    vDeltaG_As(I)= &
    & vSpc(vOrdAs(I))%G0rt - dot_product(tNuAs(I,:),vSpc(vOrdPr(:))%G0rt)
  end do
  !
  ! update deltaG's of phases (may have changed if TP changed ...)
  ! DeltaG is free energy change for G(phase) - G(prim'species)
  ! i.e. DeltaG of reaction of formatION of phase from prim'species
  do I=1,size(vFas)
    vDeltaG_Fs(I)= &
    & vFas(I)%Grt - dot_product(tNuFas(I,:),vSpc(vOrdPr(:))%G0rt)
  end do
  !----------------------------------------------------/ update deltaG's
  !
  do I=1,size(vFas)
    vAffScale(I)= sum(abs(tNuFas(I,:)))
  end do
  !
  ! retrieve data from current vSpcDat 
  vLnAct(:)=    vSpcDat(:)%LAct
  vLnGam(:)=    vSpcDat(:)%LGam
  vMolF(1:nAq)= vSpcDat(1:nAq)%Mole
  !
  ! retrieve data from current vFas
  vFasMole(:)=  vFas(:)%MolFs
  !
  ! vTotS= mole nr components in SYSTEM (fluid + min), RHS in residual
  ! == in SPC, vTotS = vTot_inFluid
  ! == in EQU, vTotS = vTot_inFluid + vTotM (amount in equilibrium minerals)
  !
  ! retrieve data from current vCpn
  vTotS(1:nCi)= vCpn(1:nCi)%Mole ! mole nrs in fluid
  !---------------------------------------------------------- REWORK !!!
  != add mole nrs from non-fluid (case of PATHEQ* calculations)
  ! vTotS= mole nr components in SYSTEM (fluid + min), RHS in residual
  do I=1,size(vFas)
    vTotS(1:nCi)= vTotS(1:nCi) + tAlfFs(1:nCi,I) *vFasMole(I)
  end do
  !
  if(iBal/=0) vTotS(iBal)= Zero
  !
  return
end subroutine Equil_Restart

subroutine Specia_InitBalance( &
!---------------------------------------------------------------NOT USED
& iBal, &
& vOrdPr,vOrdAs,vOrdMs, &
& tAlfSp,tAlfPr,tAlfAs,tAlfMs,tAlfFs, &
& vTot_)
  use M_Basis,     only: Basis_InitBalance
  !
  integer, intent(in)   :: iBal
  integer, intent(in)   :: vOrdPr(:),vOrdAs(:),vOrdMs(:)
  real(dp),intent(inout):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:),tAlfFs(:,:)
  real(dp),intent(inout):: vTot_(:)
  !
  !! if H is "INERT", mole balance on H replaced by mole balance on charges  !
  !! if(trim(vCpn(iH_)%Statut)=="INERT") iBal=iH_ !NEW_17/10/2007 17:40
  !
  if(iBal/=0) call Basis_InitBalance( & !
  & iBal,                 & !IN
  & vOrdPr,vOrdAs,vOrdMs, & !IN
  & tAlfSp,tAlfPr,tAlfAs,tAlfMs, & !INOUT
  & tAlfFs)                 !INOUT
  !
  if(iBal/=0) vTot_(iBal)= Zero
  !
end subroutine Specia_InitBalance

subroutine Equil_Save !(lEquil)
!--
!-- Update vSpc & vCpn for use in next speciations
!--
  use M_Global_Vars,only: vEle,vSpc,vFas,nAq,vSpcDat
  use M_Global_Vars,only: vMixModel,vMixFas
  use M_System_Vars,only: vCpn
  use M_Basis_Vars, only: nCi,nAx,nAs
  use M_Basis_Vars, only: tAlfPr,tAlfAs,tAlfFs,vOrdPr,vOrdAs
  use M_Basis_Vars, only: tStoikio
  use M_Equil_Vars, only: vMolF,vLnAct,vLnGam,vFasMole,nEquFas
  !
  use M_Numeric_Const,only: Ln10
  !
  ! logical,intent(in),optional:: lEquil
  !
  ! logical :: B
  integer :: iCp
  
  ! B=.false.
  ! if(present(lEquil)) B=lEquil
  
  !-------------- store results in vSpcDat to use in subsequent calcul's
  vSpcDat(1:nAq)%Mole= vMolF(1:nAq) !/vMolF(1)/MWsv
  vSpcDat(:)%LAct=     vLnAct(:) !_____id
  vSpcDat(:)%LGam=     vLnGam(:) !_____id
  
  !----------------- store results in vFas to use in subsequent calcul's
  if(nEquFas>0) call Save_EquPhase
  !----/
  !---------------- compute mole nrs components in fluid -> vCpn(:)%Mole
  do iCp=1,size(vCpn)
    vCpn(iCp)%Mole= &
    &  sum(tAlfPr(iCp,1    :nCi    ) *vMolF(vOrdPr(1:nCi)))  &
    +  sum(tAlfPr(iCp,nCi+1:nCi+nAx) *vMolF(vOrdPr(nCi+1:nCi+nAx))) &
    +  sum(tAlfAs(iCp,1    :nAs    ) *vMolF(vOrdAs(1:nAs)))
  end do
  !----/
  !
  !---------------------------------------------------------save to file
  ! if(iDebug==1) then 
    call Equil_SaveSpc_ToFile(vSpcDat)
    call Equil_SaveCpn_ToFile(nAq,tStoikio,vSpcDat)
    call Equil_SaveFas_ToFile(vMixModel,vFas,vMixFas)
  ! end if
  !--------------------------------------------------------/save to file
  !
  !----------------------------------------------------------------trace
  if(iDebug==4) then
    print '(A)',"Equil_Save, check components' mole nrs"
    do iCp=1,size(vCpn) !
      print '(2A,2G15.6)',&
      & "Equil_Save, ", vEle(vCpn(iCp)%iEle)%NamEl,vCpn(iCp)%Mole
    end do
    call Pause_
  end if
  !---------------------------------------------------------------/trace
  !
  return
end subroutine Equil_Save

subroutine Equil_SaveFas_ToFile(vMixModel,vFas,vMixFas)
  use M_IoTools,   only: GetUnit
  use M_T_Phase,   only: T_Phase
  use M_T_MixModel,only: T_MixModel
  use M_T_MixPhase,only: T_MixPhase
  !
  type(T_MixModel), intent(in):: vMixModel(:)
  type(T_Phase),    intent(in):: vFas(:)
  type(T_MixPhase), intent(in):: vMixFas(:)
  !
  integer:: fo,iFs,iModel,J 
  !
  call GetUnit(fo)
  open(fo,file="tmp_phase.tab")
  !
  do iFs=1,size(vFas)

    if(vFas(iFs)%MolFs > Zero) then

      write(fo,'(A,1X,G15.6)', advance="NO") &
      & trim(vFas(iFs)%NamFs),vFas(iFs)%MolFs

      if(vFas(iFs)%iMix>0) then
        write(fo,'(A)',advance="NO") " X(:)="
        !print *, "vFas(iFs)%iMix",vFas(iFs)%iMix
        !pause
        iModel= vMixFas(vFas(iFs)%iMix)%iModel
        !print *, "iModel, nPole=",iModel,vMixModel(iModel)%nPole
        !pause
        do J=1,vMixModel(iModel)%nPole
          write(fo,'(G15.6,1X)',advance="NO") &
          & vMixFas(vFas(iFs)%iMix)%vXPole(J)
        end do
      end if
      write(fo,*)

    end if

  end do
  !
  close(fo)
  !
  return
end subroutine Equil_SaveFas_ToFile

subroutine Equil_SaveSpc_ToFile(vSpcDat)
  use M_IoTools,  only: GetUnit
  use M_T_Species,only: T_SpcData
  !
  type(T_SpcData), intent(in) :: vSpcDat(:)
  !
  integer:: F,i
  !
  call GetUnit(F)
  !open(F,file="tmp_species_mole.dat",form="unformatted",access="stream")
  !write(F) ( vSpcDat(i)%Mole, i=1,size(vSpcDat) )
  open(F,file="tmp_species_mole.dat")
  write(F,'(*(G15.6,1X))') (vSpcDat(i)%Mole, i=1,size(vSpcDat) )
  close(F)
  call GetUnit(F)
  !open(F,file="tmp_species_lact.dat",form="unformatted",access="stream")
  !write(F) ( vSpcDat(i)%LAct, i=1,size(vSpcDat) )
  open(F,file="tmp_species_lact.dat")
  write(F,'(*(G15.6,1X))') ( vSpcDat(i)%LAct, i=1,size(vSpcDat) )
  close(F)
  call GetUnit(F)
  !open(F,file="tmp_species_lgam.dat",form="unformatted",access="stream")
  !write(F) ( vSpcDat(i)%LGam, i=1,size(vSpcDat) )
  open(F,file="tmp_species_lgam.dat")
  write(F,'(*(G15.6,1X))') ( vSpcDat(i)%LGam, i=1,size(vSpcDat) )
  close(F)
  !
  return
end subroutine Equil_SaveSpc_ToFile

subroutine Equil_SaveCpn_ToFile(nAq,tStoikio,vSpcDat)
  use M_IoTools,  only: GetUnit
  use M_T_Species,only: T_SpcData
  !
  integer,         intent(in):: nAq
  real(dp),        intent(in):: tStoikio(:,:)
  type(T_SpcData), intent(in):: vSpcDat(:)
  !
  integer :: F,iPr
  real(dp):: X
  !
  call GetUnit(F)
  ! open(F,file="tmp_cpn_fluid.dat",form="unformatted",access="stream")
  open(F,file="tmp_cpn_fluid.dat")
  do iPr=1,size(tStoikio,1)
    X= sum(tStoikio(iPr,1:nAq)*vSpcDat(1:nAq)%Mole) !-> mole nrs in fluid
    ! write(F) X
    write(F,'(G15.6,1X)') X
  end do
  close(F)
  !
  return
end subroutine Equil_SaveCpn_ToFile

subroutine Equil_Read_FromFile(vSpcDat,OK)
  use M_T_Species,only: T_SpcData
  !---------------------------------------------------------------------
  type(T_SpcData), intent(inout):: vSpcDat(:)
  logical,         intent(out)  :: OK
  !---------------------------------------------------------------------
  real(dp):: v(size(vSpcDat))
  !
  OK= .false.
  call Read_FromFile("tmp_species_mole.dat",size(vSpcDat),v,OK)
  if(OK) then
    vSpcDat(:)%Mole= v(:)
  else
    return
  end if
  call Read_FromFile("tmp_species_lact.dat",size(vSpcDat),v,OK)
  if(OK) then
    vSpcDat(:)%LAct= v(:)
  else
    return
  end if
  call Read_FromFile("tmp_species_lgam.dat",size(vSpcDat),v,OK)
  if(OK) then
    vSpcDat(:)%LGam= v(:)
  else
    return
  end if
  !
  return
end subroutine  Equil_Read_FromFile

subroutine Read_FromFile(str,N,V,OK)
  use M_IoTools,  only: GetUnit, LineToArray 
  use M_FileUtils,only: File_Exist 
  !
  character(len=*),intent(in) :: str
  integer,         intent(in) :: N
  real(dp),        intent(out):: V(N)
  logical,         intent(out):: OK
  !
  character(len=512):: Line
  integer:: F,i
  integer:: ios
  !
  OK= .false.
  v(:)= 0.D0
  !
  if(File_Exist(str)) then
    call GetUnit(F)
    !open(F,file=str,form="unformatted",access="stream")
    !do i=1,size(v)
    !  read(F) v(i)
    !end do
    open(F,file=str)
    read(F, '(A)', iostat=ios) Line
    if(ios/=0) return
    call LineToArray(trim(Line),N,V)
    close(F)
    OK= .true.
  end if
  !
  return
end subroutine Read_FromFile

subroutine Save_EquPhase
  use M_T_MixModel, only: MaxPole
  use M_T_Phase,    only: T_Phase
  use M_Equil_Vars, only: T_EquPhase
  !
  use M_Global_Vars,only: vMixModel,vMixFas,vFas
  use M_Basis_Vars, only: tAlfFs,tNuFas
  use M_Equil_Vars, only: nEquFas,vEquFas,cEquMode
  !---------------------------------------------------------------------
  type(T_Phase),allocatable:: vFas0(:)
  real(dp),     allocatable:: tAlf(:,:),tNu(:,:)
  integer:: vIPole(MaxPole) ! indexes of end-members in vFas
  integer:: nPur,nMix
  integer:: I,K,nP,C,nCp
  
  nPur= count(vFas(:)%iSpc>0)
  nCp=  size(tAlfFs,1)
  
  !--------------------------------------------------------build vMixFas
  K=0
  do I=1,nEquFas
    if(vEquFas(I)%iMix>0) K= K+1
  end do
  deallocate(vMixFas)
  allocate(vMixFas(K))
  !
  K=0
  if(size(vMixFas)>0) then
    do I=1,nEquFas
      if(vEquFas(I)%iMix>0) then
        !
        K= K+1
        nP= vMixModel(vEquFas(I)%iMix)%nPole
        vMixFas(K)%Name= trim(vEquFas(I)%NamEq)
        vMixFas(K)%iModel= 1 !vEquFas(I)%iMix
        vMixFas(K)%vXPole(1:nP)= vEquFas(I)%vXPole(1:nP)
        !
      end if
    end do
  end if
  !-------------------------------------------------------/build vMixFas
  !
  !---------------------------------------------------------rebuild vFas
  nMix= count(vEquFas(1:nEquFas)%iMix>0)
  if(nMix>0) then
    !
    ! temporary copy
    allocate(vFas0(nPur))
    vFas0(1:nPur)= vFas(1:nPur)
    !
    deallocate(vFas)
    allocate(vFas(nPur+nMix))
    vFas(1:nPur)= vFas0(1:nPur)
    !
    deallocate(vFas0)
    !
    allocate(tAlf(nCp,nPur))        ; tAlf(:,1:nPur)= tAlfFs(:,1:nPur)
    deallocate(tAlfFs)
    allocate(tAlfFs(nCp,nPur+nMix)) ; tAlfFs(:,1:nPur)= tAlf(:,1:nPur)
    deallocate(tAlf)
    !
    allocate(tNu(nPur,nCp))         ; tNu(1:nPur,:)= tNuFas(1:nPur,:)
    deallocate(tNuFas)
    allocate(tNuFas(nPur+nMix,nCp)) ; tNuFas(1:nPur,:)= tNu(1:nPur,:)
    deallocate(tNu)
    !
  end if
  !--------------------------------------------------------/rebuild vFas
  !
  vFas(:)%MolFs= Zero
  !
  K=nPur
  do I=1,nEquFas
    !
    if(vEquFas(I)%iPur>0) vFas(vEquFas(I)%iPur)%MolFs= vEquFas(I)%MolFs
    !
    if(vEquFas(I)%iMix>0) then
      !
      K= K+1
      !
      vFas(K)%MolFs= vEquFas(I)%MolFs
      vFas(K)%iSpc=  0
      vFas(K)%iSol=  0
      vFas(K)%iMix=  K -nPur
      vFas(K)%NamFs= trim(vEquFas(I)%NamEq)
      !! vFas(K)%Typ=   "MIXT"
      !
      nP= vEquFas(I)%NPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      do C=1,nCp
        tAlfFs(C,K)= sum( vEquFas(I)%vXPole(1:nP) &
        &          * tAlfFs(C,vIPole(1:nP)) )
        tNuFas(K,C)= sum( vEquFas(I)%vXPole(1:nP) &
        &          * tNuFas(vIPole(1:nP),C) )
      end do
      !
    end if
    !
  end do
  !
  return
end subroutine Save_EquPhase

subroutine Equil_CalcFasAff
!--
!-- calculate affinity of pure phases (i.e. non-aqu'species)
!--
  use M_Global_Vars,only: vSpc,vFas,vSpcDat
  use M_Basis_Vars, only: tNuFas,vOrdPr
  use M_Equil_Vars, only: vFasAff
  !---------------------------------------------------------------------
  integer  :: iFs,nFs,nCp
  
  nFs= size(vFas)
  nCp= size(vOrdPr)
  
  if (nFs>0) then
    !note:
    !vDeltaG_Ms(iMs)=  vSpc(vOrdMs(iMs))%G0rt 
    !               - dot_product(tNuMs(iMs,1:nCp),Spc(vOrdPr(1:nCp))%G0rt) 
    !
    !-> deltaG of PRECIPITATION reaction, Prim'Species -> Mineral
    !
    !affinity/RT of formatION reaction of mineral,
    != A/RT= -DotProd(Nu,Mu)= -deltaG of Prim'Sp->Mineral, A/RT=ln(K/Q)
    !
    !LnQsK(iMs)= dot_product( tNuMs(iMs,1:nCp), vSpc(vOrdPr(1:nCp))%LnAct ) 
    !          - vDeltaG_Ms(iMs)
    !
    do iFs=1,nFs
      vFasAff(iFs)= &
      & vFas(iFs)%Grt &
      - dot_product( tNuFas(iFs,1:nCp), &
      &           vSpcDat(vOrdPr(1:nCp))%LAct+vSpc(vOrdPr(1:nCp))%G0rt )
      !vFasAff(iFs)= vFasAff(iFs) /sum(abs(tNuFas(iFs,1:nCp)))
    end do !_______________________________end
  end if
  
  !expression used in Equil_1:
  !
  !  vFasAff(iFs)= vDeltaG_Fs(iFs) &
  !  &           - dot_product(tNuFas(iFs,1:nCp),vLnAct(vOrdPr(1:nCp)))
  !
  !with
  !
  !  vDeltaG_Fs(:)= vFas(:)%Grt 
  !               - matmul(tNuFas(:,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
  !
  !-> can also calculate vFasAff as
  !
  !  vFasAff(iFs)= &
  !  & vSpc(vFas(iFs)%iSpc)%G0rt &
  !  - dot_product(tNuFas(iFs,1:nCp), &
  !  &             vSpcDat(vOrdPr(1:nCp))%LAct+vSpc(vOrdPr(1:nCp))%G0rt)
  
  return
end subroutine Equil_CalcFasAff

subroutine Equil_Trace_Close
!--
!-- close special trace files for activities, Newton, etc
!--
  use M_Equil_Vars, only: fTrAct, fActiz, fTrcEq
  use M_Numeric_Tools,only: fNewtF,fNewtR
  !
  if(fTrAct>0)  then; close(fTrAct); fTrAct= 0; end if
  if(fActiZ>0)  then; close(fActiZ); fActiZ= 0; end if
  if(fTrcEq>0)  then; close(fTrcEq); fTrcEq= 0; end if
  !
  if(fNewtF>0)  then; close(fNewtF); fNewtF= 0; end if
  if(fNewtR>0)  then; close(fNewtR); fNewtR= 0; end if
  !
  return
end subroutine Equil_Trace_Close

subroutine Equil_Trace_Init
!--
!-- open special trace files for activities, Newton, etc
!--
  use M_IOTools,    only: GetUnit
  use M_Files,      only: DirOut,Files_Index_Write
  use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
  !
  use M_Global_Vars,only: vSpc
  use M_Equil_Write,only: Equil_Write_EnTete
  use M_Basis_Vars, only: vLCi,vLAs,vLAx,vPrmFw
  use M_Equil_Vars, only: DebActiv,fTrAct,fActiZ,DebNewt
  !---------------------------------------------------------------------
  integer:: I
  
  if(DebActiv) then
    !
    call GetUnit(fTrAct)
    open(fTrAct,file=trim(DirOut)//"_activ_deb.log")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_activ_deb.log",&
    & "EQUIL/LOG: activ' aqu'species")
    write(fTrAct,'(A,A1)',advance="no") "indx",T_
    do I=1,size(vSpc)
      if(vSpc(I)%Typ=="AQU") write(fTrAct,'(A,A1)',advance="no") trim(vSpc(I)%NamSp),T_
    end do
    write(fTrAct,*)
    !
    call GetUnit(fActiz); open(fActiz,file=trim(DirOut)//"_activ.log")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_activ.log",&
    & "EQUIL/LOG: log10(Activ)")
    call Equil_Write_EnTete(fActiz,vSpc,vPrmFw,vLCi.or.vLAs.or.vLAx,Str1="indx")
    !
  else
    !
    fTrAct= 0
    fActiZ= 0
    !
  end if
  
  if(DebNewt) then
    !
    call GetUnit(fNewtF)
    open(fNewtF,file=trim(DirOut)//"_newtequspc1.log")
    call Equil_Write_EnTete(fNewtF,vSpc,vPrmFw,vLCi.or.vLAs,Str1="count1",Str2="count2")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_newtequspc1.log",&
    & "EQUIL/LOG: Newton solution")
    !
    call GetUnit(fNewtR)
    open(fNewtR,file=trim(DirOut)//"_newtequspc2.log")
    call Equil_Write_EnTete(fNewtR,vSpc,vPrmFw,vLCi.or.vLAs,Str1="count1",Str2="count2")
    call Files_Index_Write(fHtm,&
    & trim(DirOut)//"_newtequspc2.log",&
    & "EQUIL/LOG: Newton residue")
    !
    fNewt_I=0
    !
  else
    !
    fNewtF=0
    fNewtR=0
    !
  end if
  
  return
end subroutine Equil_Trace_Init

subroutine Equil_FasTrace_EnTete( &
& vFas,vYesList, &
& F)
  use M_IOTools,only: GetUnit
  use M_Files,  only: DirOut
  use M_T_Phase,only: T_Phase
  !---------------------------------------------------------------------
  type(T_Phase),intent(in) :: vFas(:)
  logical,      intent(in) :: vYesList(:)
  integer,      intent(out):: F
  !---------------------------------------------------------------------
  integer:: iFs,nFs
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_equil.log")
  !
  write(F,'(3(A,A1))',advance="NO") "iCount",T_,"nIts",T_,"dXi",T_
  !
  nFs=size(vFas)
  !
  do iFs=1,nFs
    if(vYesList(iFs)) &
    & write(F,'(A,A1)',advance="NO") "Mol_"//trim(vFas(iFs)%NamFs),T_
  end do
  do iFs=1,nFs
    if(vYesList(iFs)) &
    & write(F,'(A,A1)',advance="NO") "lQsK_"//trim(vFas(iFs)%NamFs),T_
  end do
  write(F,*)
  !
  return
end subroutine Equil_FasTrace_EnTete

subroutine Equil_FasTrace_Write(F,ICount,nIts,dXi,vYesList,vFasMole,vFasAff)
  use M_Numeric_Const,  only: Ln10
  use M_IOTools,only: OutStrVec
  !---------------------------------------------------------------------
  integer, intent(in):: F
  integer, intent(in):: iCount
  integer, intent(in):: nIts
  real(dp),intent(in):: dXi
  logical, intent(in):: vYesList(:)
  real(dp),intent(in):: vFasMole(:),vFasAff(:)
  !---------------------------------------------------------------------
  write(F,'(2(I7,A1),G15.6,A1)',advance="NO") iCount,T_,nIts,T_,dXi,T_
  call OutStrVec(F,vFasMole,     vL=vYesList,CR=.false.)
  call OutStrVec(F,-vFasAff/Ln10,vL=vYesList)
  !
  return
end subroutine Equil_FasTrace_Write

end module M_Equil_Tools

