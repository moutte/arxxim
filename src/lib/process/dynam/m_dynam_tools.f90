module M_Dynam_Tools
  use M_Kinds
  use M_Trace,only: iDebug,Stop_,fTrc,fHtm,T_,Pause_
  
  implicit none
  
  private
  !
  public:: Dynam_Init
  public:: Dynam_Alloc
  public:: Dynam_TP_Update
  public:: Dynam_TimeFactor
  public:: Dynam_Zero_Box
  public:: Dynam_Zero_Time
  public:: Dynam_Zero_Numeric
  public:: Dynam_ReStart_Time
  public:: Dynam_ReStart_Box
  public:: Dynam_Save
  public:: Dynam_Save_ToFile
  public:: Dynam_Clean
  
contains

subroutine Dynam_Alloc
  use M_Global_Vars,only: vSpc,vKinFas,vKinModel
  use M_Dynam_Vars, only: vCpnBox
  use M_Dynam_Vars, only: vWeitSp,vWeitCp,vKinMod
  !
  integer :: nCp,nAq,nMk
  !
  nAq=count(vSpc%Typ=="AQU")
  if(nAq==0) call Stop_("NO Aqueous Species ???")
  !
  nCp= size(vCpnBox)
  nMk= size(vKinFas)
  !
  call Dynam_Alloc_FluidBox(nCp)
  !-> allocate vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  call Dynam_Alloc_FluidInj(nCp)
  !-> allocate vTotInj,vQTotInj
  call Dynam_Alloc_KinFas(nMk)
  !-> allocate vMolK,vMolarVol,vKinQsK,vKinMinim,vSurfK0,...
  call Dynam_KinFas_Stoik_Alloc(nMk,nCp)
  !-> allocate vDG_Kin,tNu_Kin,tAlfKin
  !
  call Dynam_TP_Alloc
  !
  allocate(vKinMod(1:size(vKinModel)))
  !
  allocate(vWeitSp(1:nAq))
  allocate(vWeitCp(1:nCp))
  !
end subroutine Dynam_Alloc

subroutine Dynam_TP_Update(TdgK,Pbar,TimeFactor)
  use M_Dynam_Vars,only: tNu_Kin,vDG_Sp,vDG_Kin,vKinMod
  !
  real(dp),intent(in):: TdgK,Pbar,TimeFactor
  !
  call Dynam_Global_TP_Update(TdgK,Pbar)
  !
  call Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)
  !
  call KinModel_T_Update(TimeFactor,TdgK,vKinMod)
  !
end subroutine Dynam_TP_Update

subroutine Dynam_Init(TdgK,Pbar,T,B)
!--
!-- initializations of variables used in "dynamic" calculations
!--
  use M_IOTools
  use M_Global_Vars,only: vSpc,SolModel,vMixModel,vEle
  use M_Global_Vars,only: vFas,vKinFas,vMixFas,vSpcDat
  use M_Basis_Vars, only: vOrdAq

  use M_Dynam_Vars,only: vCpnBox,vCpnInj
  use M_Dynam_Vars,only: T_DynTime,T_DynBox
  use M_Dynam_Vars,only: vWeitSp,vWeitCp
  use M_Dynam_Vars,only: vDG_Sp,vDG_Kin,tNu_Kin
  use M_Dynam_Vars,only: TdgK0,Pbar0
  use M_Dynam_Vars,only: vTotInj,vTotF,vMolF,vLnAct,vLnGam,vLnBuf

  !! use M_Stockvar_Kinxim,only: LSTOCK,INIT_STOCKVAR
  !---------------------------------------------------------------------
  real(dp),       intent(in)   :: TdgK,Pbar
  type(T_DynTime),intent(inout):: T
  type(T_DynBox), intent(inout):: B
  !---------------------------------------------------------------------
  integer, parameter:: MaxStockVar = 40000
  integer :: nCp,nAq,nMk
  !---------------------------------------------------------------------
  !real(dp):: &
  !& MFluid, &    !mass of fluid corresponding to a given composition
  !& MFluidBox    !mass of fluid, of known density RhoF, that can occupy PhiF*VBox
  !MFluid=    dot_product(vSpc(1:nAq)%WeitKg,vMolF(1:nAq))
  !MFluidBox= MFluidBox= RhoF*VBox*PhiF
  
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_Init"
  !
  nCp=size(vCpnBox)
  nMk=size(vKinFas)            ! ; if(nMk==0) call Stop_("NO Kinetic Species ???")
  nAq=count(vSpc%Typ=="AQU")   ; if(nAq==0) call Stop_("NO Aqueous Species ???")
  !
  !---------------- init table to store the time evolution of species --
  !! if (LSTOCK) call  init_stockvar(MaxStockVar, nAq)
  !
  !--------------------------------------------------update fluid data--
  !B%PhiFInj= B%PhiF !initial vol. fraction of fluid
  B%RhoF= 1.0D3 !! Solvent%Rho(iTP)*1.0D3 !fluid density
  if(idebug>1) write(fTrc,'(/,A,G15.3,/)') "RhoF",B%RhoF
  !
  !---------------------------------------------------update time data--
  T%Time= Zero
  !
  if(iDebug>2) call Dynam_Init_StoikioEle
  !
  call Dynam_Init_FluidBox_Compo( &   !
  & vCpnBox,vSpcDat,B%RhoF,B%VBox,B%PhiF, &   !IN
  & vTotF,vMolF,vLnAct,vLnGam,vLnBuf) !OUT
  !
  call Dynam_Init_FluidInj_Compo( & !
  & vCpnInj,B%RhoF,B%VBox, &  !IN
  & vTotInj)                  !OUT
  !
  call Dynam_Init_KinFas_Texture(nMk,B%VBox,B%PhiF) !IN
  !-> initialize/update vKinFas,vMolK,SurfMinim,vMolarVol,vKinMinim,vSurfK0
  !
  !call Dynam_KinFas_Stoik_Init(nMk,nCp)
  call Dynam_KinFas_Stoik_Alfa(nMk,nCp)
  !
  call Dynam_KinFas_Stoik_NuKin(nMk)
  if(iDebug>2) call NuTable_Sho(vCpnBox)
  !
  TdgK0= TdgK
  Pbar0= Pbar
  !
  call Dynam_Global_TP_Update(TdgK,Pbar)
  !
  call Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)
  !
  call Dynam_KinModel_Init(T%TimeFactor,TdgK) !IN
  !
  call Dynam_Init_KinFas_SatState
  !
  vWeitSp(:)= vSpc(vOrdAq(:))%WeitKg
  !-> Aqu'Species' weitkg sorted according to local Species' order
  vWeitCp(:)= vEle(vCpnBox(:)%iEle)%WeitKg
  !-> Components' weitkg, i.e. Elements' weitkg in the Components' order ...
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_Init"
  
  return
end subroutine Dynam_Init

subroutine Dynam_Init_StoikioEle
  use M_Global_Vars,only: vSpc,vEle,vKinFas,vFas,tFormula
  use M_Basis_Vars, only: vOrdAq
  use M_Dynam_Vars, only: tStoikioAqu,tStoikioKin
  !
  integer:: I,J,N,nMk
  !
  N= size(vEle)
  allocate(tStoikioAqu(size(vOrdAq),N))
  do I=1,size(vOrdAq)
    J= vOrdAq(I)
    tStoikioAqu(I,:)= tFormula(:,J)
    !vSpc(J)%vStoikio(1:N)/real(vSpc(J)%vStoikio(0))
    !! do K=1,N
    !!   write(6,'(F7.2,1X)',advance="no") tStoikioAqu(I,K)
    !! end do
    !! write(6,*)
  end do
  nMk= size(vKinFas)
  allocate(tStoikioKin(nMk,N))
  do I=1,nMk
    J= vFas(vKinFas(I)%iFas)%iSpc
    tStoikioKin(I,:)= tFormula(:,J)
    !vSpc(J)%vStoikio(1:N)/real(vSpc(J)%vStoikio(0))
    !! do K=1,N
    !!   write(6,'(F7.2,1X)',advance="no") tStoikioKin(I,K)
    !! end do
    !! write(6,*)
  end do
  !! pause
  
  return
end subroutine Dynam_Init_StoikioEle

!--------------------------------------------------allocations for fluid
subroutine Dynam_Alloc_FluidBox(nCp)
!=
!= allocate vTotF,vMolF,vLnAct,vLnGam,vLnBuf
!=
  use M_Global_Vars,only: vSpc,vEle,vSpcDat
  !! use M_Basis_Vars, only: nMx,nAx
  use M_Dynam_Vars, only: vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  !
  integer,intent(in):: nCp
  !
  integer:: nAq
  !
  allocate(vTotF(1:nCp))
  !
  nAq= count(vSpc%Typ=="AQU")
  !
  allocate(vMolF(1:nAq))
  allocate(vLnAct(1:nAq))
  allocate(vLnGam(1:nAq))
  !
  !! allocate(vLnBuf(1:nAx+nMx))
  allocate(vLnBuf(1:nAq))
  !
  return
end subroutine Dynam_Alloc_FluidBox

subroutine Dynam_Alloc_FluidInj(nCp)
  use M_Dynam_Vars, only: vTotInj,vQTotInj !,vTotInjFile
  !
  integer,intent(in):: nCp
  !
  allocate(vTotInj(1:nCp))       ; vTotInj= Zero
  allocate(vQTotInj(1:nCp))      ; vQTotInj= Zero
  !
end subroutine Dynam_Alloc_FluidInj
!-------------------------------------------------/allocations for fluid

!-----------------------------------------------------initialize vTotInj
!-- initialize mole nrs injected fluid
!-- units: mole nrs in a volume VBox for a fluid density RhoF
subroutine Dynam_Init_FluidInj_Compo( &
& vCpn,RhoF,VBox, &
& vTotInj)
  use M_T_Component,only: T_Component
  use M_Global_Vars,only: vEle !(for %WeitKg)
  !
  type(T_Component),intent(in)::  vCpn(:)
  real(dp),         intent(in)::  RhoF,VBox
  real(dp),         intent(out):: vTotInj(:)
  !
  integer:: I
  !
  vTotInj(:)= vCpn(:)%Mole
  vTotInj(:)= vTotInj(:) &
  &           /dot_product(vEle(vCpn(:)%iEle)%WeitKg,vTotInj(:))  & !-> nr.mole/Kg
  &           *RhoF & !-> nr.mole/M3
  &           *VBox   !-> in mole nrs in the box
  !
  if(iDebug>2) then
    print '(A)',"Dynam_Init_FluidInj_Compo"
    do I=1,size(vTotInj)
      print '(A,G15.6)',vEle(vCpn(I)%iEle)%NamEl,vTotInj(I)
    end do
  end if
  !
end subroutine Dynam_Init_FluidInj_Compo

!----------------------------------------------initialize fluid in box--
subroutine Dynam_Init_FluidBox_Compo( &       !
& vCpn,vSpcDat,RhoF,VBox,PhiF, &                      !
& vTotF,vMolF,vLnAct,vLnGam,vLnBuf)
  !
  use M_T_Component,only: T_Component
  use M_T_Species,  only: T_SpcData
  use M_Global_Vars,only: vSpc,vEle
  use M_Basis_Vars, only: nAx,nMx,vOrdAq,nCi
  !
  type(T_Component),intent(in) :: vCpn(:)
  type(T_SpcData),  intent(in) :: vSpcDat(:)
  real(dp),         intent(in) :: RhoF,VBox,PhiF
  real(dp),         intent(out):: vTotF(:),vMolF(:)
  real(dp),         intent(out):: vLnAct(:),vLnGam(:),vLnBuf(:)
  !
  integer:: nAq,I
  
  nAq= size(vOrdAq)
  !
  !-------------------------------------retrieve results from speciation
  vMolF(:)=      vSpcDat(vOrdAq(:))%Mole
  vLnAct(1:nAq)= vSpcDat(vOrdAq(:))%LAct
  vLnGam(1:nAq)= vSpcDat(vOrdAq(:))%LGam
  !
  ! given mole numbers of aqu.species for 1kg H2O
  ! (=results of initial speciation)
  ! calculate mole numbers of aqu.species in a box of volume VBox,
  ! assuming the fluid has a known density RhoF
  !
  vMolF(:)= vMolF(:) &
  & /dot_product(vMolF(1:nAq), vSpc(vOrdAq(:))%WeitKg) *RhoF & !-> mole/m3
  & *VBox*PhiF !-> mole in fluid in box
  !
  vTotF(:)= vCpn(:)%Mole
  !
  !--moles of Components in fluid in the box--
  vTotF(:)= vTotF(:) &
  & /dot_product(vMolF(:), vSpc(vOrdAq(:))%WeitKg) *RhoF & !-> mole/m3
  & *VBox*PhiF
  !
  !--activities of buffered species--
  vLnBuf(:)= 0.D0
  !! print *,"nAq, size(vLnBuf)",nAq, size(vLnBuf)  ; pause
  if(nAx+nMx>0) then
    do i=1,nAx+nMx
      vLnBuf(I)= vSpcDat(vCpn(nCi+I)%iSpc)%LAct
    end do
  end if
  !
  if(iDebug>2) then
    print '(A)',"Dynam_Init_FluidBox_Compo"
    do I=1,size(vCpn)
      if(vCpn(I)%Statut=="BUFFER") &
      & print '(A,A3,G15.6)', "BUFFER= ",vCpn(I)%NamCp,vSpcDat(vCpn(I)%iSpc)%LAct
      if(vCpn(I)%Statut=="INERT") &
      & print '(A,A3,G15.6)', "INERT=  ",vCpn(I)%NamCp,vCpn(I)%Mole
    end do
  end if
  
  return
end subroutine Dynam_Init_FluidBox_Compo

!--------------------------------allocations for kinetic & equil' phases
subroutine Dynam_Alloc_KinFas(N)
!--
!-- allocate tables related with kinetic minerals
!--
  use M_Dynam_Vars, only: &
  & vMolK,vSurfK,vMolarVol,vKinQsK,vKinMinim,vStatusK,&
  & vSurfK0,vMolK0,vVmQsK,vVmAct, &
  & vLKinActiv,vLEquActiv,vKinPrm
  !
  integer, intent(in):: N
  !
  allocate(vMolK(1:N))         ; vMolK=  Zero   !moles kin.species in the box
  allocate(vSurfK(1:N))        ; vSurfK= Zero   !surface kin.speices in the box
  allocate(vKinMinim(1:N))     ; vKinMinim=Zero !minimal mole numbers
  allocate(vStatusK(1:N))
  !
  allocate(vKinQsK(1:N))       ; vKinQsK= One !
  !
  allocate(vMolarVol(1:N))     ; vMolarVol= Zero !molar volumes
  !
  allocate(vMolK0(1:N))        ; vMolK0=  Zero !mole numbers at ref'time (for kin'rate)
  allocate(vSurfK0(1:N))       ; vSurfK0= Zero !surface at ref'time (for "scaled" output and for CRUNCH mode)
  !
  allocate(vVmQsK(1:N))        ; vVmQsK=  Zero !Precip /Dissol. Rate, VmF=VmAct*VmQsK
  allocate(vVmAct(1:N))        ; vVmAct=  Zero !Precip /Dissol. Rate, VmF=VmAct*VmQsK
  !
  allocate(vLKinActiv(1:N))    ; vLKinActiv=.false.
  allocate(vLEquActiv(1:N))    ; vLEquActiv=.false.
  allocate(vKinPrm   (1:N))    ; vKinPrm=   0
  !
end subroutine Dynam_Alloc_KinFas

subroutine Dynam_Init_KinFas_Texture(nMk,VBox,PhiF)
  use M_Numeric_Const,only: Pi
  use M_KinFas_Surf,  only: KinFas_Surf_Init
  !
  use M_Global_Vars,  only: vFas,vKinFas
  !
  use M_Dynam_Vars,   only:    &
  & PhiF0,vKinMinim,           &
  & vMolK,vSurfK,vMolarVol,    &
  & vMolK0,vSurfK0,            &
  & RadiusMinim,DensitMinim
  !
  integer, intent(in):: nMk
  real(dp),intent(in):: VBox,PhiF
  !
  real(dp):: VolMinim,SurfMinim
  integer :: I
  !
  !--------------------------------------------------------INIT.MINERALS
  !
  vMolarVol(1:nMk)= vFas(vKinFas(1:nMk)%iFas)%VolM3 !-> molar volumes of phases
  !
  vKinFas(:)%Dat%PhiM= (One - PhiF) *vKinFas(:)%Dat%PhiM !-> vol.fraction of box
  !
  !-------------------------------------------------------init vKinMinim
  !
  !! !test calculation radius of "cluster" of 6 molecules -> results around 4.D-9 m
  !! vKinMinim(1:nMk)= 36.D-23 != approx. 6 molecules
  !! do i=1,nMk
  !!   VolMinim=    vKinMinim(I) * vFas(vKinFas(I)%iFas)%V !-> volume of six molecules
  !!   RadiusMinim= ( 3.D0*VolMinim /4.D0/Pi) **(1.D0/3.D0)
  !!   write(*,'(A15,1X,2G15.6,1X,A)') vKinFas(i)%Name, VolMinim, RadiusMinim, vKinFas(I)%cSatur
  !! end do
  !! call Pause_
  !
  SurfMinim= DensitMinim *Pi*4._dp       *RadiusMinim**2  ! per m3
  VolMinim=  DensitMinim *Pi*4._dp/3._dp *RadiusMinim**3  ! per m3
  !
  if(iDebug>2) then
    print '(A,G15.6)', "SurfMinim=",SurfMinim
    print '(A,G15.6)', "VolMinim =",VolMinim
    call Pause_
  end if
  !
  !--- mimimal mole number in the box
  vKinMinim(1:nMk)= VolMinim*Vbox /vMolarVol(1:nMk)
  !
  if(iDebug>2) then
    do i=1,nMk
      print '(A15,2A,2(A,G15.6))', &
      & vKinFas(i)%NamKF," Sat=",vKinFas(I)%Dat%cSat,&
      & " PhiM=",vKinFas(I)%Dat%PhiM," Minim=",vKinMinim(I)
    end do
    call Pause_ !; stop
  end if
  !
  if(idebug>1) write(fTrc,'(A)') "Initial Surfaces of Minerals"
  !
  do i=1,nMk
    !
    call KinFas_Surf_Init( &
    & vFas, &
    & vKinFas(i),                & !in
    & vKinFas(i)%Dat%PhiM,       & !in,  volume_phase / volume_box
    & vKinFas(i)%Dat%Surf)         !out, surface_phase / volume_box
    !
    !! vKinFas(i)%Dat%Surf= max(vKinFas(i)%Dat%Surf, SurfMinim)
    !
  end do
  !
  !!where(trim(vKinFas%cSatur)=="2") vMolK(:)=vKinMinim(:)  !-> mole number minerals in box
  !!where(trim(vKinFas%cSatur)=="1") vMolK(:)=VBox* vKinFas(:)%PhiM /vMolarVol(:)
  !!where(vKinFas%Dat%cSat=="I")     vMolK(:)= VBox *vKinFas(:)%Dat%PhiM /vMolarVol(:)
  !
  vMolK(:)=  vKinFas(:)%Dat%PhiM *Vbox/vMolarVol(:)
  vSurfK(:)= vKinFas(:)%Dat%Surf *vBox
  !
  vSurfK0(:)= vSurfK(:) ! update vSurfK0
  PhiF0=      PhiF      ! update PhiF0
  vMolK0(:)=  vMolK(:)  ! update vMolK0
  !
  if(idebug>1) then
    write(fTrc,'(A,G15.6)') "VBox     ", VBox
    write(fTrc,'(/,A,/)') "vSurfK0=Initial Surfaces of Minerals"
    do I=1,nMk
      write(fTrc,'(A15,1X,3(A,1X,G12.5,1X))') &
      & vKinFas(i)%NamKF,&
      & "vSurfK0=",vSurfK0(i),&
      & "/Rho=",vFas(vKinFas(I)%iFas)%WeitKg /vMolarVol(I),& !mineral density
      & "/vMolK0=",vMolK0(I)
    end do
  end if
  !-------------------------------------------------------/INIT.MINERALS
  !
  return
end subroutine Dynam_Init_KinFas_Texture

subroutine Dynam_KinFas_Stoik_Alloc(nMk,nCp)
  use M_Dynam_Vars, only: tNu_Kin,tAlfKin,vDG_Kin
  !
  integer,intent(in):: nMk,nCp
  !
  allocate(vDG_Kin(1:nMk))
  allocate(tNu_Kin(1:nMk,1:nCp))
  allocate(tAlfKin(1:nCp,1:nMk))
  !
end subroutine Dynam_KinFas_Stoik_Alloc

subroutine Dynam_KinFas_Stoik_Alfa(nMk,nCp)
!--
!-- compute tAlfKin, (and tAlfPr, tAlfAs for orthogonal basis)
!--
  use M_T_MixModel, only: T_MixModel
  use M_T_MixPhase, only: T_MixPhase
  !
  use M_Global_Vars,only: vMixModel,vMixFas,vFas,vKinFas
  use M_Basis_Vars, only: tAlfSp,tAlfPr,tAlfAs
  use M_Dynam_Vars, only: tNu_Kin,tAlfKin
  !
  integer,intent(in):: nMk,nCp
  !
  type(T_MixModel):: SM
  type(T_MixPhase):: SF
  integer:: I,J,K,iCp
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_KinFas_Stoik_Alfa"
  !
  !-----------------------------------------------update stoichio tables
  do I=1,nMk
  
    J=vKinFas(I)%iFas
    
    if(vFas(J)%iSpc>0) then
      !
      K=  vFas(J)%iSpc
      tAlfKin(1:nCp,I)= tAlfSp(1:nCp,K)
      !
    else if(vFas(J)%iMix>0) then
      !
      K=  vFas(J)%iMix
      SF= vMixFas(K)
      SM= vMixModel(SF%iModel)
      !
      do iCp=1,nCp
        tAlfKin(iCp,I)= &
        & dot_product( tAlfSp(iCp,SM%vIPole(1:SM%NPole)) , SF%vXPole(1:SM%NPole) )
      end do
      !
    end if

    !! select case(trim(vFas(J)%Typ))
    !! case("PURE") !,"DISCRET")
    !!   K=  vFas(J)%iSpc
    !!   tAlfKin(1:nCp,I)= tAlfSp(1:nCp,K)
    !! case("MIXT")
    !!   K=  vFas(J)%iMix
    !!   SF= vMixFas(K)
    !!   SM= vMixModel(SF%iModel)
    !!   !
    !!   do iCp=1,nCp
    !!     tAlfKin(iCp,I)= &
    !!     & dot_product( tAlfSp(iCp,SM%vIPole(1:SM%NPole)) , SF%vXPole(1:SM%NPole) )
    !!   end do
    !! end select
    
  end do
  !----------------------------------------------/update stoichio tables
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_KinFas_Stoik_Alfa"
  !
end subroutine Dynam_KinFas_Stoik_Alfa

subroutine Dynam_KinFas_Stoik_NuKin(nMk)
!--
!-- compute tNu_Kin
!--
  use M_T_MixModel, only: T_MixModel
  use M_T_MixPhase, only: T_MixPhase
  !
  use M_Global_Vars,only: vMixModel,vMixFas,vFas,vKinFas
  use M_Basis_Vars, only: tNuSp
  use M_Dynam_Vars, only:tNu_Kin
  !
  integer,intent(in):: nMk
  !
  type(T_MixModel):: SM
  type(T_MixPhase):: SF
  integer:: I,J,K,iCp
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_KinFas_Stoik_Nu"
  !
  !-----------------------------------------------update stoichio tables
  do I=1,nMk
  
    J=vKinFas(I)%iFas

    if(vFas(J)%iSpc>0) then
      K=  vFas(J)%iSpc
      tNu_Kin(I,:)= tNuSp(K,:)
      !
    else if(vFas(J)%iMix>0) then
      K=  vFas(J)%iMix
      SF= vMixFas(K)
      SM= vMixModel(SF%iModel)
      !
      do iCp=1,size(tNu_Kin,2)
        tNu_Kin(I,iCp)= &
        & dot_product( tNuSp(SM%vIPole(1:SM%NPole),iCp) ,  SF%vXPole(1:SM%NPole) )
      end do
      !
    end if

    !! select case(trim(vFas(J)%Typ))
    !! case("PURE") !,"DISCRET")
    !!   K=  vFas(J)%iSpc
    !!   tNu_Kin(I,:)= tNuSp(K,:)
    !! case("MIXT")
    !!   K=  vFas(J)%iMix
    !!   SF= vMixFas(K)
    !!   SM= vMixModel(SF%iModel)
    !!   !
    !!   do iCp=1,size(tNu_Kin,2)
    !!     tNu_Kin(I,iCp)= &
    !!     & dot_product( tNuSp(SM%vIPole(1:SM%NPole),iCp) ,  SF%vXPole(1:SM%NPole) )
    !!   end do
    !! end select
    
  end do
  !----------------------------------------------/update stoichio tables
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_KinFas_Stoik_Nu"
  !
end subroutine Dynam_KinFas_Stoik_NuKin

subroutine Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)

  use M_Global_Vars,only: vSpc,vFas,vKinFas
  use M_Basis_Vars, only: vOrdPr

  real(dp),intent(in) :: vDG_Sp(:)
  real(dp),intent(in) :: tNu_Kin(:,:)
  real(dp),intent(out):: vDG_Kin(:)
  !
  integer:: nMk,nCp,I,J

  nCp= size(vOrdPr)
  nMk= size(vKinFas)
  !
  !-----------------------------------------------update thermo' param's
  do I=1,nMk
    J=vKinFas(I)%iFas
    !
    if(vFas(J)%iSpc>0) then
      vDG_Kin(I)= vDG_Sp(vFas(J)%iSpc)
    else if(vFas(J)%iMix>0) then
      vDG_Kin(I)= &
      & vFas(vFas(J)%iMix)%Grt &
      & - dot_product(tNu_Kin(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
    end if
    !
  end do
  !----------------------------------------------/update thermo' param's

  return
end subroutine Dynam_KinFas_DG_UpDate

subroutine Dynam_KinModel_Init(TimeFactor,TdgK)

  use M_Global_Vars,only: vKinModel,vSpc
  use M_T_Kinmodel, only: KinModel_Show
  use M_Basis_Vars, only: vPrmBk,vPrmFw
  use M_Dynam_Vars, only: vKinMod !,
  !
  real(dp),intent(in):: TimeFactor,TdgK

  call KinModel_Prm_Update( &
  & vPrmBk,vKinModel, &
  & vKinMod)
  !
  call KinModel_T_Update(TimeFactor,TdgK,vKinMod)
  !
  if(idebug>1) call KinModel_Show(fTrc,vSpc,vKinMod,vPrm=vPrmFw)

  return
end subroutine Dynam_KinModel_Init

subroutine Dynam_TP_Alloc

  use M_Global_Vars,only: vSpc
  use M_Dynam_Vars,only: vDG_Sp,vDG_As,vDLGam_As
  use M_Basis_Vars,only: nAs
  !
  if(allocated(vDG_Sp)) deallocate(vDG_Sp); allocate(vDG_Sp(1:size(vSpc)))
  !! nAs= count(vSpc%Typ(1:3)=="AQU")
  if(nAs>0) then
    if(allocated(vDG_As))    deallocate(vDG_As);    allocate(vDG_As(1:nAs))
    if(allocated(vDLGam_As)) deallocate(vDLGam_As); allocate(vDLGam_As(1:nAs))
  end if
  !
end subroutine Dynam_TP_Alloc

subroutine Dynam_Global_TP_Update(TdgK,Pbar)
!--
!-- calc' global thermo parameters of species, 
!-- sol'models and sol'phases for (TdgK,Pbar)
!-- and calc' local deltaG's (vDG_Sp,vDG_As)
!--

  use M_SolModel_Tools,only: SolModel_TP_Update
  use M_Global_Tools,  only: Global_TP_Update

  !--vars--
  use M_Global_Vars, only: vSpcDtb,vSpc,vMixModel,vDiscretModel,vDiscretParam
  use M_Global_Vars, only: vMixFas,vFas
  use M_Global_Vars, only: SolModel
  use M_System_Vars, only: System_Type
  use M_Basis_Vars,  only: nAs,tNuSp,tNuAs,vOrdPr,vOrdAs
  use M_Dynam_Vars,  only: vDG_Sp,vDG_As
  !--/

  real(dp),intent(in):: TdgK,Pbar
  !
  integer:: nCp

  call Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  !<new 200911
  if(System_Type=="AQUEOUS") call SolModel_TP_Update(TdgK,Pbar,SolModel)
  !</new 200911
  !
  nCp= size(vOrdPr)
  !
  vDG_Sp(:)= vSpc(:)%G0rt - matmul(tNuSp(:,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
  !
  if(nAs>0) &
  & vDG_As(:)= vSpc(vOrdAs(:))%G0rt &
  &          - matmul(tNuAs(:,1:nCp), vSpc(vOrdPr(1:nCp))%G0rt)

end subroutine Dynam_Global_TP_Update

subroutine Dynam_Init_KinFas_SatState
!-> vKinFas(:)%Dat%QsK, vKinFas(:)%Dat%cSatur
  use M_IOTools
  use M_Numeric_Const,only: Ln10
  use M_KinRate,    only: KinRate_SatState,KinRate_CalcQsK
  !
  use M_Global_Vars,only: vSpc,vKinFas,vFas,vSpcDat
  use M_Basis_Vars, only: vOrdPr,nCi
  !
  use M_Dynam_Vars, only: vDG_Kin,tNu_Kin,vKinMinim,Qsk_Iota
  use M_Dynam_Vars, only: vLnAct,vLnGam
  use M_Dynam_Vars, only: vMolF,vMolK
  !
  integer :: i,nCp,nAq,nCx
  real(dp):: QsK,QskIota
  character:: cStatus
  real(dp),dimension(1:size(vOrdPr)):: dQsKdLnXi
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Dynam_Init_KinFas_SatState"
  !
  QskIota= Qsk_Iota
  !
  nCp= size(vOrdPr)
  nCx= nCp-nCi
  nAq= count(vSpc%Typ=="AQU")
  !
  if(idebug>1) then
    !
    write(fTrc,'(A)') "Log10(QsK) of Mineral species"
    write(fTrc,'(/,A,/)') "iCp, Log(Act), LogK0, Prim.Species"
    do i=1,nCp
      write(fTrc,'(I2,A1,2(G15.3,A1),A)') &
      & I,                              T_, &
      & vSpcDat(vOrdPr(i))%LAct/Ln10,   T_, &
      & vSpc(vOrdPr(i))%G0rt/Ln10,      T_, &
      & trim(vSpc(vOrdPr(i))%NamSp)
    end do
    !
    write(fTrc,'(/,A,/)') "iMk, LogK0, Mineral"
    do i=1,size(vKinFas)
      write(fTrc,'(I2,A1, G15.3,A1, A)') &
      & I,T_, vFas(vKinFas(i)%iFas)%Grt/Ln10,T_,  trim(vFas(vKinFas(i)%iFas)%NamFs)
    end do
    !
  end if
  !
  !! call OutStrVec(51,vSpcDat(:)%LAct)
  !
  do I=1,size(vKinFas)
    !----------------------------calc' QsK, saturation state of minerals
    call KinRate_CalcQsK( &
    & nCi,nCx,      &
    & vDG_Kin(I),   &
    & tNu_Kin(I,:), &
    & log(vMolF(1:nAq)), & !IN:  Ln(aqu.species mole nrs)
    & vLnGam,       & !IN:  Ln(ActCoeff,aqu.species)
    & vLnAct,       & !IN:  Ln(ActivitySpecies), used for buffered species
    & vLnAct(1),    & !IN:  Ln(ActivitySolvent)
    & QsK,          & !OUT, mineral saturation
    & dQsKdLnXi)      !OUT
    !
    !------------ from QsK deduce cSatur, for branching to dissol/precip
    call KinRate_SatState( &
    & vKinFas(I),     & !IN:  mineral: for cSatur, QsKseuil
    & vMolK(I),       & !IN:  Nr Moles of mineral, used for checking nMol<=MolMinim
    & vKinMinim(I),   & !IN:
    & QsK,            & !IN:
    & QskIota,        & !IN
    & cStatus)          !OUT
    !
    !---------------------------------------------save in vKinFas(:)%Dat
    vKinFas(i)%Dat%QsK=  QsK
    vKinFas(I)%Dat%cSat= cStatus
    !
  end do
  !
  if(idebug>1) then
    write(fTrc,'(/,A,/)') "Minerals Saturation State, logQsK"
    do i=1,size(vKinFas)
      write(fTrc,'(I2,A1,G15.3,A1,2A)') &
      & I,T_, &
      & log(vKinFas(i)%Dat%QsK)/Ln10,T_, &
      & vKinFas(i)%NamKF,vKinFas(i)%Dat%cSat
    end do
  end if
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Dynam_Init_KinFas_SatState"
  !
end subroutine Dynam_Init_KinFas_SatState

subroutine KinModel_Prm_Update(vPrmBk,vKinModel,vKinMod) !
!--- update indexes of species in the current kinetic models vKinMod  --
  use M_T_Kinmodel,only: T_Kinmodel,KinModel_PrmIndex,KinModel_CalcCoeffs
  !
  integer,         intent(in) :: vPrmBk(:)
  type(T_Kinmodel),intent(in) :: vKinModel(:)
  type(T_Kinmodel),intent(out):: vKinMod(:)
  !
  integer ::I
  !
  vKinMod= vKinModel
  !vKinModel is in "global" database, vKinMod is a copy in Dynam_Vars
  !
  !calc. indexes of species in the current kinetic models vKinMod
  do I=1,size(vKinMod)
    call KinModel_PrmIndex(vPrmBk,vKinModel(I),vKinMod(I))
  end do
  !
end subroutine KinModel_Prm_Update

subroutine KinModel_T_Update(TimeFactor,TdgK,vKinMod)
!--
!--update temperature dependent param's of kinetic models
!--
  use M_T_Kinmodel,only: T_Kinmodel,KinModel_CalcCoeffs
  !
  real(dp),        intent(in)   :: TimeFactor,TdgK
  type(T_Kinmodel),intent(inout):: vKinMod(:)
  !
  integer:: I
  !
  do I=1,size(vKinMod)
    call KinModel_CalcCoeffs(vKinMod(I),TimeFactor,TdgK)
  end do
  !
  return
end subroutine KinModel_T_Update

subroutine NuTable_Sho(vCpn)
  use M_IoTools,    only: GetUnit
  use M_Files,      only: DirOut,Files_Index_Write
  use M_T_Component,only: T_Component
  use M_Global_Vars,only: vSpc,vFas,vKinFas
  use M_Dynam_Vars, only: tNu_Kin
  !
  type(T_Component),intent(in):: vCpn(:)
  !
  integer::F,iPr,I
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_stoik_nu_kin.log")
  !
  call Files_Index_Write(fHtm,&
  & trim(DirOut)//"_stoik_nu_kin.log",&
  & "Nu Table: stoikio of kin'species in terms of prim'species")
  !
  do iPr=1,size(tNu_Kin,2)
    write(F,'(A7,A1)',advance="no") trim(vSpc(vCpn(iPr)%iSpc)%NamSp),T_
  end do
  write(F,'(A)') "Name"
  !
  do I=1,size(tNu_Kin,1)
    do iPr=1,size(tNu_Kin,2)
      write(F,'(F7.2,A1)',advance="no") tNu_Kin(I,iPr), T_
    end do
    write(F,'(A)') trim(vFas(vKinFas(I)%iFas)%NamFs)
  end do
  !
  close(F)
  !
end subroutine NuTable_Sho

subroutine Dynam_ReStart_Time(T)
  use M_Dynam_Vars
  !
  type(T_DynTime),intent(in):: T
  
  TUnit=  T%TUnit
  Time=   T%Time
  dTime=  T%dTime
  TFinal= T%TFinal
  dTMin=  T%dTMin
  dTMax=  T%dTMax
  dTSav=  T%dTSav
  !
  !!--new 200912--
  if(TimeIsSeconds) then
    !
    Time=   Time   *T%TimeFactor
    TFinal= TFinal *T%TimeFactor
    dTime=  dTime  *T%TimeFactor
    dTmin=  dTmin  *T%TimeFactor
    dTMax=  dTMax  *T%TimeFactor
    dTSav=  dTSav  *T%TimeFactor
    !
    TimeFactor= One
    TimeScale=  T%TimeFactor !--------------TimeScale is used for output
    !
  else
    !
    TimeFactor= T%TimeFactor
    TimeScale=  One !-----------------------TimeScale is used for output
    !
  end if
  !!--new 200912--/
  
  if(iDebug>2) then
    print '(3(G12.3,A))', &
    & T%dTime,"=dTime /",T%dTMin,"=dTMin /",T%dTMax,"=dTMax /"
    call Pause_
  end if
  
  return
end subroutine Dynam_ReStart_Time
  !
subroutine Dynam_ReStart_Box(B,TimeFactor_)
  use M_Dynam_Vars,only: T_DynBox
  use M_Dynam_Vars,only: dX,UDarcy,VBox,FOut,PhiF,RhoF,UpdateMassFluid
  use M_Dynam_Vars,only: TimeIsSeconds
  use M_Dynam_Vars,only: Dynam_nTotalNewtIter
  !
  type(T_DynBox), intent(in):: B
  real(dp),       intent(in):: TimeFactor_
  !
  dX=      B%dX
  UDarcy=  B%UDarcy
  VBox=    B%VBox
  FOut=    B%FOut
  PhiF=    B%PhiF
  RhoF=    B%RhoF
  UpdateMassFluid= B%VFixed
  !
  !!<new 200912>
  if(TimeIsSeconds) UDarcy= UDarcy /TimeFactor_
  !!</new 200912>
  !
  !! nCell=   B%nCell
  Dynam_nTotalNewtIter = 0
end subroutine Dynam_ReStart_Box

subroutine Dynam_Zero_Time(T)
!.default dynamic parameters, in case not found in Input
  use M_Dynam_Vars,only: T_DynTime
  type(T_DynTime),intent(out):: T
  !
  !time parameters
  T%TUnit=  "DAY"  !default unit used for all time data
  T%TFinal= 5.0D+6 !duration of simulation, in TUnit
  T%dTime=  1.0D-6 !initial time step
  T%dTMin=  1.0D-16 !minimal time step
  T%dTMax=  Zero   !maximal time step
  T%dTSav=  Zero   !delta time between saves
  T%TimeFactor= Dynam_TimeFactor(T%TUnit)
  !
end subroutine Dynam_Zero_Time

real(dp) function Dynam_TimeFactor(S) result(X)
!--
!-- Time Factor --
!--
  character(len=*),intent(in):: S
  
  X=One
  select case(trim(S))
    case("SECOND")  ;  X= One
    case("MINUTE")  ;  X= 60._dp
    case("HOUR")    ;  X= 3600._dp
    case("DAY")     ;  X= 3600._dp * 24._dp
    case("YEAR")    ;  X= 3600._dp * 24._dp * 365._dp
  end select
  
end function Dynam_TimeFactor

subroutine Dynam_Zero_Box(B)
!--
!-- initialize box with default values --
!--
  use M_Dynam_Vars,only: T_DynBox, VBox0
  type(T_DynBox),intent(out):: B
  !
  !--- box parameters
  B%UDarcy=Zero   ! flux rate, metre/time
  B%VBox=  One    ! volume of box, m3 (just scaling parameter (?), useful for "outside")
  B%dX=    One    ! length of box, meter
  B%PhiF=  0.5D0  ! porosity
  B%RhoF=  1.D3   ! fluid density, kg/m^3
  B%FOut=  Zero   !
  B%VFixed=.true. ! Vfixed is .false. for free volume
  !
  return
end subroutine Dynam_Zero_Box

subroutine Dynam_Zero_Numeric
  use M_Dynam_Vars
  !
  !--debug parameters--
  DebNewt=   .false.
  DebJacob=  .false.
  TestJacob= .false.
  TestMax=   .false.
  bFinDif=   .false.  !.true. !
  !
  !--calculation mode parameters--
  UpdateMassFluid= .true.
  !
  ! provisionally, dVmAdLnX_M is not computed in this new
  ! -> we assume that this factor will not be implicited
  Implicit_ActivFactor=.false.
  !
  Implicit_Surface= .false. ! .true. !
  !
  !---------------------------------------------------NUMERIC parameters
  !
  !--cMethod= select which Newton method is used (Kelley, Walker, ...)--
  cMethod= "NEWTONPRESS"   !-> default method is Newton_Press
  cMethod= "NEWTONKELLEY"  !-> default method is Newton_Kelley
  !--/
  iCtrlTime= 1
  !
  !-----------------------------------parameters for time step selection
  !-- if(Its>NewtMaxIts)  decrease dtime, cycle Newton
  !-- if(Its>NewtIterMax) decrease dtime (outside Newton)
  !-- if(Its-NewtIterMax) increase dtime (outside Newton)
  NewtMaxIts=  40 !75 !30  !100 !
  NewtIterMax= 20 !50 !16  !60  !
  NewtIterMin= 10 !25 !8   !30  !
  !--/
  !
  NewtTolF=      1.0E-07_dp  
  !NewtTolF=     1.0E-09_dp  
  != convergence criterion on function values
  NewtTolF_Equil= 1.0E-06_dp  !
  !NewtTolF_Equil= NewtTolF
  !
  NewtTolMin=     1.0E-06_dp 
  != criterion for spurious convergence
  NewtTolX= epsilon(real(dp))
  != convergence criterion on dx -> approx. 2.E-16
  !
  ! usage in Newton:
  !   NewtTolF: if( maxval(abs(fVec(:)))<NewtTolF ) Newton_Ok
  !   NewtTolX: if( maxval( abs(vX(:)-vXOld(:))/MAX(abs(vX(:)),One)<NewtTolX ) Converge on X
  !
  !--------------------------------------------------/NUMERIC parameters
  
  return
end subroutine Dynam_Zero_Numeric

subroutine Dynam_InitCalc_Restart
  use M_Global_Vars,only: nAq,vSpc,vFas,vKinFas,vSpcDat
  use M_Basis_Vars, only: vOrdAq
  use M_Dynam_Vars, only: vMolK,vMolF
  use M_Dynam_Vars, only: vMolarVol,vKinQsk,VBox,PhiF
  use M_Dynam_Vars, only: Time, Tfinal, Dtime, Dtmin, Dtmax
  use M_Dynam_Vars, only: Dynam_nStep, Dynam_nNewtIter, Dynam_nTotalNewtIter
  !
  integer:: nMk
  
  nMk= size(vKinFas)
  !
  Dynam_nStep=     0
  Dynam_nNewtIter= 0
  Dynam_nTotalNewtIter= 0
  !
  !------------------------------------------molar volumes of kin'phases
  vMolarVol(1:nMk)= vFas(vKinFas(1:nMk)%iFas)%VolM3
  !---------------------------------------mole nr's of kin'phases in box
  vMolK(1:nMk)=     vKinFas(1:nMk)%Dat%PhiM /vMolarVol(1:nMk) *VBox
  !---------------------------------------------sat'states of kin'phases
  vKinQsk(1:nMk)=   vKinFas(1:nMk)%Dat%QsK
  !
  PhiF= One - sum(vKinFas(1:nMk)%Dat%PhiM)
  !             !-> vol'fraction of fluid
  vMolF(1:nAq)= vSpcDat(vOrdAq(1:nAq))%Mole
  !             !-> mole nr's of aqu'species in box
  !
  if (idebug>1) &
  & print '(A,3F12.5, 2G12.3)', &
  & 'Restart avec [Time, Tfinal, Dtime, Dtmin, Dtmax =]', &
  & Time, Tfinal, Dtime, Dtmin, Dtmax
  
  return
end subroutine Dynam_InitCalc_Restart

subroutine Dynam_Restart
  use M_Global_Vars,only: nAq,vSpc,vFas,vKinFas,vSpcDat
  use M_Basis_Vars, only: vOrdAq
  !
  use M_Dynam_Vars, only: vMolarVol,vMolK,vMolF
  use M_Dynam_Vars, only: vSurfK,vKinQsk,vStatusK
  use M_Dynam_Vars, only: VBox,PhiF
  use M_Dynam_Vars, only: Time, Tfinal, Dtime, Dtmin, Dtmax
  !
  integer:: nMk
  
  nMk= size(vKinFas)
  !
  vMolarVol(1:nMk)= vFas(vKinFas(1:nMk)%iFas)%VolM3
  !                 !-> molar volumes of kin'phases
  vMolK(1:nMk)=     vKinFas(1:nMk)%Dat%PhiM /vMolarVol(1:nMk) *VBox
  !                 !-> mole nr's of kin'phases in box
  vSurfK(1:nMk)=    vKinFas(1:nMk)%Dat%Surf *VBox
  !                 !-> surfaces of kin'phases in box
  vKinQsk(1:nMk)=   vKinFas(1:nMk)%Dat%QsK
  !                 !-> sat'states of kin'phases
  vStatusK(1:nMk)=  vKinFas(1:nMk)%Dat%cSat
  !
  PhiF= One - sum(vKinFas(1:nMk)%Dat%PhiM)
  !                 !-> vol'fraction of fluid
  vMolF(1:nAq)= vSpcDat(vOrdAq(1:nAq))%Mole
  !                 !-> mole nr's of aqu'species in box

  if (idebug>1) &
  & print '(A,3F12.5, 2G12.3)', &
  & 'Restart avec [Time, Tfinal, Dtime, Dtmin, Dtmax =]', &
  & Time, Tfinal, Dtime, Dtmin, Dtmax
  
  return
end subroutine Dynam_Restart

subroutine Dynam_Save

  use M_Global_Vars,only: vSpc,vFas,vKinFas,vSpcDat
  !
  use M_Basis_Vars, only: vOrdAq
  use M_Dynam_Vars, only: vCpnBox
  use M_Dynam_Vars, only: vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  use M_Dynam_Vars, only: vMolK,vSurfK,vKinQsk,vStatusK,vMolarVol,VBox
  !
  use M_Dynam_Cell
  
  type(T_Cell_State):: Cell
  !
  integer :: nCp,nMk,nAq
  real(dp):: PhiF
  
  nCp= size(vCpnBox)
  nAq= count(vSpc(:)%Typ(1:3)=="AQU")
  nMk= size(vKinFas)
  !
  vCpnBox(:)%Mole=  vTotF(:)
  vCpnBox(:)%LnAct= vSpcDat(vCpnBox(:)%iSpc)%LAct
  !
  vSpcDat(vOrdAq(1:nAq))%Mole= vMolF(1:nAq)
  vSpcDat(vOrdAq(1:nAq))%LAct= vLnAct(1:nAq)
  vSpcDat(vOrdAq(1:nAq))%LGam= vLnGam(1:nAq)
  !
  vKinFas(1:nMk)%Dat%PhiM=   vMolK(1:nMk) *vMolarVol(1:nMk) / VBox
  vKinFas(1:nMk)%Dat%Surf=   vSurfK(1:nMk) /Vbox
  vKinFas(1:nMk)%Dat%QsK=    vKinQsk(1:nMk)
  vKinFas(1:nMk)%Dat%cSat=   vStatusK(1:nMk)
  !
  ! vFas(vKinFas(1:nMk)%iFas)%MolFs= vMolK(1:nMk)
  
  PhiF= One - sum(vKinFas(1:nMk)%Dat%PhiM)
  !
  call Cell%New_(nCp,nAq,nMk)
  !
  call Cell%Set_( &  !
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
  
  call Cell%Get_( &  !
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
  
  !!print *,"vMolK(1)",vMolK(1)  ;  pause
  
  return
end subroutine Dynam_Save

subroutine Dynam_Save_ToFile
  use M_IOTools,only: GetUnit
  use M_Files,  only: DirOut
  use M_Numeric_Const,only: Ln10
  use M_T_KinFas,only: T_KinFas
  !
  use M_Global_Vars,only: vEle,vSpc,vKinFas,vKinModel,vFas
  use M_Dynam_Vars, only: vCpnBox,vCpnInj
  use M_Dynam_Vars, only: PhiF
  !
  type(T_KinFas):: K
  integer :: F, I
  real(dp):: X
  !
  call GetUnit(F)
  open(F,file=trim(DirOut)//"_restart.inn")

  !-----------------------------------------------------------SYSTEM.BOX
  write(F,'(A)') "SYSTEM.BOX"
  !! print *,vCpnBox(1)%Mole  ;  pause
  X= 55.51D0 /vCpnBox(1)%Mole
  !
  do I=1,size(vCpnBox)
    write(F,'(A)',   advance="NO") "  "
    write(F,'(A,1X)',advance="NO") trim(vEle(vCpnBox(I)%iEle)%NamEl)
    write(F,'(A,1X)',advance="NO") trim(vCpnBox(I)%Statut)
    write(F,'(A,1X)',advance="NO") trim(vSpc(vCpnBox(I)%iSpc)%NamSp)
    !
    select case(trim(vCpnBox(I)%Statut))
      case("INERT")   ;  write(F,'(G15.6)') vCpnBox(I)%Mole *X
      case("BUFFER")  ;  write(F,'(G15.6)') - vCpnBox(I)%LnAct /Ln10
    end select
    !
  end do
  write(F,'(A)') "END SYSTEM.BOX"
  !----------------------------------------------------------/SYSTEM.BOX

  !--------------------------------------------------------SYSTEM.INJECT
  write(F,'(A)') "SYSTEM.INJECT"
  X= 55.51D0 /vCpnInj(1)%Mole
  !
  do I=1,size(vCpnInj)
    write(F,'(A)',   advance="NO") "  "
    write(F,'(A,A1)',advance="NO") trim(vEle(vCpnInj(I)%iEle)%NamEl),T_
    write(F,'(A,A1)',advance="NO") trim(vCpnInj(I)%Statut),          T_
    write(F,'(A,A1)',advance="NO") trim(vSpc(vCpnInj(I)%iSpc)%NamSp),T_
    !
    select case(trim(vCpnInj(I)%Statut))
      case("INERT")   ;  write(F,'(G15.6)') vCpnInj(I)%Mole *X
      case("BUFFER")  ;  write(F,'(G15.6)') - vCpnInj(I)%LnAct /Ln10
    end select
    !
  end do
  write(F,'(A)') "END SYSTEM.INJECT"
  !-------------------------------------------------------/SYSTEM.INJECT

  !-------------------------------------------------------- DYNAMIC.ROCK
  write(F,'(A)') "DYNAMIC.ROCK"
  do I=1,size(vKinFas)
    !
    K= vKinFas(I)
    !
    write(F,'(A)',   advance="NO") "  "
    write(F,'(A,A1)',advance="NO") trim(K%NamKF),T_
    !
    write(F,'(A,A1,G15.6,A1)',advance="NO") "SURFACE",T_, K%Dat%SurfKg,T_
    write(F,'(A,A1,G15.6,A1)',advance="NO") "VOLUME", T_, K%Dat%PhiM,  T_
    !
    if(K%iKin>0) then
      if(trim(K%NamKF)/=trim(vKinModel(K%iKin)%Name)) &
      & write(F,'(2(A,A1))',advance="NO") &
      & "KINETICS", T_, &
      & trim(vKinModel(K%iKin)%Name),T_
    end if
    !
    if(trim(K%NamKF)/=trim(vFas(K%iFas)%NamFs)) &
    & write(F,'(2(A,A1))',advance="NO") &
    & "SPECIES", T_, &
    & trim(vFas(K%iFas)%NamFs),T_!
    !
    if(K%cMode/="D") write(F,'(A1,A1)',advance="NO") K%cMode, T_
    !
    write(F,*)
    !
  end do
  write(F,'(A)') "END DYNAMIC.ROCK"
  !--------------------------------------------------------/DYNAMIC.ROCK
  !
  !--------------------------------------------------------------DYNAMIC
  write(F,'(A)') "DYNAMIC"
    write(F,'(A)',   advance="NO") "  "
    write(F,'(A,A1,G15.6)',advance="NO") "POROSITY",T_, PhiF
    write(F,*)
  write(F,'(A)') "END DYNAMIC"
  !-------------------------------------------------------------/DYNAMIC
  !
  close(F)
end subroutine Dynam_Save_ToFile

subroutine Dynam_Clean
  use M_Dynam_Vars
  !
  if(allocated(vCpnInj))    deallocate(vCpnInj)
  if(allocated(vCpnBox))    deallocate(vCpnBox)
  !
  if(allocated(vTotF))      deallocate(vTotF)
  if(allocated(vMolF))      deallocate(vMolF)
  if(allocated(vTotInj))    deallocate(vTotInj)
  !
  if(allocated(vLnAct))     deallocate(vLnAct)
  if(allocated(vLnGam))     deallocate(vLnGam)
  if(allocated(vLnBuf))     deallocate(vLnBuf)
  !
  if(allocated(vMolarVol))  deallocate(vMolarVol)
  if(allocated(vMolK))      deallocate(vMolK)
  if(allocated(vSurfK))     deallocate(vSurfK)
  if(allocated(vStatusK))   deallocate(vStatusK)
  !
  if(allocated(vMolK0))      deallocate(vMolK0)
  if(allocated(vSurfK0))     deallocate(vSurfK0)
  !
  if(allocated(vKinMinim))  deallocate(vKinMinim)
  if(allocated(vLKinActiv)) deallocate(vLKinActiv)
  if(allocated(vLEquActiv)) deallocate(vLEquActiv)
  if(allocated(vKinPrm))    deallocate(vKinPrm)
  !
  if(allocated(vKinMod))    deallocate(vKinMod)
  !
  if(allocated(vQTotInj))    deallocate(vQTotInj)
  !
  if(allocated(vKinQsK))     deallocate(vKinQsK)
  if(allocated(vVmAct))      deallocate(vVmAct)
  if(allocated(vVmQsK))      deallocate(vVmQsK)
  !
  if(allocated(vDG_Sp))     deallocate(vDG_Sp)
  if(allocated(vDG_As))     deallocate(vDG_As)
  if(allocated(vDLGam_As))  deallocate(vDLGam_As)
  if(allocated(vDG_Kin))    deallocate(vDG_Kin)
  !
  if(allocated(tNu_Kin))    deallocate(tNu_Kin)
  if(allocated(tAlfKin))    deallocate(tAlfKin)
  !
  if(allocated(vLnBuf))     deallocate(vLnBuf)
  !
  if(allocated(vWeitSp))    deallocate(vWeitSp)
  if(allocated(vWeitCp))    deallocate(vWeitCp)
  !
  if(allocated(tStoikioAqu)) deallocate(tStoikioAqu)
  if(allocated(tStoikioKin)) deallocate(tStoikioKin)
  !
end subroutine Dynam_Clean

end module M_Dynam_Tools

!! subroutine SaturStateInit_Show( &
!! & Time,PhiF,vKinFas,vSurfK0,vMolK,vMolarVol)
!!   use M_T_KinFas,only: T_KinFas
!!   !
!!   real(dp),                   intent(in):: Time,PhiF
!!   type(T_KinFas),dimension(:),intent(in):: vKinFas
!!   real(dp),      dimension(:),intent(in):: vSurfK0,vMolK,vMolarVol
!!   !
!!   type(T_KinFas)::M
!!   integer::i
!!   !
!!   write(fTrc,'(/,A,/)') "!!!SaturState_Show_begin"
!!   !
!!   write(fTrc,'(8(A,A1))') &
!!   & "Name",   T_,"Satur",T_,"Time", T_,"VFluid",T_, &
!!   & "Surface",T_,"QsK",  T_,"vMolK",T_,"Rate",  T_
!!   do i=1,size(vKinFas)
!!     M=vKinFas(i)
!!     write(fTrc,&
!!     & '(   A,A1,A,A1,G15.9,A1,F12.3,A1,E12.3,A1,E12.6,A1,E12.3,A1,E7.2,A1)') &
!!     & M%Name,T_,M%Dat%cSat,T_,Time, T_,PhiF, T_,&
!!     & M%Dat%Surf,T_,M%Dat%QsK,T_,vMolK(i),T_,vMolarVol(i),T_
!!   end do
!!   !!! do i=1,size(vKinFas)
!!   !!!   M=vKinFas(i)
!!   !!!   write(fTrc,'(A,A1,G15.8, A1,G15.8, A1,G15.8, A1)') &
!!   !!!   & trim(M%Name),T_,M%Kd_A,T_,M%Kd_W,T_,M%Kd_B,T_
!!   !!! end do
!!   write(fTrc,'(/,A,/)') "!!!SaturState_Show_end"
!! end subroutine SaturStateInit_Show

