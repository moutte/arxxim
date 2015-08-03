MODULE M_Dynam_Tools
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,Stop_,fTrc,fHtm,T_,Pause_
  
  IMPLICIT NONE
  
  PRIVATE
  !
  PUBLIC:: Dynam_Init
  PUBLIC:: Dynam_Alloc
  PUBLIC:: Dynam_TP_Update
  PUBLIC:: Dynam_TimeFactor
  PUBLIC:: Dynam_Zero_Box
  PUBLIC:: Dynam_Zero_Time
  PUBLIC:: Dynam_Zero_Numeric
  PUBLIC:: Dynam_ReStart_Time
  PUBLIC:: Dynam_ReStart_Box
  PUBLIC:: Dynam_Save
  PUBLIC:: Dynam_Save_ToFile
  PUBLIC:: Dynam_Clean
  
CONTAINS

SUBROUTINE Dynam_Alloc
  USE M_Global_Vars,ONLY: vSpc,vKinFas,vKinModel
  USE M_Dynam_Vars, ONLY: vCpnBox
  USE M_Dynam_Vars, ONLY: vWeitSp,vWeitCp,vKinMod
  !
  INTEGER :: nCp,nAq,nMk
  !
  nAq=COUNT(vSpc%Typ=="AQU")
  IF(nAq==0) CALL Stop_("NO Aqueous Species ???")
  !
  nCp= SIZE(vCpnBox)
  nMk= SIZE(vKinFas)
  !
  CALL Dynam_Alloc_FluidBox(nCp)
  !-> ALLOCATE vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  CALL Dynam_Alloc_FluidInj(nCp)
  !-> ALLOCATE vTotInj,vQTotInj
  CALL Dynam_Alloc_KinFas(nMk)
  !-> ALLOCATE vMolK,vMolarVol,vKinQsK,vKinMinim,vSurfK0,...
  CALL Dynam_KinFas_Stoik_Alloc(nMk,nCp)
  !-> ALLOCATE vDG_Kin,tNu_Kin,tAlfKin
  !
  CALL Dynam_TP_Alloc
  !
  ALLOCATE(vKinMod(1:SIZE(vKinModel)))
  !
  ALLOCATE(vWeitSp(1:nAq))
  ALLOCATE(vWeitCp(1:nCp))
  !
ENDSUBROUTINE Dynam_Alloc

SUBROUTINE Dynam_TP_Update(TdgK,Pbar,TimeFactor)
  USE M_Dynam_Vars,ONLY: tNu_Kin,vDG_Sp,vDG_Kin,vKinMod
  !
  REAL(dp),INTENT(IN):: TdgK,Pbar,TimeFactor
  !
  CALL Dynam_Global_TP_Update(TdgK,Pbar)
  !
  CALL Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)
  !
  CALL KinModel_T_Update(TimeFactor,TdgK,vKinMod)
  !
ENDSUBROUTINE Dynam_TP_Update

SUBROUTINE Dynam_Init(TdgK,Pbar,T,B)
!--
!-- initializations of variables used in "dynamic" calculations
!--
  USE M_IOTools
  USE M_Global_Vars,ONLY: vSpc,vFas,vKinFas,vMixFas,vMixModel,vEle
  USE M_Global_Vars,ONLY: SolModel
  USE M_Basis_Vars, ONLY: vOrdAq

  USE M_Dynam_Vars,ONLY: vCpnBox,vCpnInj
  USE M_Dynam_Vars,ONLY: T_DynTime,T_DynBox
  USE M_Dynam_Vars,ONLY: vWeitSp,vWeitCp
  USE M_Dynam_Vars,ONLY: vDG_Sp,vDG_Kin,tNu_Kin
  USE M_Dynam_Vars,ONLY: TdgK0,Pbar0
  USE M_Dynam_Vars,ONLY: vTotInj,vTotF,vMolF,vLnAct,vLnGam,vLnBuf

  !! USE M_Stockvar_Kinxim,ONLY: LSTOCK,INIT_STOCKVAR
  !---------------------------------------------------------------------
  REAL(dp),       INTENT(IN)   :: TdgK,Pbar
  TYPE(T_DynTime),INTENT(INOUT):: T
  TYPE(T_DynBox), INTENT(INOUT):: B
  !---------------------------------------------------------------------
  INTEGER, PARAMETER:: MaxStockVar = 40000
  INTEGER :: nCp,nAq,nMk
  !---------------------------------------------------------------------
  !REAL(dp):: &
  !& MFluid, &    !mass of fluid corresponding to a given composition
  !& MFluidBox    !mass of fluid, of known density RhoF, that can occupy PhiF*VBox
  !MFluid=    DOT_PRODUCT(vSpc(1:nAq)%WeitKg,vMolF(1:nAq))
  !MFluidBox= MFluidBox= RhoF*VBox*PhiF
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_Init"
  !
  nCp=SIZE(vCpnBox)
  nMk=SIZE(vKinFas)            ! ; IF(nMk==0) CALL Stop_("NO Kinetic Species ???")
  nAq=COUNT(vSpc%Typ=="AQU")   ; IF(nAq==0) CALL Stop_("NO Aqueous Species ???")
  !
  !---------------- init table to store the time evolution of species --
  !! IF (LSTOCK) CALL  init_stockvar(MaxStockVar, nAq)
  !
  !--------------------------------------------------update fluid data--
  !B%PhiFInj= B%PhiF !initial vol. fraction of fluid
  B%RhoF= 1.0D3 !! Solvent%Rho(iTP)*1.0D3 !fluid density
  IF(iDebug>0) WRITE(fTrc,'(/,A,G15.3,/)') "RhoF",B%RhoF
  !
  !---------------------------------------------------update time data--
  T%Time= Zero
  !
  IF(iDebug>2) CALL Dynam_Init_StoikioEle
  !
  CALL Dynam_Init_FluidBox_Compo( &   !
  & vCpnBox,B%RhoF,B%VBox,B%PhiF, &   !IN
  & vTotF,vMolF,vLnAct,vLnGam,vLnBuf) !OUT
  !
  CALL Dynam_Init_FluidInj_Compo( & !
  & vCpnInj,B%RhoF,B%VBox, &  !IN
  & vTotInj)                  !OUT
  !
  CALL Dynam_Init_KinFas_Texture(nMk,B%VBox,B%PhiF) !IN
  !-> initialize/update vKinFas,vMolK,SurfMinim,vMolarVol,vKinMinim,vSurfK0
  !
  !CALL Dynam_KinFas_Stoik_Init(nMk,nCp)
  CALL Dynam_KinFas_Stoik_Alfa(nMk,nCp)
  !
  CALL Dynam_KinFas_Stoik_NuKin(nMk)
  IF(iDebug>2) CALL NuTable_Sho(vCpnBox)
  !
  TdgK0= TdgK
  Pbar0= Pbar
  !
  CALL Dynam_Global_TP_Update(TdgK,Pbar)
  !
  CALL Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)
  !
  CALL Dynam_KinModel_Init(T%TimeFactor,TdgK) !IN
  !
  CALL Dynam_Init_KinFas_SatState
  !
  vWeitSp(:)= vSpc(vOrdAq(:))%WeitKg
  !-> Aqu'Species' weitkg sorted according to local Species' order
  vWeitCp(:)= vEle(vCpnBox(:)%iEle)%WeitKg
  !-> Components' weitkg, i.e. Elements' weitkg in the Components' order ...
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_Init"
  
  RETURN
END SUBROUTINE Dynam_Init

SUBROUTINE Dynam_Init_StoikioEle
  USE M_Global_Vars,ONLY: vSpc,vEle,vKinFas,vFas,tFormula
  USE M_Basis_Vars, ONLY: vOrdAq
  USE M_Dynam_Vars, ONLY: tStoikioAqu,tStoikioKin
  !
  INTEGER:: I,J,N,nMk
  !
  N= SIZE(vEle)
  ALLOCATE(tStoikioAqu(SIZE(vOrdAq),N))
  DO I=1,SIZE(vOrdAq)
    J= vOrdAq(I)
    tStoikioAqu(I,:)= tFormula(:,J)
    !vSpc(J)%vStoikio(1:N)/REAL(vSpc(J)%vStoikio(0))
    !~ do K=1,N
      !~ write(6,'(F7.2,1X)',advance="no") tStoikioAqu(I,K)
    !~ enddo
    !~ write(6,*)
  ENDDO
  nMk= SIZE(vKinFas)
  ALLOCATE(tStoikioKin(nMk,N))
  DO I=1,nMk
    J= vFas(vKinFas(I)%iFas)%iSpc
    tStoikioKin(I,:)= tFormula(:,J)
    !vSpc(J)%vStoikio(1:N)/REAL(vSpc(J)%vStoikio(0))
    !! do K=1,N
    !!   write(6,'(F7.2,1X)',advance="no") tStoikioKin(I,K)
    !! enddo
    !! write(6,*)
  ENDDO
  !! pause
  
  RETURN
END SUBROUTINE Dynam_Init_StoikioEle

!------------------------------------------------allocations for fluid--
SUBROUTINE Dynam_Alloc_FluidBox(nCp)
!=
!= allocate vTotF,vMolF,vLnAct,vLnGam,vLnBuf
!=
  USE M_Global_Vars,ONLY: vSpc,vEle
  !! USE M_Basis_Vars, ONLY: nMx,nAx
  USE M_Dynam_Vars, ONLY: vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  !
  INTEGER,INTENT(IN):: nCp
  !
  INTEGER:: nAq
  !
  ALLOCATE(vTotF(1:nCp))
  !
  nAq= COUNT(vSpc%Typ=="AQU")
  !
  ALLOCATE(vMolF(1:nAq))
  ALLOCATE(vLnAct(1:nAq))
  ALLOCATE(vLnGam(1:nAq))
  !
  !! ALLOCATE(vLnBuf(1:nAx+nMx))
  ALLOCATE(vLnBuf(1:nAq))
  !
  RETURN
ENDSUBROUTINE Dynam_Alloc_FluidBox

SUBROUTINE Dynam_Alloc_FluidInj(nCp)
  USE M_Dynam_Vars, ONLY: vTotInj,vQTotInj !,vTotInjFile
  !
  INTEGER,INTENT(IN):: nCp
  !
  ALLOCATE(vTotInj(1:nCp))       ; vTotInj= Zero
  ALLOCATE(vQTotInj(1:nCp))      ; vQTotInj= Zero
  !
ENDSUBROUTINE Dynam_Alloc_FluidInj
!-----------------------------------------------/allocations for fluid--

!------------------------------------------------ initialize vTotInj  --
!-- initialize mole nrs injected fluid
!-- units: mole nrs in a volume VBox for a fluid density RhoF
SUBROUTINE Dynam_Init_FluidInj_Compo( &
& vCpn,RhoF,VBox, &
& vTotInj)
  USE M_T_Component,ONLY: T_Component
  USE M_Global_Vars,ONLY: vEle !(for %WeitKg)
  !
  TYPE(T_Component),INTENT(IN)::  vCpn(:)
  REAL(dp),         INTENT(IN)::  RhoF,VBox
  REAL(dp),         INTENT(OUT):: vTotInj(:)
  !
  INTEGER:: I
  !
  vTotInj(:)= vCpn(:)%Mole
  vTotInj(:)= vTotInj(:) &
  &           /DOT_PRODUCT(vEle(vCpn(:)%iEle)%WeitKg,vTotInj(:))  & !-> nr.mole/Kg
  &           *RhoF & !-> nr.mole/M3
  &           *VBox   !-> in mole nrs in the box
  !
  IF(iDebug>2) THEN
    PRINT '(A)',"Dynam_Init_FluidInj_Compo"
    DO I=1,SIZE(vTotInj)
      PRINT '(A,G15.6)',vEle(vCpn(I)%iEle)%NamEl,vTotInj(I)
    ENDDO
  ENDIF
  !
ENDSUBROUTINE Dynam_Init_FluidInj_Compo

!----------------------------------------------initialize fluid in box--
SUBROUTINE Dynam_Init_FluidBox_Compo( &       !
& vCpn,RhoF,VBox,PhiF, &                      !
& vTotF,vMolF,vLnAct,vLnGam,vLnBuf)
  !
  USE M_T_Component,ONLY: T_Component
  USE M_Global_Vars,ONLY: vSpc,vEle
  USE M_Basis_Vars, ONLY: nAx,nMx,vOrdAq,nCi
  !
  TYPE(T_Component),INTENT(IN) :: vCpn(:)
  REAL(dp),         INTENT(IN) :: RhoF,VBox,PhiF
  REAL(dp),         INTENT(OUT):: vTotF(:),vMolF(:)
  REAL(dp),         INTENT(OUT):: vLnAct(:),vLnGam(:),vLnBuf(:)
  !
  INTEGER:: nAq,I
  
  nAq= SIZE(vOrdAq)
  !
  !-----------------------------------retrieve results from speciation--
  vMolF(:)= vSpc(vOrdAq(:))%Dat%Mole
  vLnAct(1:nAq)=vSpc(vOrdAq(:))%Dat%LAct
  vLnGam(1:nAq)=vSpc(vOrdAq(:))%Dat%LGam
  !
  ! given mole numbers of aqu.species for 1kg H2O
  ! (=results of initial speciation)
  ! calculate mole numbers of aqu.species in a box of volume VBox,
  ! assuming the fluid has a known density RhoF
  !
  vMolF(:)= vMolF(:) &
  & /DOT_PRODUCT(vMolF(1:nAq), vSpc(vOrdAq(:))%WeitKg) *RhoF & !-> mole/m3
  & *VBox*PhiF !-> mole in fluid in box
  !
  vTotF(:)= vCpn(:)%Mole
  !
  !--moles of Components in fluid in the box--
  vTotF(:)= vTotF(:) &
  & /DOT_PRODUCT(vMolF(:), vSpc(vOrdAq(:))%WeitKg) *RhoF & !-> mole/m3
  & *VBox*PhiF
  !
  !--activities of buffered species--
  vLnBuf(:)= 0.D0
  !! PRINT *,"nAq, SIZE(vLnBuf)",nAq, SIZE(vLnBuf)  ; PAUSE
  IF(nAx+nMx>0) THEN
    DO i=1,nAx+nMx
      vLnBuf(I)= vSpc(vCpn(nCi+I)%iSpc)%Dat%LAct
    ENDDO
  ENDIF
  !
  IF(iDebug>2) THEN
    PRINT '(A)',"Dynam_Init_FluidBox_Compo"
    DO I=1,SIZE(vCpn)
      IF(vCpn(I)%Statut=="BUFFER") &
      & PRINT '(A,A3,G15.6)', "BUFFER= ",vCpn(I)%NamCp,vSpc(vCpn(I)%iSpc)%Dat%LAct
      IF(vCpn(I)%Statut=="INERT") &
      & PRINT '(A,A3,G15.6)', "INERT=  ",vCpn(I)%NamCp,vCpn(I)%Mole
    ENDDO
  ENDIF
  
  RETURN
ENDSUBROUTINE Dynam_Init_FluidBox_Compo

!------------------------------allocations for kinetic & equil' phases--
SUBROUTINE Dynam_Alloc_KinFas(N)
!--
!-- allocate tables related with kinetic minerals
!--
  USE M_Dynam_Vars, ONLY: &
  & vMolK,vSurfK,vMolarVol,vKinQsK,vKinMinim,vStatusK,&
  & vSurfK0,vMolK0,vVmQsK,vVmAct, &
  & vLKinActiv,vLEquActiv,vKinPrm
  !
  INTEGER, INTENT(IN):: N
  !
  ALLOCATE(vMolK(1:N))         ; vMolK=  Zero   !moles kin.species in the box
  ALLOCATE(vSurfK(1:N))        ; vSurfK= Zero   !surface kin.speices in the box
  ALLOCATE(vKinMinim(1:N))     ; vKinMinim=Zero !minimal mole numbers
  ALLOCATE(vStatusK(1:N))
  !
  ALLOCATE(vKinQsK(1:N))       ; vKinQsK= One !
  !
  ALLOCATE(vMolarVol(1:N))     ; vMolarVol= Zero !molar volumes
  !
  ALLOCATE(vMolK0(1:N))        ; vMolK0=  Zero !mole numbers at ref'time (for kin'rate)
  ALLOCATE(vSurfK0(1:N))       ; vSurfK0= Zero !surface at ref'time (for "scaled" output and for CRUNCH mode)
  !
  ALLOCATE(vVmQsK(1:N))        ; vVmQsK=  Zero !Precip /Dissol. Rate, VmF=VmAct*VmQsK
  ALLOCATE(vVmAct(1:N))        ; vVmAct=  Zero !Precip /Dissol. Rate, VmF=VmAct*VmQsK
  !!ALLOCATE(vVm0  (1:N))        ; vVm0=    Zero !Precip /Dissol. Rate, VmF=VmAct*VmQsK
  !
  ALLOCATE(vLKinActiv(1:N))    ; vLKinActiv=.FALSE.
  ALLOCATE(vLEquActiv(1:N))    ; vLEquActiv=.FALSE.
  ALLOCATE(vKinPrm   (1:N))    ; vKinPrm=   0
  !
ENDSUBROUTINE Dynam_Alloc_KinFas

SUBROUTINE Dynam_Init_KinFas_Texture(nMk,VBox,PhiF)
  USE M_Numeric_Const,ONLY: Pi
  USE M_KinFas_Surf,  ONLY: KinFas_Surf_Init
  !
  USE M_Global_Vars,  ONLY: vFas,vKinFas
  !
  USE M_Dynam_Vars,   ONLY:    &
  & PhiF0,vKinMinim,           &
  & vMolK,vSurfK,vMolarVol,    &
  & vMolK0,vSurfK0,            &
  & RadiusMinim,DensitMinim
  !
  INTEGER, INTENT(IN):: nMk
  REAL(dp),INTENT(IN):: VBox,PhiF
  !
  REAL(dp):: VolMinim,SurfMinim
  INTEGER :: I
  !
  !------------------------------------------------------INIT.MINERALS--
  !
  vMolarVol(1:nMk)= vFas(vKinFas(1:nMk)%iFas)%VolM3 !-> molar volumes of phases
  !
  vKinFas(:)%Dat%PhiM= (One - PhiF) *vKinFas(:)%Dat%PhiM !-> vol.fraction of box
  !
  !--------------------------------------------------- init vKinMinim --
  !
  !! !test calculation radius of "cluster" of 6 molecules -> results around 4.D-9 m
  !! vKinMinim(1:nMk)= 36.D-23 != approx. 6 molecules
  !! DO i=1,nMk
  !!   VolMinim=    vKinMinim(I) * vFas(vKinFas(I)%iFas)%V !-> volume of six molecules
  !!   RadiusMinim= ( 3.D0*VolMinim /4.D0/Pi) **(1.D0/3.D0)
  !!   WRITE(*,'(A15,1X,2G15.6,1X,A)') vKinFas(i)%Name, VolMinim, RadiusMinim, vKinFas(I)%cSatur
  !! ENDDO
  !! CALL Pause_
  !
  SurfMinim= DensitMinim *Pi*4._dp       *RadiusMinim**2  ! per m3
  VolMinim=  DensitMinim *Pi*4._dp/3._dp *RadiusMinim**3  ! per m3
  !
  IF(iDebug>2) THEN
    PRINT '(A,G15.6)', "SurfMinim=",SurfMinim
    PRINT '(A,G15.6)', "VolMinim =",VolMinim
    CALL Pause_
  ENDIF
  !
  !--- mimimal mole number in the box
  vKinMinim(1:nMk)= VolMinim*Vbox /vMolarVol(1:nMk)
  !
  IF(iDebug>2) THEN
    DO i=1,nMk
      PRINT '(A15,2A,2(A,G15.6))', &
      & vKinFas(i)%NamKF," Sat=",vKinFas(I)%Dat%cSat,&
      & " PhiM=",vKinFas(I)%Dat%PhiM," Minim=",vKinMinim(I)
    ENDDO
    CALL Pause_ !; STOP
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') "Initial Surfaces of Minerals"
  !
  DO i=1,nMk
    !
    CALL KinFas_Surf_Init( &
    & vFas, &
    & vKinFas(i),                & !in
    & vKinFas(i)%Dat%PhiM,       & !in,  volume_phase / volume_box
    & vKinFas(i)%Dat%Surf)         !out, surface_phase / volume_box
    !
    !! vKinFas(i)%Dat%Surf= MAX(vKinFas(i)%Dat%Surf, SurfMinim)
    !
  ENDDO
  !
  !!WHERE(TRIM(vKinFas%cSatur)=="2") vMolK(:)=vKinMinim(:)  !-> mole number minerals in box
  !!WHERE(TRIM(vKinFas%cSatur)=="1") vMolK(:)=VBox* vKinFas(:)%PhiM /vMolarVol(:)
  !!WHERE(vKinFas%Dat%cSat=="I")     vMolK(:)= VBox *vKinFas(:)%Dat%PhiM /vMolarVol(:)
  !
  vMolK(:)=  vKinFas(:)%Dat%PhiM *Vbox/vMolarVol(:)
  vSurfK(:)= vKinFas(:)%Dat%Surf *vBox
  !
  vSurfK0(:)= vSurfK(:) ! update vSurfK0
  PhiF0=      PhiF      ! update PhiF0
  vMolK0(:)=  vMolK(:)  ! update vMolK0
  !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(A,G15.6)') "VBox     ", VBox
    WRITE(fTrc,'(/,A,/)') "vSurfK0=Initial Surfaces of Minerals"
    DO I=1,nMk
      WRITE(fTrc,'(A15,1X,3(A,1X,G12.5,1X))') &
      & vKinFas(i)%NamKF,&
      & "vSurfK0=",vSurfK0(i),&
      & "/Rho=",vFas(vKinFas(I)%iFas)%WeitKg /vMolarVol(I),& !mineral density
      & "/vMolK0=",vMolK0(I)
    ENDDO
  ENDIF
  !-----------------------------------------------------/INIT.MINERALS--
  !
  RETURN
ENDSUBROUTINE Dynam_Init_KinFas_Texture

SUBROUTINE Dynam_KinFas_Stoik_Alloc(nMk,nCp)
  USE M_Dynam_Vars, ONLY: tNu_Kin,tAlfKin,vDG_Kin
  !
  INTEGER,INTENT(IN):: nMk,nCp
  !
  ALLOCATE(vDG_Kin(1:nMk))
  ALLOCATE(tNu_Kin(1:nMk,1:nCp))
  ALLOCATE(tAlfKin(1:nCp,1:nMk))
  !
ENDSUBROUTINE Dynam_KinFas_Stoik_Alloc

SUBROUTINE Dynam_KinFas_Stoik_Alfa(nMk,nCp)
!--
!-- compute tAlfKin, (and tAlfPr, tAlfAs for orthogonal basis)
!--
  USE M_T_MixModel, ONLY: T_MixModel
  USE M_T_MixPhase, ONLY: T_MixPhase
  !
  USE M_Global_Vars,ONLY: vMixModel,vMixFas,vFas,vKinFas
  USE M_Basis_Vars, ONLY: tAlfSp,tAlfPr,tAlfAs
  USE M_Dynam_Vars, ONLY: tNu_Kin,tAlfKin
  !
  INTEGER,INTENT(IN):: nMk,nCp
  !
  TYPE(T_MixModel):: SM
  TYPE(T_MixPhase):: SF
  INTEGER:: I,J,K,iCp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_KinFas_Stoik_Alfa"
  !
  !------------------------------------------- update stoichio tables --
  DO I=1,nMk
  
    J=vKinFas(I)%iFas
    
    IF(vFas(J)%iSpc>0) THEN
      !
      K=  vFas(J)%iSpc
      tAlfKin(1:nCp,I)= tAlfSp(1:nCp,K)
      !
    ELSE IF(vFas(J)%iMix>0) THEN
      !
      K=  vFas(J)%iMix
      SF= vMixFas(K)
      SM= vMixModel(SF%iModel)
      !
      DO iCp=1,nCp
        tAlfKin(iCp,I)= &
        & DOT_PRODUCT( tAlfSp(iCp,SM%vIPole(1:SM%NPole)) , SF%vXPole(1:SM%NPole) )
      ENDDO
      !
    END IF

    !! SELECT CASE(TRIM(vFas(J)%Typ))
    !! CASE("PURE") !,"DISCRET")
    !!   K=  vFas(J)%iSpc
    !!   tAlfKin(1:nCp,I)= tAlfSp(1:nCp,K)
    !! CASE("MIXT")
    !!   K=  vFas(J)%iMix
    !!   SF= vMixFas(K)
    !!   SM= vMixModel(SF%iModel)
    !!   !
    !!   DO iCp=1,nCp
    !!     tAlfKin(iCp,I)= &
    !!     & DOT_PRODUCT( tAlfSp(iCp,SM%vIPole(1:SM%NPole)) , SF%vXPole(1:SM%NPole) )
    !!   ENDDO
    !! END SELECT
    
  ENDDO
  !------------------------------------------/ update stoichio tables --
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_KinFas_Stoik_Alfa"
  !
ENDSUBROUTINE Dynam_KinFas_Stoik_Alfa

SUBROUTINE Dynam_KinFas_Stoik_NuKin(nMk)
!--
!-- compute tNu_Kin
!--
  USE M_T_MixModel, ONLY: T_MixModel
  USE M_T_MixPhase, ONLY: T_MixPhase
  !
  USE M_Global_Vars,ONLY: vMixModel,vMixFas,vFas,vKinFas
  USE M_Basis_Vars, ONLY: tNuSp
  USE M_Dynam_Vars, ONLY:tNu_Kin
  !
  INTEGER,INTENT(IN):: nMk
  !
  TYPE(T_MixModel):: SM
  TYPE(T_MixPhase):: SF
  INTEGER:: I,J,K,iCp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_KinFas_Stoik_Nu"
  !
  !------------------------------------------- update stoichio tables --
  DO I=1,nMk
  
    J=vKinFas(I)%iFas

    IF(vFas(J)%iSpc>0) THEN
      K=  vFas(J)%iSpc
      tNu_Kin(I,:)= tNuSp(K,:)
      !
    ELSE IF(vFas(J)%iMix>0) THEN
      K=  vFas(J)%iMix
      SF= vMixFas(K)
      SM= vMixModel(SF%iModel)
      !
      DO iCp=1,SIZE(tNu_Kin,2)
        tNu_Kin(I,iCp)= &
        & DOT_PRODUCT( tNuSp(SM%vIPole(1:SM%NPole),iCp) ,  SF%vXPole(1:SM%NPole) )
      ENDDO
      !
    END IF

    !! SELECT CASE(TRIM(vFas(J)%Typ))
    !! CASE("PURE") !,"DISCRET")
    !!   K=  vFas(J)%iSpc
    !!   tNu_Kin(I,:)= tNuSp(K,:)
    !! CASE("MIXT")
    !!   K=  vFas(J)%iMix
    !!   SF= vMixFas(K)
    !!   SM= vMixModel(SF%iModel)
    !!   !
    !!   DO iCp=1,SIZE(tNu_Kin,2)
    !!     tNu_Kin(I,iCp)= &
    !!     & DOT_PRODUCT( tNuSp(SM%vIPole(1:SM%NPole),iCp) ,  SF%vXPole(1:SM%NPole) )
    !!   ENDDO
    !! END SELECT
    
  ENDDO
  !------------------------------------------/ update stoichio tables --
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_KinFas_Stoik_Nu"
  !
ENDSUBROUTINE Dynam_KinFas_Stoik_NuKin

SUBROUTINE Dynam_KinFas_DG_UpDate(vDG_Sp,tNu_Kin,vDG_Kin)

  USE M_Global_Vars,ONLY: vSpc,vFas,vKinFas
  USE M_Basis_Vars, ONLY: vOrdPr

  REAL(dp),INTENT(IN) :: vDG_Sp(:)
  REAL(dp),INTENT(IN) :: tNu_Kin(:,:)
  REAL(dp),INTENT(OUT):: vDG_Kin(:)
  !
  INTEGER:: nMk,nCp,I,J

  nCp= SIZE(vOrdPr)
  nMk= SIZE(vKinFas)
  !
  !------------------------------------------- update thermo' param's --
  DO I=1,nMk
    J=vKinFas(I)%iFas
    !
    IF(vFas(J)%iSpc>0) THEN
      vDG_Kin(I)= vDG_Sp(vFas(J)%iSpc)
    ELSE IF(vFas(J)%iMix>0) THEN
      vDG_Kin(I)= &
      & vFas(vFas(J)%iMix)%Grt &
      & - DOT_PRODUCT(tNu_Kin(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
    END IF
    !
  ENDDO
  !------------------------------------------/ update thermo' param's --

  RETURN
ENDSUBROUTINE Dynam_KinFas_DG_UpDate

SUBROUTINE Dynam_KinModel_Init(TimeFactor,TdgK)

  USE M_Global_Vars,ONLY: vKinModel,vSpc
  USE M_T_Kinmodel, ONLY: KinModel_Show
  USE M_Basis_Vars, ONLY: vPrmBk,vPrmFw
  USE M_Dynam_Vars, ONLY: vKinMod !,
  !
  REAL(dp),INTENT(IN):: TimeFactor,TdgK

  CALL KinModel_Prm_Update( &
  & vPrmBk,vKinModel, &
  & vKinMod)
  !
  CALL KinModel_T_Update(TimeFactor,TdgK,vKinMod)
  !
  IF(iDebug>0) CALL KinModel_Show(fTrc,vSpc,vKinMod,vPrm=vPrmFw)

  RETURN
ENDSUBROUTINE Dynam_KinModel_Init

SUBROUTINE Dynam_TP_Alloc

  USE M_Global_Vars,ONLY: vSpc
  USE M_Dynam_Vars,ONLY: vDG_Sp,vDG_As,vDLGam_As
  USE M_Basis_Vars,ONLY: nAs
  !
  IF(ALLOCATED(vDG_Sp)) DEALLOCATE(vDG_Sp); ALLOCATE(vDG_Sp(1:SIZE(vSpc)))
  !! nAs= COUNT(vSpc%Typ(1:3)=="AQU")
  IF(nAs>0) THEN
    IF(ALLOCATED(vDG_As))    DEALLOCATE(vDG_As);    ALLOCATE(vDG_As(1:nAs))
    IF(ALLOCATED(vDLGam_As)) DEALLOCATE(vDLGam_As); ALLOCATE(vDLGam_As(1:nAs))
  ENDIF
  !
ENDSUBROUTINE Dynam_TP_Alloc

SUBROUTINE Dynam_Global_TP_Update(TdgK,Pbar)
!--
!-- calc' global thermo parameters of species, sol'models and sol'phases for (TdgK,Pbar)
!-- and calc' local deltaG's (vDG_Sp,vDG_As)
!--

  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  USE M_Global_Tools,  ONLY: Global_TP_Update

  !--vars--
  USE M_Global_Vars, ONLY: vSpcDtb,vSpc,vMixModel,vDiscretModel,vDiscretParam
  USE M_Global_Vars, ONLY: vMixFas,vFas
  USE M_Global_Vars, ONLY: SolModel
  USE M_System_Vars, ONLY: System_Type
  USE M_Basis_Vars,  ONLY: nAs,tNuSp,tNuAs,vOrdPr,vOrdAs
  USE M_Dynam_Vars,  ONLY: vDG_Sp,vDG_As
  !--/

  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  INTEGER:: nCp

  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  !<new 200911
  IF(System_Type=="AQUEOUS") CALL SolModel_TP_Update(TdgK,Pbar,SolModel)
  !</new 200911
  !
  nCp= SIZE(vOrdPr)
  !
  vDG_Sp(:)= vSpc(:)%G0rt - MATMUL(tNuSp(:,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
  !
  IF(nAs>0) &
  & vDG_As(:)= vSpc(vOrdAs(:))%G0rt &
  &          - MATMUL(tNuAs(:,1:nCp), vSpc(vOrdPr(1:nCp))%G0rt)

ENDSUBROUTINE Dynam_Global_TP_Update

SUBROUTINE Dynam_Init_KinFas_SatState
!-> vKinFas(:)%Dat%QsK, vKinFas(:)%Dat%cSatur
  USE M_IOTools
  USE M_Numeric_Const,ONLY: Ln10
  USE M_KinRate,    ONLY: KinRate_SatState,KinRate_CalcQsK
  !
  USE M_Global_Vars,ONLY: vSpc,vKinFas,vFas
  USE M_Basis_Vars, ONLY: vOrdPr,nCi
  !
  USE M_Dynam_Vars, ONLY: vDG_Kin,tNu_Kin,vKinMinim,Qsk_Iota
  USE M_Dynam_Vars, ONLY: vLnAct,vLnGam
  USE M_Dynam_Vars, ONLY: vMolF,vMolK
  !
  INTEGER :: i,nCp,nAq,nCx
  REAL(dp):: QsK,QskIota
  CHARACTER:: cStatus
  REAL(dp),DIMENSION(1:SIZE(vOrdPr)):: dQsKdLnXi
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_Init_KinFas_SatState"
  !
  QskIota= Qsk_Iota
  !
  nCp= SIZE(vOrdPr)
  nCx= nCp-nCi
  nAq= COUNT(vSpc%Typ=="AQU")
  !
  IF(iDebug>0) THEN
    !
    WRITE(fTrc,'(A)') "Log10(QsK) of Mineral species"
    WRITE(fTrc,'(/,A,/)') "iCp, Log(Act), LogK0, Prim.Species"
    DO i=1,nCp
      WRITE(fTrc,'(I2,A1,2(G15.3,A1),A)') &
      & I,                              T_, &
      & vSpc(vOrdPr(i))%Dat%LAct/Ln10,  T_, &
      & vSpc(vOrdPr(i))%G0rt/Ln10,      T_, &
      & TRIM(vSpc(vOrdPr(i))%NamSp)
    ENDDO
    !
    WRITE(fTrc,'(/,A,/)') "iMk, LogK0, Mineral"
    DO i=1,SIZE(vKinFas)
      WRITE(fTrc,'(I2,A1, G15.3,A1, A)') &
      & I,T_, vFas(vKinFas(i)%iFas)%Grt/Ln10,T_,  TRIM(vFas(vKinFas(i)%iFas)%NamFs)
    ENDDO
    !
  ENDIF
  !
  !! CALL OutStrVec(51,vSpc(:)%Dat%LAct)
  !
  DO I=1,SIZE(vKinFas)
    !------------------------ calc' QsK, saturation state of minerals --
    CALL KinRate_CalcQsK( &
    & nCi,nCx,      &
    & vDG_Kin(I),   &
    & tNu_Kin(I,:), &
    & LOG(vMolF(1:nAq)), & !IN:  Ln(aqu.species mole nrs)
    & vLnGam,       & !IN:  Ln(ActCoeff,aqu.species)
    & vLnAct,       & !IN:  Ln(ActivitySpecies), USEd for buffered species
    & vLnAct(1),    & !IN:  Ln(ActivitySolvent)
    & QsK,          & !OUT, mineral saturation
    & dQsKdLnXi)      !OUT
    !
    !--------- from QsK deduce cSatur, for branching to dissol/precip --
    CALL KinRate_SatState( &
    & vKinFas(I),     & !IN:  mineral: for cSatur, QsKseuil
    & vMolK(I),       & !IN:  Nr Moles of mineral, USEd for checking nMol<=MolMinim
    & vKinMinim(I),   & !IN:
    & QsK,            & !IN:
    & QskIota,        & !IN
    & cStatus)          !OUT
    !
    !----------------------------------------- save in vKinFas(:)%Dat --
    vKinFas(i)%Dat%QsK=  QsK
    vKinFas(I)%Dat%cSat= cStatus
    !
  ENDDO
  !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(/,A,/)') "Minerals Saturation State, logQsK"
    DO i=1,SIZE(vKinFas)
      WRITE(fTrc,'(I2,A1,G15.3,A1,2A)') &
      & I,T_, &
      & LOG(vKinFas(i)%Dat%QsK)/Ln10,T_, &
      & vKinFas(i)%NamKF,vKinFas(i)%Dat%cSat
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_Init_KinFas_SatState"
  !
ENDSUBROUTINE Dynam_Init_KinFas_SatState

SUBROUTINE KinModel_Prm_Update(vPrmBk,vKinModel,vKinMod) !
!--- update indexes of species in the current kinetic models vKinMod  --
  USE M_T_Kinmodel,ONLY: T_Kinmodel,KinModel_PrmIndex,KinModel_CalcCoeffs
  !
  INTEGER,         INTENT(IN) :: vPrmBk(:)
  TYPE(T_Kinmodel),INTENT(IN) :: vKinModel(:)
  TYPE(T_Kinmodel),INTENT(OUT):: vKinMod(:)
  !
  INTEGER ::I
  !
  vKinMod= vKinModel
  !vKinModel is in "global" database, vKinMod is a copy in Dynam_Vars
  !
  !calc. indexes of species in the current kinetic models vKinMod
  DO I=1,SIZE(vKinMod)
    CALL KinModel_PrmIndex(vPrmBk,vKinModel(I),vKinMod(I))
  ENDDO
  !
ENDSUBROUTINE KinModel_Prm_Update

SUBROUTINE KinModel_T_Update(TimeFactor,TdgK,vKinMod)
!------------ update temperature dependent param's of kinetic models  --
  USE M_T_Kinmodel,ONLY: T_Kinmodel,KinModel_CalcCoeffs
  !
  REAL(dp),        INTENT(IN)   :: TimeFactor,TdgK
  TYPE(T_Kinmodel),INTENT(INOUT):: vKinMod(:)
  !
  INTEGER:: I
  !
  DO I=1,SIZE(vKinMod)
    CALL KinModel_CalcCoeffs(vKinMod(I),TimeFactor,TdgK)
  ENDDO
  !
  RETURN
ENDSUBROUTINE KinModel_T_Update

SUBROUTINE NuTable_Sho(vCpn)
  USE M_IoTools,    ONLY: GetUnit
  USE M_Files,      ONLY: DirOut,Files_Index_Write
  USE M_T_Component,ONLY: T_Component
  USE M_Global_Vars,ONLY: vSpc,vFas,vKinFas
  USE M_Dynam_Vars, ONLY: tNu_Kin
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  !
  INTEGER::F,iPr,I
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_stoik_nu_kin.log")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_stoik_nu_kin.log",&
  & "Nu Table: stoikio of kin'species in terms of prim'species")
  !
  DO iPr=1,SIZE(tNu_Kin,2)
    WRITE(F,'(A7,A1)',ADVANCE='NO') TRIM(vSpc(vCpn(iPr)%iSpc)%NamSp),T_
  ENDDO
  WRITE(F,'(A)') "Name"
  !
  DO I=1,SIZE(tNu_Kin,1)
    DO iPr=1,SIZE(tNu_Kin,2)
      WRITE(F,'(F7.2,A1)',ADVANCE='NO') tNu_Kin(I,iPr), T_
    ENDDO
    WRITE(F,'(A)') TRIM(vFas(vKinFas(I)%iFas)%NamFs)
  ENDDO
  !
  CLOSE(F)
  !
ENDSUBROUTINE NuTable_Sho

SUBROUTINE Dynam_ReStart_Time(T)
  USE M_Dynam_Vars
  !
  TYPE(T_DynTime),INTENT(IN):: T
  
  TUnit=  T%TUnit
  Time=   T%Time
  dTime=  T%dTime
  TFinal= T%TFinal
  dTMin=  T%dTMin
  dTMax=  T%dTMax
  dTSav=  T%dTSav
  !
  !!--new 200912--
  IF(TimeIsSeconds) THEN
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
  ELSE
    !
    TimeFactor= T%TimeFactor
    TimeScale=  One !-----------------------TimeScale is used for output
    !
  ENDIF
  !!--new 200912--/
  
  IF(iDebug>2) THEN
    PRINT '(3(G12.3,A))', &
    & T%dTime,"=dTime /",T%dTMin,"=dTMin /",T%dTMax,"=dTMax /"
    CALL Pause_
  ENDIF
  
  RETURN
ENDSUBROUTINE Dynam_ReStart_Time
  !
SUBROUTINE Dynam_ReStart_Box(B,TimeFactor_)
  USE M_Dynam_Vars,ONLY: T_DynBox
  USE M_Dynam_Vars,ONLY: dX,UDarcy,VBox,FOut,PhiF,RhoF,UpdateMassFluid
  USE M_Dynam_Vars,ONLY: TimeIsSeconds
  USE M_Dynam_Vars,ONLY: Dynam_nTotalNewtIter
  !
  TYPE(T_DynBox), INTENT(IN):: B
  REAL(dp),       INTENT(IN):: TimeFactor_
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
  IF(TimeIsSeconds) UDarcy= UDarcy /TimeFactor_
  !!</new 200912>
  !
  !! nCell=   B%nCell
  Dynam_nTotalNewtIter = 0
ENDSUBROUTINE Dynam_ReStart_Box

SUBROUTINE Dynam_Zero_Time(T)
!.default dynamic parameters, in case not found in Input
  USE M_Dynam_Vars,ONLY: T_DynTime
  TYPE(T_DynTime),INTENT(OUT):: T
  !
  !time parameters
  T%TUnit=  "DAY"  !default unit used for all time data
  T%TFinal= 5.0D+6 !duration of simulation, in TUnit
  T%dTime=  1.0D-6 !initial time step
  T%dTMin=  1.0D-16 !minimal time step
  T%dTMax=  Zero   !maximal time step
  T%dTSav=  Zero   !delta time between SAVEs
  T%TimeFactor= Dynam_TimeFactor(T%TUnit)
  !
ENDSUBROUTINE Dynam_Zero_Time

REAL(dp) FUNCTION Dynam_TimeFactor(S) RESULT(X)
!--
!-- Time Factor --
!--
  CHARACTER(LEN=*),INTENT(IN):: S
  
  X=One
  SELECT CASE(TRIM(S))
    CASE("SECOND")  ;  X= One
    CASE("MINUTE")  ;  X= 60._dp
    CASE("HOUR")    ;  X= 3600._dp
    CASE("DAY")     ;  X= 3600._dp * 24._dp
    CASE("YEAR")    ;  X= 3600._dp * 24._dp * 365._dp
  END SELECT
  
END FUNCTION Dynam_TimeFactor

SUBROUTINE Dynam_Zero_Box(B)
!--
!-- initialize box with default values --
!--
  USE M_Dynam_Vars,ONLY: T_DynBox, VBox0
  TYPE(T_DynBox),INTENT(OUT):: B
  !
  !--- box parameters
  B%UDarcy=Zero   ! flux rate, metre/time
  B%VBox=  One    ! volume of box, m3 (just scaling parameter (?), useful for "outside")
  B%dX=    One    ! length of box, meter
  B%PhiF=  0.5D0  ! porosity
  B%RhoF=  1.D3   ! fluid density, kg/m^3
  B%FOut=  Zero   !
  B%nCell= 1      !
  B%VFixed=.TRUE. ! Vfixed is .false. for free volume
  
  RETURN
ENDSUBROUTINE Dynam_Zero_Box

SUBROUTINE Dynam_Zero_Numeric
  USE M_Dynam_Vars
  !
  !--debug parameters--
  DebNewt=   .FALSE.
  DebJacob=  .FALSE.
  TestJacob= .FALSE.
  TestMax=   .FALSE.
  bFinDIF=   .FALSE.  !.true. !
  !
  !--calculation mode parameters--
  UpdateMassFluid= .TRUE.
  !
  ! provisionally, dVmAdLnX_M is not computed in this new
  ! -> we assume that this factor will not be implicited
  Implicit_ActivFactor=.FALSE.
  !
  Implicit_Surface= .FALSE. ! .TRUE. !
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
  !-- IF(Its>NewtMaxIts)  decrease dtime, cycle Newton
  !-- IF(Its>NewtIterMax) decrease dtime (outside Newton)
  !-- IF(Its-NewtIterMax) increase dtime (outside Newton)
  NewtMaxIts=  40 !75 !30  !100 !
  NewtIterMax= 20 !50 !16  !60  !
  NewtIterMin= 10 !25 !8   !30  !
  !--/
  !
  NewtTolF=       1.0E-09_dp  ! convergence criterion on function values
  !NewtTolF=      1.0E-07_dp  ! convergence criterion on function values
  !NewtTolF_Equil= NewtTolF
  NewtTolF_Equil= 1.0E-06_dp  !
  !
  NewtTolMin=     1.0E-06_dp  !criterion for spurious convergence
  !
  NewtTolX= EPSILON(REAL(dp)) !convergence criterion on dx -> approx. 2.E-16
  !NewtTolX=   1.0E-09_dp  !convergence criterion on dx
  !
  ! usage in Newton:
  !   NewtTolF: IF( MAXVAL(ABS(fVec(:)))<NewtTolF ) Newton_Ok
  !   NewtTolX: IF( MAXVAL( ABS(vX(:)-vXOld(:))/MAX(ABS(vX(:)),One)<NewtTolX ) Converge on X
  !
  !--------------------------------------------------/NUMERIC parameters
  
  RETURN
ENDSUBROUTINE Dynam_Zero_Numeric

SUBROUTINE Dynam_InitCalc_Restart
  USE M_Global_Vars,ONLY: nAq,vSpc,vFas,vKinFas
  USE M_Basis_Vars, ONLY: vOrdAq
  USE M_Dynam_Vars, ONLY: vMolK,vMolF
  USE M_Dynam_Vars, ONLY: vMolarVol,vKinQsk,VBox,PhiF
  USE M_Dynam_Vars, ONLY: Time, Tfinal, Dtime, Dtmin, Dtmax
  USE M_Dynam_Vars, ONLY: Dynam_nStep, Dynam_nNewtIter, Dynam_nTotalNewtIter
  !
  INTEGER:: nMk
  
  nMk= SIZE(vKinFas)
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
  vKinQsk(1:nMk)=  vKinFas(1:nMk)%Dat%QsK
  !
  PhiF= One - SUM(vKinFas(1:nMk)%Dat%PhiM)
  !             !-> vol'fraction of fluid
  vMolF(1:nAq)= vSpc(vOrdAq(1:nAq))%Dat%Mole
  !             !-> mole nr's of aqu'species in box
  !
  IF (iDebug>0) &
  & PRINT '(A,3F12.5, 2G12.3)', &
  & 'Restart avec [Time, Tfinal, Dtime, Dtmin, Dtmax =]', &
  & Time, Tfinal, Dtime, Dtmin, Dtmax
  
  RETURN
ENDSUBROUTINE Dynam_InitCalc_Restart

SUBROUTINE Dynam_Restart
  USE M_Global_Vars,ONLY: nAq,vSpc,vFas,vKinFas
  USE M_Basis_Vars, ONLY: vOrdAq
  !
  USE M_Dynam_Vars, ONLY: vMolarVol,vMolK,vSurfK,vKinQsk,vStatusK
  USE M_Dynam_Vars, ONLY: VBox,PhiF,vMolF
  USE M_Dynam_Vars, ONLY: Time, Tfinal, Dtime, Dtmin, Dtmax
  !
  INTEGER:: nMk
  
  nMk= SIZE(vKinFas)
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
  PhiF= One - SUM(vKinFas(1:nMk)%Dat%PhiM)
  !                 !-> vol'fraction of fluid
  vMolF(1:nAq)= vSpc(vOrdAq(1:nAq))%Dat%Mole
  !                 !-> mole nr's of aqu'species in box

  IF (iDebug>0) &
  & PRINT '(A,3F12.5, 2G12.3)', &
  & 'Restart avec [Time, Tfinal, Dtime, Dtmin, Dtmax =]', &
  & Time, Tfinal, Dtime, Dtmin, Dtmax
  
  RETURN
ENDSUBROUTINE Dynam_Restart

SUBROUTINE Dynam_Save

  USE M_Global_Vars,ONLY: vSpc,vFas,vKinFas
  !
  USE M_Basis_Vars, ONLY: vOrdAq
  USE M_Dynam_Vars, ONLY: vCpnBox
  USE M_Dynam_Vars, ONLY: vTotF,vMolF,vLnAct,vLnGam,vLnBuf
  USE M_Dynam_Vars, ONLY: vMolK,vSurfK,vKinQsk,vStatusK,vMolarVol,VBox
  !
  USE M_Dynam_Cell
  
  TYPE(T_Cell_State):: Cell
  !
  INTEGER :: nCp,nMk,nAq
  REAL(dp):: PhiF
  
  nCp= SIZE(vCpnBox)
  nAq= COUNT(vSpc(:)%Typ(1:3)=="AQU")
  nMk= SIZE(vKinFas)
  !
  vCpnBox(:)%Mole=  vTotF(:)
  vCpnBox(:)%LnAct= vSpc(vCpnBox(:)%iSpc)%Dat%LAct
  !
  vSpc(vOrdAq(1:nAq))%Dat%Mole= vMolF(1:nAq)
  vSpc(vOrdAq(1:nAq))%Dat%LAct= vLnAct(1:nAq)
  vSpc(vOrdAq(1:nAq))%Dat%LGam= vLnGam(1:nAq)
  !
  vKinFas(1:nMk)%Dat%PhiM=   vMolK(1:nMk) *vMolarVol(1:nMk) / VBox
  vKinFas(1:nMk)%Dat%Surf=   vSurfK(1:nMk) /Vbox
  vKinFas(1:nMk)%Dat%QsK=    vKinQsk(1:nMk)
  vKinFas(1:nMk)%Dat%cSat=   vStatusK(1:nMk)
  !
  vFas(vKinFas(1:nMk)%iFas)%Mole= vMolK(1:nMk)
  
  PhiF= One - SUM(vKinFas(1:nMk)%Dat%PhiM)
  
  !! CALL Cell_State_New(Cell, nAq,nMk)
  CALL Cell%New_(nCp,nAq,nMk)
  !! CALL Cell_State_Write( &  !
  CALL Cell%Set_( &  !
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
  
  CALL Cell%Get_( &  !
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
  
  !!PRINT *,"vMolK(1)",vMolK(1)  ;  PAUSE
  
  RETURN
ENDSUBROUTINE Dynam_Save

SUBROUTINE Dynam_Save_ToFile
  USE M_IOTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_KinFas,ONLY: T_KinFas
  !
  USE M_Global_Vars,ONLY: vEle,vSpc,vKinFas,vKinModel,vFas
  USE M_Dynam_Vars, ONLY: vCpnBox,vCpnInj
  USE M_Dynam_Vars, ONLY: PhiF
  !
  TYPE(T_KinFas):: K
  INTEGER :: F, I
  REAL(dp):: X
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_restart.inn")

  !-----------------------------------------------------------SYSTEM.BOX
  WRITE(F,'(A)') "SYSTEM.BOX"
  !! print *,vCpnBox(1)%Mole  ;  pause
  X= 55.51D0 /vCpnBox(1)%Mole
  !
  DO I=1,SIZE(vCpnBox)
    WRITE(F,'(A)',   ADVANCE="NO") "  "
    WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vEle(vCpnBox(I)%iEle)%NamEl)
    WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vCpnBox(I)%Statut)
    WRITE(F,'(A,1X)',ADVANCE="NO") TRIM(vSpc(vCpnBox(I)%iSpc)%NamSp)
    !
    SELECT CASE(TRIM(vCpnBox(I)%Statut))
      CASE("INERT")   ;  WRITE(F,'(G15.6)') vCpnBox(I)%Mole *X
      CASE("BUFFER")  ;  WRITE(F,'(G15.6)') - vCpnBox(I)%LnAct /Ln10
    ENDSELECT
    !
  ENDDO
  WRITE(F,'(A)') "END SYSTEM.BOX"
  !----------------------------------------------------------/SYSTEM.BOX

  !--------------------------------------------------------SYSTEM.INJECT
  WRITE(F,'(A)') "SYSTEM.INJECT"
  X= 55.51D0 /vCpnInj(1)%Mole
  !
  DO I=1,SIZE(vCpnInj)
    WRITE(F,'(A)',   ADVANCE="NO") "  "
    WRITE(F,'(A,A1)',ADVANCE="NO") TRIM(vEle(vCpnInj(I)%iEle)%NamEl),T_
    WRITE(F,'(A,A1)',ADVANCE="NO") TRIM(vCpnInj(I)%Statut),          T_
    WRITE(F,'(A,A1)',ADVANCE="NO") TRIM(vSpc(vCpnInj(I)%iSpc)%NamSp),T_
    !
    SELECT CASE(TRIM(vCpnInj(I)%Statut))
      CASE("INERT")   ;  WRITE(F,'(G15.6)') vCpnInj(I)%Mole *X
      CASE("BUFFER")  ;  WRITE(F,'(G15.6)') - vCpnInj(I)%LnAct /Ln10
    ENDSELECT
    !
  ENDDO
  WRITE(F,'(A)') "END SYSTEM.INJECT"
  !-------------------------------------------------------/SYSTEM.INJECT

  !-------------------------------------------------------- DYNAMIC.ROCK
  WRITE(F,'(A)') "DYNAMIC.ROCK"
  DO I=1,SIZE(vKinFas)
    !
    K= vKinFas(I)
    !
    WRITE(F,'(A)',   ADVANCE="NO") "  "
    WRITE(F,'(A,A1)',ADVANCE="NO") TRIM(K%NamKF),T_
    !
    WRITE(F,'(A,A1,G15.6,A1)',ADVANCE="NO") "SURFACE",T_, K%Dat%SurfKg,T_
    WRITE(F,'(A,A1,G15.6,A1)',ADVANCE="NO") "VOLUME", T_, K%Dat%PhiM,  T_
    !
    IF(K%iKin>0) THEN
      IF(TRIM(K%NamKF)/=TRIM(vKinModel(K%iKin)%Name)) &
      & WRITE(F,'(2(A,A1))',ADVANCE="NO") &
      & "KINETICS", T_, &
      & TRIM(vKinModel(K%iKin)%Name),T_
    ENDIF
    !
    IF(TRIM(K%NamKF)/=TRIM(vFas(K%iFas)%NamFs)) &
    & WRITE(F,'(2(A,A1))',ADVANCE="NO") &
    & "SPECIES", T_, &
    & TRIM(vFas(K%iFas)%NamFs),T_!
    !
    IF(K%cMode/="D") WRITE(F,'(A1,A1)',ADVANCE="NO") K%cMode, T_
    !
    WRITE(F,*)
    !
  ENDDO
  WRITE(F,'(A)') "END DYNAMIC.ROCK"
  !--------------------------------------------------------/DYNAMIC.ROCK
  !
  !--------------------------------------------------------------DYNAMIC
  WRITE(F,'(A)') "DYNAMIC"
    WRITE(F,'(A)',   ADVANCE="NO") "  "
    WRITE(F,'(A,A1,G15.6)',ADVANCE="NO") "POROSITY",T_, PhiF
    WRITE(F,*)
  WRITE(F,'(A)') "END DYNAMIC"
  !-------------------------------------------------------------/DYNAMIC
  !
  CLOSE(F)
ENDSUBROUTINE Dynam_Save_ToFile

SUBROUTINE Dynam_Clean
  USE M_Dynam_Vars
  !
  IF(ALLOCATED(vCpnInj))    DEALLOCATE(vCpnInj)
  IF(ALLOCATED(vCpnBox))    DEALLOCATE(vCpnBox)
  !
  IF(ALLOCATED(vTotF))      DEALLOCATE(vTotF)
  IF(ALLOCATED(vMolF))      DEALLOCATE(vMolF)
  IF(ALLOCATED(vTotInj))    DEALLOCATE(vTotInj)
  !
  IF(ALLOCATED(vLnAct))     DEALLOCATE(vLnAct)
  IF(ALLOCATED(vLnGam))     DEALLOCATE(vLnGam)
  IF(ALLOCATED(vLnBuf))     DEALLOCATE(vLnBuf)
  !
  IF(ALLOCATED(vMolarVol))  DEALLOCATE(vMolarVol)
  IF(ALLOCATED(vMolK))      DEALLOCATE(vMolK)
  IF(ALLOCATED(vSurfK))     DEALLOCATE(vSurfK)
  IF(ALLOCATED(vStatusK))   DEALLOCATE(vStatusK)
  !
  IF(ALLOCATED(vMolK0))      DEALLOCATE(vMolK0)
  IF(ALLOCATED(vSurfK0))     DEALLOCATE(vSurfK0)
  !
  IF(ALLOCATED(vKinMinim))  DEALLOCATE(vKinMinim)
  IF(ALLOCATED(vLKinActiv)) DEALLOCATE(vLKinActiv)
  IF(ALLOCATED(vLEquActiv)) DEALLOCATE(vLEquActiv)
  IF(ALLOCATED(vKinPrm))    DEALLOCATE(vKinPrm)
  !
  IF(ALLOCATED(vKinMod))    DEALLOCATE(vKinMod)
  !
  !IF(ALLOCATED(vTooLow))    DEALLOCATE(vTooLow)
  IF(ALLOCATED(vQTotInj))    DEALLOCATE(vQTotInj)
  !
  IF(ALLOCATED(vKinQsK))     DEALLOCATE(vKinQsK)
  IF(ALLOCATED(vVmAct))      DEALLOCATE(vVmAct)
  IF(ALLOCATED(vVmQsK))      DEALLOCATE(vVmQsK)
  !!IF(ALLOCATED(vVm0)  )     DEALLOCATE(vVm0)
  !
  IF(ALLOCATED(vDG_Sp))     DEALLOCATE(vDG_Sp)
  IF(ALLOCATED(vDG_As))     DEALLOCATE(vDG_As)
  IF(ALLOCATED(vDLGam_As))  DEALLOCATE(vDLGam_As)
  IF(ALLOCATED(vDG_Kin))    DEALLOCATE(vDG_Kin)
  !
  IF(ALLOCATED(tNu_Kin))    DEALLOCATE(tNu_Kin)
  IF(ALLOCATED(tAlfKin))    DEALLOCATE(tAlfKin)
  !
  IF(ALLOCATED(vLnBuf))     DEALLOCATE(vLnBuf)
  !
  IF(ALLOCATED(vWeitSp))    DEALLOCATE(vWeitSp)
  IF(ALLOCATED(vWeitCp))    DEALLOCATE(vWeitCp)
  !
  IF(ALLOCATED(tStoikioAqu)) DEALLOCATE(tStoikioAqu)
  IF(ALLOCATED(tStoikioKin)) DEALLOCATE(tStoikioKin)
  !
ENDSUBROUTINE Dynam_Clean

ENDMODULE M_Dynam_Tools

!! SUBROUTINE SaturStateInit_Show( &
!! & Time,PhiF,vKinFas,vSurfK0,vMolK,vMolarVol)
!!   USE M_T_KinFas,ONLY: T_KinFas
!!   !
!!   REAL(dp),                   INTENT(IN):: Time,PhiF
!!   TYPE(T_KinFas),DIMENSION(:),INTENT(IN):: vKinFas
!!   REAL(dp),      DIMENSION(:),INTENT(IN):: vSurfK0,vMolK,vMolarVol
!!   !
!!   TYPE(T_KinFas)::M
!!   INTEGER::i
!!   !
!!   WRITE(fTrc,'(/,A,/)') "!!!SaturState_Show_begin"
!!   !
!!   WRITE(fTrc,'(8(A,A1))') &
!!   & "Name",   T_,"Satur",T_,"Time", T_,"VFluid",T_, &
!!   & "Surface",T_,"QsK",  T_,"vMolK",T_,"Rate",  T_
!!   DO i=1,SIZE(vKinFas)
!!     M=vKinFas(i)
!!     WRITE(fTrc,&
!!     & '(   A,A1,A,A1,G15.9,A1,F12.3,A1,E12.3,A1,E12.6,A1,E12.3,A1,E7.2,A1)') &
!!     & M%Name,T_,M%Dat%cSat,T_,Time, T_,PhiF, T_,&
!!     & M%Dat%Surf,T_,M%Dat%QsK,T_,vMolK(i),T_,vMolarVol(i),T_
!!   ENDDO
!!   !!! DO i=1,SIZE(vKinFas)
!!   !!!   M=vKinFas(i)
!!   !!!   WRITE(fTrc,'(A,A1,G15.8, A1,G15.8, A1,G15.8, A1)') &
!!   !!!   & TRIM(M%Name),T_,M%Kd_A,T_,M%Kd_W,T_,M%Kd_B,T_
!!   !!! ENDDO
!!   WRITE(fTrc,'(/,A,/)') "!!!SaturState_Show_END"
!! ENDSUBROUTINE SaturStateInit_Show

