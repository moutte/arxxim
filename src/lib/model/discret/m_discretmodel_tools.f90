MODULE M_DiscretModel_Tools
!--
!-- tools for "discretization" of a mixture phase
!-- to an array of pure phases
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DiscretParam_Init
  PUBLIC:: DiscretSpecies_TP_Update
  PUBLIC:: DiscretSpecies_Stoikio_Calc
  !
CONTAINS

SUBROUTINE DiscretModel_Init( &
& vEle,         & !IN
& vSpc,         & !IN
& MixModel,     & !IN
& DiscretModel, & !IN
& vParam,       & !OUT
& vSpcOut)        !OUT properties of species resulting from discretization
!--
!--!!! called by DiscretParam_Init !!!
!--
!-- discretize one mixture phase, following the mixing model MixModel,
!-- into nDiscret pure species (for a binary solution);
!-- -> for each species generated,
!--    build the stoichiometric formula, e.g. CA(10)MG(01)FE(09)SI(20)O(60)/(10)
!--    build the species name from the names of the end members, e.g. DIO01HED09
!--
  USE M_T_Element, ONLY: T_Element
  USE M_T_Species, ONLY: T_Species, Species_Zero
  USE M_T_MixPhase,ONLY: T_MixPhase, MixPhase_NormCompo
  USE M_T_MixModel,ONLY: T_MixModel, MixModel_Index
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  !
  TYPE(T_Element),    INTENT(IN) :: vEle(:)
  TYPE(T_Species),    INTENT(IN) :: vSpc(:)
  TYPE(T_MixModel),   INTENT(IN) :: MixModel
  TYPE(T_DiscretModel),INTENT(IN):: DiscretModel
  !
  TYPE(T_DiscretParam),INTENT(OUT):: vParam(:)
  TYPE(T_Species),     INTENT(OUT):: vSpcOut(:)
  !
  TYPE(T_MixPhase):: P
  TYPE(T_Species) :: S
  ! CHARACTER(LEN=7):: cName
  ! INTEGER      :: nDiscret
  INTEGER :: nEl,N,I,J,K,i0,Dims
  ! INTEGER :: Z,nDiv1,nDiv2 !,nDiv3
  REAL(dp):: ntot
  LOGICAL :: Ternary
  CHARACTER(LEN=23):: DiscretName
  CHARACTER(LEN=71):: Formul
  INTEGER,ALLOCATABLE:: vStoik1(:),vStoik2(:),vStoik3(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_Init"
  !
  Ternary= (MixModel%NPole==3)
  !! MixModel= vMixModel(DiscretModel%iMix)
  Dims= DiscretModel%Dim1 +1
  !
  P%Name=      TRIM(DiscretModel%Name)
  P%iModel=    DiscretModel%iMix !!MixModel_Index(vMixModel,MixModel%Name)
  P%vLPole=    .FALSE.
  P%vXPole=    Zero
  !
  ! compute stoikio vectors of end members
  IF(iDebug>0) WRITE(fTrc,'(A)') TRIM(vSpc(MixModel%vIPole(1))%Formula)
  IF(iDebug>0) WRITE(fTrc,'(A)') TRIM(vSpc(MixModel%vIPole(2))%Formula)
  !
  IF(MixModel%NPole==1) RETURN
  !
  nEl= SIZE(vEle)
  !
  ALLOCATE(vStoik1(0:nEl))  ;  vStoik1(:)= 0
  ALLOCATE(vStoik2(0:nEl))  ;  vStoik2(:)= 0
  ALLOCATE(vStoik3(0:nEl))  ;  vStoik3(:)= 0
  !
  vStoik1(0:nEl)= vSpc(MixModel%vIPole(1))%vStoikio(0:nEl)
  vStoik2(0:nEl)= vSpc(MixModel%vIPole(2))%vStoikio(0:nEl)
  IF(Ternary) &
  & vStoik3(0:nEl)= vSpc(MixModel%vIPole(3))%vStoikio(0:nEl)
  !
  IF(Ternary) THEN  ;  i0= 2
  ELSE              ;  i0= 1
  ENDIF
  !
  N=0
  i=0
  DO
    i= i+1
    k= 0
    DO
      IF(Ternary) k= k +1
      j= Dims -i -k
      N= N+1
      !----------------------------------------------------------------!
      !
      vParam(N)%I= I
      vParam(N)%J= J
      vParam(N)%K= K
      !
      ntot= REAL(I+J+K)
      P%vLPole(1)=I>0   ; P%vXPole(1)= I /ntot
      P%vLPole(2)=J>0   ; P%vXPole(2)= J /ntot
      P%vLPole(3)=K>0   ; P%vXPole(3)= K /ntot
      !
      CALL MixPhase_NormCompo(P)
      !
      CALL Species_Zero(S)
      !
      S%Typ= vSpc(MixModel%vIPole(1))%Typ
      S%WeitKg= DOT_PRODUCT(P%vXPole(1:MixModel%NPole), &
      &         vSpc(MixModel%vIPole(1:MixModel%NPole))%WeitKg)
      !
      IF(Ternary) THEN
        CALL MixtureName( &
        & vSpc(MixModel%vIPole(1))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(2))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(3))%NamSp(1:3), &
        & I,J,K, &
        & DiscretName)
      ELSE
        CALL MixtureName( &
        & vSpc(MixModel%vIPole(1))%NamSp(1:3), &
        & vSpc(MixModel%vIPole(2))%NamSp(1:3), &
        & "xxx", &
        & I,J,K, &
        & DiscretName)
      ENDIF
      !
      S%NamSp= TRIM(DiscretName)
      !
      S%vStoikio(0:nEl)= I*vStoik1(0:nEl) + J*vStoik2(0:nEl) + K*vStoik3(0:nEl)
      !
      CALL DiscretModel_Formula_Build( &
      & vEle,S, &
      & Formul)
      S%Formula= TRIM(Formul)
      !
      vSpcOut(N)= S
      !
      IF(iDebug>0) WRITE(fTrc,'(2(A,A1))') &
      & TRIM(DiscretName),T_,TRIM(Formul),T_
      !
      !----------------------------------------------------------------!
      IF(.not. Ternary) EXIT
      IF(k==Dims-1 -i) EXIT
    ENDDO
    !pause_
    IF(i==Dims -i0) EXIT
  ENDDO
  !
  DEALLOCATE(vStoik1,vStoik2,vStoik3)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretModel_Init"
  !
  RETURN
ENDSUBROUTINE DiscretModel_Init

SUBROUTINE MixtureName(S1,S2,S3,I,J,K,S)
  CHARACTER(LEN=3),INTENT(IN):: S1,S2,S3
  INTEGER,         INTENT(IN) ::I,J,K
  CHARACTER(LEN=*),INTENT(OUT)::S
  !
  CHARACTER(LEN=3):: Str
  !
  S=""
  WRITE(Str,'(I2)') I  ; IF(I<10) Str(1:1)='0'
  S=TRIM(S)//TRIM(S1)//"-"//TRIM(ADJUSTL(Str))
  !
  S=TRIM(S)//"_"
  WRITE(Str,'(I2)') J  ; IF(J<10) Str(1:1)='0'
  S=TRIM(S)//TRIM(S2)//"-"//TRIM(ADJUSTL(Str))
  !
  IF(K==0) RETURN
  !
  S=TRIM(S)//"_"
  WRITE(Str,'(I2)') K  ; IF(K<10) Str(1:1)='0'
  S=TRIM(S)//TRIM(S3)//"-"//TRIM(ADJUSTL(Str))
  !
  RETURN
ENDSUBROUTINE MixtureName

SUBROUTINE DiscretSpecies_Stoikio_Calc( &
& vEle,          & !IN
& vMixModel,     & !IN
& vDiscretModel, & !IN
& vDiscretParam, & !IN
& vSpc)            !INOUT with updated stoichio's for discretized species
!--
!-- update the stoichiometries
!-- of all discretized mixture phases
!-- according to models read from vDiscretModel,
!--
  USE M_T_Element, ONLY: T_Element
  USE M_T_Species, ONLY: T_Species
  USE M_T_MixPhase,ONLY: T_MixPhase
  USE M_T_MixModel,ONLY: T_MixModel, MixModel_Index
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  !
  TYPE(T_Element),     INTENT(IN):: vEle(:)
  TYPE(T_MixModel),    INTENT(IN):: vMixModel(:)
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  TYPE(T_DiscretParam),INTENT(IN):: vDiscretParam(:)
  !
  TYPE(T_Species),  INTENT(INOUT):: vSpc(:)
  !
  TYPE(T_DiscretModel):: DsModl
  TYPE(T_MixModel):: MxModl
  ! TYPE(T_MixPhase):: MxPhas
  !
  LOGICAL :: Ternary
  INTEGER :: I,J,K
  INTEGER :: iSp,iMix0,iMix,iDis,nEl
  !
  INTEGER,ALLOCATABLE:: vStoik1(:),vStoik2(:),vStoik3(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretSpecies_Stoikio_Calc"
  !
  nEl= SIZE(vEle)
  !
  iMix0= 0
  !
  ALLOCATE(vStoik1(0:nEl))  ;  vStoik1(:)= 0
  ALLOCATE(vStoik2(0:nEl))  ;  vStoik2(:)= 0
  ALLOCATE(vStoik3(0:nEl))  ;  vStoik3(:)= 0
  !
  DO iSp=1,SIZE(vSpc)
    !
    IF(vSpc(iSp)%iDiscret==0) CYCLE
    !
    iDis=    vSpc(iSp)%iDiscret
    DsModl=  vDiscretModel(vDiscretParam(iDis)%iModel)
    iMix=    DsModl%iMix
    MxModl=  vMixModel(iMix)
    !nDiscret= DsModl%Dim1
    !
    IF(iMix/=iMix0) THEN !> change parameters only when needed
      !
      Ternary= (MxModl%NPole==3)
      !
      vStoik1(0:nEl)= vSpc(MxModl%vIPole(1))%vStoikio(0:nEl)
      vStoik2(0:nEl)= vSpc(MxModl%vIPole(2))%vStoikio(0:nEl)
      IF(Ternary) &
      & vStoik3(0:nEl)= vSpc(MxModl%vIPole(3))%vStoikio(0:nEl)
      !
      IF(iDebug>0) &
      WRITE(fTrc,'(4(2A,/))') &
      & "disc.model=", DsModl%Name, &
      & "mixt.model=", MxModl%Name, &
      & "pole 1=    ", vSpc(MxModl%vIPole(1))%Formula, &
      & "pole 2=    ", vSpc(MxModl%vIPole(2))%Formula
      !
      iMix0= iMix
    ENDIF
    !
    I= vDiscretParam(iDis)%I
    J= vDiscretParam(iDis)%J
    K= vDiscretParam(iDis)%K
    !
    vSpc(iSp)%vStoikio(0:nEl)= &
    & I*vStoik1(0:nEl) + J*vStoik2(0:nEl) + K*vStoik3(0:nEl)
    !
    !~ CALL DiscretModel_Formula_Build( &
    !~ & vEle,Dims-2,vStoik1(0),vStoik2(0),I,J, &
    !~ & vStoik1(1:nEl),vStoik2(1:nEl), &
    !~ & Formul)
    !~ S%Formula= TRIM(Formul)
    !
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretSpecies_Stoikio_Calc"
  !
  RETURN
ENDSUBROUTINE DiscretSpecies_Stoikio_Calc

SUBROUTINE DiscretParam_Init( &
& vEle,vSpc,vMixModel,vDiscretModel, &
& vDiscretParam,vSpcOut)
!--
!--initialize vDiscretParam
!--
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  USE M_T_Element, ONLY: T_Element
  USE M_T_Species, ONLY: T_Species
  USE M_T_MixModel,ONLY: T_MixModel
  !
  TYPE(T_Element),     INTENT(IN):: vEle(:)
  TYPE(T_Species),     INTENT(IN):: vSpc(:)
  TYPE(T_MixModel),    INTENT(IN):: vMixModel(:)
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  !
  TYPE(T_DiscretParam),INTENT(OUT):: vDiscretParam(:)
  TYPE(T_Species),     INTENT(OUT):: vSpcOut(:)
  !
  TYPE(T_Species),     ALLOCATABLE:: vSpcTmp(:)
  TYPE(T_DiscretParam),ALLOCATABLE:: vParam(:)
  INTEGER:: I, iM, iMix, dimDis, N
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretParam_Init"
  !
  N= 0
  DO iM=1,SIZE(vDiscretModel)
    !
    dimDis= vDiscretModel(iM)%Dim1
    dimDis= vDiscretModel(iM)%DimTot
    !
    ALLOCATE(vSpcTmp(1:dimDis))
    ALLOCATE(vParam(1:dimDis))
    !
    iMix= vDiscretModel(iM)%iMix
    !
    CALL DiscretModel_Init( &
    & vEle,vSpc,vMixModel(iMix),vDiscretModel(iM), &
    & vParam,vSpcTmp)
    !
    DO I=1,dimDis
      N= N+1
      !
      vSpcOut(N)= vSpcTmp(i)
      vSpcOut(N)%iDiscret= N
      !
      vDiscretParam(N)= vParam(i)  !-> for I,J,K
      vDiscretParam(N)%iModel= iM  !must come after
      !
      IF(iDebug>0) THEN
        WRITE(fTrc,'(2A)') vSpcOut(N)%NamSp, vSpcOut(N)%Typ
      ENDIF
      !
    ENDDO
    !
    DEALLOCATE(vSpcTmp)
    DEALLOCATE(vParam)
    !
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretParam_Init"
  !
ENDSUBROUTINE DiscretParam_Init

SUBROUTINE DiscretSpecies_TP_Update( &
!--
!--.update the T,P dependent parameters
!--.of all mixture phases discretized according
!--.to models read from vDiscretModel,
!--
& vMixModel,     & !IN
& TdgK,Pbar,     & !IN
& vDiscretModel, & !IN
& vDiscretParam, & !IN
& vSpc)            !INOUT with updated parameters for discretized species
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  USE M_Dtb_Const, ONLY: T_CK,R_jk
  USE M_T_Species, ONLY: T_Species
  USE M_T_MixPhase,ONLY: T_MixPhase,MixPhase_CalcMixing,MixPhase_NormCompo
  USE M_T_MixModel,ONLY: T_MixModel, MixModel_Margul_Wg_Calc
  USE M_T_MixModel,ONLY: MixModel_Index,MixModel_XPoleToXSite
  !
  TYPE(T_MixModel),    INTENT(IN):: vMixModel(:)
  REAL(dp),            INTENT(IN):: TdgK,Pbar
  TYPE(T_DiscretModel),INTENT(IN):: vDiscretModel(:)
  TYPE(T_DiscretParam),INTENT(IN):: vDiscretParam(:)
  !
  TYPE(T_Species),  INTENT(INOUT):: vSpc(:)
  !
  TYPE(T_DiscretModel):: DsModl
  TYPE(T_MixModel):: MxModl
  TYPE(T_MixPhase):: MxPhas
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  INTEGER :: nDiscret
  REAL(dp):: G_IdMix,GMix,G_XsMix,GMecaRT,ntot !,GTotRT !,V0,WeitKg
  INTEGER :: I,J,K
  INTEGER :: iSp,iMix0,iMix,iDis
  !
  !CHARACTER(LEN=30):: cFormat
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretSpecies_TP_Update"
  !
  iMix0= 0
  !
  DO iSp=1,SIZE(vSpc)
    !
    IF(vSpc(iSp)%iDiscret==0) CYCLE
    !
    iDis=    vSpc(iSp)%iDiscret
    DsModl=  vDiscretModel(vDiscretParam(iDis)%iModel)
    !
    iMix=    DsModl%iMix
    MxModl=  vMixModel(iMix)
    nDiscret= DsModl%Dim1
    !
    !! PRINT *,vSpc(iSp)%NamSp,iDis," ",DsModl%Name,MxModl%Name
    !
    MxPhas%iModel=    iMix
    MxPhas%vLPole=    .FALSE.
    MxPhas%vXPole=    Zero
    !
    !---------------------check that parameters are not already uptodate
    IF(iMix/=iMix0) THEN
      !--------------update Margules parameters for the current T,MxPhas
      IF(MxModl%NMarg>0) &
      & CALL MixModel_Margul_Wg_Calc(TdgK,Pbar,MxModl)
      !----------------------------------------------------------------/
      iMix0= iMix
      !
      IF(iDebug>2) THEN
        WRITE(fTrc,'(2(2A,/))') &
        & "disc.model= ", TRIM(DsModl%Name), &
        & "mixt.model= ", TRIM(MxModl%Name)
        WRITE(fTrc,'(2(2A,1X,G15.6,/))') &
        & "pole 1    = ", TRIM(vSpc(MxModl%vIPole(1))%Formula),vSpc(MxModl%vIPole(1))%G0rt, &
        & "pole 2    = ", TRIM(vSpc(MxModl%vIPole(2))%Formula),vSpc(MxModl%vIPole(2))%G0rt
      ENDIF
      !
    ENDIF
    !------------------------------------------------------------------/
    !
    I= vDiscretParam(iDis)%I
    J= vDiscretParam(iDis)%J
    K= vDiscretParam(iDis)%K
    !
    ntot= REAL(I+J+K)
    MxPhas%vLPole(1)=I>0 ; MxPhas%vXPole(1)= I /ntot
    MxPhas%vLPole(2)=J>0 ; MxPhas%vXPole(2)= J /ntot
    MxPhas%vLPole(3)=K>0 ; MxPhas%vXPole(3)= K /ntot
    !
    IF(TRIM(MxModl%Model)=="SITE") THEN
      ALLOCATE(vXAtom(MxModl%NAtom))
      CALL MixModel_XPoleToXSite(  &
      & MxModl,MxPhas%vXPole(1:MxModl%NPole), &
      & vXAtom)
    ENDIF
    !
    !WRITE(cFormat,'(A,I3,A)') '(A,1X',MxModl%NPole,'(G12.3,1X))'
    !IF(iDebug>2) &
    !& WRITE(fTrc,cFormat) "MxModl%vIPole=",(MxPhas%vXPole(I),I=1,MxModl%NPole)
    !
    !~ M= 4 !SIZE(vSpc(iSp)%vStoikio)
    !~ WRITE(cFormat,'(A,I3,A)') '(A,1X,I6,1X,',M+1,'(I6,1X))'
    !~ IF(iDebug>2) THEN
      !~ WRITE(91,cFormat) "Pol1%vStoik=",&
      !~ & vSpc(MxModl%vIPole(1))%Div, &
      !~ & (vSpc(MxModl%vIPole(1))%vStoikio(I),I=0,M)

      !~ WRITE(91,cFormat) "Pol2%vStoik=",&
      !~ & vSpc(MxModl%vIPole(2))%Div, &
      !~ & (vSpc(MxModl%vIPole(2))%vStoikio(I),I=0,M)

      !~ WRITE(91,cFormat) "Spc%vStoik= ", &
      !~ & vSpc(iSp)%Div,(vSpc(iSp)%vStoikio(I),I=0,M)

      !~ WRITE(91,*)
    !~ ENDIF
    !
    IF(iDebug>2) WRITE(fTrc,'(A,3(F7.3,1X))',ADVANCE="NO") &
    & "I,J,K=", I/ntot, J/ntot, K/ntot
    !
    !!CALL MixPhase_NormCompo(MxPhas)
    !
    vSpc(iSp)%V0= DOT_PRODUCT( &
    & MxPhas%vXPole(1:MxModl%NPole), &
    & vSpc(MxModl%vIPole(1:MxModl%NPole))%V0)
    !
    CALL MixPhase_CalcMixing( &
    & TdgK,Pbar,MxModl,MxPhas, & !IN
    & GMix,G_IdMix,G_XsMix) !OUT
    ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
    !
    GMecaRT= DOT_PRODUCT( &
    & MxPhas%vXPole(1:MxModl%NPole), &
    & vSpc(MxModl%vIPole(1:MxModl%NPole))%G0rt)
    !
    vSpc(iSp)%G0rt= GMecaRT + GMix/R_jk/TdgK
    !
    IF(iDebug>2) WRITE(fTrc,'(4(G15.6,1X))') &
    & G_XsMix/R_jk/TdgK, &
    & GMix/R_jk/TdgK, &
    & GMecaRT, &
    & vSpc(iSp)%G0rt
    !
    IF(ALLOCATED(vXAtom)) DEALLOCATE(vXAtom)
    !
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretSpecies_TP_Update"
  !
ENDSUBROUTINE DiscretSpecies_TP_Update

SUBROUTINE DiscretModel_Formula_Build( &
& vEle,Spc, &
& S)
!--
!-- build the stoikio formula  of "pseudophases" in terms of elements
!--
  USE M_T_Element,ONLY: T_Element
  USE M_T_Species,ONLY: T_Species
  !
  TYPE(T_Element), INTENT(IN) :: vEle(:)
  TYPE(T_Species), INTENT(IN) :: Spc
  CHARACTER(LEN=*),INTENT(OUT):: S
  !
  CHARACTER(LEN=5):: Str
  INTEGER:: iEl, I, nDiv
  !
  S=""
  nDiv= Spc%vStoikio(0)
  !
  DO iEl=1,SIZE(vEle)

    I= Spc%vStoikio(iEl)
    IF(I>0) THEN
      WRITE(Str,'(I5)') I
      S= TRIM(S)//vEle(iEl)%NamEl(1:2)//"("//TRIM(ADJUSTL(Str))//")"
    ENDIF

  ENDDO
  !
  WRITE(Str,'(I3)') nDiv; IF(nDiv<10) Str(1:1)='0'
  !
  S= TRIM(S)//"/("//TRIM(ADJUSTL(Str))//")"
  !
ENDSUBROUTINE DiscretModel_Formula_Build

SUBROUTINE FelsparName(I,J,K,S)
  INTEGER,         INTENT(IN) ::I,J,K
  CHARACTER(LEN=*),INTENT(OUT)::S
  !
  INTEGER         ::I_,J_,K_
  CHARACTER(LEN=2)::Str
  !
  S=""
  I_=I; J_=J; K_=K
  WRITE(Str,'(I2)') I_; IF(I_<10) Str(1:1)='0'; S=TRIM(S)//"AB"//Str//"_" !IF(I_>0) I_=3*I_+1;
  WRITE(Str,'(I2)') J_; IF(J_<10) Str(1:1)='0'; S=TRIM(S)//"OR"//Str//"_" !IF(J_>0) J_=3*J_+1;
  WRITE(Str,'(I2)') K_; IF(K_<10) Str(1:1)='0'; S=TRIM(S)//"AN"//Str//"_" !IF(K_>0) K_=3*K_+1;
  !
END SUBROUTINE FelsparName

ENDMODULE M_DiscretModel_Tools

!~ SUBROUTINE DiscretModel_Formula_Build_( &
!~ & vEle,nDiscret,nDiv1,nDiv2,I,J, &
!~ & vStoik1,vStoik2, &
!~ & S)
!~ !--
!~ !-- build the stoikio formula  of "pseudophases" in terms of elements
!~ !--
  !~ USE M_T_Element,ONLY: T_Element
  !~ !
  !~ TYPE(T_Element), INTENT(IN) :: vEle(:)
  !~ INTEGER,         INTENT(IN) :: nDiscret,nDiv1,nDiv2,I,J
  !~ INTEGER,         INTENT(IN) :: vStoik1(:),vStoik2(:)
  !~ CHARACTER(LEN=*),INTENT(OUT):: S
  !~ !
  !~ INTEGER,ALLOCATABLE:: vStoik(:)
  !~ CHARACTER(LEN=3):: Str
  !~ INTEGER:: iEl
  !~ !
  !~ ALLOCATE(vStoik(1:SIZE(vEle)))
  !~ !
  !~ S=""
  !~ DO iEl=1,SIZE(vEle)
    !~ !
    !~ vStoik(iEl)= I*vStoik1(iEl)*nDiv1 + J*vStoik2(iEl)*nDiv2
    !~ IF(vStoik(iEl)>0) THEN
      !~ WRITE(Str,'(I3)') vStoik(iEl)
      !~ S= TRIM(S)//vEle(iEl)%NamEl(1:2)//"("//TRIM(ADJUSTL(Str))//")"
      !~ !!IF(iDebug>0) WRITE(fTrc,'(A)') TRIM(S)
    !~ ENDIF
    !~ !
  !~ ENDDO
  !~ !
  !~ WRITE(Str,'(I2)') nDiscret+2; IF(nDiscret<10) Str(1:1)='0'
  !~ S= TRIM(S)//"/("//TRIM(ADJUSTL(Str))//")"
  !~ !
  !~ DEALLOCATE(vStoik)
  !~ !
!~ ENDSUBROUTINE DiscretModel_Formula_Build_

!~ SUBROUTINE FelsparDiscretize
  !~ USE M_Dtb_Const,   ONLY: T_CK,R_jk
  !~ USE M_T_Species,   ONLY: T_Species
  !~ USE M_Global_Vars, ONLY: vEle,vSpc,vMixModel,tFormula
  !~ USE M_Path_Vars,   ONLY: tGrt,vTPpath
  !~ USE M_T_MixModel,  ONLY: T_MixModel,MixModel_Margul_Wg_Calc
  !~ USE M_T_MixPhase
  !~ !
  !~ INTEGER::I,J,K,N
  !~ !! TYPE(T_Species),DIMENSION(:),ALLOCATABLE::vSpc0 !store vSpc
  !~ !
  !~ TYPE(T_MixModel)::S
  !~ TYPE(T_MixPhase)::P
  !~ REAL(dp):: Pbar,TdgK
  !~ REAL(dp):: GIdMix,GMix,GXsMix,GMeca,GTot
  !~ INTEGER::nEl,iTP,iEl
  !~ CHARACTER(LEN=15)::SolName
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A)') "FelsparDiscretize"
  !~ !
  !~ N=0
  !~ DO I=0,Discret
    !~ DO J=0,Discret-I
      !~ N=N+1
      !~ K=Discret-I-J
    !~ ENDDO
  !~ ENDDO
  !~ WRITE(*,*) "N_tot=",N  !!; CALL Pause_
  !~ !
  !~ iTP=1
  !~ nEl=SIZE(vEle)
  !~ P%Name= "FELSPAR"
  !~ P%iModel= 1
  !~ !
  !~ S=vMixModel(1)
  !~ P%vLPole=.FALSE.; P%vXPole=0.0D0
  !~ !
  !~ DO iTP=1,SIZE(vTPpath)
    !~ Pbar= vTPpath(iTP)%Pbar
    !~ TdgK= vTPpath(iTP)%TdgC+T_CK
    !~ IF(iDebug>0) WRITE(fTrc,'(I3,2(1X,G12.3))') iTP,Pbar,TdgK
    !~ !
    !~ !update Margules PARAMETERs for the current T,P
    !~ IF(S%NMarg>0) CALL MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
    !~ vMixModel(1)=S
    !~ !
    !~ N=0
    !~ DO I=0,Discret
      !~ DO J=0,Discret-I
        !~ N=N+1
        !~ K=Discret-I-J
        !~ P%vLPole(1)=I>0; P%vXPole(1)= REAL(I)/REAL(Discret)
        !~ P%vLPole(2)=J>0; P%vXPole(2)= REAL(J)/REAL(Discret)
        !~ P%vLPole(3)=K>0; P%vXPole(3)= REAL(K)/REAL(Discret)
        !~ !
        !~ CALL MixPhase_NormCompo(P)
        !~ CALL MixPhase_CalcMixing(TdgK,Pbar,S,P,GMix,GIdMix,GXsMix) !all OUT in joules
        !~ !GMeca=DOT_PRODUCT(P%vXPole(1:S%NPole),vSpc(S%vIPole(1:S%NPole))%tG0(iTP))*TdgK*R_jk
        !~ GMeca= DOT_PRODUCT(P%vXPole(1:S%NPole),tGrt(S%vIPole(1:S%NPole),iTP)) &
        !~ &    * TdgK*R_jk
        !~ GTot=GMeca+GMix
        !~ !
        !~ !! vGibFas(nFasPur+N,iTP)=GTot*REAL(Discret)
        !~ !
        !~ IF(iTP==1) THEN !DO ONLY for first TP value
          !~ IF(iDebug>0) WRITE(fTrc,'(A,A1,I3,A1,4(G15.6,A1))') &
          !~ & Solname,T_,N,T_,P%vXPole(1),T_,P%vXPole(2),T_,P%vXPole(3),T_,GTot,T_
          !~ !
          !~ CALL FelsparName(I,J,K,SolName)
          !~ !! v%namSol(nFasPur+N)=TRIM(SolName)
          !~ !
          !~ !calculate stoikio of "pseuDOphases" in terms of elements
          !~ DO iEl=1,nEl
            !~ tFormula(nFasPur+N,iEl)= &
            !~ & I*tFormula(S%vIPole(1),iEl) &
            !~ + J*tFormula(S%vIPole(2),iEl) &
            !~ + K*tFormula(S%vIPole(3),iEl)
            !~ !! WRITE(*,'(A3,I6)') vNamEle(iEl),tFormula(nFasPur+N,iEl)
          !~ ENDDO
          !~ !CALL Pause_
        !~ ENDIF
      !~ ENDDO
    !~ ENDDO
  !~ ENDDO
  !~ WRITE(fTrc,'(/,A,I3,/)') "N_tot=",N
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A)') "FelsparDiscretize"
  !~ !! nFasSol=N
!~ END SUBROUTINE FelsparDiscretize

!~ SUBROUTINE DiscretModel_TP_Update( &
!~ !--
!~ !-- §§§!!! OBSOLETE !!!§§§!!! OBSOLETE !!!§§§!!! OBSOLETE !!!§§§!!!
!~ !-- for a mixture phase discretized according to model DiscretModel,
!~ !-- update the T,P dependent parameters
!~ !--
!~ & vSpc,        & !IN
!~ & vMixModel,   & !IN
!~ & TdgK,Pbar,   & !IN
!~ & DiscretModel, & !IN
!~ & vSpcDiscret)    !INOUT discretized species with updated PARAMETERs
  !~ USE M_T_DiscretModel,ONLY: T_DiscretModel !,T_DiscretParam
  !~ USE M_Dtb_Const, ONLY: T_CK,R_jk
  !~ USE M_T_Species, ONLY: T_Species
  !~ USE M_T_MixPhase,ONLY: T_MixPhase,MixPhase_CalcMixing,MixPhase_NormCompo
  !~ USE M_T_MixModel,ONLY: MixModel_Margul_Wg_Calc, T_MixModel
  !~ USE M_T_MixModel,ONLY: MixModel_Index,MixModel_XPoleToXSite
  !~ !
  !~ TYPE(T_Species),     INTENT(IN):: vSpc(:)
  !~ TYPE(T_MixModel),    INTENT(IN):: vMixModel(:)
  !~ REAL(dp),            INTENT(IN):: TdgK,Pbar
  !~ TYPE(T_DiscretModel),INTENT(IN):: DiscretModel
  !~ !
  !~ TYPE(T_Species),  INTENT(INOUT):: vSpcDiscret(:)
  !~ !
  !~ TYPE(T_MixModel):: MxModl
  !~ TYPE(T_MixPhase):: MxPhas
  !~ INTEGER :: nDiscret
  !~ REAL(dp):: G_IdMix,GMix,G_XsMix,GMecaRT !,GTotRT !,V0,WeitKg
  !~ INTEGER :: N,I,J,K
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_TP_Update"
  !~ !
  !~ MxModl=   vMixModel(DiscretModel%iMix)
  !~ nDiscret= DiscretModel%DimTot
  !~ !
  !~ MxPhas%iModel=    DiscretModel%iMix
  !~ MxPhas%vLPole=    .FALSE.
  !~ MxPhas%vXPole=    Zero
  !~ !
  !~ !--- compute stoikio vectors of end members
  !~ IF(iDebug>0) WRITE(fTrc,'(A)') vSpc(MxModl%vIPole(1))%Formula !debug_
  !~ IF(iDebug>0) WRITE(fTrc,'(A)') vSpc(MxModl%vIPole(2))%Formula
  !~ !
  !~ !--- update Margules parameters for the current T,MxPhas
  !~ IF(MxModl%NMarg>0) &
  !~ & CALL MixModel_Margul_Wg_Calc(TdgK,Pbar,MxModl)
  !~ !
  !~ N= 0
  !~ K= 0
  !~ !!DO I=0,nDiscret-1
  !~ DO I=1,nDiscret
    !~ N=N+1 !-> for binary solution
    !~ !DO J=0,Discret-I !-> for ternary solution
      !~ !!J=nDiscret-1 -I
      !~ J= nDiscret+1 -I
      !~ !!N=N+1 !-> for ternary solution
      !~ !!K=Discret-I-J !-> for ternary solution
      !~ !
      !~ MxPhas%vLPole(1)= I>0  ; MxPhas%vXPole(1)= REAL(I) /REAL(nDiscret+1)
      !~ MxPhas%vLPole(2)= J>0  ; MxPhas%vXPole(2)= REAL(J) /REAL(nDiscret+1)
      !~ MxPhas%vLPole(3)= K>0  ; MxPhas%vXPole(3)= REAL(K) /REAL(nDiscret+1)
      !~ !
      !~ CALL MixPhase_NormCompo(MxPhas)
      !~ !
      !~ !IF(TRIM(MxModl%Model)=="SITE") THEN
      !~ !  ALLOCATE(vXAtom(MxModl%NAtom))
      !~ !  CALL MixModel_XPoleToXSite(  &
      !~ !  & MxModl,MxPhas%vXPole(1:MxModl%NPole), &
      !~ !  & vXAtom)
      !~ !ENDIF
      !~ !
      !~ vSpcDiscret(N)%V0= &
      !~ & DOT_PRODUCT(MxPhas%vXPole(1:MxModl%NPole),vSpc(MxModl%vIPole(1:MxModl%NPole))%V0)
      !~ !
      !~ CALL MixPhase_CalcMixing( &
      !~ & TdgK,Pbar,MxModl,MxPhas, &  !IN
      !~ & GMix,G_IdMix,G_XsMix) !OUT
      !~ ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
      !~ !
      !~ GMecaRT= DOT_PRODUCT( &
      !~ & MxPhas%vXPole(1:MxModl%NPole), &
      !~ & vSpc(MxModl%vIPole(1:MxModl%NPole))%G0rt)
      !~ !
      !~ vSpcDiscret(N)%G0rt= GMecaRT + GMix/R_jk/TdgK
      !~ !
    !~ !ENDDO
  !~ ENDDO
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretModel_TP_Update"
  !~ !
!~ ENDSUBROUTINE DiscretModel_TP_Update
