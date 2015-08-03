MODULE M_MixPhase_Test
!--
!--- tools for testing mixture phases
!--- (mixture- phase of variable composition, following a symmetric model)
!--
  USE M_Kinds
  USE M_Trace
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: MixPhase_Minimize
  PUBLIC:: MixPhase_Test
  PUBLIC:: MixPhase_Show
  
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tGMix,tGtot,tGMarg
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tActiv1,tActiv2
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tGamma1,tGamma2
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tAtom
  !
  INTEGER,PARAMETER:: nPt=99

CONTAINS

SUBROUTINE MixPhase_Minimize
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_T_MixModel, ONLY: T_MixModel
  USE M_IoTools,    ONLY: GetUnit
  !
  USE M_Global_Vars,ONLY: vSpc,vMixModel
  !
  USE M_Optimsolver_Theriak
  USE M_MixModel_Optim
  !
  TYPE(T_MixModel):: MM
  REAL(dp),ALLOCATABLE:: vX(:)
  REAL(dp),ALLOCATABLE:: vMu(:) !
  REAL(dp):: G
  INTEGER :: N
  INTEGER :: i,j,k
  REAL(dp):: TdgK,Pbar
  !
  REAL(dp):: TolX,DeltaInit
  INTEGER :: its,nCallG
  INTEGER :: f,ff
  LOGICAL :: Converge
  CHARACTER(LEN=80):: sFMT
  !
  REAL(dp),PARAMETER:: &
  & Tmin= 600.0D0,   &
  & Tmax= 1500.0D0,  &
  & Tstp= 50.0D0
  !
  IF(SIZE(vMixModel)<1) THEN
    CALL Warning_("NO MIXING MODEL TO BE TESTED")
    RETURN
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixPhase_Minimize"
  !
  IF(iDebug>2) THEN
    CALL GetUnit(ff)
    OPEN(ff,FILE= 'optimsolver_theriak.log')
  ENDIF
  !
  CALL GetUnit(f)
  OPEN(f,FILE= 'out_mixminim.log')
  !
  TolX= 1.0D-3
  DeltaInit= 0.05D0
  !
  MM= vMixModel(1)
  N= MM%NPole
  !
  WRITE(sFMT,'(a,i3,a,i3,a)') &
  & '(3(G15.6,1X),', N,'(G15.6,1X),', N,'(G15.6,1X))'
  !WRITE(sFMT,'(a,i3,a)') '(2G15.6,',N,'(G15.6,1X))'
  !print *,sFMT  ;  pause
  !
  ALLOCATE(vX(1:N))
  ALLOCATE(vMu(1:N))
  !
  TdgK= Tmin +T_CK
  Pbar= 1.0D3
  !
  ALLOCATE(Mixmodel_Optim_vMu0rt(N))
  ALLOCATE(Mixmodel_Optim_vLPole(N))
  
  DoTP: DO
  
    CALL MixModel_Optim_SetParams(TdgK,Pbar,MM)
    DO i=1,N
      Mixmodel_Optim_vMu0rt(i)= vSpc(MM%vIPole(i))%G0rt
      Mixmodel_Optim_vLPole(i)= .TRUE.  !! MM%vHasPole(i)   !!
    END DO
    
    DO i=1,N
      !
      vX(i)= One - 1.0D-3
      DO j=1,N
        IF(j/=i) vX(j)= 1.0D-3/REAL(N-1)
      ENDDO
      !
      CALL Optimsolver_Theriak( & !
      !& MixModel_Optim_GetGibbs,      & !
      !& MixModel_Optim_GetPotentials, & !
      & MixModel_Optim_GetMu, & !
      & ff,            & ! IN
      & DeltaInit,     & ! IN
      & TolX,          & ! IN
      & vX,            & ! INOUT
      & vMu,           & ! OUT
      & G,             & ! OUT
      & Converge,      & ! OUT
      & its,nCallG)      ! OUT
      !
      WRITE(f,sFMT) TdgK-T_CK,Pbar,G,(vX(k),k=1,N),(vMu(k),k=1,N)
      !WRITE(f,sFMT) TdgK-T_CK,Pbar,(vX(k),k=1,N)
      !
    ENDDO
    
    TdgK= TdgK +Tstp
    IF(TdgK -T_CK > Tmax) EXIT DoTP
    
    !~ PAUSE
  
  ENDDO DoTP
  !
  DEALLOCATE(vX)
  DEALLOCATE(vMu)
  DEALLOCATE(Mixmodel_Optim_vMu0rt)
  DEALLOCATE(Mixmodel_Optim_vLPole)
  !
  IF(iDebug>2) CLOSE(ff)
  CLOSE(f)
  !
  WRITE(6,'(A)') 'RESULTS IN out_mixminim.log'
  CALL Pause_
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "</ MixPhase_Minimize"
  !
END SUBROUTINE MixPhase_Minimize

SUBROUTINE MixPhase_Show
!--
!-- check the mixture model for the phase vMixFas(1) of the current system
!--
  USE M_IoTools,ONLY: GetUnit,OutStrVec
  USE M_Files,  ONLY: DirOut
  !
  USE M_T_MixModel,ONLY: T_MixModel,T_Margul !, MixModel_Margul_Wg_Calc
  USE M_T_MixPhase !, ONLY: T_MixPhase, MixPhase_NormCompo
  !
  !--- global variables
  USE M_Global_Vars,ONLY: vMixFas,vMixModel
  !---/
  !
  TYPE(T_MixModel):: S
  TYPE(T_MixPhase):: P
  TYPE(T_Margul)  :: M
  INTEGER:: I,J,K
  INTEGER:: f1
  INTEGER:: iP,iEl
  !
  IF(SIZE(vMixFas)<1) THEN
    CALL Warning_("NO MIXTURE PHASE TO BE TESTED")
    RETURN
  ENDIF
  !
  P= vMixFas(1)
  S= vMixModel(P%iModel)
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirOut)//"_check.res")
  !
  WRITE(f1,'(2A,/)') "Model= ",TRIM(S%Name)
  !
  IF(TRIM(S%Model)=="SITE") THEN
    !
    WRITE(f1,'(A)') "==tPoleAtom=="
    !
    WRITE(f1,'(A12,1X)',ADVANCE="NO") "_"
    DO iEl=1,S%NAtom
      WRITE(f1,'(1X,A6,1X)',ADVANCE="NO") S%vNamAtom(iEl)
    ENDDO
    WRITE(f1,*)
    !
    DO iP=1,S%NPole
      WRITE(f1,'(A12,1X)',ADVANCE="NO") S%vNamPole(iP)
      DO iEl=1,S%NAtom
        WRITE(f1,'(I7,1X)',ADVANCE="NO") S%tPoleAtom(iP,iEl)
      ENDDO
      WRITE(f1,*)
    ENDDO
    !
    WRITE(f1,'(A)')  "==vIAtomSite=="
    WRITE(f1,'(A12,1X)',ADVANCE="NO") "_"
    DO iEl=1,S%NAtom
      WRITE(f1,'(I7,1X)',ADVANCE="NO") S%vIAtomSite(iEl)
    ENDDO
    WRITE(f1,*)
    !
    WRITE(f1,'(A)')  "==vAtomMulti=="
    WRITE(f1,'(A12,1X)',ADVANCE="NO") "_"
    DO iEl=1,S%NAtom
      WRITE(f1,'(I7,1X)',ADVANCE="NO") S%vAtomMulti(iEl)
    ENDDO
    WRITE(f1,*)
    !
    WRITE(f1,'(A)') "==tPoleAtom (fractions)=="
    DO iP=1,S%NPole
      WRITE(f1,'(A12,1X)',ADVANCE="NO") S%vNamPole(iP)
      DO iEl=1,S%NAtom
        WRITE(f1,'(F7.2,1X)',ADVANCE="NO") &
        & S%tPoleAtom(iP,iEl) /REAL(S%vAtomMulti(iEl))
      ENDDO
      WRITE(f1,*)
    ENDDO
    !
    WRITE(f1,'(/,A)') "==Normalization=="
    DO iP=1,S%NPole
      WRITE(f1,'(A12,1X,F7.2)') S%vNamPole(iP), S%vPoleCoeff(iP)
    ENDDO
    WRITE(f1,*)
    !
  ENDIF
  !
  IF(S%NMarg>0) THEN
    WRITE(f1,'(A)') "==Margules=="
    DO I=1,S%NMarg
      M= S%vMarg(I)
      J= S%vIMargSite(I)
      WRITE(f1,'(A,I1,2A)',ADVANCE="NO") "Site= ",J,"= ",S%vNamSite(J)
      WRITE(f1,'(A)',ADVANCE="NO") ",  degrees= "
      DO K=1,S%NAtom
        WRITE(f1,'(1X,I2)',ADVANCE="NO") M%vDegree(K)
      ENDDO
      WRITE(f1,*)
    ENDDO
  ENDIF
  !
END SUBROUTINE MixPhase_Show

SUBROUTINE MixPhase_Test(Cod)
!--
!-- test the mixture model for the phase vMixFas(1) of the current system
!-- using data from theriak database for end members
!-- !!! must have initialized vDtbMinThr before calling MixPhase_Test_1 !!!
!--
  !
  USE M_Dtb_Const,  ONLY: R_jk,T_CK,Tref,Pref
  USE M_Dtb_Calc,   ONLY: SpeciesMin_TP_Update_fromDtb
  USE M_T_Species,  ONLY: Species_Index
  USE M_TPcond_Read,ONLY: TPpath_Read
  USE M_T_MixModel, ONLY: T_MixModel, MixModel_Margul_Wg_Calc,MixModel_XPoleToXSite
  USE M_T_MixPhase, ONLY: T_MixPhase, MixPhase_NormCompo
  USE M_T_MixPhase, ONLY: MixPhase_CalcActivs,MixPhase_CalcMixing
  !
  !--- global variables --
  USE M_Global_Vars,ONLY: vEle,vSpcDtb,vSpc,vMixFas,vMixModel
  USE M_Path_Vars,  ONLY: vTPpath
  !---/
  !
  CHARACTER(LEN=*):: Cod
  !
  TYPE(T_MixModel):: S
  TYPE(T_MixPhase):: P
  REAL(dp):: TdgK,Pbar
  REAL(dp):: GTotRT,GMecaRT,GIdMix,GMix,GXsMix
  REAL(dp):: GMix_FromAct,GXsMix_FromAct
  INTEGER :: iP,jTP,iX,N
  INTEGER :: I,I1,I2
  !
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixPhase_Test "//TRIM(Cod)
  !
  TdgK= Tref ! default
  Pbar= Pref ! values
  !
  CALL TPpath_Read(TdgK,Pbar)
  N= SIZE(vTPpath)
  !
  ALLOCATE(tGMix(N,nPt))
  ALLOCATE(tGTot(N,nPt))
  ALLOCATE(tGMarg(N,nPt))
  !
  ALLOCATE(tActiv1(N,nPt))
  ALLOCATE(tActiv2(N,nPt))
  !  
  ALLOCATE(tGamma1(N,nPt))
  ALLOCATE(tGamma2(N,nPt))
  !
  IF(SIZE(vMixFas)<1) THEN
    CALL Warning_("No Mixture phase in current run")
    RETURN
  ENDIF
  !
  P= vMixFas(1)
  !
  P%vLPole(:)= .FALSE.
  P%vXPole(:)= Zero
  !
  S= vMixModel(P%iModel)
  !
  IF(S%NPole==1) RETURN
  !
  !--------------------- update the vector S%vIPole according to vSpc --
  DO iP=1,S%NPole
    I= Species_Index(S%vNamPole(iP),vSpc)
    IF(I<1) THEN
      CALL Stop_("Species "//TRIM(S%vNamPole(iP))//" NOT in this database !!!")
    ELSE
      S%vIPole(iP)= I !J
    ENDIF
  ENDDO
  vMixModel(P%iModel)= S
  !
  IF(TRIM(S%Model)=="SITE") ALLOCATE(tAtom(nPt,S%NAtom))
  IF(TRIM(S%Model)=="SITE") ALLOCATE(vXAtom(S%NAtom))
  !
  I1=1  ; I2=3    !-> test for feldspar_ss, albite(=1)-anorthite(=3) join
  I1=2  ; I2=1    !-> test for binary solution
  P%vLPole(I1)=.TRUE.
  P%vLPole(I2)=.TRUE. !; S%vLPole(3)=.TRUE.
  !
  LoopTP: DO jTP=1,N
    !
    Pbar= vTPpath(jTP)%Pbar
    TdgK= vTPpath(jTP)%TdgC +T_CK
    !
    IF(iDebug>0) WRITE(*,'(I3,1X,G12.3,1X,G12.3)') jTP,Pbar,TdgK
    !
    !----------------- update values of GR for the end-members at T,P --
    DO I=1,S%NPole
      CALL SpeciesMin_TP_Update_fromDtb(TdgK,Pbar,vSpcDtb,vSpc(S%vIPole(I)))
    ENDDO
    !-----------------/
    !
    !----------------- update Margules parameters for the current T,P --
    IF(S%NMarg>0) CALL MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
    !-----------------/
    !
    !-------------------- test calculation of Gibbs energy / Activity --
    !S%vXPole(3)=0.3D0
    !! P%vXPole(1)= One !-S%vXPole(3)
    !! P%vXPole(2)= 0.0D0
    LoopX: DO iX=1,nPt !increment 0.01 is for nPt=99 !!! 
      !
      P%vXPole(I1)= One - iX*0.01D0 !P%vXPole(1)-0.01D0
      P%vXPole(I2)=       iX*0.01D0 !P%vXPole(2)+0.01D0
      !
      CALL MixPhase_NormCompo(P)
      !
      IF(TRIM(S%Model)=="SITE") THEN
        CALL MixModel_XPoleToXSite( &
        & S, &
        & P%vXPole(1:S%NPole), & 
        & vXAtom(1:S%NAtom))
        IF(jTP==1) tAtom(iX,:)= vXAtom(:)
        !DO I=1,S%NAtom; WRITE(fTrc,'(F12.3,A1)',ADVANCE='NO') P%vXAtom(I),T_; ENDDO
        !WRITE(fTrc,*)
      ENDIF
      !
      !
      CALL MixPhase_CalcActivs(TdgK,Pbar,S,P)
      !
      !--- save in tables for printing
      tActiv1(jTP,iX)= EXP(P%vLnAct(I1))
      tActiv2(jTP,iX)= EXP(P%vLnAct(I2))
      tGamma1(jTP,iX)= EXP(P%vLGam(I1))
      tGamma2(jTP,iX)= EXP(P%vLGam(I2))
      !---/
      !
      CALL MixPhase_CalcMixing( &
      & TdgK,Pbar,S,P, &    !IN
      & GMix,GIdMix,GXsMix) !OUT
      ! NB: in MixPhase_CalcMixing, results are NOT divided by RT !!
      !
      GMecaRT= Zero
      DO I=1,S%NPole
        IF(P%vLPole(I)) &
        & GMecaRT= GMecaRT &
        &        + P%vXPole(I) *vSpc(S%vIPole(I))%G0rt
      ENDDO
      GTotRT= GMecaRT + GMix/TdgK/R_jk
      !
      !--- save in tables for printing
      tGMix(jTP,iX)=  GMix   /TdgK/R_jk
      tGMarg(jTP,iX)= GXsMix /TdgK/R_jk
      tGTot(jTP,iX)=  GTotRT
      !---/
      !
      GMix_FromAct=   P%vXPole(I1)*P%vLnAct(I1) + P%vXPole(I2)*P%vLnAct(I2)
      GXsMix_FromAct= P%vXPole(I1)*P%vLGam(I1)  + P%vXPole(I2)*P%vLGam(I2)
      !PRINT '(A,3G15.6)',"CHECK XsMix=", GXsMix_FromAct, GXsMix/TdgK/R_jk
      !PRINT '(A,3G15.6)',"CHECK Mix=  ", GMix_FromAct,   GMix/TdgK/R_jk
      !
    ENDDO LoopX
      !PAUSE
    !------------------ / test calculation of Gibbs energy / Activity --
    !
  ENDDO LoopTP
  !
  CALL WriteResults(S)
  !
  DEALLOCATE(tGMix,tGtot,tGMarg,tActiv1,tActiv2,tGamma1,tGamma2)
  !
  IF(TRIM(S%Model)=="SITE") DEALLOCATE(vXAtom,tAtom)
  !
  DEALLOCATE(vTPpath)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ MixPhase_Test_"//TRIM(Cod)
  !
  RETURN
ENDSUBROUTINE MixPhase_Test

SUBROUTINE WriteResults(S)
  USE M_IoTools,   ONLY: GetUnit,OutStrVec
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Files,     ONLY: DirOut
  USE M_Path_Vars, ONLY: vTPpath
  USE M_T_MixModel,ONLY: T_MixModel
  !
  TYPE(T_MixModel),INTENT(IN):: S
  !
  INTEGER:: f1,f1x,f21,f22,f31,f32,fAtom
  INTEGER:: N,I,J,iX
  CHARACTER(LEN=30):: sFormat
  !
  ! line= TP-path
  ! column= composition (from 1% to 99% end-member 2)
  !
  ! file out_mix_gmix.restab  = G_mixing and G_mixing_xs (on different lines)
  ! file out_mix_activ.restab = activities (e-m 1, then e-m 2, on same mine)
  !
  CALL GetUnit(f1);  OPEN(f1, FILE=TRIM(DirOut)//"_mix_gmix.restab")
  CALL GetUnit(f1x); OPEN(f1x,FILE=TRIM(DirOut)//"_mix_gxs.restab")
  CALL GetUnit(f21); OPEN(f21,FILE=TRIM(DirOut)//"_mix_activ1.restab")
  CALL GetUnit(f22); OPEN(f22,FILE=TRIM(DirOut)//"_mix_activ2.restab")
  CALL GetUnit(f31); OPEN(f31,FILE=TRIM(DirOut)//"_mix_gamma1.restab")
  CALL GetUnit(f32); OPEN(f32,FILE=TRIM(DirOut)//"_mix_gamma2.restab")
  !
  N= SIZE(vTPpath)
  !
  IF(TRIM(S%Model)=="SITE") THEN
    CALL GetUnit(fAtom)
    OPEN(fAtom,FILE=TRIM(DirOut)//"_mix_atom.restab")
    DO I=1,S%NAtom; WRITE(fAtom,'(A,A1)',ADVANCE='NO') S%vNamAtom(I),T_; ENDDO
    WRITE(fAtom,*)
  ENDIF
  !
  WRITE(sFormat,'(A,I3,A)') '(2(A,A1),',N,'(A15,A1))'
  !
  !------------------------------------------------ header Gibbs file --
  WRITE(f1,sFormat) ".WHAT",T_,"index",T_, (vTPpath(J)%Name,T_,J=1,N)
  !---------------------------------------------- header GibbsXs file --
  WRITE(f1x,sFormat) ".WHAT",T_,"index",T_,(vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Activ file --
  WRITE(f21,sFormat) ".WHAT",T_,"iX",T_,   (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Activ file --
  WRITE(f22,sFormat) ".WHAT",T_,"iX",T_,   (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Gamma file --
  WRITE(f31,sFormat) ".Typ",T_,"iX",T_,    (vTPpath(J)%Name,T_,J=1,N)
  !------------------------------------------------ header Gamma file --
  WRITE(f32,sFormat) ".Typ",T_,"iX",T_,    (vTPpath(J)%Name,T_,J=1,N)
  !
  !---------------------------------------------------- WRITE RESULTS --
  WRITE(sFormat,'(A,I3,A)') '(A,A1,I3,A1,',N,'(G15.8,A1))'
  !print *,sFormat  ;  pause
  !
  LoopX_: DO iX=1,nPt
    ! 
    !------------------------------- total Gibbs energy of mixing /RT --
    WRITE(f1,sFormat)  "GMix/T",T_,  iX,T_, (tGMix(J,iX),T_,J=1,N)
    !------------------------------ excess Gibbs energy of mixing /RT --
    WRITE(f1x,sFormat) "GXsMix/T",T_,iX,T_, (tGMarg(J,iX),T_,J=1,N)
    !
    !------------------------------------------------------- activity --
    ! results for component 1
    WRITE(f21,sFormat) "ACTIV1",T_,  iX,T_, (tActiv1(J,iX),T_,J=1,N)
    ! results for component 2
    WRITE(f22,sFormat) "ACTIV2",T_,  iX,T_, (tActiv2(J,iX),T_,J=1,N)
    !
    !---------------------------------------------------------- gamma --
    ! results for component 1
    WRITE(f31,sFormat) "GAMMA1",T_,  iX,T_, (tGamma1(J,iX),T_,J=1,N)
    ! results for component 2
    WRITE(f32,sFormat) "GAMMA2",T_,  iX,T_, (tGamma2(J,iX),T_,J=1,N)
    !
    IF(TRIM(S%Model)=="SITE") THEN
      DO I=1,S%NAtom; WRITE(fAtom,'(F12.3,A1)',ADVANCE='NO') tAtom(iX,I),T_; ENDDO
      WRITE(fAtom,*)
    ENDIF
  ENDDO LoopX_
  !
  CLOSE(f1)   ;  CLOSE(f1x)
  CLOSE(f21)  ;  CLOSE(f22)
  CLOSE(f31)  ;  CLOSE(f32)
  IF(TRIM(S%Model)=="SITE") CLOSE(fAtom)
  !-------------------------------------------------- / WRITE RESULTS --
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') "results in files "//TRIM(DirOut)//"_mix_*.restab"
  IF(iDebug>0) PRINT      '(A)', "results in files "//TRIM(DirOut)//"_mix_*.restab"
  IF(iDebug>2) CALL Pause_
  !
ENDSUBROUTINE WriteResults

ENDMODULE M_MixPhase_Test
