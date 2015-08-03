MODULE M_Solmodel_Calc_Pitzer
!--
!-- Module routines associées à Pitzer
!-- Livraison 10/2006: Nicolas Ferrando
!-- Remoduling : Anthony Michel
!-- Modified : J.Moutte, 11/2006, 01/2014
!-- -> M_Solmodel_Pitzer_Dtb-  database reading, build vPitz, vIPitz
!-- -> M_Solmodel_Pitzer_Calc- calculations
!--

  USE M_Kinds
  USE M_Trace,ONLY: fTrc,iDebug,T_
  
  IMPLICIT NONE
  
  PRIVATE
  !
  PUBLIC:: Solmodel_Calc_Pitzer
  PUBLIC:: EThetaCalc
  
  !! PUBLIC:: Pitzer_Calc_APhiMonnin
  
  INTEGER:: fTrcPitz
  
CONTAINS

SUBROUTINE Solmodel_Calc_Pitzer( &
& vSpc,                     & !IN
& MWsv,Rho,Eps,TdgK,vMolal, & !IN
& vLnGam,LnActSv,Osmo)        !OUT

  USE M_Numeric_Const,ONLY: Pi
  USE M_T_Species,    ONLY: T_Species
  !
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  REAL(dp),INTENT(IN) :: MWsv,Rho,Eps,TdgK
  REAL(dp),INTENT(IN) :: vMolal(:)
  !
  REAL(dp),INTENT(OUT):: vLnGam(:)
  REAL(dp),INTENT(OUT):: LnActSv   !solvent activity
  REAL(dp),INTENT(OUT):: Osmo
  !
  TYPE(T_Species),ALLOCATABLE:: vSolut(:)
  REAL(dp),       ALLOCATABLE:: vMolalSolut(:)
  INTEGER :: I,J,N
  REAL(dp):: APhi
  
  CALL Pitzer_Calc_APhi( &
  & Rho,Eps,TdgK, &
  & APhi)
  ! result at 25C/1ATM = 0.392521
  !
  ! CALL Pitzer_Calc_APhiMonnin(TdgK,1.013D0,Aphi)
  ! result at 25C/1ATM = 0.3914566
  !  
  ! APhi= 2.303D0 *0.5092D0 / 3.D0 !! Aphi
  ! = 0.380896
  !
  !! comparison of values of APhi at 298K/1atm
  !! 1- according to Pitzer formula    = 0.392521
  !! 2- according to Monnin regression = 0.391457
  !! 3- according to Bethke            = 0.380896
  
  N= COUNT(vSpc(:)%Typ=="AQU") !- 1 !mod 19/06/2008 09:38
  ALLOCATE(vSolut(N))
  ALLOCATE(vMolalSolut(N))
  
  J= 0
  DO I=1,SIZE(vSpc)
    IF(vSpc(I)%Typ=="AQU" .AND. vSpc(I)%NamSp/="H2O") THEN
      J= J+1
      vSolut(J)= vSpc(I)
      vMolalSolut(J)= vMolal(I)
    ENDIF
  ENDDO
  
  CALL Pitzer_Calc_Gamma( &
  & vSolut,Aphi,vMolalSolut, & !IN
  & Osmo,vLnGam) !OUT
  !
  LnActSv= -Osmo *MWsv * (SUM(vMolalSolut))
  
  DEALLOCATE(vSolut)
  DEALLOCATE(vMolalSolut)
  
  RETURN
ENDSUBROUTINE Solmodel_Calc_Pitzer

SUBROUTINE Pitzer_Calc_Gamma( &
& vSpc,Aphi,vMolal, &
& Osmo,vLnGam) !,ERR)
!--
!-- Calcule les coefficients d'activité des solute's
!-- et le coefficient osmotique de l'eau par le mod'ele de Pitzer
!--
!-- vSpc is the list of solute species only
!--
!--IN
!--  vSpc:   solute's
!--  Aphi:   pente du terme de Debye-Huckel
!--  vMolal: molalite's des solute's
!--
!--OUT
!--  Osmo: Coefficient osmotique 
!--  vLnGam: log neperien des coefficient d'activite des constituants
!--  ERR: indicateur d'erreur (non utilise -> a completer si besoin)
!--
  USE M_IoTools,ONLY:GetUnit
  USE M_T_Species,ONLY: T_Species
  USE M_Solmodel_Pitzer_Dtb,ONLY: tBeta0,tBeta1,tBeta2,tAlfa1,tAlfa2
  USE M_Solmodel_Pitzer_Dtb,ONLY: tTheta,tPsi,tCPhi,tZeta,tLamda
  
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  REAL(dp),       INTENT(IN) :: Aphi
  REAL(dp),       INTENT(IN) :: vMolal(:)
  !
  REAL(dp),       INTENT(OUT):: Osmo
  REAL(dp),       INTENT(OUT):: vLnGam(:)
  !
  !------------------------------------------------ Variables locales --
  INTEGER :: iC,jC,kC,iA,jA,kA,jN,A1,A2,C1,C2,K1,K2
  INTEGER :: I,J,K,nSp
  INTEGER :: FF !file index
  REAL(dp):: Ionic, M, BigZ, BigF
  REAL(dp):: RI, Fi0, Fi1, fOsmo
  REAL(dp):: AC, Tmp, X
  REAL(dp),DIMENSION(1:SIZE(vSpc),1:SIZE(vSpc)):: &
  & Bmx1, & ! B_MX  of Pitzer equations, p117, Bethke'96
  & Bmx2, & ! B'_MX of Pitzer equations, p117, Bethke'96
  & Cmx,  & ! C_MX  of Pitzer equations, p118, Bethke'96
  & tETheta, tETheta1
  !
  INTEGER,ALLOCATABLE:: vCharge(:),vIonType(:)
  !---------------------------------------------------------------------
  
  IF(iDebug>0) WRITE(fTrc,'(A)') "< Pitzer_Calc_Gamma"
  !
  nSp= SIZE(vSpc)
  !
  ALLOCATE(vCharge(1:nSp))
  vCharge(:)=vSpc(:)%Z
  !
  ALLOCATE(vIonType(1:nSp))
  vIonType(:)= 0
  WHERE(vCharge(:)>0) vIonType(:)=  1
  WHERE(vCharge(:)<0) vIonType(:)= -1
  !
  IF(iDebug>2) THEN
    !
    CALL GetUnit(ff)
    OPEN(ff,FILE="debug_pitzer.log")
    !
    WRITE(ff,'(A,/)') "trace of routine Pitzer_Calc_Gamma"
    !
    CALL CheckCoeffs
    !
  ENDIF
  
  !---------------------------------------------------- FORCE IONIQUE --
  M=    Zero  !! total molality
  BigZ= Zero  !! total molal charge, sum(|z_i|.m_i)
  Ionic=Zero  !! ionic strength
  DO I=1,nSp
    !print *,I
    M=     M     + vMolal(I)
    BigZ=  BigZ  + vMolal(I)*vCharge(I)*vIonType(I)
    Ionic= Ionic + vMolal(I)*vCharge(I)*vCharge(I)
    !print *,Ionic
  ENDDO
  !pause
  Ionic= Ionic/2.0D0
  RI=  SQRT(Ionic)
  !
  IF(iDebug>2) WRITE(ff,'(A,G15.6)') "Ionic= ", Ionic
  !---------------------------------------------------/ FORCE IONIQUE --
  
  !--------------------------------------- COEFFICIENTS LONGUE PORTEE --
  Fi0= -4.0d0*Aphi*Ionic/1.2D0*LOG(One + 1.2D0*RI)
  Fi1= -4.0d0*Aphi      /1.2D0*LOG(One + 1.2D0*RI) &
  &    -2.0d0*Aphi*RI/(One + 1.2D0*RI)
  IF(iDebug>2) WRITE(ff,'(A,G15.6)') "Fi0=   ", Fi0
  IF(iDebug>2) WRITE(ff,'(A,G15.6)') "Fi1=   ", Fi1
  !--------------------------------------/ COEFFICIENTS LONGUE PORTEE --
  
  !------------------------------------------------ tETheta, tETheta1 --
  tETheta(:,:)= Zero
  tETheta1(:,:)= Zero
  DO I=1,nSp
    DO J=1,nSp
      !
      IF(I==J) CYCLE
      IF(vIonType(I)*vIonType(J)==1) THEN ! ions same sign
        !
        CALL EThetaCalc( &
        & vCharge(I),vCharge(J), RI, Aphi, &
        & tETheta(I,J), tETheta1(I,J))
        !
        IF(iDebug>2) WRITE(ff,'(A,2(G15.6,1X),2(A,A1))') &
        & "ETheta, ETheta1",&
        & tETheta(I,J), tETheta1(I,J), &
        & TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_
        !
      END IF
      !
    END DO
  END DO
  !------------------------------------------------/tETheta, tETheta1 --
  
  Bmx1= Zero
  Bmx2= Zero
  Cmx=  Zero
  !
  !-------------------------------------- PARAMETRES DE PITZER B et C --
  !------------------------------------------------------ B_MX, B'_MX --
  !------------------------------ using tBeta0, tBeta1, tBeta2, tCPhi --
  IF (Ionic /= Zero) THEN
    DO I=1,nSp
      DO J=1,nSp
        IF(I==J) CYCLE
        !
        !X1=tAlfa1(I,J)*RI
        !X2=tAlfa2(I,J)*RI
        !
        Bmx1(I,J)= tBeta0(I,J) &
        &        + tBeta1(I,J)*FONC1(tAlfa1(I,J)*RI) & 
        &        + tBeta2(I,J)*FONC1(tAlfa2(I,J)*RI)
        !
        Bmx2(I,J)=(tBeta1(I,J)*FONC2(tAlfa1(I,J)*RI)    &
        &        + tBeta2(I,J)*FONC2(tAlfa2(I,J)*RI)  ) &
        &        /Ionic
        !
        IF (vCharge(I)*vCharge(J)/=0) & !C_cation-anion eq24, Moeller'98
        & Cmx(I,J)=  0.5D0*tCPhi(I,J) / SQRT(ABS(REAL(vCharge(I)*vCharge(J))))
        !
      END DO
    END DO
    !
    !------------------------------------------------------------trace--
    IF(iDebug>2) THEN
      DO I=1,nSp
        DO J=1,nSp
          IF(I==J) CYCLE
          IF(Bmx1(I,J)/=Zero) &
          & WRITE(ff,'(A,3(G15.6,1X),2(A,A1))') &
          & "Beta-0-1-2   ", &
          & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), &
          & TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_
          !
        END DO
      END DO
      DO I=1,nSp
        DO J=1,nSp
          IF(I==J) CYCLE
          IF(Bmx1(I,J)/=Zero) &
          & WRITE(ff,'(A,3(G15.6,1X),2(A,A1))') &
          & "Bmx1,Bmx2,Cmx", &
          & Bmx1(I,J), Bmx2(I,J), Cmx(I,J), &
          & TRIM(vSpc(I)%NamSp),t_,TRIM(vSpc(J)%NamSp),t_
        END DO
      END DO
    END IF
  END IF
  !---/
  !--------------------------------------/PARAMETRES DE PITZER B et C --
  
  !---------------------------------------------------------------------
  !-------------------------------------------- COEFFICIENT OSMOTIQUE --
  IF (Ionic>Zero) THEN
  
    fOsmo= -Aphi*RI**3 /(One + 1.2D0*RI)
    
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') I," fOsmo=", fOsmo
    !fOsmo=  Ionic*Fi1 - Fi0
    !
    !------------------------------------- somme croisee anion/cation --
    tmp= Zero
    DO I=1,nSp
      DO J=1,nSp
        IF(I==J) CYCLE
        IF(vIonType(I)*vIonType(J)==-1) THEN ! anion/cation
          tmp= tmp &
          &  + vMolal(I)*vMolal(J) &
          &    * ( Bmx1(I,J) + Ionic*Bmx2(I,J) + BigZ*Cmx(I,J) )
        ENDIF
      ENDDO
    ENDDO
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') I," fOsmo, anion/cation=", tmp
    fOsmo= fOsmo + tmp
    !CALL Pause_
    
    !---------------------------------- somme croisee Cation1/Cation2 --
    !--------------------------------------------- using tPsi, tTheta --
    tmp= Zero
    DO C1=1,nSp-1
      IF(vIonType(C1)/=1) CYCLE     ! C1 is CATION
      DO C2=C1+1,nSp
        IF(vIonType(C2)/=1) CYCLE   ! C2 is CATION
        !
        ! SUM_ON_ANIONS PSI(C1,C2,ANI) * vMolal(ANI)
        X=  Zero
        DO K=1,nSp
          IF(vIonType(K)==-1) &     !K is ANION
          & X= X + vMolal(K) *tPsi(C1,C2,K)
        ENDDO
        !X=  Zero
        !
        tmp= tmp &
        &  + vMolal(C1)*vMolal(C2) &
        &    *(tTheta(C1,C2) +tETheta(C1,C2) +RI*RI*tETheta1(C1,C2) +X)
        !
      ENDDO
      !
    ENDDO
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') C1," fOsmo, cation1/cation2=", tmp
    fOsmo= fOsmo + tmp
    !-----------------------------------------------/ Cation1/Cation2 --
    
    !------------------------------------ Somme Croisee Anion1/Anion2 --
    !--------------------------------------------- using tPsi, tTheta --
    tmp= Zero
    DO A1=1,nSp-1
      IF(vIonType(A1)/=-1) CYCLE     ! A1 is ANION
      DO A2=A1+1,nSp
        IF(vIonType(A2)/=-1) CYCLE   ! A2 is ANION
        !
        ! SUM_ON_CATIONS PSI(A1,A2,CAT) * vMolal(CAT)
        X=  Zero
        DO K=1,nSp
          IF(vIonType(K)==1) & !K is cation
          & X= X + vMolal(K) *tPsi(A1,A2,K)
        ENDDO
        ! X=  Zero
        !
        ! BigPhi^phi_ij= theta_ij + ETheta_ij + IStrength * EThetaPrim_ij
        tmp= tmp &
        &    + vMolal(A1) *vMolal(A2) &
        &      *(tTheta(A1,A2) +tETheta(A1,A2) +RI*RI*tETheta1(A1,A2) +X)
        !
      ENDDO
      !
    ENDDO
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') A1," fOsmo, Anion1/Anion2=", tmp
    fOsmo= fOsmo + tmp
    !--------------------------------------------------/Anion1/Anion2 --
    
    !---------------------------------- Somme Sur Les Especes Neutres --
    !--------------------------------------------------- using tLamda --
    tmp= Zero
    DO jN=1,nSp
      IF(vIonType(jN)/=0) CYCLE !jN is Neutral
      DO J=1,nSp
        !
        X= Zero
        DO K=J,nSp
          X= X +vMolal(K)*tZeta(J,K,jN)
        ENDDO
        !
        tmp= tmp + vMolal(jN)*vMolal(J) *(tLamda(J,jN) +X)
        !
      ENDDO
    END DO
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') jN," fOsmo, Neutres=", tmp
    fOsmo= fOsmo + tmp
    !--------------------------------------------------------/neutres --
    !
    fOsmo= fOsmo *2.0D0
    IF(iDebug>2) WRITE(ff,'(I3,A,G15.6)') jN," fOsmo=", fOsmo
    !
    Osmo= One + 2.0D0*fOsmo/M
    
  ELSE ! Force ionique=  0
  
    IF (M==Zero) THEN
      Osmo=  One
    ELSE
      Osmo=  One + fOsmo/M
    ENDIF
  
  ENDIF
  !
  !--------------------------------------------/COEFFICIENT OSMOTIQUE --
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  !------------------------------------------ COEFFICIENTS D'ACTIVITE --
  !
  tmp= 2.303D0 *0.5092D0 / 3.D0
  Fi1= -4.0d0*tmp      /1.2D0*LOG(One + 1.2D0*RI) &
  &    -2.0d0*tmp*RI/(One + 1.2D0*RI)
  BigF= 0.5d0 *Fi1
  ! BigF= BigF
  ! + SUM( m_c  m_a  B'_ca    ) 
  ! + SUM( m_c1 m_c2 PHI'_c1c2) 
  ! + SUM( m_a1 m_a2 PHI'_a1a2) 
  DO I=1,nSp
    DO J=1,nSp
      IF(I==J) CYCLE
      tmp= vMolal(I)*vMolal(J)
      IF(vIonType(I)*vIonType(J)== 1) BigF= BigF + Tmp *tETheta1(I,J)
      IF(vIonType(I)*vIonType(J)==-1) BigF= BigF + Tmp *Bmx2(I,J)
    END DO
  END DO
  IF(iDebug>2) WRITE(ff,'(A,G15.6)') " BigF=", BigF
  
  !-------------------------------------------------- CALCUL LN GAMMA --
  DO iC=1,nSp
    !
    !---------------------------------------- CALCUL POUR LES CATIONS --
    IF(vIonType(iC)==1) THEN
      !
      AC= 0.5d0 *vCharge(iC)**2 *Fi1
      AC= vCharge(iC)**2 *BigF
      !
      IF(iDebug>2) WRITE(ff,'(/,A,/)') TRIM(vSpc(iC)%NamSp)
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " F-term=", AC
      !
      !----------------------------------------- Somme Sur Les Anions --
      Tmp= Zero
      DO jA=1,nSp
        IF(vIonType(jA)==-1) THEN !-> anion
          K1= MIN(iC,jA)
          K2= MAX(iC,jA)
          Tmp= Tmp + vMolal(jA) &
          & *( 2.0d0*Bmx1(K1,K2) +BigZ*Cmx(K1,K2) )
        END IF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 1=anions_=", Tmp
      AC= AC + Tmp
      !
      !---------------------------------------- Somme Sur Les Cations --
      Tmp= Zero
      DO jC=1,nSp
        IF(iC==jC) CYCLE
        IF(vIonType(jC)==1) THEN !-> cation
          !
          K1= MIN(iC,jC)
          K2= MAX(iC,jC)
          !
          X= Zero
          DO jA=1,nSp
            IF(vIonType(jA)==-1) & ! iC-jC-jA
            & X= X &
            &  + vMolal(jA) *tPsi(K1,K2,jA)
          END DO
          X= Zero
          !
          Tmp= Tmp &
          &  + vMolal(jC) &
          &    *( 2.0D0*tTheta(K1,K2) + 2.0D0*tETheta(K1,K2) + X)
          
          IF(iDebug>2) WRITE(FF,'(2A12,G12.3)') &
          & TRIM(vSpc(K1)%NamSp),TRIM(vSpc(K2)%NamSp), &
          & tTheta(K1,K2)+tETheta(K1,K2)
          
          !! Tmp= Tmp &
          !! &  + vMolal(jC) &
          !! &    *( 2.0D0*tTheta(jC,iC) + 2.0D0*tETheta(jC,iC) )
          !
        ENDIF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 2=cations=", Tmp
      AC= AC + Tmp
      !
      !------------------------------------------ Somme Anion - Anion --
      Tmp= Zero
      DO jA=1,nSp-1
        IF(vIonType(jA)==-1) THEN !ANI
          DO kA=jA+1,nSp
            IF(vIonType(jA)==-1) THEN !ANI
              Tmp= Tmp &
              & + vMolal(jA) *vMolal(kA) *tPsi(jA,kA,iC) !ANI-ANI-CAT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 3=ani/ani=", Tmp
      AC=  AC + Tmp
      !
      !------------------ Somme Croisée Sur Les Cations Et Les Anions --
      Tmp=  Zero
      DO I=1,nSp
        DO J=1,nSp
          IF(vIonType(I)*vIonType(J)==-1) &
          & Tmp= Tmp &
          &    + vMolal(I)*vMolal(J) *vCharge(iC) *Cmx(I,J)
        ENDDO
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 4=ani/cat=", Tmp
      AC= AC + Tmp
      !
      !-------------------------------- Somme Sur Les Especes Neutres --
      Tmp= Zero
      DO jN=1,nSp
        IF(vIonType(jN)==0) &
        & Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iC,jN)
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 5=neutral=", Tmp
      AC= AC +Tmp
      !---/
      !
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iC)=AC
    ENDIF
    !---------------------------------------/ CALCUL POUR LES CATIONS --
  END DO
    
  DO iA=1,nSp
    !
    !----------------------------------------- CALCUL POUR LES ANIONS --
    IF(vIonType(iA)==-1) THEN
      !
      IF(iDebug>2) WRITE(ff,'(/,A,/)') TRIM(vSpc(iA)%NamSp)
      !
      AC= 0.5d0 *vCharge(iA)**2 *Fi1
      AC= vCharge(iA)**2 *BigF
      !
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " LnGamma=", AC
      !
      !-------------------------------------------- somme sur cations --
      Tmp= Zero
      DO jC=1,nSp
        IF(vIonType(jC)==1) THEN
          K1= MIN(iA,jC)
          K2= MAX(iA,jC)
          Tmp= Tmp + vMolal(jC) &
          &         *(2.0d0 *Bmx1(K1,K2) &
          &            + BigZ*Cmx(K1,K2))
        END IF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 1=cations=", Tmp
      AC= AC + Tmp
      !
      !--------------------------------------------- somme sur anions --
      Tmp= Zero
      DO jA=1,nSp
        IF(iA==jA) CYCLE
        IF(vIonType(jA)==-1) THEN
          !
          K1= MIN(iA,jA)
          K2= MAX(iA,jA)
          !
          X= Zero
          DO jC=1,nSp
            IF(vIonType(jC)==1) &
            & X= X &
            &  + vMolal(jC) *tPsi(K1,K2,jC)
          END DO
          !X= Zero
          Tmp=  Tmp &
          &  + 2.0d0*vMolal(jA)*(tTheta(K1,K2) +tETheta(K1,K2)) 
          !
          IF(iDebug>2) WRITE(FF,'(2A12,G12.3)') &
          & TRIM(vSpc(K1)%NamSp),TRIM(vSpc(K2)%NamSp), &
          & tTheta(K1,K2)+tETheta(K1,K2)
          !
        ENDIf
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 2=anions =", Tmp
      AC= AC + Tmp
      !
      !-------------------------------------- somme cations - cations --
      Tmp= Zero
      DO I=1,nSp-1
        IF(vIonType(I)==1) THEN !CAT
          DO J=I+1,nSp
            IF(vIonType(J)==1) THEN !CAT
              !
              Tmp= Tmp &
              &  + vMolal(I)*vMolal(J) *tPsi(I,J,iA) !CAT-CAT-ANI
              !
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 3=cat/cat=", Tmp
      AC= AC + Tmp
      !
      !------------------ somme croisée sur les cations et les anions --
      Tmp= Zero
      DO I=1,nSp
        DO J=1,nSp
          IF(vIonType(I)*vIonType(J)==-1) &
          & Tmp= Tmp &
          &    - vMolal(I)*vMolal(J) *vCharge(iA) *Cmx (I,J)
        ENDDO
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 4=ani/cat=", Tmp
      AC=  AC + Tmp
      !
      !------------------------------------------ Somme Anion - Anion --
      !! Tmp= Zero
      !! DO jA=1,nSp-1
      !!   IF(vCharge(jA)<0) THEN
      !!     !! DO kA=jA+2,nSp
      !!     DO kA=jA+1,nSp
      !!       IF(vCharge(kA)<0) &
      !!       & Tmp= Tmp + vMolal(jA)*vMolal(kA)*tETheta1(jA,kA)
      !!     ENDDO
      !!   ENDIF
      !! ENDDO
      !! Tmp= Tmp *vCharge(iA)**2
      !! IF(iDebug>2) WRITE(ff,'(A,G15.6)') " ani/ani=", Tmp
      !! AC= AC + Tmp
      
      !-------------------------------- Somme Sur Les Especes Neutres --
      Tmp= Zero
      DO jN=1,nSp
        IF(vIonType(jN)==0) &
        & Tmp= Tmp +2.0d0*vMolal(jN)*tLamda(iA,jN)
      ENDDO
      !! Tmp= SUM(vMolal(:)*tLamda(iA,:),MASK=vIonType(:)==0)
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " 5=neutral=", Tmp
      AC= AC + Tmp
      !---/
      !
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " LnGamma=", AC
      !
      vLnGam(iA)=AC
    ENDIF
    !----------------------------------------/ CALCUL POUR LES ANIONS --
  END DO
  
  DO jN=1,nSp
    !
    !-------------------------------- CALCUL POUR LES ESPECES NEUTRES --
    IF(vIonType(jN)==0) THEN
      !
      IF(iDebug>2) WRITE(ff,'(/,A,/)') TRIM(vSpc(jN)%NamSp)
      !
      AC= Zero
      !
      Tmp= Zero
      DO J=1,nSp
        Tmp= Tmp + 2.0*vMolal(J)*tLamda(J,jN)
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " tLamda_=", Tmp
      AC= AC + Tmp
      !
      !----------------------------------------------- terme ternaire --
      Tmp= Zero
      DO J=1,nSp
        IF(vCharge(J)/=0) THEN
          DO K=J+1,nSp
            IF(vCharge(K)/=0) THEN
              Tmp= Tmp + vMolal(J)*vMolal(K)*tZeta(J,K,jN)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(iDebug>2) WRITE(ff,'(A,G15.6)') " ternary=", Tmp
      !-----------------------------------------------/terme ternaire --
      !
      AC= AC + Tmp
      !
      vLnGam(jN)=AC
    ENDIF
    !-------------------------------/ CALCUL POUR LES ESPECES NEUTRES --
    !
  ENDDO
  !
  DO I=1,nSp
    IF(iDebug>2) &
    & WRITE(ff,'(I3,1X,A15,G15.6)') I,vSpc(I)%NamSp,EXP(vLnGam(I))
  ENDDO
  !
  IF(iDebug>2) CLOSE(ff)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Pitzer_Calc_Gamma"
  
  RETURN

CONTAINS

SUBROUTINE CheckCoeffs

  WRITE(FF,'(A)') "==Beta0-1-2 Cphi=="
  DO I=1,nSp
    DO J=1,nSp
      IF(J>I) &
      & WRITE(FF,'(2A12,4G12.3)') &
      & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp), &
      & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), tCPhi(I,J)
    END DO
  END DO
  
  WRITE(FF,'(A)') "==Psi, CAT-CAT-ANI=="
  DO K=1,nSp
    DO I=1,nSp-1
      DO J=I+1,nSp
        IF( vIonType(I)== 1 .AND. &
        &   vIonType(J)== 1 .AND. &
        &   vIonType(K)==-1 )     &
        & WRITE(FF,'(3A12,G12.3)') &
        & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp),TRIM(vSpc(K)%NamSp), &
        & tPsi(I,J,K)
      END DO
    END DO
  END DO
  
  WRITE(FF,'(A)') "==Psi, ANI-ANI-CAT=="
  DO K=1,nSp
    DO I=1,nSp-1
      DO J=I+1,nSp
        IF( vIonType(I)==-1 .AND. &
        &   vIonType(J)==-1 .AND. &
        &   vIonType(K)== 1 )     &
        & WRITE(FF,'(3A12,G12.3)') &
        & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp),TRIM(vSpc(K)%NamSp), &
        & tPsi(I,J,K)
      END DO
    END DO
  END DO
  
  WRITE(FF,'(A)') "==Lambda, ION-NEU=="
  DO J=1,nSp
    IF( vIonType(J)/=0 ) CYCLE
    DO I=1,nSp
      IF( vIonType(I)/=0 )     &
      & WRITE(FF,'(2A12,G12.3)') &
      & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp), &
      & tLamda(I,J)
    END DO
  END DO
  
  !! WRITE(FF,'(A)') "==Beta0-1-2 Cphi=="
  !! DO I=1,nSp
  !!   DO J=1,nSp
  !!     IF(J<I) &
  !!     & WRITE(FF,'(2A12,4G12.3)') &
  !!     & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp), &
  !!     & tBeta0(I,J), tBeta1(I,J), tBeta2(I,J), tCPhi(I,J)
  !!   END DO
  !! END DO
  !! 
  !! WRITE(FF,'(A)') "==Psi=="
  !! DO K=1,nSp
  !!   DO I=1,nSp
  !!     DO J=1,nSp
  !!       IF( J<I &
  !!       .AND. vIonType(I)== 1  &
  !!       .AND. vIonType(J)== 1  &
  !!       .AND. vIonType(K)==-1) &
  !!       & WRITE(FF,'(3A12,G12.3)') &
  !!       & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp),TRIM(vSpc(K)%NamSp), &
  !!       & tPsi(I,J,K)
  !!     END DO
  !!   END DO
  !! END DO
  !! DO K=1,nSp
  !!   DO I=1,nSp
  !!     DO J=1,nSp
  !!       IF( J<I &
  !!       .AND. vIonType(I)==-1  &
  !!       .AND. vIonType(J)==-1  &
  !!       .AND. vIonType(K)== 1) &
  !!       & WRITE(FF,'(3A12,G12.3)') &
  !!       & TRIM(vSpc(I)%NamSp),TRIM(vSpc(J)%NamSp),TRIM(vSpc(K)%NamSp), &
  !!       & tPsi(I,J,K)
  !!     END DO
  !!   END DO
  !! END DO
  
  RETURN
END SUBROUTINE CheckCoeffs

ENDSUBROUTINE Pitzer_Calc_Gamma

REAL(dp) FUNCTION FONC1(X)
!--
!-- function g, used for computing B_MX and B_MX'
!-- from Beta0, Beta1, Beta2
!--
  REAL(dp)::X
  IF (X==Zero) THEN  ;  FONC1= Zero
  ELSE               ;  FONC1= 2.0d0*(One - (One + X)*EXP(-X)) /X/X
  END IF
END FUNCTION FONC1

REAL(dp) FUNCTION FONC2(X)
!--
!-- function g', used for computing B_MX and B_MX'
!-- from Beta0, Beta1, Beta2
!--
  REAL(dp)::X
  IF (X==Zero) THEN  ;  FONC2= Zero
  ELSE               ;  FONC2=-2.0d0*(One - (One + X + X*X/2.0d0)*EXP(-X)) /X/X
  END IF
END FUNCTION FONC2

SUBROUTINE EThetaCalc ( &
& Zi, Zj, SqrtI, Aphi,  &
& Eth, Eth1)
!--
!-- Calcul de la fonction ETheta(I) --
!-- (effet electrostatique asymetrique) --
!-- par la theorie de Pitzer (1975) --
!-- IN --
!  Zi: charge du premier ion
!  Zj: charge du deuxieme ion
!  SqrtI: racine carree de la force ionique
!  Aphi:  terme de Debye-Huckel
!-- OUT --
!  Eth : effet electrostatique asymetrique
!  Eth1 : derivee de Eth par rapport a la temperature
!--
  INTEGER, INTENT(IN) :: ZI, ZJ
  REAL(dp),INTENT(IN) :: SqrtI, APHI
  REAL(dp),INTENT(OUT):: eth, eth1
  !
  REAL(dp):: J0mn,J1mn,J0mm,J1mm,J0nn,J1nn
  REAL(dp):: XMN,XMM,XNN,RI2
  !
  eth=  Zero
  eth1= Zero
  RI2=  SqrtI**2
  
  IF (SqrtI==Zero .OR. Zi==Zj) RETURN
  
  Xmn= 6.0d0 *Zi*Zj *Aphi *SqrtI
  Xmm= 6.0d0 *Zi*Zi *Aphi *SqrtI
  Xnn= 6.0d0 *Zj*Zj *Aphi *SqrtI
  
  CALL EThetaParam (Xmn, J0mn, J1mn)
  CALL EThetaParam (Xmm, J0mm, J1mm)
  CALL EThetaParam (Xnn, J0nn, J1nn)
  
  Eth=  Zi*Zj *(J0mn     - 0.5d0*(    J0mm +     J0nn)) /4.D0/RI2
  Eth1= Zi*Zj *(Xmn*J1mn - 0.5d0*(Xmm*J1mm + Xnn*J1nn)) /8.D0/RI2/RI2 &
  &   - eth/RI2
  
  RETURN
END SUBROUTINE EThetaCalc

SUBROUTINE EThetaParam (X, JX, J1X)
!--
!-- Element pour le calcul de la fonction ETheta(I)
!-- (effet electrostatique asymetrique)
!-- par la theorie de Pitzer (1975)
!-- IN-
!--   X : charge de l'ion considere
!-- OUT-
!--   JX : intermediaire de calcul
!--   J1X : intermediaire de calcul

  REAL(dp),INTENT(IN) :: X
  REAL(dp),INTENT(OUT):: JX, J1X
  
  REAL(dp):: Z, DZ
  INTEGER :: I
  REAL(dp),DIMENSION(0:22):: aik, aiik, bk, dk
  !
  DATA (aik(i), aiik(i), i=0,22)                &
  & /1.925154014814667D0,  0.628023320520852D0, &
  & -0.060076477753119D0,  0.462762985338493D0, &
  & -0.029779077456514D0,  0.150044637187895D0, &
  & -0.007299499690937D0, -0.028796057604906D0, &
  &  0.000388260636404D0, -0.036552745910311D0, &
  &  0.000636874599598D0, -0.001668087945272D0, &
  &  0.000036583601823D0,  0.006519840398744D0, &
  & -0.000045036975204D0,  0.001130378079086D0, &
  & -0.000004537895710D0, -0.000887171310131D0, &
  &  0.000002937706971D0, -0.000242107641309D0, &
  &  0.000000396566462D0,  0.000087294451594D0, &
  & -0.000000202099617D0,  0.000034682122751D0, &
  & -0.000000025267769D0, -0.000004583768938D0, &
  &  0.000000013522610D0, -0.000003548684306D0, &
  &  0.000000001229405D0, -0.000000250453880D0, &
  & -0.000000000821969D0,  0.000000216991779D0, &
  & -0.000000000050847D0,  0.000000080779570D0, &
  &  0.000000000046333D0,  0.000000004558555D0, &
  &  0.000000000001943D0, -0.000000006944757D0, &
  & -0.000000000002563D0, -0.000000002849257D0, &
  & -0.000000000010991D0,  0.000000000237816D0, &
  &  Zero,                Zero,               &
  &  Zero,                Zero                /
  
  DO i=20,0,-1
  
    IF (x <= One) THEN
      Z=     4.0D0 * x**0.2D0    - 2.0D0
      dZ=    0.8D0 * x**(-0.8D0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aik(i)
    ELSE
      Z=     40.0D0/9.0D0 *x**(-0.1D0) - 22.0D0/9.0D0
      dZ=   -40.0D0/90.0D0*x**(-1.1D0)
      bk(i)= Z*bk(i+1) - bk(i+2) + aiik(i)
    ENDIF
    dk(i)=  bk(i+1) + Z*dk(i+1) - dk(i+2)
  
  ENDDO
  
  JX=  0.25D0*x - One + 0.5D0*(bk(0) - bk(2))
  J1X= 0.25D0   +    0.5D0*dZ*(dk(0) - dk(2))
  
  RETURN
END SUBROUTINE EThetaParam

SUBROUTINE Pitzer_Calc_APhi(RhoW,EpsW,TdgK,APhi)
!--
!-- Aphi calculation is slightly different from that used in DH
!-- refs: Anderson-Crerar, p449; or Lin_2003_FPE ...
!--
  ! USE M_Numeric_Const,ONLY: Pi
  ! USE M_Dtb_Const,ONLY: Electron,Boltzmann,Avogadro,EpsVacuum
  !
  REAL(dp),INTENT(IN) :: RhoW,EpsW,TdgK
  REAL(dp),INTENT(OUT):: APhi
  
  !! REAL(dp):: X
  
  !! X= SQRT( (2.0D0*pi*Avogadro *1.0D-3)  &
  !! &          /(4.0D0*pi *EpsVacuum *Boltzmann)**3 ) &
  !! &   * Electron**3 &
  !! &   *1.0D+3 /3.0D0
  !! print '(F24.6)',X
  !! pause
  
  !! Aphi= SQRT( (2.0D0*pi*Avogadro *RhoW*1.0D-3)  &
  !! &          /(4.0D0*pi *EpsVacuum*EpsW *Boltzmann*TdgK)**3 ) &
  !! &   * Electron**3 &
  !! &   *1.0D+3 /3.0D0
  
  Aphi= 1400684.314D0 *SQRT( RhoW / (EpsW*TdgK)**3 )
  
  !~ APhi_DH= 1.824829238D6 *SQRT( RhoW / (EpsW*TdgK)**3)
  !~ call Pitzer_Calc_APhiMonnin(TdgK,1.0D0,Aphi_)
  !~ print '(A,3G15.6)',"Aphi, AphiMonnin,APhi_DH:", Aphi,Aphi_,APhi_DH
  !~ pause
  
  RETURN
ENDSUBROUTINE Pitzer_Calc_APhi

SUBROUTINE Pitzer_Calc_APhiMonnin(TdgK,Pbar,Aphi_)
!--
!-- Module de calcul du terme APhi de Debye-Huckel
!-- (version spéciale Pitzer)
!-- Correlations de Monnin, Chemical Geology, 153(1999) 187-209
!--
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(OUT):: Aphi_
  !
  REAL(dp)::Psat,T,A0,AV
  !
  A0=Zero
  AV=Zero
  APhi_=Zero
  T=TdgK
  !
  A0=  3.36901531D-1                 &
  &  - 6.32100430D-4 * T             &
  &  + 9.14252359    / T             & 
  &  - 1.35143986D-2 * LOG(T)        & 
  &  + 2.26089488D-3 / (T-263.0D0)   &
  &  + 1.92118597D-6 * T*T           &
  &  + 4.52586464D1  / (680.D0-T)     
  !
  IF (T<=373.15D0) THEN
    !
    AV=  8.106377         &
    & - 1.256008D-1 *T    &
    & + 7.760276D-4 *T**2 &
    & - 2.098163D-6 *T**3 &
    & + 2.25777D-9  *T**4 !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0D-06 !conversion en m3
    !
  ELSE
    !
    AV= 3.849971D2          &
    & - 6.982754     * T    &
    & + 3.877068D-2  * T**2 &
    & - 1.11381D-4   * T**3 &
    & + 1.589736D-7  * T**4 &
    & - 9.395266D-11 * T**5 &
    & + 6.469241D4 / (680.D0-T) !cm3.kg^0.5.mol(-1.5), Monnin,1989
    !
    AV=  AV * 1.0D-06 !conversion en m3
    !
  ENDIF
  !  
  Psat= EXP( 7.3649D1 -7.2582D3/T -7.3037D0*LOG(T) +4.1653D-6*T*T )
  ! in Pascal
  !
  Aphi_ =  A0 + (Pbar*1.0D+5 - Psat) * (-AV/(4.0*8.314*T))  !Av=-1/4RT * dAphi/dP_, Monnin, 1989
  !  
  RETURN
ENDSUBROUTINE Pitzer_Calc_APhiMonnin

ENDMODULE M_Solmodel_Calc_Pitzer
