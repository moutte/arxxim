MODULE M_GEM_Write
  !
  USE M_Kinds
  USE M_Trace, ONLY: fTrc,iDebug,T_,Stop_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: GEM_Write_Phases
  PUBLIC:: GEM_Write_Mixtures
  PUBLIC:: GEM_Write_Log_Entete
  !
CONTAINS

SUBROUTINE GEM_Write_Phases( &
& DimPath, &
& vFasIsPresent, &
& vSimplex_Ok, &
& tResult, &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  USE M_Files,        ONLY: DirOut
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_IoTools,      ONLY: GetUnit
  USE M_T_MixModel,   ONLY: T_MixModel
  USE M_GEM_Vars,     ONLY: T_SavPhase
  !
  USE M_Global_Vars,ONLY: vSpc,vFas,vMixModel
  USE M_GEM_Vars,   ONLY: vCpnGEM
  !---------------------------------------------------------------------
  INTEGER,       INTENT(IN):: DimPath
  LOGICAL,       INTENT(IN):: vFasIsPresent(:)
  LOGICAL,       INTENT(IN):: vSimplex_Ok(DimPath)
  REAL(dp),      INTENT(IN) :: tResult(:,:)
  TYPE(T_SavPhase),INTENT(IN):: tResultMix(:,:)
  !---------------------------------------------------------------------
  INTEGER :: iPath,iFs,I,J,K
  INTEGER :: nC,nF,nP,nFpur,nFmix
  INTEGER :: F,F1,F2
  !REAL(dp):: vX(MaxPole)
  REAL(dp):: Tot,X,Y
  TYPE(T_SavPhase):: Fas0
  TYPE(T_MixModel):: MM
  !
  LOGICAL,PARAMETER:: ComputeXMean= .TRUE. !.FALSE.
  !---------------------------------------------------------------------
  nC= SIZE(vCpnGEM)
  nF= SIZE(vFasIsPresent)
  nFpur= SIZE(vFas)
  nFmix= SIZE(vMixModel)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_phase_mole.restab")
  CALL  Write_Title(F)
  !
  CALL GetUnit(F1)
  OPEN(F1,FILE=TRIM(DirOut)//"_phase_cm3.restab")
  CALL  Write_Title(F1)
  !
  CALL GetUnit(F2)
  OPEN(F2,FILE=TRIM(DirOut)//"_phase_grams.restab")
  CALL  Write_Title(F2)
  !
  K=0
  DO iPath=1,DimPath

    IF(vSimplex_Ok(iPath)) THEN
      !
      K=K+1
      !
      WRITE(F, '(I3,A1)',      ADVANCE='NO') K,T_  != count
      WRITE(F1,'(I3,A1)',      ADVANCE='NO') K,T_  != count
      WRITE(F2,'(I3,A1)',      ADVANCE='NO') K,T_  != count
      !
      WRITE(F,'(2(F12.3,A1))',ADVANCE='NO') &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      DO iFs=1,nC != (3:nC+2)                     != components
        WRITE(F,'(G15.6,A1)',ADVANCE='NO') tResult(2+iFs,iPath),T_
      ENDDO
      !
      WRITE(F1,'(2(F12.3,A1))',ADVANCE='NO') &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      DO iFs=1,nC != (3:nC+2)                     != components
        WRITE(F1,'(G15.6,A1)',ADVANCE='NO') tResult(2+iFs,iPath),T_
      ENDDO
      !
      WRITE(F2,'(2(F12.3,A1))',ADVANCE='NO') &     !
      & tResult(1,iPath), T_, &                   != TdgC
      & tResult(2,iPath), T_                      != Pbar
      !
      DO iFs=1,nC != (3:nC+2)                     != components
        WRITE(F2,'(G15.6,A1)',ADVANCE='NO') tResult(2+iFs,iPath),T_
      ENDDO
      !
      !~ Tot=SUM(tResult(1:nC,iPath))
      !~ DO iFs=1,nC
        !~ WRITE(F,'(G15.6,A1)',ADVANCE='NO') tResult(iFs,iPath)/Tot,T_
      !~ ENDDO
      
      !-----------------------------------------------------mole numbers
      DO iFs=1,nF != (nC+3:nC+nF+2)
        IF(vFasIsPresent(iFs)) &
        & WRITE(F,'(G15.3,A1)',ADVANCE='NO') &
        ! nr'moles phase          *nr'oxygen in formula
        ! & tResult(nC+iFs,iPath)*tFormula(iFs,1),T_
        & tResult(2+nC+iFs,iPath),T_
      ENDDO
      WRITE(F,*)
      !----------------------------------------------------/mole numbers

      !---------------------------------------------volume/cm3, weight/g
      DO iFs=1,nFpur != (nC+3:nC+nF+2)
        IF(vFasIsPresent(iFs)) &
        & WRITE(F1,'(G15.3,A1)',ADVANCE='NO') &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%VolM3*1.D6,T_
        IF(vFasIsPresent(iFs)) &
        & WRITE(F2,'(G15.3,A1)',ADVANCE='NO') &
        & tResult(2+nC+iFs,iPath)*vFas(iFs)%WeitKg*1.D3,T_
      END DO
      DO iFs=1,nFmix
        IF(vFasIsPresent(nFpur +iFs)) THEN
        
          Fas0= tResultMix(iFs,iPath)
          MM= vMixModel(Fas0%iModel)
          X= 0.D0
          Y= 0.D0
          DO J=1,MM%nPole
            DO I=1,MM%nPole
              X= X + vSpc(MM%vIPole(I))%V0     * Fas0%tXPole(J,I)
              Y= Y + vSpc(MM%vIPole(I))%WeitKg * Fas0%tXPole(J,I)
            END DO
          END DO
          
          !X= DOT_PRODUCT( &
          !& tResultMix(iFs,iPath)%tXPole(1,:), &
          !& tResultMix(iFs,iPath)%vVol0(:))
          
          WRITE(F1,'(G15.3,A1)',ADVANCE='NO') &
          & tResult(2+nC+nFpur+iFs,iPath)*X*1.D6,T_
          WRITE(F2,'(G15.3,A1)',ADVANCE='NO') &
          & tResult(2+nC+nFpur+iFs,iPath)*Y*1.D3,T_
          
        END IF
      END DO
      WRITE(F1,*)
      WRITE(F2,*)
      !--------------------------------------------/volume cm3, weight/g
    ENDIF !!IF(vSimplex_Ok(iPath))

  ENDDO !!DO iPath
  !
  WRITE(F,'(A1)') "_"
  !
  CLOSE(F)
  CLOSE(F1)
  CLOSE(F2)
  !
  IF(iDebug>2) &
  & PRINT '(2A)', "Phases: Results in ",TRIM(DirOut)//"_phase_xx.restab"
  !
CONTAINS

  SUBROUTINE Write_Title(FF)
    INTEGER,INTENT(in):: FF
    !--------------------------------------------------------title lines
    WRITE(FF,'(3(A,A1))',ADVANCE='NO') "count",T_, "TdgC",T_, "Pbar",T_
    !
    !--component names
    DO i=1,nC
      WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vCpnGEM(i)%NamCp)//"_cpn", T_
    ENDDO
    !
    !--pure phase names
    DO iFs=1,nFpur
      IF(vFasIsPresent(iFs)) &
      & WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vFas(iFs)%NamFs), T_
    ENDDO
    !--mixture phase names
    DO iFs=1,nFmix
      IF(vFasIsPresent(nFpur +iFs)) &
      & WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vMixModel(iFs)%Name), T_
    ENDDO
    WRITE(FF,*)
    !-------------------------------------------------------/title lines
  END SUBROUTINE Write_Title

END SUBROUTINE GEM_Write_Phases

SUBROUTINE GEM_Write_Mixtures( &
& DimPath, &
& vFasIsPresent, &
& vSimplex_Ok,&
& TolX, &
& tResult, &
& tResultMix)
!--
!-- write tabulated results for all phases found stable at any step
!--
  USE M_Files,        ONLY: DirOut
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_IoTools,      ONLY: GetUnit
  USE M_Numeric_Tools,ONLY: iMaxLoc_R
  !
  USE M_T_Phase,      ONLY: T_Phase
  USE M_T_MixModel,   ONLY: T_MixModel,MaxPole
  !
  USE M_Global_Vars,  ONLY: vSpc,vFas,vMixModel
  USE M_GEM_Vars,     ONLY: T_SavPhase
  USE M_GEM_Vars,     ONLY: vCpnGEM
  !---------------------------------------------------------------------
  INTEGER, INTENT(IN):: DimPath
  LOGICAL, INTENT(IN):: vFasIsPresent(:)
  LOGICAL, INTENT(IN):: vSimplex_Ok(DimPath)
  REAL(dp),INTENT(IN):: TolX ! <-MixMinim_TolX
  REAL(dp),INTENT(IN) :: tResult(:,:) ! for T,P values
  TYPE(T_SavPhase),INTENT(INOUT):: tResultMix(:,:)
  !---------------------------------------------------------------------
  INTEGER :: iPath,iFs,I,J,K,P,Q
  INTEGER :: nP,nFmix,nFPur,nPP
  INTEGER :: FMIX
  REAL(dp):: vX(MaxPole)
  REAL(dp):: Tot
  TYPE(T_SavPhase):: Fas0,Fas1,SavPhaseZero
  TYPE(T_MixModel):: MM
  !
  LOGICAL,PARAMETER:: ComputeXMean= .TRUE. !.FALSE.
  !---------------------------------------------------------------------
  CALL GetUnit(FMIX)
  OPEN(FMIX,FILE=TRIM(DirOut)//"_mixtures.restab")
  !
  nFpur= SIZE(vFas)
  nFmix= SIZE(vMixModel)
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  !-----------------------------------------------------------title line
  WRITE(FMIX,'(3(A,A1))',ADVANCE='NO') "count",T_, "TdgC",T_, "Pbar",T_
  !
  DO iFs=1,nFmix
    !
    IF(vFasIsPresent(nFpur +iFs)) THEN
      !
      !Fas= tResultMix(iFs,iPath)
      MM= vMixModel(iFs)
      nP= MM%nPole
      !IF(Fas0%iModel /= 0) THEN
      IF( TRIM(MM%Model)=="MOLECULAR" .AND. MM%NMarg==0 ) THEN
        nPP= 1
      ELSE
        nPP= MM%nPole
      END IF
      DO I=1,nPP
        WRITE(FMIX,'(A,A1)',ADVANCE="NO") TRIM(MM%Name),T_
        DO P=1,nP
          WRITE(FMIX,'(A,A1)',ADVANCE="NO") &
          & TRIM(vSpc(MM%vIPole(P))%NamSp),T_
        ENDDO
      ENDDO
      !
    END IF
  
  ENDDO
  WRITE(FMIX,*)
  !----------------------------------------------------------/title line
  
  K=0
  DO iPath=1,DimPath

    IF(vSimplex_Ok(iPath)) THEN
      !
      K=K+1
      !------------------------------compute mean comp'n of each mixture
      !----------------------------------------------and its free energy
      IF(ComputeXMean) THEN
        !
        DO J=1,nFmix
          !
          Fas0=        tResultMix(J,iPath)
          Fas1=        SavPhaseZero
          Fas1%iModel= Fas0%iModel
          nP=          vMixModel(Fas0%iModel)%nPole
          Fas1%nFas=   nP
          !
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          !IF(Fas0%iModel /= 0) THEN
          DO P=1,Fas0%nFas
            DO Q=1,Fas1%nFas

              vX(1:nP)= ABS( Fas0%tXPole(P,1:nP) - Fas0%tXPole(Q,1:nP) )

              IF(MAXVAL(vX(1:nP)) < TolX *1.0D2) THEN
              
                Tot= Fas0%vMole(P) +Fas1%vMole(Q)
                vX(1:nP)= Fas0%vMole(P) *Fas0%tXPole(P,1:nP) &
                &       + Fas1%vMole(Q) *Fas0%tXPole(Q,1:nP)
                Fas1%tXPole(Q,1:nP)= vX(1:nP) /Tot
                Fas1%vMole(Q)= Tot
                
                EXIT
              
              ENDIF

            ENDDO
          ENDDO
          !
          tResultMix(J,iPath)= Fas1
          !END IF
          !
        ENDDO
        !
      ENDIF
      !----------------------------------/compute mean comp'n of mixture

      !---------------------------------------write mixture compositions
      !-----------------------------------------------------in file FMIX
      WRITE(FMIX,'(I3,A1)',     ADVANCE='NO') K,T_   != count
      !
      WRITE(FMIX,'(2(F7.2,A1))',ADVANCE='NO') &      !
      & tResult(1,iPath), T_, &                      != TdgC
      & tResult(2,iPath), T_                         != Pbar

      DO J=1,nFmix
        IF(vFasIsPresent(nFpur +J)) THEN
          Fas0= tResultMix(J,iPath)
          MM= vMixModel(Fas0%iModel)
          nP= MM%nPole
          !IF(Fas0%iModel /= 0) THEN
          IF( TRIM(MM%Model)=="MOLECULAR" .AND. MM%NMarg==0 ) THEN
            nPP= 1
          ELSE
            nPP= MM%nPole
          END IF
          ! print *,"Fas0%iModel=",Fas0%iModel
          ! print *,"nP=",nP ;  pause
          DO I=1,nPP
            WRITE(FMIX,'(G15.6,A1)',ADVANCE="NO") Fas0%vMole(I),T_
            DO P=1,nP
              WRITE(FMIX,'(F7.3,A1)',ADVANCE="NO") &
              & Fas0%tXPole(I,P),T_
            END DO
          END DO
          !END IF
        END IF
      ENDDO
      !
      WRITE(FMIX,*)
      !--------------------------------------/write mixture compositions
      !
    ENDIF !!IF(vSimplex_Ok(iPath))

  ENDDO !!DO iPath
  !
  IF(iDebug>2) &
  & PRINT '(2A)', "Mixtures: Results in ",TRIM(DirOut)//"_mixtures.restab"
  !
END SUBROUTINE GEM_Write_Mixtures

SUBROUTINE GEM_Write_Log_Entete(vFas)
  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut
  USE M_T_Phase,ONLY: T_Phase
  !
  TYPE(T_Phase),INTENT(IN):: vFas(:)
  !
  INTEGER:: F1,F2
  !
  CALL GetUnit(F1)
  OPEN(F1,FILE=TRIM(DirOut)//"_gibbs.log")
  WRITE(F1,'(A)') "g of all phases, relative to current assemblage"
  !
  CALL GetUnit(F2)
  OPEN(F2,FILE=TRIM(DirOut)//"_mixcomp.log")
  WRITE(F2,'(A)')   "properties of mixtures with min(g)<0 at each iteration"
  WRITE(F2,'(A,/)') "mix.model.name, Xi, Gmin"
  !
END SUBROUTINE GEM_Write_Log_Entete

SUBROUTINE Simplex_WriteTable(iFil,M,N)
  USE M_IoTools,     ONLY: OutStrVec 
  USE M_Simplex_Vars,ONLY: tSimplex 
  INTEGER,INTENT(IN)::iFil,M,N
  !
  REAL(dp),DIMENSION(1:N+1)::V
  INTEGER:: iLin !, iCol
  !
  DO iLin=0,M+1
    V(1:N+1)=tSimplex(iLin,0:N)
    CALL OutStrVec(iFil,V,Opt_C="7")
  ENDDO
  !
  WRITE(iFil,'(A1)') "_"
  !
ENDSUBROUTINE Simplex_WriteTable

ENDMODULE M_GEM_Write

!Simplex, notes from NumRecp:
!
!linear programming= 
!maximize the linear FUNCTION 
!  Z=  a(0,1).x(1) + ... + a(0,n).x(n)
!subject to 
!  primary constraints: a(0,1:n)>=0
!and simultaneously subject to 
!  M=  m1 + m2 + m3 additional constraints,
!  m1 of them of the form: a(i,1).x1 + ... + a(i,N).x(N) >= b(i)   //i=1..m1
!  m2 of them of the form: a(j,1).x1 + ... + a(j,N).x(N) <= b(j)   //j=m1+1..m1+m2
!  m3 of them of the form: a(k,1).x1 + ... + a(k,N).x(N) <= b(k)   //k=m1+m2+1..m1+m2+m3

!On output, the tableau a is indexed by two RETURNed arrays of INTEGERs
!iposv(j) CONTAINS, for j= 1...M, 
!the number i whose original variable x(i) is now represented by row j+1 of a
!These are thus the left-hand variables in the solution
!(The first row of a is of course the z-row.)
!A value i > N indicates that the variable is a y(i) rather than an x(i), x(N+j)==y(j)
!Likewise, izrov(j) contains, for j= 1..N,
!the number i whose original variable x(i) is now a right-hand variable, represented by column j+1 of a. 
!These variables are all zero in the solution. 
!The meaning ofi > N is the same as above, 
!except that i >N +m1 +m2 denotes an artIFicial or slack variable 
!which was used only internally and should now be entirely ignored.


