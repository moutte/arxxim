MODULE M_T_Kinmodel
!--
!-- data structure for kinetic model --
!--
  USE M_Kinds
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: MaxKinTerm
  PUBLIC:: T_KinModel
  PUBLIC:: KinModel_CalcCoeffs
  PUBLIC:: KinModel_Index
  PUBLIC:: KinModel_PrmIndex
  PUBLIC:: KinModel_Show
  
  INTEGER,PARAMETER:: MaxKinTerm=4 !max nr of terms in the rate law
  
  TYPE:: T_KinModel
    !implementation of a 'general' form of rate law, a la USGS (?)
    !rate/surf= SUM_i ( k_i * exp(-Ea_i /RT) * a_i^n_i * finh )  !with k_i=10^-pk_i
    !(finh=  inhibitor=  e.g. Al-species <-not implemented yet
    CHARACTER(23):: Name !name of kinetic model
    !CHARACTER(7):: Special !code for special case; not used now
    INTEGER      :: NTermP,NTermD !number of terms, in Precip. and Dissol. laws, 
    !                             !up to MaxKinTerm each
    INTEGER,      DIMENSION(1:MaxKinTerm):: &
    & ISpcPrec, ISpcDiss, &    !index of "activator" species in current species array
    & JSpcPrec, JSpcDiss       !idem, for "DOuble depENDency", e.g. sulfides
    !
    REAL(dp),     DIMENSION(1:MaxKinTerm)::&
    & Kd,pK_d,E_d,N_d,NJ_d,& !dissolution rate parameters !Kd_A=pKdA*FUNCTION(-EdA/RT)
    & Kp,pK_p,E_p,N_p,NJ_p   !precipitation rate parameters !Kp_A=pKpA*FUNCTION(-EdA/RT)
    !
    REAL(dp)::&
    & AlfaD,BetaD,& !VmQsK=(QsK**M%AlfaD - One)**M%BetaD, etc
    & AlfaP,BetaP
    !
  ENDTYPE T_KinModel
  !
CONTAINS
  !
SUBROUTINE KinModel_CalcCoeffs(M,TimeFactor,TdgK) !-> Kd*EXP(-Ea/RT)
!--
!-- temperature dependency of kinetic parameters --
!--
  USE M_Dtb_Const,ONLY: R_jK
  !
  TYPE(T_KinModel),INTENT(INOUT)::M
  REAL(dp),        INTENT(IN)   ::TimeFactor,TdgK
  !
  REAL(dp)::xRT
  !
  M%Kd=Zero
  M%Kp=Zero
  !
  ! PalandriKharaka, USGS_OpenFileReport, 2004-1068
  ! k(T)=k_298 * exp( - Ea_298 (1/T - 1/298) / R )
  ! M%Kd_A= TimeFactor * 10**(-M%pKdA) *EXP(M%EdA*xRT)
  !
  xRT= (One/298.15D0 - One/TdgK) /R_jk
  !
  M%Kd(1:M%NTermD)= TimeFactor * 10**(-M%pK_d(1:M%NTermD)) *EXP(M%E_d(1:M%NTermD)*xRT)
  M%Kp(1:M%NTermP)= TimeFactor * 10**(-M%pK_p(1:M%NTermP)) *EXP(M%E_p(1:M%NTermP)*xRT)
  !
ENDSUBROUTINE KinModel_CalcCoeffs

INTEGER FUNCTION KinModel_Index(Str,V)
!--
!-- -> position of T_KinModel with %Name--Str in V
!--
  CHARACTER(LEN=*),INTENT(IN)::Str
  TYPE(T_KinModel),INTENT(IN)::V(:)
  !
  INTEGER::I
  
  KinModel_Index=0
  !
  IF(SIZE(V)==0) RETURN
  !
  I=0
  DO
    I=I+1 
    IF(TRIM(Str)==TRIM(V(I)%Name)) THEN
      KinModel_Index=I
      EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO ! if Str not found, KinModel_Index=0
  
  RETURN
ENDFUNCTION KinModel_Index

SUBROUTINE KinModel_PrmIndex(vPrmBk,M_In,M_Out)
!--
!-- update indexes of species in the kinetic model of phase M
!-- according to current species permutation
!--
  !
  INTEGER,DIMENSION(:),INTENT(IN)   :: vPrmBk
  TYPE(T_KinModel),    INTENT(IN)   :: M_in
  TYPE(T_KinModel),    INTENT(INOUT):: M_out
  !
  !-- precipitation law
  M_Out%ISpcPrec(1:M_Out%NTermP)= vPrmBk(M_In%ISpcPrec(1:M_In%NTermP))
  ! IF another "activator" species
  WHERE(M_Out%JSpcPrec(1:M_Out%NTermP) >0) &
  & M_Out%JSpcPrec(1:M_Out%NTermP)=vPrmBk(M_In%JSpcPrec(1:M_In%NTermP))
  !
  !-- dissolution law
  M_Out%ISpcDiss(1:M_Out%NTermD)=vPrmBk(M_In%ISpcDiss(1:M_In%NTermD))
  !-- when two activator species
  WHERE(M_Out%JSpcDiss(1:M_Out%NTermD) >0) &
  & M_Out%JSpcDiss(1:M_Out%NTermD)=vPrmBk(M_In%JSpcDiss(1:M_In%NTermD))
  !
ENDSUBROUTINE KinModel_PrmIndex

SUBROUTINE KinModel_Show(f,vSpc,vKinMod,vPrm)
!--
  USE M_T_Species,ONLY: T_Species
  !
  INTEGER,         INTENT(IN):: F
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  TYPE(T_KinModel),INTENT(IN):: vKinMod(:)
  INTEGER,         INTENT(IN),OPTIONAL:: vPrm(:)
  !
  CHARACTER:: T_=CHAR(9)
  INTEGER::I,J,K
  TYPE(T_KinModel)::M
  !
  !CALL GetUnit(f)
  !OPEN(f,FILE=TRIM(DirLog)//"kinmodel.log")
  !WRITE(f,'(A,/)') "!.-> DATA on all minerals found in vSpc"
  !
  WRITE(F,'(/,A)') "< KinModel_Show"
  WRITE(F,'(A,/)')   "species, exponent, pK, ActivEnergy"
  
  DO I=1,SIZE(vKinMod)
    
    WRITE(F,'(A)') TRIM(vKinMod(I)%Name)
    
    M=vKinMod(I)
    
    !-- dissolution terms --
    DO J=1,M%NTermD
      
      K=M%ISpcDiss(J)
      
      IF(PRESENT(vPrm)) K=vPrm(K) !for vSpc
      
      WRITE(F,'(I3,2A,A1,3(F12.3,A1))') &
      & J,"_Dissol ",TRIM(vSpc(K)%NamSp),T_,M%N_D(J),T_,M%pK_D(J),T_,M%E_D(J),T_
      
      K=M%JSpcDiss(J)
      IF(K>0) THEN
        IF(PRESENT(vPrm)) K=vPrm(K) !for vSpc
        WRITE(F,'(I3,2A,A1,F12.3)') &
        & J,"_Dissol ",TRIM(vSpc(K)%NamSp),T_,M%NJ_D(J) !,T_,M%pK_D(J),T_,M%E_D(J),T_
      ENDIF
      
    ENDDO
    
    !-- precipitation terms--
    DO J=1,M%NTermP
    
      K=M%ISpcPrec(J)
      
      IF(PRESENT(vPrm)) K=vPrm(K) !for vSpc
      
      WRITE(F,'(I3,2A,A1,3(F12.3,A1))') &
      & J,"_Precip ",TRIM(vSpc(K)%NamSp),T_,M%N_P(J),T_,M%pK_P(J),T_,M%E_P(J),T_
      
      K=M%JSpcPrec(J)
      IF(K>0) THEN
        IF(PRESENT(vPrm)) K=vPrm(K) !for vSpc
        WRITE(F,'(I3,2A,A1,3(F12.3,A1))') &
        & J,"_Precip ",TRIM(vSpc(K)%NamSp),T_,M%NJ_P(J) !,T_,M%pK_D(J),T_,M%E_D(J),T_
      ENDIF
      
    ENDDO
    
  ENDDO
  WRITE(F,'(A,/)') "</ KinModel_Show"
  
ENDSUBROUTINE KinModel_Show

ENDMODULE M_T_Kinmodel 

