MODULE M_Simplex_Calc
  !
  USE M_Kinds
  USE M_Trace,ONLY:Stop_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Simplex_Calc
  !
  !!!  using Simplex_Calc,
  !!!  for a system of (nCn Components, nFs Fases)
  !!!  where
  !!!    vMolCpn(1:nCpn) is the components' abundance
  !!!    vFasGibbs(1:nFs) is the array of the fases' free energy
  !!!    tStoikSpl is the stoikiometry table of the Fases in terms of the components
  !!!
  !!!  CALL Simplex_Vars_Alloc(nCpn,nFs)
  !!!  !
  !!!  !first row= Gibbs energy of phases 1:nFs
  !!!  tSimplex(0,1:nFs)= vFasGibbs(1:nFs)
  !!!  !first column= bulk compos'n  
  !!!  tSimplex(1:nCpn,0)= vMolCpn(1:nCpn)       
  !!!  !main = stoikiometry matrix              
  !!!  tSimplex(1:nCpn,1:nFs)=-TRANSPOSE(tStoikSpl(1:nFs,1:nCpn))
  !!!  !                      
  !!!  CALL Simplex_Calc(iError)
  !!!  !
  !!!  iError is  0 =  Ok
  !!!  iError is  1 = "Unbounded Objective Function"
  !!!  iError is -1 = "No Solutions Satisfy Constraints Given"
  !
CONTAINS

SUBROUTINE Simplex_Calc(iError) ! &
!! & tSimplex,IZROV,IPOSV, &
!! & iError,n1,n2)
!--
!-- simplex routine, 
!-- modified from Press et al, Numerical Recipes in Fortran
!--
  USE M_Simplex_Vars, ONLY: tSimplex,IZROV,IPOSV
  USE M_Numeric_Tools,ONLY: iFirstLoc,iMaxLoc_R,OuterProd_R, Swap_I
  
  ! M=nCpn
  ! N=nFs
  
  !! REAL(dp),INTENT(INOUT):: tSimplex(:,:)
  !! INTEGER, INTENT(OUT)  :: IZROV(:),IPOSV(:)
  INTEGER, INTENT(OUT)  :: iError !,n1,n2
  !
  INTEGER:: M,N
  INTEGER:: n1,n2
  INTEGER,DIMENSION(:),ALLOCATABLE:: L1,L2
  !
  REAL(dp),PARAMETER:: EPS= 1.0E-12_dp  !1.0E-6_dp 
  !
  INTEGER :: NL1,NL2 !IP,K,KH,KP,
  INTEGER :: IP,KP,K,KH !,NL1,NL2,M,N
  REAL(dp):: bMax
  LOGICAL :: INIT,PROCEED
  
  iError= 0
  !
  M= SIZE(tSimplex,1) -2 !or M= SIZE(IPOSV)
  N= SIZE(tSimplex,2) -1 !or N= SIZE(IZROV)
  !
  ALLOCATE(L1(1:N+1)) !SIZE(tSimplex,2)
  ALLOCATE(L2(1:M+2)) !SIZE(tSimplex,1)
  !
  IF (ANY(tSimplex(1:M,0) < Zero)) THEN
    !! CALL Stop_('BAD INPUT TABLEAU IN SIMPLX')
    iError= -2
    RETURN
  END IF
  !
  NL1=      N
  L1(1:N)=  Arth_I(1,1,N)
  IZROV(:)= L1(1:N)
  !
  NL2=      M
  L2(1:M)=  Arth_I(1,1,M)
  iPosV(:)= N+L2(1:M)
  !
  INIT=    .TRUE.
  n1= 0
  n2= 0
  
  DoMain: DO
    !------------------------------------------------------------phase 1
    Phase1: DO
      
      IF (INIT) THEN
        INIT=.FALSE.
        PROCEED=.FALSE.
        tSimplex(M+1,0:N)= -SUM(tSimplex(1:M,0:N),DIM=1)
      END IF
      !
      KP=L1(iMaxLoc_R(tSimplex(M+1,L1(1:NL1))))
      bMax=tSimplex(M+1,KP)
      
      !---------------------------------------------------------phase 1A
      Phase1A: DO
        
        IF (bMax<=EPS.AND.tSimplex(M+1,0)<-EPS) THEN
          
          iError=-1
          RETURN !------------------------------------------------RETURN
        
        ELSE IF (bMax <= EPS .AND. tSimplex(M+1,0) <= EPS) THEN
        
          !M12=1
          DO IP=1,M
            IF (iPosV(IP)==IP+N) THEN
              KP=   L1(iMaxLoc_R(ABS(tSimplex(IP,L1(1:NL1)))))
              bMax= tSimplex(IP,KP)
              IF (bMax>Zero) EXIT Phase1A
            END IF
          END DO
          PROCEED=.TRUE.
          
          EXIT Phase1
          
        END IF
        
        CALL SIMP1(IP) ; n1= n1 +1
        
        IF (IP==0) THEN
          iError=-1
          RETURN !------------------------------------------------RETURN
        END IF
        
        EXIT Phase1A
        
      ENDDO Phase1A
      !--------------------------------------------------------/phase 1A
      
      !---------------------------------------------------------phase 1B
      Phase1B: DO
        !
        CALL SIMP2(M+1,N) ; n2= n2 +1
        IF (iPosV(IP) >= N+1) THEN
          K=iFirstLoc(L1(1:NL1)==  KP)
          NL1=NL1-1
          L1(K:NL1)=L1(K+1:NL1+1)
        ELSE
          IF (iPosV(IP) < N+1) EXIT Phase1B
          KH=iPosV(IP)-N
          !IF (L3(KH)==0) EXIT Phase1B
          !L3(KH)=0
        END IF
        tSimplex(M+1,KP)=   tSimplex(M+1,KP)  +One
        tSimplex(0:M+1,KP)=-tSimplex(0:M+1,KP)
        EXIT Phase1B
      ENDDO Phase1B
      !--------------------------------------------------------/phase 1B
      
      CALL Swap_I(IZROV(KP),iPosV(IP))
      
      IF (PROCEED) EXIT Phase1
      
    ENDDO Phase1
    !-----------------------------------------------------------/phase 1
    !
    !------------------------------------------------------------phase 2
    PHASE2: DO
      !
      KP=L1(iMaxLoc_R(tSimplex(0,L1(1:NL1))))
      bMax=tSimplex(0,KP)
      
      IF (bMax<=Zero) THEN
        iError=0
        RETURN !--------------------------------------------------RETURN
      END IF
      
      CALL SIMP1(IP)
      n1= n1 +1
      IF (IP==0) THEN
        iError=1
        RETURN !--------------------------------------------------RETURN
      ENDIF
      
      CALL SIMP2(M,N) ; n2= n2 +1
      CALL Swap_I(IZROV(KP),iPosV(IP))
      
      IF (.NOT. PROCEED) CYCLE DoMain
      
    END DO PHASE2
    !-----------------------------------------------------------/phase 2
    
    EXIT DoMain
    
  ENDDO DoMain
  !
  DEALLOCATE(L1,L2) !,L3)
  !
CONTAINS
  !
  SUBROUTINE SIMP1(IP_)
    USE M_Numeric_Tools, ONLY: iFirstLoc
    
    INTEGER, INTENT(OUT)::IP_
    
    INTEGER :: I,J,K
    REAL(dp):: Q,Q0,Q1,QP
    
    IP_=0
    I=iFirstLoc(tSimplex(L2(1:NL2),KP) < -EPS)
    IF (I > NL2) RETURN
    
    Q1=-tSimplex(L2(I),0)/tSimplex(L2(I),KP)
    IP_=L2(I)
    DO I=I+1,NL2
      
      J=L2(I)
      IF (tSimplex(J,KP) < -EPS) THEN
        Q= -tSimplex(J,0) /tSimplex(J,KP)
        IF (Q < Q1) THEN
          IP_= J
          Q1=  Q
        ELSEIF (Q==Q1) THEN
          DO K=1,N
            QP= -tSimplex(IP_,K) /tSimplex(IP_,KP)
            Q0= -tSimplex(J,  K) /tSimplex(J,  KP)
            IF (Q0/=QP) EXIT
          ENDDO
          IF (Q0<QP) IP_=J
        ENDIF
      ENDIF
      
    ENDDO
    
    RETURN
  END SUBROUTINE SIMP1

  SUBROUTINE SIMP2(I1,K1)
    INTEGER, INTENT(IN):: I1,K1
    !
    INTEGER            :: IP1,KP1
    REAL(dp)           :: PIV
    INTEGER, DIMENSION(K1) :: iCol
    INTEGER, DIMENSION(I1) :: iRow
    INTEGER, DIMENSION(1:MAX(I1,K1)+1):: iTmp
    
    IP1= IP+1
    KP1= KP+1
    !
    PIV= One /tSimplex(IP1-1,KP1-1)
    !
    iTmp(1:K1+1)= Arth_I(1,1,K1+1)
    iCol= PACK(iTmp(1:K1+1),iTmp(1:K1+1) /= KP1)
    !
    iTmp(1:I1+1)=  Arth_I(1,1,I1+1)
    iRow= PACK(iTmp(1:I1+1),iTmp(1:I1+1) /= IP1)
    !
    tSimplex(iRow-1,KP1-1)=  tSimplex(iRow-1,KP1-1)*PIV
    tSimplex(iRow-1,iCol-1)= tSimplex(iRow-1,iCol-1) &
    &        - OuterProd_R(tSimplex(iRow-1,KP1-1),tSimplex(IP1-1,iCol-1))
    tSimplex(IP1-1,iCol-1)= -tSimplex(IP1-1,iCol-1)*PIV
    tSimplex(IP1-1,KP1-1)=   PIV
    !
    RETURN
  END SUBROUTINE SIMP2

END SUBROUTINE Simplex_Calc

FUNCTION Arth_I(First,Increment,dimN) !-> array(1:dimN)
!--
!-- build array(dimN) of INTEGER (First,First+Increment,First+2*Increment,...)
!-- ex: Arth_I(1,1,dimN) gives 1,2,3,4,...dimN
!--
  INTEGER, INTENT(IN)  :: First,Increment,dimN
  INTEGER, DIMENSION(dimN):: Arth_I
  INTEGER :: K,K2,TEMP
  !
  IF(dimN>0) Arth_I(1)=First
  !
  !IF (dimN<=nPar_Arth) THEN !INTEGER, PARAMETER :: nPar_Arth=16,NPAR2_ARTH=8
  !
  IF (dimN<=16) THEN
    DO K=2,dimN; Arth_I(K)=Arth_I(K-1)+Increment; END DO
  ELSE
  !DO K=2,NPAR2_ARTH; Arth_I(K)=Arth_I(K-1)+Increment; END DO
    DO K=2,8; Arth_I(K)=Arth_I(K-1)+Increment; END DO
    TEMP=Increment*8 !NPAR2_ARTH
    K=8 !NPAR2_ARTH
    DO
      IF(K>=dimN) EXIT
      K2=   K+K
      Arth_I(K+1:MIN(K2,dimN))= TEMP+Arth_I(1:MIN(K,dimN-K))
      Temp= Temp +Temp
      K=  K2
    ENDDO
  ENDIF
  RETURN
END FUNCTION Arth_I
  
END MODULE M_Simplex_Calc

!Simplex, notes from NumRecp:
!
! linear programming= 
! maximize the linear FUNCTION 
!  Z=  a(0,1).x(1) + ... + a(0,n).x(n)
! subject to 
!  primary constraints: a(0,1:n)>=0
! and simultaneously subject to 
!  M=  m1 + m2 + m3 additional constraints,
!  m1 of them of the form: a(i,1).x1 + ... + a(i,N).x(N) >= b(i)   //i=1..m1
!  m2 of them of the form: a(j,1).x1 + ... + a(j,N).x(N) <= b(j)   //j=m1+1..m1+m2
!  m3 of them of the form: a(k,1).x1 + ... + a(k,N).x(N) <= b(k)   //k=m1+m2+1..m1+m2+m3

! On output, the tableau a is indexed by two RETURNed arrays of INTEGERs
! iposv(j) CONTAINS, for j= 1...M, 
! the number i whose original variable x(i) is now represented by row j+1 of a
! These are thus the left-hand variables in the solution
! (The first row of a is of course the z-row.)
! A value i > N indicates that the variable is a y(i) rather than an x(i), x(N+j)==y(j)
! Likewise, izrov(j) CONTAINS, for j= 1..N,
! the number i whose original variable x(i) is now a right-hand variable, represented by column j+1 of a. 
! These variables are all zero in the solution. 
! The meaning of i > N is the same as above, 
! except that i >N +m1 +m2 denotes an artificial or slack variable 
! which was used only internally and should now be entirely ignored.

