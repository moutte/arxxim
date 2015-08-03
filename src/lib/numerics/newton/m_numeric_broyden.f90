MODULE M_Numeric_Broyden
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_
  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC:: Broyden
  !
  LOGICAL,PUBLIC:: TestJacob=.FALSE., ShoResult=.FALSE.
  !
CONTAINS

SUBROUTINE Broyden( &
& vX,      & !inout= the initial guess, and the root RETURNed
& Residual, Jacobian, &
& TolF,    & !in=    convergence criterion on FUNCTION values 
& TolX,    & !in=    convergence criterion on dx
!!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
& bFinDIF, & !in=    USE numeric Jacobian
& MaxIts,  & !in=    maximum number of iterations
& Error_F, & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Gradient, & !out=   
& Nits,    & !out=   number of iterations
& Check,   & !out=   IF Check, should check convergence
& iErr)      !out=   error code
!.NumericalRecipes, Broyden’s method embedded in a globally convergent strategy. 
!The LENgth N vector of FUNCTIONs to be zeroed, called fvec in the routine below,
!is RETURNed by a USEr-supplied routine that must be called funcv and have the declaration FUNCTION funcv(x).
!USE SUBROUTINE Jacobian_Numeric and FUNCTION fmin from newt
!
!check=
!-- false on a normal RETURN,
!-- true IF the routine has converged to a local minimum of the FUNCTION fmin
!-- or IF Broyden’s method can make no further progress. 
!In this CASE try restarting from a dIFferent initial guess.
!
  USE M_Numeric_Const,      ONLY: Ln10
  USE M_Numeric_Tools,ONLY: OuterProd_R,VAbs,Unit_Matrix,Jacobian_Numeric
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR,fNewt_I
  USE M_IoTools,    ONLY: OutStrVec
  !USE M_NRBroyden,ONLY: &
  !& Get_Diag_Rv,Put_Diag_Rv,Lower_Triangle, &
  !& QRDCMP,QRUPDT,RSOLV
  !!USE M_IOTools,ONLY: OutStrVec !for debugging
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  REAL(dp),INTENT(IN)   :: TolF,TolX !!,TolMin
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X,Gradient
  LOGICAL, INTENT(OUT)  :: Check
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))    :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
  ENDINTERFACE
  !
  REAL(dp),PARAMETER:: STPMX=  100._dp, TolMin= 1.0E-9_dp
  INTEGER  :: I,ITS,K,N
  REAL(dp) :: F,FOLD,STPMAX,EPS
  !
  REAL(dp),DIMENSION(SIZE(vX)):: vFunc,vFuncOld,C,D,G,P,S,T,W,XOLD
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: QT,R
  LOGICAL  :: RESTRT,SING
  !
  EPS= TolX
  N=SIZE(vX)
  vFunc= Residual(vX)
  F=    DOT_PRODUCT(vFunc,vFunc)/2.0D0
  !
  !IF (MAXVAL(ABS(vFunc(:))) < 0.01_dp*NewtTolF) THEN
  !  iErr=1; RETURN !__________________________________________________________RETURN
  !ENDIF
  !
  STPMAX= STPMX*MAX(VABS(vX(:)),REAL(N,dp))
  RESTRT= .TRUE.
  !
  iErr= -1
  DO0: DO ITS=1,MAXITS
    nIts=Its
    !
  fNewt_I= fNewt_I +1
  IF(fNewtF>0) CALL OutStrVec(fNewtF,vX(:)/Ln10,Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  IF(fNewtR>0) CALL OutStrVec(fNewtR,vFunc(:),  Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
    !IF(DebNewt .AND. bOpenNewt) THEN
    !  CALL OutStrVec(fNewt1,vX(1:SIZE(vX))/Ln10,Opt_I=Its,Opt_C="F")
    !  CALL OutStrVec(fNewt2,vFunc(1:SIZE(vX)),  Opt_I=Its,Opt_C="F")
    !ENDIF
    !
    IF (RESTRT) THEN
      !CALL Jacobian_Numeric(vX,vFunc,R)
      IF(bFinDIF) THEN
        CALL Jacobian_Numeric(Residual,vX,vFunc,R)
      ELSE
        CALL Jacobian(vX,R)
      ENDIF
      !
      CALL QRDCMP(R,C,D,SING)
      !
      IF (SING) THEN !CALL Stop_('SINGULAR JACOBIAN IN Residual')
        iErr=-2; EXIT DO0
      ENDIF
      !
      CALL UNIT_MATRIX(QT)
      !
      DO K=1,N-1
        IF(C(K)/=Zero) &
        & QT(K:N,:)= QT(K:N,:) &
        &          - OUTERPROD_R(R(K:N,K),MATMUL(R(K:N,K),QT(K:N,:)))/C(K)
      ENDDO
      WHERE (LOWER_TRIANGLE(N,N)) R(:,:)=Zero
      CALL Put_Diag_Rv(D(:),R(:,:))
    ELSE
      S(:)=vX(:)-XOLD(:)
      DO I=1,N; T(I)=DOT_PRODUCT(R(I,I:N),S(I:N)); ENDDO
      !
      W(:)=vFunc(:)-vFuncOld(:) -MATMUL(T(:),QT(:,:))
      !
      WHERE (ABS(W(:)) < EPS*(ABS(vFunc(:))+ABS(vFuncOld(:)))) W(:)=Zero
      IF (ANY(W(:) /= Zero)) THEN
        T(:)=MATMUL(QT(:,:),W(:))
        S(:)=S(:)/DOT_PRODUCT(S,S)
        CALL QRUPDT(R,QT,T,S)
        D(:)=GET_DIAG_RV(R(:,:)) !DO J=1,SIZE(MAT,1) GET_DIAG_RV(J)=MAT(J,J); ENDDO
        IF (ANY(D(:)==Zero)) THEN !CALL Stop_('R SINGULAR IN Residual')
          iErr=-2; EXIT DO0
        ENDIF
      ENDIF
    ENDIF
    P(:)=-MATMUL(QT(:,:),vFunc(:))
    DO I=1,N; G(I)=-DOT_PRODUCT(R(1:I,I),P(1:I)); ENDDO
    XOLD(:)=vX(:)
    vFuncOld(:)=vFunc(:)
    FOLD=F
    CALL RSOLV(R,D,P)
    !
    CALL LNSRCH(XOLD,FOLD,G,STPMAX,P,vX,F,CHECK) !-> this updates also vFunc !!
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Gradient= MAXVAL(ABS(G(:))*MAX(ABS(vX(:)),One)/MAX(F,0.5_dp*N))
    Delta_X= MAXVAL((ABS(vX(:)-XOLD(:)))/MAX(ABS(vX(:)),One))
    !
    IF (Error_F < TolF) THEN
      iErr=0; RETURN !__________________________________________________________RETURN
    ENDIF
    IF (CHECK) THEN
      IF(RESTRT .OR. Gradient < TolMin) THEN
        iErr=-3; RETURN !________________________________________________________RETURN
      ENDIF
      RESTRT=.TRUE.
    ELSE
      RESTRT=.FALSE.
      IF(Delta_X < TolX) THEN
        iErr=-4; RETURN !________________________________________________________RETURN
      ENDIF
    ENDIF
  ENDDO DO0
  RETURN

CONTAINS

SUBROUTINE LNSRCH( & !from NR
& vXOld, & !IN= an N-DIMENSIONal point
& rFOld, & !IN= the value of the FUNCTION at point vXOld
& vGrad, & !IN= the gradient of the FUNCTION at that point
& rStpMax, & !IN=limits the LENgth of the steps
& vDX,   & !INOUT= a direction
& vXnew, & !OUT= a new point vXnew(1:N) along the direction vDX from vXOld
& rFnew, & !OUT= the new FUNCTION value
& CHECK)   !FALSE= normal EXIT, TRUE= vXnew is too CLOSE to vXOld 
!Given 
!  vXOld(1:N),IN= an N-DIMENSIONal point, 
!  rFOld(1:N),IN= the value of the FUNCTION at vXOld,
!  vGrad(1:N),IN= the gradient of the FUNCTION at vXOld,
!  vDX(1:N),INOUT= a direction vDX
!finds
!  a new point vXnew(1:N) along the direction vDX from vXOld
!  WHERE the FUNCTION FUNC has decreased "sufficiently".
!rFnew(1),OUT= the new FUNCTION value
!rStpMax,IN=
!  an input quantity that limits the LENgth of the steps 
!  so that you DO not try to evaluate the FUNCTION in regions WHERE it is undefined 
!  or subject to overflow. 
!vDX,INOUT= usually the Newton direction.
!Check,OUT=  
!  FALSE on a normal EXIT, 
!  TRUE when vXnew is too CLOSE to vXOld,
!  in a minimization algorithm, this usually signals convergence and can be ignored
!  in a zero-finding algorithm, the calling program should check whether the convergence is spurious.
!Parameters
!  ALF=  ensures sufficient decrease in FUNCTION value;
!  TolX= the convergence criterion on delta(vXnew)
  !
  REAL(dp),DIMENSION(:), INTENT(IN)    :: vXOld,vGrad
  REAL(dp),DIMENSION(:), INTENT(INOUT) :: vDX
  REAL(dp),              INTENT(IN)    :: rFOld,rStpMax
  REAL(dp),DIMENSION(:), INTENT(OUT)   :: vXnew
  REAL(dp),              INTENT(OUT)   :: rFnew
  LOGICAL,               INTENT(OUT)   :: CHECK
  !
  REAL(dp),PARAMETER ::&
  & TolX=EPSILON(vXnew), &
  & ALF= EPSILON(vXnew) !1.0E-4_dp  ! !1.0E-6_dp
  REAL(dp)::&
  & A,Alam,Alam2,AlaMin,B,DISC,F2,PAbs,&
  & RHS1,RHS2,Slope,tmpLam,x
  !
  CHECK=.FALSE.
  PAbs=vAbs(vDX(:))
  IF (PAbs>rStpMax) vDX(:)=vDX(:)*rStpMax/PAbs !Scale IF attempted step is too big.
  !
  Slope=DOT_PRODUCT(vGrad,vDX) !______________Compute lambda_min
  IF (slope>=Zero) CALL Stop_('rounDOff problem in lnsrch') !____________________________STOP
  AlaMin=TolX / MAXVAL(ABS(vDX(:))/MAX(ABS(vXOld(:)),One))
  Alam=One !____________________________Always try full Newton step first
  DO !__________________________________Start of iteration loop
    vXnew(:)= vXOld(:) +Alam *vDX(:)
    !rFnew=FUNC(vXnew) !rFnew=scalar,vXnew=Vector
    vFunc= Residual(vXnew)
    rFnew= DOT_PRODUCT(vFunc,vFunc)/2.0D0
    IF (Alam<AlaMin) THEN !_____________Convergence on vXnew
      vXnew(:)=vXOld(:) !For zero finding, the calling program should verIFy the convergence.
      CHECK=.TRUE.
      RETURN !__________________________________________________RETURN!!!
    ELSEIF (rFnew <= rFOld+ALF*Alam*Slope) THEN
      RETURN !sufficient FUNCTION decrease______________________RETURN!!!
    ELSE                        !Backtrack.
      IF (Alam==One) THEN; tmpLam=-Slope/(2.0_dp*(rFnew-rFOld-Slope)) !First time.
      ELSE !____________________________Subsequent backtracks
        RHS1=rFnew -rFOld -Alam *Slope
        RHS2=F2-rFOld -Alam2*Slope
        A= (       RHS1/Alam/Alam -      RHS2/Alam2/Alam2) / (Alam-Alam2)
        B= (-Alam2*RHS1/Alam/Alam + Alam*RHS2/Alam2/Alam2) / (Alam-Alam2)
        IF (A==Zero) THEN; tmpLam=-Slope/(2.0_dp*B)
        ELSE
          DISC=B*B-3.0_dp*A*Slope
          IF (DISC<Zero)   THEN; tmpLam=0.5_dp*Alam
          ELSEIF (B<=Zero) THEN; tmpLam=(-B+SQRT(DISC))/(3.0_dp*A)
          ELSE                 ; tmpLam=-Slope/(B+SQRT(DISC))
          ENDIF
        ENDIF
        IF (tmpLam>0.5_dp*Alam) tmpLam=0.5_dp*Alam
      ENDIF
    ENDIF
    Alam2=Alam
    Alam=MAX(tmpLam,0.1_dp*Alam)
    F2=rFnew
  ENDDO
ENDSUBROUTINE LNSRCH

ENDSUBROUTINE Broyden

!!MODULE M_NRBroyden
FUNCTION Get_Diag_Rv(MAT) !->vector of diagonal elements, USEd in Broyden method
  USE M_Trace,ONLY:Stop_
  REAL(dp),DIMENSION(:,:),INTENT(IN)::MAT
  REAL(dp),DIMENSION(SIZE(MAT,1))::Get_Diag_Rv
  !
  INTEGER::J
  !
  IF(SIZE(MAT,1)/=SIZE(MAT,2)) CALL Stop_("SIZE !!! in Get_Diag_Rv") !check that MAT is Square Matrix
  DO J=1,SIZE(MAT,1); Get_Diag_Rv(J)=MAT(J,J); ENDDO
ENDFUNCTION Get_Diag_Rv

SUBROUTINE Put_Diag_Rv(DIAGV,MAT) !USEd in Broyden method
  USE M_Trace,ONLY:Stop_
  REAL(dp),DIMENSION(:),  INTENT(IN)   :: DIAGV
  REAL(dp),DIMENSION(:,:),INTENT(INOUT):: MAT
  INTEGER::J,N
  !
  IF(SIZE(DIAGV) /= MIN(SIZE(MAT,1),SIZE(MAT,2))) CALL Stop_("SIZE !!! PUT_DIAG_RV")
  N=SIZE(DIAGV)
  DO J=1,N; MAT(J,J)=DIAGV(J); ENDDO
ENDSUBROUTINE Put_Diag_Rv

FUNCTION ARTH_I(First,Increment,N) !USEd by Broyden
  INTEGER,INTENT(IN)  :: First,Increment,N
  INTEGER,DIMENSION(N):: ARTH_I
  INTEGER             :: K,K2,TEMP
  IF(N>0) ARTH_I(1)=FIRST
  IF(N<=16) THEN !(N<=NPAR_ARTH)
    DO K=2,N; Arth_I(K)=Arth_I(K-1)+Increment; ENDDO
  ELSE
    DO K=2,8; Arth_I(K)=Arth_I(K-1)+Increment; ENDDO !K=2,NPAR2_ARTH
    Temp=Increment*8 !NPAR2_ARTH
    K=8 !NPAR2_ARTH
    DO
      IF (K >= N) EXIT
      K2=K+K
      Arth_I(K+1:MIN(K2,N))=Temp+Arth_I(1:MIN(K,N-K))
      Temp=Temp+Temp
      K=K2
    ENDDO
  ENDIF
ENDFUNCTION ARTH_I

FUNCTION OuterDIFf_R(A,B)
  REAL(dp),DIMENSION(:), INTENT(IN)  :: A,B
  REAL(dp),DIMENSION(SIZE(A),SIZE(B)):: OuterDIFf_R
  OuterDIFf_R= SPREAD(A,DIM=2,NCOPIES=SIZE(B)) &
  &          - SPREAD(B,DIM=1,NCOPIES=SIZE(A))
ENDFUNCTION OuterDIFf_R

FUNCTION OuterdIFf_I(A,B)
  INTEGER,DIMENSION(:),INTENT(IN) :: A,B
  INTEGER,DIMENSION(SIZE(A),SIZE(B)):: OuterdIFf_I
  !
  OUTERDIFF_I= SPREAD(A,DIM=2,NCOPIES=SIZE(B)) &
  &          - SPREAD(B,DIM=1,NCOPIES=SIZE(A))
ENDFUNCTION OuterdIFf_I

FUNCTION Lower_Triangle(J,K,EXTRA) !matricial mask (J,K) of LOGICAL
  INTEGER,          INTENT(IN):: J,K
  INTEGER,OPTIONAL, INTENT(IN):: EXTRA
  LOGICAL,DIMENSION(J,K)      :: Lower_Triangle
  INTEGER :: N
  N=0; IF (PRESENT(EXTRA)) N=EXTRA
  Lower_Triangle=(OuterdIFf_I(Arth_I(1,1,J),Arth_I(1,1,K)) > -N)
ENDFUNCTION Lower_Triangle

FUNCTION Assert_Eq4(N1,N2,N3,N4,Str) !"NR"
  USE M_Trace,ONLY:Stop_
  CHARACTER(LEN=*), INTENT(IN) :: Str
  INTEGER, INTENT(IN)          :: N1,N2,N3,N4
  INTEGER                      :: Assert_Eq4
  IF (N1==N2 .AND. N2==N3 .AND. N3==N4) THEN
    Assert_Eq4=N1
  ELSE
    CALL Stop_(Str)
  ENDIF
END FUNCTION Assert_Eq4

FUNCTION Assert_EqN(NN,Str) !"NR"
  USE M_Trace,ONLY:Stop_
  CHARACTER(LEN=*),    INTENT(IN):: Str
  INTEGER,DIMENSION(:),INTENT(IN):: NN
  INTEGER                        :: Assert_Eqn
  IF (ALL(NN(2:)==NN(1))) THEN; Assert_EqN=NN(1)
  ELSE;                         CALL Stop_(Str)
  ENDIF
ENDFUNCTION Assert_EqN

SUBROUTINE RSolv(a,d,b) !right triangular backsubstitution
  REAL(dp),DIMENSION(:,:),INTENT(IN)   :: a
  REAL(dp),DIMENSION(:),  INTENT(IN)   :: d
  REAL(dp),DIMENSION(:),  INTENT(INOUT):: b
  INTEGER:: i,n
  n=assert_eq4(SIZE(a,1),SIZE(a,2),SIZE(b),SIZE(d),'rsolv')
  b(n)=b(n)/d(n)
  DO i=n-1,1,-1; b(i)=(b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/d(i); ENDDO
ENDSUBROUTINE RSolv

SUBROUTINE QRUpdt(r,qt,u,v) !from "NR" 
  USE M_Numeric_Tools,ONLY:IFirstloc
  REAL(dp),DIMENSION(:,:),INTENT(INOUT):: r,qt
  REAL(dp),DIMENSION(:),  INTENT(INOUT):: u
  REAL(dp),DIMENSION(:),  INTENT(IN)   :: v
  INTEGER::i,k,n
  n=Assert_EqN((/SIZE(r,1),SIZE(r,2),SIZE(qt,1),SIZE(qt,2),SIZE(u),SIZE(v)/),'qrupdt')
  k= n+1 -IFirstloc(u(n:1:-1) /= Zero)
  IF (k < 1) k=1
  DO i=k-1,1,-1
    CALL Rotate(r,qt,i,u(i),-u(i+1))
    u(i)=pythag(u(i),u(i+1))
  ENDDO
  r(1,:)=r(1,:)+u(1)*v
  DO i=1,k-1; CALL Rotate(r,qt,i,r(i,i),-r(i+1,i)); ENDDO
ENDSUBROUTINE QRUpdt

FUNCTION Pythag(a,b) !from "NR"
  USE M_Kinds
  REAL(dp),INTENT(IN):: a,b
  REAL(dp)           :: pythag
  REAL(dp):: absa,absb
  absa=ABS(a)
  absb=ABS(b)
  IF (absa > absb) THEN
    pythag=absa*SQRT(One+(absb/absa)**2)
  ELSE
    IF (absb==Zero) THEN; pythag=Zero
    ELSE                ; pythag=absb*SQRT(One+(absa/absb)**2)
    ENDIF
  ENDIF
ENDFUNCTION Pythag

SUBROUTINE Rotate(r,qt,i,a,b) !from "NR"
  REAL(dp),DIMENSION(:,:),TARGET,INTENT(INOUT):: r,qt
  INTEGER,                        INTENT(IN)   :: i
  REAL(dp),                     INTENT(IN)   :: a,b
  REAL(dp),DIMENSION(SIZE(r,1)):: temp
  INTEGER:: n
  REAL:: c,fact,s
  n=assert_eq4(SIZE(r,1),SIZE(r,2),SIZE(qt,1),SIZE(qt,2),'rotate')
  IF(a== Zero) THEN
    c=   Zero
    s=   SIGN(One,b)
  ELSEIF (ABS(a)>ABS(b)) THEN
    fact=b/a
    c=   sign(One/SQRT(One+fact**2),a)
    s=   fact*c
  ELSE
    fact=a/b
    s=   sign(One/SQRT(One+fact**2),b)
    c=   fact*s
  ENDIF
  temp(i:n)= r(i,i:n)
  r(i,  i:n)= c *temp(i:n) -s *r(i+1,i:n)
  r(i+1,i:n)= s *temp(i:n) +c *r(i+1,i:n)
  temp=       qt(i,:)
  qt(i,  :)= c *temp -s *qt(i+1,:)
  qt(i+1,:)= s *temp +c *qt(i+1,:)
ENDSUBROUTINE Rotate

SUBROUTINE QRDcmp(a,c,d,sing) !QR decomposition, from "NR"
  USE M_Numeric_Tools,ONLY:OuterProd_R,VAbs
  REAL(dp),DIMENSION(:,:),INTENT(INOUT):: a
  REAL(dp),DIMENSION(:),  INTENT(OUT)  :: c,d
  LOGICAL,                 INTENT(OUT)  :: sing
  INTEGER :: k,n
  REAL(dp):: scale,sigma
  n=Assert_Eq4(SIZE(a,1),SIZE(a,2),SIZE(c),SIZE(d),'qrdcmp')
  sing=.false.
  DO k=1,n-1
    scale=MAXVAL(ABS(a(k:n,k)))
    IF (scale==Zero) THEN
      sing= .true.
      c(k)= Zero
      d(k)= Zero
    ELSE
      a(k:n,k)=a(k:n,k)/scale
      sigma=   SIGN(vabs(a(k:n,k)),a(k,k))
      a(k,k)=  a(k,k)+sigma
      c(k)=    sigma*a(k,k)
      d(k)=   -scale*sigma
      a(k:n,k+1:n)= a(k:n,k+1:n) &
      &           - OuterProd_R(a(k:n,k),MATMUL(a(k:n,k),a(k:n,k+1:n)))/c(k)
    ENDIF
  ENDDO
  d(n)=a(n,n)
  IF (d(n)==Zero) sing=.true.
ENDSUBROUTINE QRDcmp
!!ENDMODULE M_NRBroyden

ENDMODULE M_Numeric_Broyden

