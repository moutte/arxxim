MODULE M_Numeric_Newton
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_
  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC:: NewtLnsrch
  PUBLIC:: Newton_Press
  PUBLIC:: Newton_Walker
  PUBLIC:: Newton_Kelley
  PUBLIC:: Newton
  PUBLIC:: NewtonChess
  !
  LOGICAL,PUBLIC:: TestJacob=.FALSE.
  LOGICAL,PUBLIC:: ShoResult=.FALSE.
  !
  LOGICAL:: DeltaX_Ok= .TRUE.
  !
CONTAINS

SUBROUTINE NewtLnsrch( & !from Press et al, Numerical Recipes
& vX,       & !inout= the initial guess, and the root returned
!!& vLPos,  & !in=    enforce vX(vLPos=TRUE)>0
& Residual, & !
& Jacobian, & !
& Converge, & !
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
!!& TolMin, & !in=    whether spurious convergence to a minimum of fmin has occurred
& bFinDIF,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   MAXVAL(ABS(vFunc(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Gradient, & !out=
& Nits,     & !out=   number of iterations
& Check,    & !out=   IF Check, should check convergence
& iErr)      !out=   error code
!
!------------------------------------------------------------error codes
! iErr= 0 : convergence reached -> OK
! iErr=-1 : MaxIts reached without convergence
! iErr=-2 : singular jacobian
! iErr=-3 : roundoff problem in linesearch
! iErr=-4 : no convergence on vFunc, problem with gradient ?
! iErr=-5 : no convergence on vFunc, vX very close to previous vX -> stationary point ??
!-----------------------------------------------------------/error codes
!
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  USE M_Numeric_Const,ONLY: Ln10 !for debugging
  USE M_IoTools,      ONLY: GetUnit, OutStrVec !for debugging
  !
  !!LOGICAL, DIMENSION(:), INTENT(IN)   :: vLPos
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
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),              INTENT(IN) :: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(vF,vTolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: vF(:),vTolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  !INTEGER, PARAMETER :: MAXITS=NewtMaxIts
  !NewtTolF=   convergence criterion on FUNCTION values
  !NewtTolMin= criterion for spurious convergence
  !NewtTolX=   convergence criterion on dx
  !
  REAL(dp),PARAMETER:: STPMX= 100.0
  !scaled maximum step LENgth allowed in line searches
  REAL(dp),PARAMETER:: Alf= Epsilon(vX) !1.0D-12
  !
  LOGICAL :: bSingul
  INTEGER :: ITS
  REAL(dp):: D,F,rFOld,rStpMax,rSlope,Norm_vDX !,FMin
  !
  REAL(dp),DIMENSION(SIZE(vX)):: vFunc,vGrad,vDX,vXOld
  INTEGER, DIMENSION(SIZE(vX)):: vIndex
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJacob,tJacobNum
  !
  REAL(dp):: A,Alam,Alam2,AlaMin,B,DISC,F2
  REAL(dp):: RHS1,RHS2,tmpLam
  !
  vFunc= Residual(vX)
  F= 0.5_dp*DOT_PRODUCT(vFunc,vFunc)
  !
  iErr=-1
  !Error_F= MAXVAL(ABS(vFunc(:))) !-> error on FUNCTION
  !IF (Error_F < 0.01_dp*TolF) THEN !Test for initial guess being a root
  !  nIts=1
  !  iErr=1 !-----------Use more stringent test than simply NewtTolF ???
  !  RETURN
  !ENDIF
  !
  IF(ShoResult) PRINT '(/,A,/)',"iTs,Error_F,Gradient,Delta_X"
  !
  rStpMax=STPMX*MAX(vAbs(vX(:)),REAL(SIZE(vX),dp)) !Calculate StpMax for line searches
  !
  DO its=1,MaxIts
    !--------------------------------------------start of iteration loop
    !
    nIts=its
    !
    CALL Newton_Sho(vX,vFunc,Its)
    !
    IF(bFinDIF) THEN
      CALL Jacobian_Numeric(Residual,vX,vFunc,tJacob)
    ELSE
      CALL Jacobian(vX,tJacob)
      IF(TestJacob) THEN
      !for test, calc. numeric Jacobian and compare with analytical one
        CALL Jacobian_Numeric(Residual,vX,vFunc,tJacobNum)
        CALL TestJacob_Sho(tJacob,tJacobNum)
        TestJacob=.FALSE. !calculate ONLY on first CALL to NewtLnsrch
      ENDIF
    ENDIF
    !
    vGrad(:)= MATMUL(vFunc(:),tJacob(:,:)) ! compute grad(f) for the line search.
    vXOld(:)= vX(:)                        ! store vX
    rFOld=    F                            ! store F
    vDX(:)=  -vFunc(:)                     ! right-hand side for linear equations.
    !
    !--- solve linear equations by LU decomposition
    CALL LU_Decomp(tJacob,vIndex,D,bSingul)
    !
    IF(bSingul) THEN
      iErr= -2 !-> "SINGULAR JACOBIAN, cf Log File"
      CALL Jacobian_Sho(tJacob)
      RETURN !System not solved --------------------------------- RETURN
    ENDIF
    !
    CALL LU_BakSub(tJacob,vIndex,vDX)
    !---/ solve linear equations by LU decomposition
    !
    Norm_vDX= vAbs(vDX(:))
    !-- Scale if attempted step is too big --
    IF (Norm_vDX>rStpMax) vDX(:)= rStpMax *vDX(:) /Norm_vDX
    !
    rSlope= DOT_PRODUCT(vGrad,vDX)
    !
    IF (rSlope>=Zero) THEN
      IF(iDebug>0) THEN
        WRITE(fTrc,'(/,A,/)') "ROUNDOFF PROBLEM -> vGrad, vDX ="
        CALL OutStrVec(fTrc,vGrad(1:SIZE(vX)),OPt_C="G")
        CALL OutStrVec(fTrc,vDX(1:SIZE(vX)),  OPt_C="G")
      ENDIF
      iErr= -3 !-> 'roundoff problem in linesearch, cf log file'
      RETURN ! System not solved -------------------------------- RETURN
    ENDIF
    !
    !----------------------------------------linesearch and backtracking
    !
    CHECK= .FALSE.
    AlaMin= TolX /MAXVAL(ABS(vDX(:))/MAX(ABS(vXOld(:)),One))
    Alam=   One !always try full Newton step first
    !
    DoLineSearch: DO !Start of iteration loop
      !
      !!write(51,'(I3,G15.6)') Its, Alam
      !
      vX(:)= vXOld(:) +Alam*vDX(:)
      !
      vFunc= Residual(vX)
      F= DOT_PRODUCT(vFunc,vFunc)/2.0D0
      !
      IF (Alam<AlaMin) THEN !Convergence on vX
        !
        vX(:)=vXOld(:)
        !-- For zero finding, the calling program should verify the convergence
        CHECK=.TRUE.
        EXIT DoLineSearch !-----------------------------------------EXIT
        !
      ELSEIF (F <= rFOld + ALF*Alam*rSlope) THEN
        !
        EXIT DoLineSearch !sufficient FUNCTION decrease ----------- EXIT
        !
      ELSE !==< Backtrack.
        !
        IF (Alam==One) THEN !---------------------------------First time
          tmpLam= -rSlope /(2.0_dp*(F-rFOld-rSlope))
          !
        ELSE !-------------------------------------Subsequent backtracks
          !PRINT *,Alam2   ; PAUSE_
          RHS1= F -rFOld -Alam *rSlope
          RHS2= F2-rFOld -Alam2*rSlope
          A= (       RHS1/Alam**2 -      RHS2/Alam2**2) /(Alam-Alam2)
          B= (-Alam2*RHS1/Alam**2 + Alam*RHS2/Alam2**2) /(Alam-Alam2)
          !
          IF (A==Zero) THEN
            tmpLam= -rSlope/(2.0_dp*B)
            !
          ELSE
            IF(B<-1.0D100) B=-1.0D100
            IF(B> 1.0D100) B= 1.0D100
            !PRINT *,A     ; PAUSE_
            !PRINT *,B     ; PAUSE_
            !PRINT *,"B*B",B*B   ; PAUSE_
            DISC= B*B -3.0_dp*A*rSlope
            !PRINT *,DISC  ; PAUSE_
            IF     (DISC<Zero) THEN  ; tmpLam= 0.5_dp*Alam
            ELSEIF (B<=Zero)   THEN  ; tmpLam= (-B+SQRT(DISC))/3.0_dp/A
            ELSE                     ; tmpLam= -rSlope/(B+SQRT(DISC))
            ENDIF
          ENDIF
          !
          IF (tmpLam>0.5_dp*Alam) tmpLam=0.5_dp*Alam
          !
        ENDIF
        !
      ENDIF
      !
      Alam2= Alam
      Alam=  MAX(tmpLam,0.1_dp*Alam)
      !
      F2=F
      !
    ENDDO DoLineSearch
    !
    !---------------------------------------/linesearch and backtracking
    !
    Error_F=  MAXVAL( ABS(vFunc(:)) ) !-> error on FUNCTION value
    !
    !-- Press-> max variation on vX
    !Delta_X= MAXVAL( ABS(vX(:)-vXOld(:))/MAX(ABS(vX(:)),One) )
    !
    !-- Archim-> max relative variation on vX
    Delta_X=  MAXVAL( ABS( (vX(:)-vXOld(:)) /vX(:)) )
    !
    Gradient= MAXVAL( ABS(vGrad(:)) *MAX(ABS(vX(:)),One) /MAX(F,0.5_dp*SIZE(vX)) )
    !
    IF(ShoResult) PRINT '(I3,3(G15.6,A))',iTs,Error_F,"=F ",Gradient,"=G ",Delta_X,"=X"
    !
    !-------------------------------- Test for convergence on FUNCTION values --
    !! IF(Converge(vFunc)) THEN
    IF(Error_F < TolF) THEN
      iErr= 0
      CALL Newton_Sho(vX,vFunc,Its) !mod add 19/06/2008 08:22
      RETURN !----------------------------------------------------RETURN
    ENDIF
    !! IF(Check) THEN !Check for gradient of f being zero
    !!   Check= (Gradient < TolMin ) !convergence on dx
    !!   iErr= -4
    !!   RETURN !-------------------------------------------------RETURN
    !! ENDIF
    IF(DeltaX_Ok .AND. Delta_X < TolX) THEN
    !-> while Error_F still high, current vX very close to previous vX
    !-> arrived to a stationary point ??
      iErr= -5 !
      RETURN !----------------------------------------------------RETURN
    ENDIF
    !
    !----------------------------------------------end of iteration loop
  ENDDO
ENDSUBROUTINE NewtLnsrch

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE Newton_Press( & !from "NR"
& vX,       & !inout= the initial guess, and the root returned
& vIsPlus,  & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
& Residual, & !INTERFACE
& Jacobian, & !INTERFACE
& Converge, & !INTERFACE
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& bFinDIF,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   MAXVAL(ABS(vFunc(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)      !out=   error code
!
!iErr= 0 : convergence reached -> OK
!iErr=-1 : MaxIts reached without convergence
!iErr=-2 : singular jacobian
!iErr=-3 : roundoff problem in linesearch
!iErr=-4 : no convergence on vFunc, problem with gradient ?
!iErr=-5 : no convergence on vFunc, vX very close to previous vX -> stationary point ??
!
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Tools,ONLY: fNewtJ
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  LOGICAL, INTENT(IN)   :: vIsPlus(:)
  REAL(dp),INTENT(IN)   :: TolF,TolX !!,TolMin
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X !,Gradient
  !! LOGICAL, INTENT(OUT)  :: Check
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),              INTENT(IN) :: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(vF,vTolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: vF(:),vTolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  !INTEGER, PARAMETER :: MAXITS=NewtMaxIts
  ! NewtTolF=   convergence criterion on FUNCTION values
  ! NewtTolMin= criterion for spurious convergence
  ! NewtTolX=   convergence criterion on dx
  !
  ! scaled maximum step length allowed in line searches
  REAL(dp),PARAMETER:: STPMX= 100.0
  !
  REAL(dp),PARAMETER:: Alf= Epsilon(vX) !1.0D-12
  !
  LOGICAL :: Check
  LOGICAL :: bSingul
  LOGICAL :: ForcePositive
  INTEGER :: ITS
  REAL(dp):: D,F,rFOld,rStpMax,rSlope,Norm_vDX,Gradient !,FMin
  !
  REAL(dp),DIMENSION(SIZE(vX)):: vFunc,vDX,vX0,vGrad
  INTEGER, DIMENSION(SIZE(vX)):: vIndex
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJacob
  REAL(dp),DIMENSION(SIZE(vX)):: vTolF
  !
  REAL(dp):: A,Alam,Alam2,AlaMin,B,DISC,F2,X
  REAL(dp):: RHS1,RHS2,tmpLam
  !
  ForcePositive= COUNT(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  vFunc= Residual(vX)
  F= DOT_PRODUCT(vFunc,vFunc) *0.5_dp
  !
  iErr=-1
  !
  IF(ShoResult) PRINT '(/,A,/)',"iTs,Error_F,Gradient,Delta_X"
  !
  ! Calculate StpMax for line searches
  rStpMax= STPMX*MAX(vAbs(vX(:)),REAL(SIZE(vX),dp))
  !
  DO its=1,MaxIts
    !
    nIts=its
    !
    CALL Newton_Sho(vX,vFunc,Its)
    !
    IF(bFinDIF) THEN
      CALL Jacobian_Numeric(Residual,vX,vFunc,tJacob)
    ELSE
      CALL Jacobian(vX,tJacob)
    ENDIF
    !
    vGrad(:)= MATMUL(vFunc(:),tJacob(:,:)) ! compute grad(f) for the line search.
    vX0(:)= vX(:)                          ! store vX
    rFOld=    F                            ! store F
    vDX(:)=  -vFunc(:)                     ! right-hand side for linear equations.
    !
    !-------------------------solve linear equations by LU decomposition
    CALL LU_Decomp(tJacob,vIndex,D,bSingul)
    !
    IF(bSingul) THEN
      iErr= -2 !-> "SINGULAR JACOBIAN, cf Log File"
      !CALL Jacobian_Sho(fNewtJ,tJacob)
      RETURN !-------------------------------System not solved -> RETURN
    ENDIF
    !
    CALL LU_BakSub(tJacob,vIndex,vDX)
    !------------------------/solve linear equations by LU decomposition
    !
    !~ IF(fNewtJ>0) THEN
      !~ CALL Jacobian_Sho(fNewtJ,tJacob)
      !~ fNewtJ= 0  ;  CLOSE(fNewtJ)
    !~ ENDIF
    !
    Norm_vDX= vAbs(vDX(:))
    !-- Scale if attempted step is too big --
    IF (Norm_vDX>rStpMax) vDX(:)= rStpMax *vDX(:) /Norm_vDX
    !
    rSlope= DOT_PRODUCT(vGrad,vDX)
    !
    IF (rSlope>=Zero) THEN
      iErr= -3 !-> 'roundoff problem in linesearch, cf log file'
      RETURN !-------------------------------System not solved -> RETURN
    ENDIF
    !
    !-------------------------------------------------enforce positivity
    IF(ForcePositive) THEN
      !--- under relaxation technique from Bethke
      !--- if (Xi/2 + dXi)-0
      X= MAXVAL(-vDX(:)/vX0(:)*2.0D0,MASK=vIsPlus(:))
      if(iDebug>2 .AND. X>1.0D0) write(69,*) " Press-UR, ", X
      IF(X>1.0D0) vDX(:)= vDX(:) /X
      !DO
      !  X= MINVAL(vX(:),MASK=vIsPlus(:))
      !  IF(X>Zero) EXIT
      !  vDX(:)= vDX(:)*0.5D0
      !  vX(:)= vX0(:) + vDX(:)
      !ENDDO
    ENDIF
    !------------------------------------------------/enforce positivity
    !
    vX(:)= vX0(:) +vDX(:)
    !
    !----------------------------------------linesearch and backtracking
    !
    CHECK= .FALSE.
    AlaMin= TolX /MAXVAL(ABS(vDX(:))/MAX(ABS(vX0(:)),One))
    Alam=   One !always try full Newton step first
    !
    DoLineSearch: DO !Start of iteration loop
      !
      !!write(51,'(I3,G15.6)') Its, Alam
      !
      vX(:)= vX0(:) +Alam*vDX(:)
      !
      vFunc= Residual(vX)
      F= DOT_PRODUCT(vFunc,vFunc)/2.0D0
      !
      IF (Alam<AlaMin) THEN !Convergence on vX
        !
        vX(:)=vX0(:)
        !-- For zero finding, the calling program should verify the convergence
        CHECK=.TRUE.
        EXIT DoLineSearch !-----------------------------------------EXIT
        !
      ELSEIF (F <= rFOld + ALF*Alam*rSlope) THEN
        !
        EXIT DoLineSearch !sufficient FUNCTION decrease ------------EXIT
        !
      ELSE !--Backtrack.
        !
        IF (Alam==One) THEN !---------------------------------First time
          tmpLam= -rSlope /(2.0_dp*(F-rFOld-rSlope))
          !
        ELSE !-------------------------------------Subsequent backtracks
          RHS1= F -rFOld -Alam *rSlope
          RHS2= F2-rFOld -Alam2*rSlope
          A= (       RHS1/Alam**2 -      RHS2/Alam2**2) /(Alam-Alam2)
          B= (-Alam2*RHS1/Alam**2 + Alam*RHS2/Alam2**2) /(Alam-Alam2)
          !
          IF (A==Zero) THEN
            tmpLam= -rSlope/(2.0_dp*B)
            !
          ELSE
            B= MIN(B, 1.0D100)
            B= MAX(B,-1.0D100)
            DISC= B*B -3.0_dp*A*rSlope
            IF     (DISC<Zero) THEN  ; tmpLam= 0.5_dp*Alam
            ELSEIF (B<=Zero)   THEN  ; tmpLam= (-B+SQRT(DISC))/3.0_dp/A
            ELSE                     ; tmpLam= -rSlope/(B+SQRT(DISC))
            ENDIF
          ENDIF
          !
          tmpLam= MIN(tmpLam,0.5_dp*Alam)
          !
        ENDIF
        !
      ENDIF
      !
      Alam2= Alam
      Alam=  MAX(tmpLam,0.1_dp*Alam)
      !
      F2=F
      !
    ENDDO DoLineSearch
    !
    !---------------------------------------/linesearch and backtracking
    !
    Error_F=  MAXVAL( ABS(vFunc(:)) ) !-> error on FUNCTION value
    !
    !-- Press-> max relative variation on vX
    Delta_X= MAXVAL( ABS(vX(:)-vX0(:))/MAX(ABS(vX(:)),One) )
    !
    !-- Archim-> max relative variation on vX
    ! Delta_X=  MAXVAL( ABS( (vX(:)-vX0(:)) /vX(:)) )
    !
    Gradient= MAXVAL( ABS(vGrad(:)) *MAX(ABS(vX(:)),One) /MAX(F,0.5_dp*SIZE(vX)) )
    !
    IF(ShoResult) PRINT '(I3,3(G15.6,A))',iTs,Error_F,"=F ",Gradient,"=G ",Delta_X,"=X"
    !
    !----------------------------test for convergence on function values
    IF(Converge(vFunc,vTolF)) THEN
    !IF(Error_F < TolF) THEN
      iErr= 0
      CALL Newton_Sho(vX,vFunc,Its) !mod add 19/06/2008 08:22
      RETURN !----------------------------------------------------RETURN
    ENDIF
    !---/
    !! IF(Check) THEN !Check for gradient of f being zero
    !!   Check= (Gradient < TolMin ) !convergence on dx
    !!   iErr= -4
    !!   RETURN !-------------------------------------------------RETURN
    !! ENDIF
    !------------------------------------------test for stationary point
    IF(DeltaX_Ok .AND. Delta_X < TolX) THEN
    !-> while Error_F still high, current vX very close to previous vX
    !-> arrived to a stationary point ??
      iErr= -5 !
      RETURN !----------------------------------------------------RETURN
    ENDIF
    !---/
    !
    !----------------------------------------------end of iteration loop
  ENDDO
ENDSUBROUTINE Newton_Press

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE Newton_Walker( & !
!--
!-- version adapted for direct substitution
!-- and for solving with moles -> enforce vX(:)>0
!--
& vX,       & !inout=  the initial guess, and the root returned
& vIsPlus,  & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
& Residual, & !
& Jacobian, & !
& Converge, & !
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& bFinDIF,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!--
!-- translated from matlab code by Homer Walekr
!-- http://users.wpi.edu/~walker/MA590
!-- cf http://users.wpi.edu/~walker/MA590/HANDOUTS/newton-backtracking.pdf
!--
!-- Given an initial guess x for a root in n dimensions,
!-- take Nits Newton-Raphson steps to improve the root
!-- Stop if the root converges,
!-- in either summed absolute variable increments NewtTolX
!-- or summed absolute function values NewtTolF
!--
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Tools,ONLY: fNewtJ
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  !
  !USE M_Numeric_Const,  ONLY: Ln10
  !USE M_IOTools,ONLY: OutStrVec
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  LOGICAL, INTENT(IN)   :: vIsPlus(:)
  REAL(dp),INTENT(IN)   :: TolF
  REAL(dp),INTENT(IN)   :: TolX
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(vF,vTolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: vF(:),vTolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vFunc0,vDX,vX0
  REAL(dp),DIMENSION(SIZE(vX))         :: vTolF
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  !
  INTEGER :: Its
  REAL(dp):: D
  LOGICAL :: bSingul
  LOGICAL :: ForcePositive
  !
  REAL(dp),PARAMETER:: &
  & SigmaMax= 0.5_dp,  &
  & SigmaMin= 0.1_dp,  &
  & Tau=      1.0D-4
  INTEGER, PARAMETER:: iArmMax= 12
  INTEGER :: iArm
  REAL(dp):: Norm_vF0,Norm_vF,delta,lambda,Sigma
  REAL(dp):: X
  !
  ForcePositive= COUNT(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  iErr=-1
  !
  DO Its=1,MaxIts
    !
    vX0(:)= vX(:)
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(SUM(vFunc0(:)*vFunc0(:)))
    !
    CALL Newton_Sho(vX0,vFunc0,Its)
    !
    IF(bFinDIF) THEN
      CALL Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
    ELSE
      CALL Jacobian(vX0,tJac)
    ENDIF
    !
    !~ IF(fNewtJ>0) THEN
      !~ CALL Jacobian_Sho(fNewtJ,tJac)
      !~ fNewtJ= 0  ;  CLOSE(fNewtJ)
    !~ ENDIF
    !
    !----------------------solve linear equations using LU decomposition
    CALL LU_Decomp(tJac, vIndex, D, bSingul)
    IF(bSingul) THEN
      iErr=-2  ;  EXIT !--------------------------------------------EXIT
    ENDIF
    vDX(:)= -vFunc0(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    !---------------------/solve linear equations using LU decomposition
    !
    IF(ForcePositive) THEN
      X= MAXVAL(-vDX(:)/vX0(:)*2.0D0,MASK=vIsPlus(:))
      if(iDebug>2 .AND. X>1.0D0) write(69,*) "Walker-UR, ", X
      IF(X>1.0D0) vDX(:)= vDX(:) /X
      !DO
      !  X= MINVAL(vX(:),MASK=vIsPlus(:))
      !  IF(X>Zero) EXIT
      !  vDX(:)= vDX(:)*0.5D0
      !  vX(:)= vX0(:) + vDX(:)
      !ENDDO
    ENDIF
    !
    vX(:)= vX0(:) + vDX(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
    !
    lambda= 1.0d0
    iArm=   0
    !
    !---------------------------Test the step and backtrack as necessary
    DO WHILE ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      iArm= iArm +1
      !
      IF(iArm > iArmMax) THEN !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  EXIT
      ENDIF
      !
      delta= (Norm_vF /Norm_vF0)**2 -1.0d0 +2.0d0*lambda
      !
      IF (delta > Zero) THEN
        Sigma=  lambda/delta
        Sigma=  MIN(Sigma,SigmaMax)
        Sigma=  MAX(Sigma,SigmaMin)
      ELSE
        Sigma=  SigmaMax
      END IF
      !
      vDX(:)=  Sigma*vDX(:)
      lambda=  Sigma*lambda
      vX(:)=   vX0(:) + vDX(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
      !
    END DO
    !--------------------------/Test the step and backtrack as necessary
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vDX(:)))
    !
    !-- convergence --
    !IF (Error_F < TolF) THEN
    IF(Converge(vFunc,vTolF)) THEN  ;  iErr= 0  ;  EXIT !-----------EXIT
    ENDIF
    !
    !-- stationary --
    IF (Delta_X < TolX) THEN        ;  iErr=-5  ;  EXIT !-----------EXIT
    ENDIF
    !
  ENDDO
  !
  nIts=Its
  !
END SUBROUTINE Newton_Walker

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE Newton_Kelley( & !
& vX,        & !inout= the initial guess, and the root returned
& vIsPlus,   & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
& Residual,  & !
& Jacobian,  & !
& Converge,  & !
& TolF,      & !in=    convergence criterion on function values
& TolX,      & !in=    convergence criterion on dx
& bFinDIF,   & !in=    use numeric Jacobian
& MaxIts,    & !in=    maximum number of iterations
& Error_F,   & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,   & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,      & !out=   number of iterations
& iErr)        !out=   error code
!--
!-- translated from MATLAB code nsold.m,
!-- http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsold.m
!-- cf CT Kelley, Solving nonlinear equations with Newton's method, SIAM, 2003
!--
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Tools,ONLY: fNewtJ
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  LOGICAL, INTENT(IN)   :: vIsPlus(:)
  REAL(dp),INTENT(IN)   :: TolF,TolX
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(F,tolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: F(:),tolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vFunc0,vDX,vX0,vStep,vTolF
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  !
  INTEGER :: Its
  REAL(dp):: D,X
  LOGICAL :: bSingul
  LOGICAL :: ForcePositive
  !
  !-- Armijo parameters --
  REAL(dp),PARAMETER:: SigmaMax= 0.5_dp
  REAL(dp),PARAMETER:: SigmaMin= 0.1_dp
  REAL(dp),PARAMETER:: Tau=      1.0D-4
  INTEGER, PARAMETER:: iArmMax= 10
  !--
  !
  INTEGER :: iArm
  REAL(dp):: Norm_vF0,Norm_vF,Norm_Ratio
  REAL(dp):: lambda,lam_m,lam_c
  REAL(dp):: ff_0, ff_m, ff_c
  !
  ForcePositive= COUNT(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  iErr=-1
  Norm_Ratio= One
  !
  DO Its=1,MaxIts
    !
    vX0(:)= vX(:)
    !
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(DOT_PRODUCT(vFunc0(:),vFunc0(:)))
    !
    CALL Newton_Sho(vX0,vFunc0,Its)
    !
    !----------------------solve linear equations using LU decomposition
    !!IF(Its==1 .OR. Norm_Ratio > 0.5_dp) THEN
      !-- evaluate and factorize Jacobian
      !-- only when norm(Residual) strongly decreases
      !
      IF(bFinDIF) THEN
        CALL Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
      ELSE
        CALL Jacobian(vX0,tJac)
      ENDIF
      !
      CALL LU_Decomp(tJac, vIndex, D, bSingul)
      IF(bSingul) THEN
        iErr=-2  ;  EXIT
      ENDIF
      !
    !!ENDIF
    !
    vDX(:)= -vFunc0(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    !-------------------------------------------------------------/solve
    !
    !~ IF(fNewtJ>0) THEN
      !~ CALL Jacobian_Sho(fNewtJ,tJac)
      !~ fNewtJ= 0  ;  CLOSE(fNewtJ)
    !~ ENDIF
    !
    IF(ForcePositive) THEN
      X= MAXVAL(-vDX(:)/vX0(:)*2.0D0,MASK=vIsPlus(:))
      if(iDebug>2 .AND. X>1.0D0) write(69,*) "Kelley-UR, ", X
      IF(X>1.0D0) vDX(:)= vDX(:) /X
      !DO
      !  X= MINVAL(vX(:),MASK=vIsPlus(:))
      !  IF(X>Zero) EXIT
      !  vStep(:)= vStep(:)*0.5D0
      !  vX(:)= vX0(:) + vStep(:)
      !ENDDO
    ENDIF
    !
    lambda= 1.0d0
    lam_m=  lambda
    lam_c=  lambda
    iArm=   0
    !
    vStep(:)= lambda *vDX(:)
    vX(:)= vX0(:) + vStep(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(DOT_PRODUCT(vFunc(:),vFunc(:)))
    !
    ff_0= Norm_vF0**2  ! initial value
    ff_c= Norm_vF**2   ! current value
    ff_m= ff_c         ! previous value in linesearch
    !
    !---------------------------test the step and backtrack as necessary
    DO WHILE ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      IF(iArm==0) THEN
        lambda=  SigmaMax*lambda
      ELSE
        CALL Parab3p(SigmaMin,SigmaMax, lam_c,lam_m, ff_0,ff_c,ff_m, lambda)
      ENDIF
      !
      iArm= iArm + 1
      IF(iArm > iArmmax) THEN !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  EXIT
      ENDIF
      !
      lam_m= lam_c    ! lam_c=  current steplength
      lam_c= lambda   ! lam_m=  previous steplength
      !
      vStep(:)= lambda *vDX(:)
      vX(:)= vX0(:) + vStep(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(DOT_PRODUCT(vFunc(:),vFunc(:)))
      !
      ff_m= ff_c
      ff_c= Norm_vF**2
      !
    END DO
    !--------------------------/test the step and backtrack as necessary
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vStep(:)))
    !
    Norm_Ratio= Norm_vF /Norm_vF0
    !! write(71,'(A,G16.6)') "Kelley, Norm_Ratio ", Norm_Ratio
    !
    !IF (Error_F < TolF) THEN
    IF(Converge(vFunc,vTolF)) THEN  ;  iErr= 0   ;  EXIT
    ENDIF
    !
    IF (Delta_X<TolX) THEN          ;  iErr=-5   ;  EXIT
    ENDIF
    !
  ENDDO
  !
  !! write(71,*)
  nIts=Its
  !
END SUBROUTINE Newton_Kelley

SUBROUTINE parab3p(sigma0,sigma1, lambdac,lambdam, ff0,ffc,ffm, lambdap)
! Apply three-point safeguarded parabolic model for a line search.
! C. T. Kelley, April 1, 2003
! input:
!       sigma0,sigma1= parameters, safeguarding bounds for the linesearch
!       lambdac=  current steplength
!       lambdam=  previous steplength
!       ff0=  |F(x_c)|^2
!       ffc=  |F(x_c + lambdac)|^2
!       ffm=  |F(x_c + lambdam)|^2
! output:
!       lambdap=  new value of lambda given by parabolic model
  REAL(dp),INTENT(IN) :: sigma0,sigma1
  REAL(dp),INTENT(IN) :: lambdac, lambdam, ff0, ffc, ffm
  REAL(dp),INTENT(OUT):: lambdap
  !
  REAL(dp):: c1, c2
  !
  !-------------------- Compute coefficients of interpolation polynomial
  ! p(lambda)=  ff0 + (-c1 lambda + c2 lambda^2)/d1
  ! d1=  (lambdac - lambdam)*lambdac*lambdam < 0
  ! so, if c2 > 0 we have negative curvature and default to
  ! lambdap=  sigma1 * lambda.
  !
  c2= lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
  !
  IF(c2 >=0.0d0) THEN
    lambdap=  sigma1*lambdac
  ELSE
    c1= lambdam *lambdam *(ffc-ff0) - lambdac *lambdac *(ffm-ff0)
    lambdap= c1 *0.5d0 /c2
    lambdap= MAX(lambdap, sigma0*lambdac)
    lambdap= MIN(lambdap, sigma1*lambdac)
  END IF
  !
  RETURN
END SUBROUTINE parab3p

SUBROUTINE Newton_Walker_old( & !
& vX,       & !inout=  the initial guess, and the root returned
& Residual, & !
& Jacobian, & !
& Converge, & !
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& bFinDIF,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!--
!-- translated from matlab code by Homer Walekr
!-- http://users.wpi.edu/~walker/MA590
!-- cf http://users.wpi.edu/~walker/MA590/HANDOUTS/newton-backtracking.pdf
!--
!-- Given an initial guess x for a root in n dimensions,
!-- take Nits Newton-Raphson steps to improve the root
!-- Stop if the root converges,
!-- in either summed absolute variable increments NewtTolX
!-- or summed absolute function values NewtTolF
!--
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  !
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_IOTools,ONLY: OutStrVec
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  REAL(dp),INTENT(IN)   :: TolF,TolX
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(vF,vTolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: vF(:),vTolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vFunc0,vDX,vX0,vTolF
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  !
  INTEGER :: Its
  REAL(dp):: D
  LOGICAL :: bSingul
  !
  REAL(dp),PARAMETER:: &
  & SigmaMax= 0.5_dp,  &
  & SigmaMin= 0.1_dp,  &
  & Tau=      1.0D-4
  INTEGER, PARAMETER:: iArmMax= 12
  INTEGER :: iArm
  REAL(dp):: Norm_vF0,Norm_vF,delta,lambda,Sigma
  !
  vTolF(:)= TolF
  !
  iErr=-1
  !
  DO Its=1,MaxIts
    !
    vX0(:)= vX(:)
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(SUM(vFunc0(:)*vFunc0(:)))
    !
    CALL Newton_Sho(vX0,vFunc0,Its)
    !
    IF(bFinDIF) THEN
      CALL Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
    ELSE
      CALL Jacobian(vX0,tJac)
    ENDIF
    !
    !---solve linear equations using LU decomposition --
    CALL LU_Decomp(tJac, vIndex, D, bSingul)
    IF(bSingul) THEN
      iErr=-2  ;  EXIT
    ENDIF
    vDX(:)= -vFunc0(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    !---/solve
    !
    vX(:)= vX0(:) + vDX(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
    !
    lambda= 1.0d0
    iArm=   0
    !
    !---------------------------Test the step and backtrack as necessary
    DO WHILE ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      iArm= iArm +1
      !
      IF(iArm > iArmMax) THEN !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  EXIT
      ENDIF
      !
      delta= (Norm_vF /Norm_vF0)**2 -1.0d0 +2.0d0*lambda
      !
      IF (delta > Zero) THEN
        Sigma=  lambda/delta
        Sigma=  MIN(Sigma,SigmaMax)
        Sigma=  MAX(Sigma,SigmaMin)
      ELSE
        Sigma=  SigmaMax
      END IF
      !
      vDX(:)=  Sigma*vDX(:)
      lambda=  Sigma*lambda
      vX(:)=   vX0(:) + vDX(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(SUM(vFunc(:)*vFunc(:)))
      !
    END DO
    !--------------------------/Test the step and backtrack as necessary
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vDX(:)))
    !
    !IF (Error_F < TolF) THEN
    IF(Converge(vFunc,vTolF)) THEN
      iErr= 0  ;  EXIT
    ENDIF
    IF (Delta_X < TolX) THEN !== stationary ==
      iErr=-5  ;  EXIT
    ENDIF
    !
  ENDDO
  !
  nIts=Its
  !
END SUBROUTINE Newton_Walker_old

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE Newton_Kelley_old( & !
& vX,        & !inout= the initial guess, and the root returned
& Residual,  & !
& Jacobian,  & !
& Converge,  & !
& TolF,      & !in=    convergence criterion on function values
& TolX,      & !in=    convergence criterion on dx
& bFinDif,   & !in=    use numeric Jacobian
& MaxIts,    & !in=    maximum number of iterations
& Error_F,   & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,   & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,      & !out=   number of iterations
& iErr)        !out=   error code
!--
!-- translated from MATLAB code nsold.m,
!-- http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsold.m
!-- cf CT Kelley, Solving nonlinear equations with Newton's method, SIAM, 2003
!--
  USE M_Numeric_Tools,ONLY: vAbs,Jacobian_Numeric
  USE M_Numeric_Mat,  ONLY: LU_BakSub,LU_Decomp
  !
  USE M_Numeric_Const,ONLY: Ln10
  USE M_IOTools,      ONLY: OutStrVec
  !
  REAL(dp),INTENT(INOUT):: vX(:)
  REAL(dp),INTENT(IN)   :: TolF,TolX
  LOGICAL, INTENT(IN)   :: bFinDIF
  INTEGER, INTENT(IN)   :: MaxIts
  REAL(dp),INTENT(OUT)  :: Error_F,Delta_X
  INTEGER, INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
    LOGICAL FUNCTION Converge(F,tolF)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),INTENT(IN):: F(:),tolF(:)
    ENDFUNCTION Converge
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vFunc0,vDX,vX0,vStep,vTolF
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  !
  INTEGER :: Its
  REAL(dp):: D
  LOGICAL :: bSingul
  !
  !--- Armijo parameters --
  REAL(dp),PARAMETER:: SigmaMax= 0.5_dp
  REAL(dp),PARAMETER:: SigmaMin= 0.1_dp
  REAL(dp),PARAMETER:: Tau=      1.0D-4
  INTEGER, PARAMETER:: iArmMax= 10
  !---/
  !
  INTEGER :: iArm
  REAL(dp):: Norm_vF0,Norm_vF,Norm_Ratio
  REAL(dp):: lambda,lam_m,lam_c
  REAL(dp):: ff_0, ff_m, ff_c
  !
  vTolF(:)= TolF
  !
  iErr=-1
  Norm_Ratio= One
  !
  DO Its=1,MaxIts
    !
    vX0(:)= vX(:)
    !
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(DOT_PRODUCT(vFunc0(:),vFunc0(:)))
    !
    CALL Newton_Sho(vX0,vFunc0,Its)
    !
    !IF(DebNewt .AND. bOpenNewt) THEN
    !  CALL OutStrVec(fNewt1,X(1:SIZE(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  CALL OutStrVec(fNewt2,vFunc0(1:SIZE(X)),  Opt_I=Its,Opt_C="F")
    !ENDIF
    !
    !-------------------------- solve linear equations using LU decomposition --
    !!IF(Its==1 .OR. Norm_Ratio > 0.5_dp) THEN
      !-- evaluate and factorize Jacobian
      !-- only when norm(Residual) strongly decreases
      !
      IF(bFinDIF) THEN
        CALL Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
      ELSE
        CALL Jacobian(vX0,tJac)
      ENDIF
      !
      CALL LU_Decomp(tJac, vIndex, D, bSingul)
      IF(bSingul) THEN
        iErr=-2  ;  EXIT
      ENDIF
      !
    !!ENDIF
    !
    vDX(:)= -vFunc0(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    !-------------------------------------------------------------/solve
    !
    lambda= 1.0d0
    lam_m=  lambda
    lam_c=  lambda
    iArm=   0
    !
    vStep(:)= lambda *vDX(:)
    vX(:)= vX0(:) + vStep(:)
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(DOT_PRODUCT(vFunc(:),vFunc(:)))
    !
    ff_0= Norm_vF0**2
    ff_c= Norm_vF**2
    ff_m= Norm_vF**2
    !
    !---------------------------test the step and backtrack as necessary
    DO WHILE ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      IF(iArm==0) THEN
        lambda=  SigmaMax*lambda
      ELSE
        CALL Parab3p(SigmaMin,SigmaMax, lam_c,lam_m, ff_0,ff_c,ff_m, lambda)
      ENDIF
      !
      iArm=  iArm + 1
      IF(iArm > iArmmax) THEN !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  EXIT
      ENDIF
      !
      ! lam_c=  current steplength
      ! lam_m=  previous steplength
      lam_m= lam_c
      lam_c= lambda
      !
      vStep(:)= lambda *vDX(:)
      vX(:)= vX0(:) + vStep(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(DOT_PRODUCT(vFunc(:),vFunc(:)))
      !
      ff_m= ff_c
      ff_c= Norm_vF**2
      !
    END DO
    !------------------------------/ test the step and backtrack as necessary --
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vStep(:)))
    !
    Norm_Ratio= Norm_vF /Norm_vF0
    !! write(71,'(G16.6)') Norm_Ratio
    !
    !IF (Error_F < TolF) THEN
    IF(Converge(vFunc,vTolF)) THEN
      iErr= 0; EXIT
    ENDIF
    IF (Delta_X<TolX) THEN
      iErr=-5; EXIT
    ENDIF
    !
  ENDDO
  !
  nIts=Its
  !
END SUBROUTINE Newton_Kelley_old

SUBROUTINE parab3p_old(sigma0,sigma1, lambdac,lambdam, ff0,ffc,ffm, lambdap)
! Apply three-point safeguarded parabolic model for a line search.
! C. T. Kelley, April 1, 2003
! input:
!       sigma0,sigma1= parameters, safeguarding bounds for the linesearch
!       lambdac=  current steplength
!       lambdam=  previous steplength
!       ff0=  |F(x_c)|^2
!       ffc=  |F(x_c + lambdac)|^2
!       ffm=  |F(x_c + lambdam)|^2
! output:
!       lambdap=  new value of lambda given by parabolic model
  REAL(dp),INTENT(IN) :: sigma0,sigma1
  REAL(dp),INTENT(IN) :: lambdac, lambdam, ff0, ffc, ffm
  REAL(dp),INTENT(OUT):: lambdap
  !
  REAL(dp):: c1, c2
  !
  !---------------------Compute coefficients of interpolation polynomial
  ! p(lambda)=  ff0 + (-c1 lambda + c2 lambda^2)/d1
  ! d1=  (lambdac - lambdam)*lambdac*lambdam < 0
  ! so, if c2 > 0 we have negative curvature and default to
  ! lambdap=  sigma1 * lambda.
  !
  c2= lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
  !
  IF(c2 >=0.0d0) THEN
    lambdap=  sigma1*lambdac
  ELSE
    c1= lambdam *lambdam *(ffc-ff0) - lambdac *lambdac *(ffm-ff0)
    lambdap= c1 *0.5d0 /c2
    lambdap= MAX(lambdap, sigma0*lambdac)
    lambdap= MIN(lambdap, sigma1*lambdac)
  END IF
  !
  RETURN
END SUBROUTINE parab3p_old

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE NewtonChess( &
& vX, &      !initial guess x for a root in N DIMENSIONs
& Residual, Jacobian, &
& TolF,    & !in=    convergence criterion on FUNCTION values
& TolX,    & !in=    convergence criterion on dx
!!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
& MaxIts,  & !in=    maximum number of iterations
& Error_F, & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
!& Gradient, & !out=
& Nits,    & !out=   number of iterations
!& Check,   & !out=   IF Check, should check convergence
& iErr)      !out=   error code
!Given an initial guess x for a root in N DIMENSIONs,
!take Nits Newton-Raphson steps to improve the root
!Stop IF the root converges,
!in either summed absolute variable increments NewtTolX
!or summed absolute FUNCTION values NewtTolF.
  USE M_Numeric_Tools,   ONLY: vAbs, iMaxLoc_R
  USE M_Numeric_Mat,ONLY: LU_BakSub,LU_Decomp
  !
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_IOTools,ONLY: OutStrVec
  !
  REAL(dp),DIMENSION(:), INTENT(INOUT):: vX
  REAL(dp),              INTENT(IN)   :: TolF,TolX !!,TolMin
  INTEGER,               INTENT(IN)   :: MaxIts
  REAL(dp),              INTENT(OUT)  :: Error_F,Delta_X !,Gradient
  !LOGICAL,               INTENT(OUT)  :: Check
  INTEGER,               INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vDX
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  INTEGER :: Its
  REAL(dp):: D
  LOGICAL ::bSingul
  !
  !Chess/Concepts, p71, Newton with "polishing factor"
  INTEGER :: I
  REAL(dp):: R,Alfa
  REAL(dp),PARAMETER::  A=0.5D0, B=3.0D0, C=-0.9D0
  !iErr=-1: MaxIts reached
  !iErr=-2: singular Jacobian
  !iErr=-5:
  !
  iErr=-1
  !
  DO Its=1,MaxIts
    !
    vFunc(:)= Residual(vX)
    CALL Jacobian(vX,tJac)
    !
    !IF(DebNewt .AND. bOpenNewt) THEN
    !  CALL OutStrVec(fNewt1,X(1:SIZE(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  CALL OutStrVec(fNewt2,vFunc(1:SIZE(X)),  Opt_I=Its,Opt_C="F")
    !ENDIF
    !
    !--- Solve linear equations using LU decomposition --
    CALL LU_Decomp(tJac,vIndex,D,bSingul)
    IF(bSingul) THEN
      iErr=-2; EXIT
    ENDIF
    vDX(:)= -vFunc(:)
    CALL LU_BakSub(tJac,vIndex,vDX) !solve vDX= -vFunc.inv(tJac)
    !---/ Solve
    !
    R= One
    !Chess/Concepts, p71, Newton with "polishing factor"
    Alfa= MAXVAL(ABS(vDX(:)/vX(:)))  !alfa_i= |dX_i|/X_i
    IF(Alfa > A) THEN !-> far from root
      I= iMaxLoc_R(ABS(vDX(:)/vX(:)))
      IF(vDX(I)>Zero) R=   (Alfa*B - A*A)/(B + Alfa -2.0D0*A)*vX(I)/vDX(I)
      IF(vDX(I)<Zero) R= C*(Alfa   - A*A)/(B + Alfa -2.0D0*A)*vX(I)/vDX(I)
    ENDIF
    !
    vX= vX + R *vDX
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vDX(:)))
    !
    IF (Error_F<TolF) THEN
      iErr= 0; EXIT
    ENDIF
    !
    IF (Delta_X<TolX) THEN
      iErr=-5; EXIT
    ENDIF
  ENDDO
  !
  nIts=Its
ENDSUBROUTINE NewtonChess

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE Newton( &
& vX, &      ! inout= the initial guess, and the root returned
& Residual, Jacobian, &
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!-----------------------------------------------------------------------
! Given an initial guess x for a root in n dimensions,
! take Nits Newton-Raphson steps to improve the root
! Stop if the root converges,
! in either summed absolute variable increments NewtTolX
! or summed absolute function values NewtTolF
!-----------------------------------------------------------------------
  USE M_Numeric_Tools,   ONLY: vAbs
  USE M_Numeric_Mat,ONLY: LU_BakSub,LU_Decomp
  !
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_IOTools,ONLY: OutStrVec
  !
  REAL(dp),DIMENSION(:), INTENT(INOUT):: vX
  REAL(dp),              INTENT(IN)   :: TolF,TolX !!,TolMin
  INTEGER,               INTENT(IN)   :: MaxIts
  REAL(dp),              INTENT(OUT)  :: Error_F,Delta_X !,Gradient
  !LOGICAL,               INTENT(OUT)  :: Check
  INTEGER,               INTENT(OUT)  :: nIts,iErr
  INTERFACE
    FUNCTION Residual(v)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:),INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v))     :: Residual
    ENDFUNCTION Residual
    SUBROUTINE Jacobian(v,t)
      USE M_Kinds
      IMPLICIT NONE
      REAL(dp),DIMENSION(:), INTENT(IN):: v
      REAL(dp),DIMENSION(SIZE(v),SIZE(v)),INTENT(OUT):: t
    ENDSUBROUTINE Jacobian
  ENDINTERFACE
  !
  REAL(dp),DIMENSION(SIZE(vX))         :: vFunc,vDX
  REAL(dp),DIMENSION(SIZE(vX),SIZE(vX)):: tJac
  INTEGER, DIMENSION(SIZE(vX))         :: vIndex
  INTEGER :: Its
  REAL(dp):: D
  LOGICAL :: bSingul
  REAL(dp),PARAMETER:: Alfa= 0.5_dp
  !
  iErr=-1
  DO Its=1,MaxIts
    !
    vFunc(:)= Residual(vX)
    CALL Jacobian(vX,tJac)
    !
    !IF(DebNewt .AND. bOpenNewt) THEN
    !  CALL OutStrVec(fNewt1,X(1:SIZE(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  CALL OutStrVec(fNewt2,vFunc(1:SIZE(X)),  Opt_I=Its,Opt_C="F")
    !ENDIF
    !
    !--- solve linear equations using LU decomposition --
    CALL LU_Decomp(tJac,vIndex,D,bSingul)
    IF(bSingul) THEN
      iErr=-2  ;  EXIT
    ENDIF
    vDX(:)= -vFunc(:)
    CALL LU_BakSub(tJac,vIndex,vDX)
    !---/ solve
    !
    vX= vX + vDX *Alfa
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Delta_X= MAXVAL(ABS(vDX(:)))
    !
    IF (Error_F < TolF) THEN
      iErr= 0; EXIT
    ENDIF
    IF (Delta_X<TolX) THEN
      iErr=-5; EXIT
    ENDIF
  ENDDO
  !
  nIts=Its
  !
  RETURN
END SUBROUTINE Newton

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

SUBROUTINE TestJacob_Sho(tJacob,tJacobNum)
  USE M_IoTools,ONLY: GetUnit, OutStrVec
  !
  REAL(dp),DIMENSION(:,:),INTENT(IN):: tJacob,tJacobNum
  !
  REAL(dp),DIMENSION(SIZE(tJacob,1))  :: xx,xr
  INTEGER:: ff,I
  !
  WRITE(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN"
  WRITE(fTrc,'(A)') "results in zzDynam_TestJacob.log"
  CALL GetUnit(ff)
  OPEN(ff,FILE="debug_dynam_jacob.log")
  DO I=1,SIZE(tJacob,1)
    CALL OutStrVec(ff,tJacob(I,:),   Opt_I=I,Opt_S="JacobAna",Opt_C="G")
    CALL OutStrVec(ff,tJacobNum(I,:),Opt_I=I,Opt_S="JacobDIF",Opt_C="G")
    !
    xx(:)= tJacob(I,:) -tJacobNum(I,:)
    CALL OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltJac",   Opt_C="G")
    !
    xr(:)= tJacob(I,:) +tJacobNum(I,:)
    WHERE(xr(:)>1.D-16) ; xx(:)= xx(:) / xr(:)
    ELSEWHERE           ; xx(:)= Zero
    ENDWHERE
    CALL OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltRel",   Opt_C="G")
  ENDDO
  WRITE(ff,*)
  CLOSE(ff)
  WRITE(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN END"
ENDSUBROUTINE TestJacob_Sho

SUBROUTINE Newton_Sho(vX,vFunc,Its) !,nAq)
  USE M_Numeric_Const,        ONLY: Ln10
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR,fNewt_I
  USE M_IoTools,      ONLY: OutStrVec
  !
  REAL(dp),INTENT(IN):: vX(:),vFunc(:)
  INTEGER, INTENT(IN):: Its
  !
  fNewt_I= fNewt_I +1
  IF(fNewtF>0) CALL OutStrVec(fNewtF,vX(:)/Ln10,   Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  IF(fNewtR>0) CALL OutStrVec(fNewtR,ABS(vFunc(:)),Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  !
ENDSUBROUTINE Newton_Sho

SUBROUTINE Jacobian_Sho(t)
  REAL(dp),INTENT(IN):: t(:,:)
  INTEGER:: I,J
  DO I=1,SIZE(t,1)
    DO J=1,SIZE(t,2)
      IF(t(I,J)/=0) THEN  ;  WRITE(fTrc,'(F7.1,A1)',ADVANCE='NO') t(I,J), T_
      ELSE                ;  WRITE(fTrc,'(A1,  A1)',ADVANCE='NO') "0",    T_
      ENDIF
    ENDDO
    WRITE(fTrc,*)
  ENDDO
  RETURN
ENDSUBROUTINE Jacobian_Sho

ENDMODULE M_Numeric_Newton
