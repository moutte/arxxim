MODULE M_CMM_SPLINE

  IMPLICIT NONE
  PRIVATE

  INTERFACE CMM_SPLINE_COMPUTE
    MODULE PROCEDURE SPLINE
  END INTERFACE

  INTERFACE CMM_SPLINE_EVAL
    MODULE PROCEDURE SEVAL
  END INTERFACE

  PUBLIC :: CMM_SPLINE_COMPUTE
  PUBLIC :: CMM_SPLINE_EVAL

CONTAINS

  SUBROUTINE SPLINE (N,X,Y,B,C,D)
    !-------------------------------------------------------------------
    !     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
    !     SPLINE TO BEST APPROXIMATE A DISCREET FONCTION GIVEN BY N POINTS
    !
    !     INPUTS:
    !     N       NUMBER OF GIVEN POINTS
    !     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
    !             OF FUNCTION F(X)
    !
    !     OUTPUTS:
    !     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
    !             OF THE CUBIC SPLINE
    !
    !     REFERENCE:
    !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
    !     COMPUTATIONS. PRENTICE-HALL,INC.
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(kind=8), INTENT(IN)  :: X(N),Y(N)
    REAL(kind=8), INTENT(OUT) :: B(N),C(N),D(N)
    !----
    INTEGER :: NM1, I, L
    REAL(kind=8) :: T
    !----
    NM1 = N-1
    IF (N.LT.2) RETURN

    IF (N.LT.3) THEN

      !     CAS N = 2

      B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.D0
      D(1) = 0.D0
      B(2) = B(1)
      C(2) = 0.D0
      D(2) = 0.D0
      RETURN

    ELSE

      !     BUILD THE TRIDIAGONAL SYSTEM
      !     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO I = 2,NM1
        D(I) = X(I+1)-X(I)
        B(I) = 2.D0*(D(I-1)+D(I))
        C(I+1) = (Y(I+1)-Y(I))/D(I)
        C(I) = C(I+1)-C(I)
      END DO

      !     CONDITIONS AT LIMITS
      !     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0

      IF (N.EQ.3) THEN
        C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
        C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
        C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
        C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
      END IF

      !     FORWARD ELIMINATION

      DO I = 2,N
        T = D(I-1)/B(I-1)
        B(I) = B(I)-T*D(I-1)
        C(I) = C(I)-T*C(I-1)
      END DO

      !     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO L = 1,NM1
        I = N-L
        C(I) = (C(I)-D(I)*C(I+1))/B(I)
      END DO

      !     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      DO I = 1,NM1
        B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
        D(I) = (C(I+1)-C(I))/D(I)
        C(I) = 3.D0*C(I)
      END DO

      C(N) = 3.D0*C(N)
      D(N) = D(NM1)

      RETURN

    END IF

  END SUBROUTINE SPLINE

  SUBROUTINE SEVAL (N, X,Y,B,C,D,U,S)
    !------------------------------------------------------------------------
    !     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCREET FUNCTION F(X),
    !     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
    !     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
    !     BY THE SPLINE SUBROUTINE.
    !
    !     INPUTS:
    !     N       NUMBER OF POINTS OF CURVE Y = F(X)
    !     U       ABSCISSA OF POINT TO BE INTERPOLATED
    !     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
    !             OF CURVE F(X)
    !     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
    !             CUBIC SPLINE
    !
    !     OUTPUTS:
    !     S   INTERPOLATED VALUE
    !             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
    !             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
    !
    !     REFERENCE :
    !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
    !     COMPUTATIONS. PRENTICE-HALL,INC.
    !------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: N
    REAL (KIND=8) :: X(N),Y(N)
    REAL (KIND=8) :: B(N),C(N),D(N)
    REAL (KIND=8) :: U
    REAL (KIND=8) :: S

    REAL (KIND=8) :: DX
    INTEGER :: J, K
    INTEGER, SAVE :: I = 1

    IF(N<2) THEN
      S= Y(1)
      RETURN
    ENDIF

    !     BINARY SEARCH
    IF (I.GE.N) I = 1
    IF (U.LT.X(I)) GO TO 10
    IF (U.LE.X(I+1)) GO TO 30
10  I = 1
    J = N+1
20  K = (I+J)/2
    IF (U.LT.X(K)) J = K
    IF (U.GE.X(K)) I = K
    IF (J.GT.I+1) GO TO 20

    !     SPLINE EVALUATION

30  DX = U-X(I)
    S = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
    RETURN
  END SUBROUTINE SEVAL

END MODULE M_CMM_SPLINE
