module M_CMM_SPLINE
  use M_Kinds
  implicit none
  private
  !
  interface CMM_SPLINE_COMPUTE
    module procedure SPLINE
  end interface
  !
  interface CMM_SPLINE_EVAL
    module procedure SEVAL
  end interface
  !
  public :: CMM_SPLINE_COMPUTE
  public :: CMM_SPLINE_EVAL

contains

  subroutine SPLINE(N,X,Y,B,C,D)
    !-------------------------------------------------------------------
    !THIS subroutine CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
    !SPLINE TO BEST APPROXIMATE A DISCREET FONCTION GIVEN BY N POINTS
    !
    !   INPUTS:
    !     N       NUMBER OF GIVEN POINTS
    !     X,Y     VECTORS OF dimension N, STORING THE COORDINATES
    !             OF function F(X)
    !
    !   OUTPUTS:
    !     A,B,C   VECTORS OF dimension N, STORING THE COEFFICIENTS
    !             OF THE CUBIC SPLINE
    !
    !  REFERENCE:
    !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
    !     COMPUTATIONS. PRENTICE-HALL,INC.
    !-------------------------------------------------------------------
    integer, intent(in) :: N
    real(dp),intent(in) :: X(N),Y(N)
    real(dp),intent(out):: B(N),C(N),D(N)
    !----
    integer :: NM1, I, L
    real(dp):: T
    !----
    NM1= N-1
    if (N<2) return
    !
    if (N<3) then
      !--CAS N= 2
      B(1)= (Y(2)-Y(1))/(X(2)-X(1))
      C(1)= 0.D0
      D(1)= 0.D0
      B(2)= B(1)
      C(2)= 0.D0
      D(2)= 0.D0
      return
    else
      !-----------------------------------BUILD THE TRIDIAGONAL SYSTEM--
      !------------B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)--
      D(1)= X(2)-X(1)
      C(2)= (Y(2)-Y(1))/D(1)
      do I= 2,NM1
        D(I)= X(I+1)-X(I)
        B(I)= 2.D0*(D(I-1)+D(I))
        C(I+1)= (Y(I+1)-Y(I))/D(I)
        C(I)= C(I+1)-C(I)
      end do
      !-----------------------------------------CONDITIONS AT LIMITS----
      !------------THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES----
      B(1)= -D(1)
      B(N)= -D(N-1)
      C(1)= 0.D0
      C(N)= 0.D0
      if (N==3) then
        C(1)= C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
        C(N)= C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
        C(1)= C(1)*D(1)*D(1)/(X(4)-X(1))
        C(N)= -C(N)*D(N-1)**2/(X(N)-X(N-3))
      end if
      !------------------------------------------FORWARD ELIMINATION----
      do I= 2,N
        T= D(I-1)/B(I-1)
        B(I)= B(I)-T*D(I-1)
        C(I)= C(I)-T*C(I-1)
      end do
      !--------------------------------------------BACK SUBSTITUTION----
      C(N)= C(N)/B(N)
      do L= 1,NM1
        I= N-L
        C(I)= (C(I)-D(I)*C(I+1))/B(I)
      end do
      !     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL
      B(N)= (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      do I= 1,NM1
        B(I)= (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
        D(I)= (C(I+1)-C(I))/D(I)
        C(I)= 3.D0*C(I)
      end do
      !
      C(N)= 3.D0*C(N)
      D(N)= D(NM1)
      !
      return
    end if
    !
  end subroutine SPLINE

  subroutine SEVAL(N,X,Y,B,C,D,U,S)
  !---------------------------------------------------------------------
  !  EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCREET function F(X),
  !GIVEN IN N POINTS X(I), Y(I).
  !  THE B, C AND D COEFFICIENTS DEFINING THE BEST CUBIC SPLINE 
  !FOR THE GIVEN POINTS, ARE CALCULATED BEFORE BY THE SPLINE subroutine.
  !
  !   INPUTS:
  !     N       NUMBER OF POINTS OF CURVE Y= F(X)
  !     U       ABSCISSA OF POINT TO BE INTERPOLATED
  !     X,Y     TABLES OF dimension N, STORING THE COORDINATES
  !             OF CURVE F(X)
  !     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
  !             CUBIC SPLINE
  !
  !   OUTPUTS:
  !     S   INTERPOLATED value
  !            = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
  !             WITH DX= U-X(I), U BETWEEN X(I) AND X(I+1)
  !
  !   REFERENCE :
  !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
  !     COMPUTATIONS. PRENTICE-HALL,INC.
  !---------------------------------------------------------------------
    integer :: N
    real(dp):: X(N),Y(N)
    real(dp):: B(N),C(N),D(N)
    real(dp):: U
    real(dp):: S

    real(dp):: DX
    integer :: J, K
    integer, save:: I= 1

    if(N<2) then
      S= Y(1)
      return
    end if

    !--BINARY SEARCH
    if (I>=N) I= 1
    if (U<X(I)) GO TO 10
    if (U<=X(I+1)) GO TO 30
10  I= 1
    J= N+1
20  K= (I+J)/2
    if (U<X(K)) J= K
    if (U>=X(K)) I= K
    if (J>I+1) GO TO 20

    !--SPLINE EVALUATION
30  DX= U-X(I)
    S= Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
    return
  end subroutine SEVAL

end module M_CMM_SPLINE
