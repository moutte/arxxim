MODULE M_Eos_Utils
!! M_ThermoUtils

  !----------------------------------------------------------
  ! Outils for Thermodynamics Computation
  !----------------------------------------------------------
  !
  ! Power functions : 
  !  - f_Square
  !  - f_Cube
  !  - f_Quart
  !  - f_Quint
  !
  ! Solvers for cubic equations : 
  !  - CUBIC
  !  - KUBIK
  !
  !----------------------------------------------------------
  ! Last Modification :  A.Michel 
  !----------------------------------------------------------

  USE M_Kinds
  IMPLICIT NONE
  PRIVATE
  REAL(dp), PARAMETER :: Pi= 3.14159265359979_dp 
 
  PUBLIC :: f_Square
  PUBLIC :: f_Cube
  PUBLIC :: f_Quart
  PUBLIC :: f_Quint

  PUBLIC :: CUBIC
  PUBLIC :: KUBIK

CONTAINS

  !----------------------------------------------------

  REAL(dp) FUNCTION f_Square(x)
    REAL(dp):: x
    f_Square=  x*x
    RETURN
  END FUNCTION f_Square

  REAL(dp) FUNCTION f_Cube(x)
    REAL(dp):: x
    f_Cube=  x*x*x
    RETURN
  ENDFUNCTION f_Cube

  REAL(dp) FUNCTION f_Quart(x)
    REAL(dp):: x
    f_Quart=  x*x*x*x
    RETURN
  END FUNCTION f_Quart

  REAL(dp) FUNCTION f_Quint(x)
    REAL(dp):: x
    f_Quint=  x*x*x*x*x
    RETURN
  END FUNCTION f_Quint

!----------------------------------------------------

SUBROUTINE CUBIC(B,C,D,X1,X2,X2I,X3) 
!----------------------------------------
! Compute Solution of a Cubic Equation 
! Equation  = X**3 + B.X**2 + C.X + D
! Solution 1 = X1
! Solution 2 = X2 + X2I * j
! Solution 3 = X3
!----------------------------------------
  IMPLICIT NONE
  REAL(dp),INTENT(IN) :: B,C,D
  REAL(dp),INTENT(OUT):: X1,X2,X2I,X3
  !
  REAL(dp):: Q,P,R,PHI3,FF
  !
  X1=  Zero
  X2=  Zero
  X2I= Zero
  X3=  Zero
  IF (C==Zero .AND. D==Zero) THEN
    X1=-B
  ELSE
    Q= (B*B*B*Two/27.D0 -B*C/3.0D0 +D) /Two
    P= (C    *3.0_dp    -B*B         ) /9.0D0
    FF=ABS(P)
    R=SQRT(FF)
    FF=R*Q
    IF (FF < Zero) R=-R
    FF=Q/(R*R*R)
    IF (P > Zero) THEN
      PHI3= LOG( FF+SQRT(FF*FF+One) ) /3.0D0
      X1=  -R*( EXP(PHI3) -EXP(-PHI3) ) -B/3.0D0
      X2I=  One
    ELSE
      IF (Q*Q+P*P*P > Zero) THEN
        PHI3= LOG( FF+SQRT(FF*FF-One) ) /3.0D0
        X1=  -R*( EXP(PHI3)+EXP(-PHI3) ) -B/3.0D0
        X2I=  One
      ELSE
        PHI3=ATAN(SQRT(One-FF*FF)/FF)/3.0D0
        X1=-R*COS( PHI3        )*Two -B/3.0D0
        X2= R*COS(PI/3.0D0-PHI3)*Two -B/3.0D0
        X3= R*COS(PI/3.0D0+PHI3)*Two -B/3.0D0
        X2I=Zero
      ENDIF
    ENDIF
  ENDIF
  
  RETURN
ENDSUBROUTINE CUBIC

!-----------------------------------------------------------------

SUBROUTINE KUBIK(b,c,d,x1,x2,x2i,x3)
!----------------------------------------
! Compute Solution of a Cubic Equation 
! Another Method ?? 
! Equation  = X**3 + B.X**2 + C.X + D
! Solution 1 = x1
! Solution 2 = x2 + x2I * j
! Solution 3 = x3
!----------------------------------------
  REAL(dp),INTENT(in) :: b,c,d
  REAL(dp),INTENT(out):: x1,x2,x2i,x3
  REAL(dp):: ff,p,phi3,q,r
  !
  x2=  Zero
  x2i= Zero
  x3=  Zero
  IF ((c==Zero).AND.(d==Zero)) THEN
    x1=-b
    RETURN
  ENDIF
  q=  (2.0*f_Cube(b)/27.0 - b*c/3.0 + d)/2.0
  p=  (3.0*c - b*b)/9.0
  ff=  ABS(p)
  r=  SQRT(ff)
  ff=  r*q
  IF (ff<Zero) r=  -r
  ff=  q/f_Cube(r)
  IF (p>Zero) THEN
    phi3=  LOG(ff + SQRT(ff*ff+One))/3.0
    x1=   -r*(EXP(phi3) - EXP(-phi3)) - b/3.0
    x2i=   One
  ELSE
    IF (q*q + p*p*p > Zero) THEN
      phi3= LOG(ff + SQRT(ff*ff-One))/3.0
      x1=  -r*(EXP(phi3) + EXP(-phi3)) - b/3.0
      x2i=  One
    ELSE 
      phi3=  ATAN(SQRT(One-ff*ff)/ff)/3.0
      x1= -2.0*r*COS(phi3) - b/3.0
      x2=  2.0*r*COS(pi/3.0-phi3) - b/3.0
      x2i= Zero
      x3=  2.0*r*COS(pi/3.0+phi3) - b/3.0
    ENDIF
  ENDIF
  !
  RETURN 
END SUBROUTINE Kubik

ENDMODULE M_Eos_Utils
