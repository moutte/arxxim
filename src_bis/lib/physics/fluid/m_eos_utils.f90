module M_Eos_Utils
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

  use M_Kinds
  implicit none
  private
  real(dp), parameter :: Pi= 3.14159265359979_dp 
 
  public :: f_Square
  public :: f_Cube
  public :: f_Quart
  public :: f_Quint

  public :: CUBIC
  public :: KUBIK

contains

  !----------------------------------------------------

  real(dp) function f_Square(x)
    real(dp):: x
    f_Square=  x*x
    return
  end function f_Square

  real(dp) function f_Cube(x)
    real(dp):: x
    f_Cube=  x*x*x
    return
  end function f_Cube

  real(dp) function f_Quart(x)
    real(dp):: x
    f_Quart=  x*x*x*x
    return
  end function f_Quart

  real(dp) function f_Quint(x)
    real(dp):: x
    f_Quint=  x*x*x*x*x
    return
  end function f_Quint

!----------------------------------------------------

subroutine CUBIC(B,C,D,X1,X2,X2I,X3) 
!----------------------------------------
! Compute Solution of a Cubic Equation 
! Equation  = X**3 + B.X**2 + C.X + D
! Solution 1 = X1
! Solution 2 = X2 + X2I * j
! Solution 3 = X3
!----------------------------------------
  implicit none
  real(dp),intent(in) :: B,C,D
  real(dp),intent(out):: X1,X2,X2I,X3
  !
  real(dp):: Q,P,R,PHI3,FF
  !
  X1=  Zero
  X2=  Zero
  X2I= Zero
  X3=  Zero
  if (C==Zero .and. D==Zero) then
    X1=-B
  else
    Q= (B*B*B*Two/27.D0 -B*C/3.0D0 +D) /Two
    P= (C    *3.0_dp    -B*B         ) /9.0D0
    FF=ABS(P)
    R=SQRT(FF)
    FF=R*Q
    if (FF < Zero) R=-R
    FF=Q/(R*R*R)
    if (P > Zero) then
      PHI3= log( FF+SQRT(FF*FF+One) ) /3.0D0
      X1=  -R*( exp(PHI3) -exp(-PHI3) ) -B/3.0D0
      X2I=  One
    else
      if (Q*Q+P*P*P > Zero) then
        PHI3= log( FF+SQRT(FF*FF-One) ) /3.0D0
        X1=  -R*( exp(PHI3)+exp(-PHI3) ) -B/3.0D0
        X2I=  One
      else
        PHI3=ATAN(SQRT(One-FF*FF)/FF)/3.0D0
        X1=-R*COS( PHI3        )*Two -B/3.0D0
        X2= R*COS(PI/3.0D0-PHI3)*Two -B/3.0D0
        X3= R*COS(PI/3.0D0+PHI3)*Two -B/3.0D0
        X2I=Zero
      end if
    end if
  end if
  
  return
end subroutine CUBIC

!-----------------------------------------------------------------

subroutine KUBIK(b,c,d,x1,x2,x2i,x3)
!----------------------------------------
! Compute Solution of a Cubic Equation 
! Another Method ?? 
! Equation  = X**3 + B.X**2 + C.X + D
! Solution 1 = x1
! Solution 2 = x2 + x2I * j
! Solution 3 = x3
!----------------------------------------
  real(dp),intent(in) :: b,c,d
  real(dp),intent(out):: x1,x2,x2i,x3
  real(dp):: ff,p,phi3,q,r
  !
  x2=  Zero
  x2i= Zero
  x3=  Zero
  if ((c==Zero).and.(d==Zero)) then
    x1=-b
    return
  end if
  q=  (2.0*f_Cube(b)/27.0 - b*c/3.0 + d)/2.0
  p=  (3.0*c - b*b)/9.0
  ff=  ABS(p)
  r=  SQRT(ff)
  ff=  r*q
  if (ff<Zero) r=  -r
  ff=  q/f_Cube(r)
  if (p>Zero) then
    phi3=  log(ff + SQRT(ff*ff+One))/3.0
    x1=   -r*(exp(phi3) - exp(-phi3)) - b/3.0
    x2i=   One
  else
    if (q*q + p*p*p > Zero) then
      phi3= log(ff + SQRT(ff*ff-One))/3.0
      x1=  -r*(exp(phi3) + exp(-phi3)) - b/3.0
      x2i=  One
    else 
      phi3=  ATAN(SQRT(One-ff*ff)/ff)/3.0
      x1= -2.0*r*COS(phi3) - b/3.0
      x2=  2.0*r*COS(pi/3.0-phi3) - b/3.0
      x2i= Zero
      x3=  2.0*r*COS(pi/3.0+phi3) - b/3.0
    end if
  end if
  !
  return 
end subroutine Kubik

end module M_Eos_Utils
