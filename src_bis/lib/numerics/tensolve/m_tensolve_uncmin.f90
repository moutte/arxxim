module M_TenSolve_UncMin
  ! UNCMIN (R. B. Schnabel, J. E. Koontz, and B. E. Weiss,
  ! "A Modular System of Algorithms of Unconstrained Minimization",
  ! ACM Trans. Math. Softw., 11 (1985), 419-440).
  ! Translated to free-format Fortran 90 style by Alan Miller
  ! Alan.Miller @ vic.cmis.csiro.au
  ! Latest revision - 10 November 2003
  ! 24 May 2001: Automatic array wk(n) added to routine lltslv.
  ! 10 Nov. 2003: Removed the private declaration from dp
  use M_Kinds
  !
  use M_TenSolve_BlasPartial
  implicit none
  !
contains
  !
subroutine bakslv(n, a, x, b)
  ! PURPOSE
  ! -------
  ! SOLVE  AX=B  where A IS UPPER TRIANGULAR MATRIX.
  ! NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND
  ! THAT THIS ROUTINE TAKES ITS TRANSPOSE implicitLY.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED)
  ! X(N)        <--  SOLUTION VECTOR
  ! B(N)         --> RIGHT-HAND SIDE VECTOR
  ! NOTE
  ! ----
  ! if B IS NO LONGER REQUIRED BY callING ROUTINE,
  ! then VECTORS B AND X MAY SHARE THE SAME STORAGE.
  integer, intent(in)     :: n
  real(dp),intent(in)   :: a(:,:)
  real(dp),intent(out)  :: x(:)
  real(dp),intent(in)   :: b(:)
  integer   :: i,ip1
  real(dp) :: sum
  ! SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE)
  i = n
  x(i) = b(i) / a(i,i)
  if(n == 1) return
  do
    ip1=i
    i=i-1
    sum = dot_product( a(ip1:n,i), x(ip1:n) )
    x(i) = (b(i) - sum)/a(i,i)
    if(i <= 1) exit
  end do
  return
end subroutine bakslv

subroutine choldc(n, a, diagmx, tol, addmax)
  ! PURPOSE
  ! -------
  ! FIND THE PERTURBED L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION
  ! OF A+D, where D IS A NON-NEGATIVE DIAGONAL MATRIX ADDED TO A if
  ! NECESSARY TO ALLOW THE CHOLESKY DECOMPOSITION TO continue.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! A(N,N)      <--> ON entry: MATRIX FOR WHICH TO FIND PERTURBED
  !                       CHOLESKY DECOMPOSITION
  !                  ON exit:  contains L OF LL+ DECOMPOSITION
  !                  IN LOWER TRIANGULAR PART AND DIAGONAL OF "A"
  ! DIAGMX       --> MAXIMUM DIAGONAL ELEMENT OF "A"
  ! TOL          --> TOLERANCE
  ! ADDMAX      <--  MAXIMUM AMOUNT implicitLY ADDED TO DIAGONAL OF "A"
  !                  IN FORMING THE CHOLESKY DECOMPOSITION OF A+D
  ! INTERNAL VARIABLES
  ! ------------------
  ! AMINL    SMALLEST ELEMENT ALLOWED ON DIAGONAL OF L
  ! AMNLSQ   = AMINL**2
  ! OFFMAX   MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN OF A
  ! DESCRIPTION
  ! -----------
  ! THE NORMAL CHOLESKY DECOMPOSITION IS PERFORMED.  HOWEVER, if AT ANY
  ! POINT THE ALGORITHM WOULD ATTEMPT TO SET L(I,I)=SQRT(TEMP)
  ! WITH TEMP < TOL*DIAGMX, then L(I,I) IS SET TO SQRT(TOL*DIAGMX)
  ! INSTEAD.  THIS IS EQUIVAlenT TO ADDING TOL*DIAGMX-TEMP TO A(I,I)
  integer, intent(in)   :: n
  real(dp),intent(inout):: a(:,:)
  real(dp),intent(in)   :: diagmx
  real(dp),intent(in)   :: tol
  real(dp),intent(out), optional :: addmax
  integer   :: j, i
  real(dp) :: aminl, total, temp, amnlsq, offmax
  if (present(addmax)) addmax = zero
  aminl = SQRT(diagmx*tol)
  amnlsq = aminl*aminl
  ! FORM COLUMN J OF L
  do j = 1,n
    ! FIND DIAGONAL ELEMENTS OF L
    temp = a(j,j) - SUM( a(j,1:j-1)**2 )
    if(temp > amnlsq) then
      a(j,j) = SQRT(temp)
    else
    ! FIND MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN
      offmax = zero
      do i = j+1,n
        if(ABS(a(i,j)) > offmax) offmax = ABS(a(i,j))
      end do
      if(offmax <= amnlsq) offmax = amnlsq
      ! ADD TO DIAGONAL ELEMENT  TO ALLOW CHOLESKY DECOMPOSITION TO continue
      a(j,j) = SQRT(offmax)
      if (present(addmax)) addmax = MAX(addmax,offmax-temp)
    end if
    ! FIND I,J ELEMENT OF LOWER TRIANGULAR MATRIX
    do i = j+1,n
      total = dot_product( a(i,1:j-1), a(j,1:j-1) )
      a(i,j) = (a(i,j) - total) / a(j,j)
    end do
  end do
  return
end subroutine choldc

subroutine chlhsn(n, a, epsm, sx, udiag)
  ! PURPOSE
  ! -------
  ! FIND THE L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION OF THE PERTURBED
  ! MODEL HESSIAN MATRIX A + MU*I (where MU\=0 AND I IS THE IDENTITY MATRIX)
  ! WHICH IS SAFELY POSITIVE DEFINITE.  if A IS SAFELY POSITIVE DEFINITE
  ! UPON entry, then MU = 0.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! A(N,N)      <--> ON entry; "A" IS MODEL HESSIAN (only LOWER
  !                  TRIANGULAR PART AND DIAGONAL STORED)
  !                  ON exit:  A contains L OF LL+ DECOMPOSITION OF PERTURBED
  !                  MODEL HESSIAN IN LOWER TRIANGULAR PART AND DIAGONAL
  !                  AND contains HESSIAN IN UPPER TRIANGULAR PART AND UDIAG
  ! EPSM         --> MACHINE EPSILON
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! UDIAG(N)    <--  ON exit: contains DIAGONAL OF HESSIAN
  ! INTERNAL VARIABLES
  ! ------------------
  ! TOL              TOLERANCE
  ! DIAGMN           MINIMUM ELEMENT ON DIAGONAL OF A
  ! DIAGMX           MAXIMUM ELEMENT ON DIAGONAL OF A
  ! OFFMAX           MAXIMUM OFF-DIAGONAL ELEMENT OF A
  ! OFFROW           SUM OF OFF-DIAGONAL ELEMENTS IN A ROW OF A
  ! EVMIN            MINIMUM EIGENvalue OF A
  ! EVMAX            MAXIMUM EIGENvalue OF A
  ! DESCRIPTION
  ! -----------
  ! 1. if "A" HAS ANY NEGATIVE DIAGONAL ELEMENTS, then CHOOSE MU>0
  ! SUCH THAT THE DIAGONAL OF A:=A+MU*I IS ALL POSITIVE
  ! WITH THE RATIO OF ITS SMALLEST TO LARGEST ELEMENT ON THE ORDER OF SQRT(EPSM).
  ! 2. "A" UNDERGOES A PERTURBED CHOLESKY DECOMPOSITION WHICH resultS IN AN LL+
  ! DECOMPOSITION OF A+D, where D IS A NON-NEGATIVE DIAGONAL MATRIX WHICH IS
  ! implicitLY ADDED TO "A" DURING THE DECOMPOSITION if "A" IS NOT POSITIVE
  ! DEFINITE.
  ! "A" IS RETAINED AND NOT CHANGED DURING THIS PROCESS BY COPYING L INTO THE
  ! UPPER TRIANGULAR PART OF "A" AND THE DIAGONAL INTO UDIAG.  then THE CHOLESKY
  ! DECOMPOSITION ROUTINE IS callED.  ON return, ADDMAX contains MAXIMUM ELEMENT
  ! OF D.
  ! 3. if ADDMAX = 0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2 AND return
  ! IS MADE TO callING program.  OTHERWISE, THE MINIMUM NUMBER SDD WHICH MUST
  ! BE ADDED TO THE DIAGONAL OF A TO MAKE IT SAFELY STRICTLY DIAGONALLY doMINANT
  ! IS CALCULATED.  SINCE A+ADDMAX*I AND A+SDD*I ARE SAFELY POSITIVE DEFINITE,
  ! CHOOSE MU = MIN(ADDMAX,SDD) AND DECOMPOSE A+MU*I TO OBTAIN L.
  integer, intent(in)        :: n
  real(dp),intent(inout)  :: a(:,:)
  real(dp),intent(in)      :: epsm
  real(dp),intent(in)      :: sx(:)
  real(dp),intent(out)     :: udiag(:)
  integer   :: i, j
  real(dp) :: tol, diagmx, diagmn, posmax, amu, offmax
  real(dp) :: addmax, evmin, evmax, offrow, sdd
  ! SCALE HESSIAN
  ! PRE- AND POST- MULTIPLY "A" BY INV(SX)
  do j = 1,n
    a(j:n,j) = a(j:n,j) / (sx(j:n)*sx(j))
  end do
  ! STEP1
  ! -----
  ! NOTE:  if A DifFERENT TOLERANCE IS DESIRED THROUGHOUT THIS
  ! ALGORITHM, CHANGE TOLERANCE HERE:
  tol = SQRT(epsm)
  diagmx = a(1,1)
  diagmn = a(1,1)
  do i = 2,n
    if(a(i,i) < diagmn) diagmn = a(i,i)
    if(a(i,i) > diagmx) diagmx = a(i,i)
  end do
  posmax = MAX(diagmx, zero)
  ! DIAGMN <= 0
  if(diagmn <= posmax*tol) then
    amu = tol*(posmax - diagmn) - diagmn
    if(amu == zero) then
    ! FIND LARGEST OFF-DIAGONAL ELEMENT OF A
      offmax = zero
      do i = 2,n
        do j = 1,i-1
          if(ABS(a(i,j)) > offmax) offmax = ABS(a(i,j))
        end do
      end do
      amu = offmax
      if(amu == zero) then
        amu = one
      else
        amu = amu*(one + tol)
      end if
    end if
    ! A = A + MU*I
    do i = 1,n
      a(i,i) = a(i,i) + amu
    end do
    diagmx = diagmx + amu
  end if
  ! STEP2
  ! -----
  ! COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART
  ! AND DIAGONAL OF "A" TO UDIAG
  do j = 1,n
    udiag(j) = a(j,j)
    if(j == n) cycle
    a(j,j+1:n) = a(j+1:n,j)
  end do
  call choldc(n, a, diagmx, tol, addmax)
  ! STEP3
  ! -----
  ! if ADDMAX = 0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2,
  ! THE LL+ DECOMPOSITION HAS BEEN doNE, AND WE return.
  ! OTHERWISE, ADDMAX >0.  PERTURB "A" SO THAT IT IS SAFELY
  ! DIAGONALLY doMINANT AND FIND LL+ DECOMPOSITION
  if(addmax > zero) then
    ! RESTORE ORIGINAL "A" (LOWER TRIANGULAR PART AND DIAGONAL)
    do j = 1,n
      a(j,j) = udiag(j)
      do i = j+1,n
        a(i,j) = a(j,i)
      end do
    end do
    ! FIND SDD SUCH THAT A+SDD*I IS SAFELY POSITIVE DEFINITE
    ! NOTE:  EVMIN<0 SINCE A IS NOT POSITIVE DEFINITE;
    evmin = zero
    evmax = a(1,1)
    do i = 1,n
      offrow = zero
      do j = 1,i-1
        offrow = offrow + ABS(a(i,j))
      end do
      do j = i+1,n
        offrow = offrow + ABS(a(j,i))
      end do
      evmin = MIN(evmin, a(i,i) - offrow)
      evmax = MAX(evmax, a(i,i) + offrow)
    end do
    sdd = tol*(evmax-evmin) - evmin
    ! PERTURB "A" AND DECOMPOSE AGAIN
    amu = MIN(sdd,addmax)
    do i = 1,n
      a(i,i) = a(i,i) + amu
      udiag(i) = a(i,i)
    end do
    ! "A" NOW GUARANTEED SAFELY POSITIVE DEFINITE
    call choldc(n, a, zero, tol, addmax)
  end if
  ! UNSCALE HESSIAN AND CHOLESKY DECOMPOSITION MATRIX
  do j = 1,n
    a(j:n,j) = sx(j:n)*a(j:n,j)
    a(1:j-1,j) = sx(1:j-1)*sx(j)*a(1:j-1,j)
    udiag(j) = udiag(j)*sx(j)*sx(j)
  end do
  return
end subroutine chlhsn

subroutine dfaut( &
! SET default valueS FOR EACH INPUT VARIABLE TO MINIMIZATION ALGORITHM.
& n, typsiz, fscale, method, iexp, msg, ndigit,  &
& itnlim, iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl)
  ! parameterS
  ! N            --> dimension OF PROBLEM
  ! TYPSIZ(N)   <--  TYPICAL size FOR EACH COMPONENT OF X
  ! FSCALE      <--  ESTIMATE OF SCALE OF MINIMIZATION function
  ! METHOD      <--  ALGORITHM TO use TO SOLVE MINIMIZATION PROBLEM
  ! IEXP        <--  = 0 if MINIMIZATION function NOT EXPENSIVE TO EVALUATE
  ! MSG         <--  MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  ! NDIGIT      <--  NUMBER OF GOOD DIGITS IN MINIMIZATION function
  ! ITNLIM      <--  MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  ! IAGFLG      <--  = 0 if ANALYTIC GRADIENT NOT SUPPLIED
  ! IAHFLG      <--  = 0 if ANALYTIC HESSIAN NOT SUPPLIED
  ! IPR         <--  DEVICE TO WHICH TO Send OUTPUT
  ! DLT         <--  TRUST REGION RADIUS
  ! GRADTL      <--  TOLERANCE AT WHICH GRADIENT CONSIDERED close ENOUGH
  !                  TO ZERO TO TERMINATE ALGORITHM
  ! STEPMX      <--  value OF ZERO TO TRIP default MAXIMUM IN OPTCHK
  ! STEPTL      <--  TOLERANCE AT WHICH SUCCESSIVE ITERATES CONSIDERED
  !                  close ENOUGH TO TERMINATE ALGORITHM
  integer, intent(in) :: n
  real(dp),intent(out):: typsiz(:)
  real(dp),intent(out):: fscale
  integer, intent(out):: method
  integer, intent(out):: iexp
  integer, intent(out):: msg
  integer, intent(out):: ndigit
  integer, intent(out):: itnlim
  integer, intent(out):: iagflg
  integer, intent(out):: iahflg
  integer, intent(out):: ipr
  real(dp),intent(out):: dlt
  real(dp),intent(out):: gradtl
  real(dp),intent(out):: stepmx
  real(dp),intent(out):: steptl
  real(dp), parameter :: epsm = EPSILON( 1.0_dp ), three = 3.0_dp
  ! SET TYPICAL size OF X AND MINIMIZATION function
  typsiz(1:n) = one
  fscale = one
  ! SET TOLERANCES
  dlt = -one
  gradtl = epsm**(one/three)
  stepmx = zero
  steptl = SQRT(epsm)
  ! SET FLAGS
  method = 1
  iexp = 1
  msg = 0
  ndigit = -1
  itnlim = 150
  iagflg = 0
  iahflg = 0
  ipr = 6
  return
  end subroutine dfaut
  subroutine dogdrv(n, x, f, g, a, p, xpls, fpls, fcn, sx, stepmx,  &
                    steptl, dlt, iretcd, mxtake)
  ! PURPOSE
  ! -------
  ! FIND A NEXT NEWTON ITERATE (XPLS) BY THE doUBLE doGLEG METHOD
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! X(N)         --> OLD ITERATE X[K-1]
  ! F            --> function value AT OLD ITERATE, F(X)
  ! G(N)         --> GRADIENT  AT OLD ITERATE, G(X), OR APPROXIMATE
  ! A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN
  !                  IN LOWER TRIANGULAR PART AND DIAGONAL
  ! P(N)         --> NEWTON STEP
  ! XPLS(N)     <--  NEW ITERATE X[K]
  ! FPLS        <--  function value AT NEW ITERATE, F(XPLS)
  ! FCN          --> NAME OF subroutine TO EVALUATE function
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! DLT         <--> TRUST REGION RADIUS
  !                  [RETAIN value BETWEEN SUCCESSIVE callS]
  ! IRETCD      <--  return CODE
  !                    = 0 SATISFACTORY XPLS FOUND
  !                    = 1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
  !                        DISTINCT FROM X
  ! MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  integer, intent(in) :: n
  real(dp),intent(in) :: x(:)
  real(dp),intent(in) :: f
  real(dp),intent(in) :: g(:)
  real(dp),intent(in) :: a(:,:)
  real(dp),intent(in) :: p(:)
  real(dp),intent(out):: xpls(:)
  real(dp),intent(out):: fpls
  real(dp),intent(in) :: sx(:)
  real(dp),intent(in) :: stepmx
  real(dp),intent(in) :: steptl
  real(dp),intent(inout):: dlt
  integer, intent(out):: iretcd
  logical, intent(out):: mxtake
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)    :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f
    end subroutine fcn
  end interface
  ! Workspace
  real(dp) :: wrk1(n), wrk2(n), wrk3(n), sc(n)
  real(dp) :: rnwtln, cln
  real(dp) :: eta, fplsp
  logical   :: fstdog, nwtake
  iretcd = 4
  fstdog = .true.
  rnwtln = dnrm2(n, sx(1:n)*p(1:n), 1)
  ! FIND NEW STEP BY doUBLE doGLEG ALGORITHM
  100 call dogstp(n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog,  &
                  wrk1, wrk2, cln, eta, sc, stepmx)
  ! CHECK NEW POINT AND UPDATE TRUST REGION
  call tregup(n, x, f, g, a, fcn, sc, sx, nwtake, stepmx, steptl, dlt,  &
              iretcd, wrk3, fplsp, xpls, fpls, mxtake, 2, wrk1)
  if(iretcd <= 1) return
  GO TO 100
end subroutine dogdrv

subroutine dogstp(n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog,  &
                    ssd, v, cln, eta, sc, stepmx)
  ! PURPOSE
  ! -------
  ! FIND NEW STEP BY doUBLE doGLEG ALGORITHM
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
  ! A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
  !                  LOWER PART AND DIAGONAL
  ! P(N)         --> NEWTON STEP
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! RNWTLN       --> NEWTON STEP lenGTH
  ! DLT         <--> TRUST REGION RADIUS
  ! NWTAKE      <--> BOOLEAN, =.true. if NEWTON STEP TAKEN
  ! FSTdoG      <--> BOOLEAN, =.true. if ON FIRST LEG OF doGLEG
  ! SSD(N)      <--> WORKSPACE [CAUCHY STEP TO THE MINIMUM OF THE
  !                  QUADRATIC MODEL IN THE SCALED STEEPEST DESCENT
  !                  DIRECTION] [RETAIN value BETWEEN SUCCESSIVE callS]
  ! V(N)        <--> WORKSPACE  [RETAIN value BETWEEN SUCCESSIVE callS]
  ! CLN         <--> CAUCHY lenGTH
  !                  [RETAIN value BETWEEN SUCCESSIVE callS]
  ! ETA              [RETAIN value BETWEEN SUCCESSIVE callS]
  ! SC(N)       <--  CURRENT STEP
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! INTERNAL VARIABLES
  ! ------------------
  ! CLN              lenGTH OF CAUCHY STEP
  integer, intent(in):: n
  real(dp),intent(in):: g(:)
  real(dp),intent(in):: a(:,:)
  real(dp),intent(in):: p(:)
  real(dp),intent(in):: sx(:)
  real(dp),intent(in):: rnwtln
  real(dp),intent(inout):: dlt
  logical, intent(inout):: nwtake
  logical, intent(inout):: fstdog
  real(dp),intent(inout):: ssd(:)
  real(dp),intent(inout):: v(:)
  real(dp),intent(inout):: cln
  real(dp),intent(inout):: eta
  real(dp),intent(out):: sc(:)
  real(dp),intent(in):: stepmx
  
  integer   :: i, j
  real(dp) :: alpha, beta, tmp
  real(dp) :: dot1, dot2, alam
  
  ! CAN WE TAKE NEWTON STEP
  if(rnwtln > dlt) GO TO 100
  
  nwtake = .true.
  sc(1:n) = p(1:n)
  dlt = rnwtln
  GO TO 700
  
  ! NEWTON STEP TOO LONG
  ! CAUCHY STEP IS ON doUBLE doGLEG CURVE
  100 nwtake = .false.
  if(.not.fstdog) GO TO 200
  !         CALCULATE doUBLE doGLEG CURVE (SSD)
  fstdog= .false.
  ssd(1:n)= g(1:n)/sx(1:n)
  alpha= SUM( ssd(1:n)**2 )
  beta= zero
  do i= 1,n
    tmp= zero
    do j= i,n
      tmp= tmp + (a(j,i)/sx(j))*ssd(j)
    end do
    beta= beta + tmp*tmp
  end do
  ssd(1:n) = -(alpha/beta)*ssd(1:n)
  cln = alpha*SQRT(alpha)/beta
  eta = .2 + (.8*alpha*alpha) / (-beta * dot_product( g(1:n), p(1:n) ))
  v(1:n) = eta*sx(1:n)*p(1:n) - ssd(1:n)
  if (dlt == -one) dlt = MIN(cln, stepmx)
  200 if(eta*rnwtln > dlt) GO TO 220
  !         TAKE PARTIAL STEP IN NEWTON DIRECTION
  sc(1:n) = (dlt/rnwtln)*p(1:n)
  GO TO 700
  220 if(cln < dlt) GO TO 240
  !         if(CLN >= DLT) then TAKE STEP IN STEEPEST DESCENT DIRECTION
  sc(1:n) = (dlt/cln)*ssd(1:n)/sx(1:n)
  GO TO 700
  !           CALCULATE CONVEX COMBINATION OF SSD AND ETA*P
  !           WHICH HAS SCALED lenGTH DLT
  240 dot1 = dot_product( v(1:n), ssd(1:n) )
  dot2 = SUM( v(1:n)**2 )
  alam = (-dot1 + SQRT((dot1*dot1) - dot2*(cln*cln - dlt*dlt)))/dot2
  sc(1:n) = (ssd(1:n) + alam*v(1:n))/sx(1:n)
  700 return
end subroutine dogstp
  
subroutine forslv (n, a, x, b)
  ! PURPOSE
  ! --------
  ! SOLVE  AX = B  where A  IS LOWER TRIANGULAR  MATRIX
  ! parameterS
  ! ---------
  ! N         -----> dimension OF PROBLEM
  ! A(N,N)    -----> LOWER TRIANGULAR MATRIX (PRESERVED)
  ! X(N)      <----  SOLUTION VECTOR
  ! B(N)       ----> RIGHT-HAND SIDE VECTOR
  ! NOTE
  !-----
  ! then VECTORS B AND X MAY SHARE THE SAME STORAGE
  integer, intent(in)     :: n
  real(dp),intent(in)   :: a(:,:)
  real(dp),intent(out)  :: x(:)
  real(dp),intent(in)   :: b(:)
  integer   :: i
  real(dp) :: sum
  ! SOLVE LX = B.  (FORWARD  SOLVE)
  x(1) = b(1)/a(1,1)
  do i = 2,n
    sum = dot_product( a(i,1:i-1), x(1:i-1) )
    x(i) = (b(i) - sum)/a(i,i)
  end do
  return
end subroutine forslv

subroutine fstocd (n, x, fcn, sx, rnoise, g)
  ! PURPOSE
  ! -------
  ! FIND CENTRAL DifFERENCE APPROXIMATION G TO THE FIRST DERIVATIVE
  ! (GRADIENT) OF THE function DEFINED BY FCN AT THE POINT X.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! X            --> POINT AT WHICH GRADIENT IS TO BE APPROXIMATED.
  ! FCN          --> NAME OF subroutine TO EVALUATE function.
  ! SX           --> DIAGONAL SCALING MATRIX FOR X.
  ! RNOISE       --> RELATIVE NOISE IN FCN [F(X)].
  ! G           <--  CENTRAL DifFERENCE APPROXIMATION TO GRADIENT.
  integer, intent(in)        :: n
  real(dp),intent(inout)  :: x(:)
  real(dp),intent(in)      :: sx(:)
  real(dp),intent(in)      :: rnoise
  real(dp),intent(out)     :: g(:)
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)    :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f
    end subroutine fcn
  end interface
  real(dp) :: third, xtempi, fplus, fminus, stepi
  integer   :: i
  ! FIND I TH  STEPsize, EVALUATE TWO NEIGHBORS IN DIRECTION OF I TH
  ! UNIT VECTOR, AND EVALUATE I TH  COMPONENT OF GRADIENT.
  third = one/3.0
  do i = 1, n
    stepi = rnoise**third * MAX(ABS(x(i)), one/sx(i))
    xtempi = x(i)
    x(i) = xtempi + stepi
    call fcn (n, x, fplus)
    x(i) = xtempi - stepi
    call fcn (n, x, fminus)
    x(i) = xtempi
    g(i) = (fplus - fminus)/(2.0*stepi)
  end do
  return
end subroutine fstocd

subroutine fstofd(m, n, xpls, fcn, fpls, a, sx, rnoise, icase)
  ! PURPOSE
  ! -------
  ! FIND FIRST ORDER FORWARD FINITE DifFERENCE APPROXIMATION "A" TO THE
  ! FIRST DERIVATIVE OF THE function DEFINED BY THE SUBprogram "FNAME"
  ! EVALUATED AT THE NEW ITERATE "XPLS".
  ! FOR OPTIMIZATION use THIS ROUTINE TO ESTIMATE:
  ! 1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION function "FCN
  !    ANALYTIC useR ROUTINE HAS BEEN SUPPLIED;
  ! 2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION function
  !    if NO ANALYTIC useR ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT
  !    ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND if THE
  !    OPTIMIZATION function IS INEXPENSIVE TO EVALUATE
  ! NOTE
  ! ----
  ! _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE function
  !      (FCN).   FCN(X) # F: R(N)-->R(1)
  ! _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE function
  !      FCN(X) # F: R(N)-->R(N).
  ! _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO
  !      function, where THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN"
  ! parameterS
  ! ----------
  ! M            --> NUMBER OF ROWS IN A
  ! N            --> NUMBER OF COLUMNS IN A; dimension OF PROBLEM
  ! XPLS(N)      --> NEW ITERATE:  X[K]
  ! FCN          --> NAME OF subroutine TO EVALUATE function
  ! FPLS(M)      --> _M=1 (OPTIMIZATION) function value AT NEW ITERATE:
  !                       FCN(XPLS)
  !                  _M=N (OPTIMIZATION) value OF FIRST DERIVATIVE
  !                       (GRADIENT) GIVEN BY useR function FCN
  !                  _M=N (SYSTEMS)  function value OF associateD
  !                       MINIMIZATION function
  ! A(NR,N)     <--  FINITE DifFERENCE APPROXIMATION (SEE NOTE).  only
  !                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE returnED
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! RNOISE       --> RELATIVE NOISE IN FCN [F(X)]
  ! Icase        --> =1 OPTIMIZATION (GRADIENT)
  !                  =2 SYSTEMS
  !                  =3 OPTIMIZATION (HESSIAN)
  ! INTERNAL VARIABLES
  ! ------------------
  ! STEPSZ - STEPsize IN THE J-TH VARIABLE DIRECTION
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(inout)  :: xpls(:)
  real(dp),intent(in)      :: fpls(:)
  real(dp),intent(out)     :: a(:,:)
  real(dp),intent(in)      :: sx(:)
  real(dp),intent(in)      :: rnoise
  integer, intent(in)        :: icase
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in) :: n
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: f
    end subroutine fcn
  end interface
  ! Workspace
  real(dp):: fhat(m)
  integer :: i, j, nm1, jp1
  real(dp):: stepsz, xtmpj, sqrtr, rstepsz
  real(dp), parameter :: half = 0.5_dp
  ! FIND J-TH COLUMN OF A
  ! J-TH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J)
  sqrtr = SQRT(rnoise)
  do j = 1,n
    stepsz = sqrtr * MAX(ABS(xpls(j)), one/sx(j))
    xtmpj = xpls(j)
    xpls(j) = xtmpj + stepsz
    call fcn(n, xpls, fhat(j))
    xpls(j) = xtmpj
    rstepsz = one/stepsz
    do i = 1,m
      a(i,j) = (fhat(i) - fpls(i))*rstepsz
    end do
  end do
  if(icase /= 3) return
  ! if COMPUTING HESSIAN, A MUST BE SYMMETRIC
  if(n == 1) return
  nm1 = n-1
  do j = 1,nm1
    jp1 = j+1
    do i = jp1,m
      a(i,j) = (a(i,j) + a(j,i))*half
    end do
  end do
  return
end subroutine fstofd

subroutine estimate_Hessian(m, n, xpls, d1fcn, fpls, a, sx, rnoise, icase)
  ! PURPOSE
  ! -------
  ! FIND FIRST ORDER FORWARD FINITE DifFERENCE APPROXIMATION "A" TO THE
  ! HESSIAN OF THE function DEFINED BY THE SUBprogram "FNAME"
  ! EVALUATED AT THE NEW ITERATE "XPLS".
  ! FOR OPTIMIZATION use THIS ROUTINE TO ESTIMATE:
  !   THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION function
  !   if NO ANALYTIC useR ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT
  !   ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND if THE
  !   OPTIMIZATION function IS INEXPENSIVE TO EVALUATE
  ! NOTE
  ! ----
  ! _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE function
  !      (FCN).   FCN(X) # F: R(N)-->R(1)
  ! _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE function
  !      FCN(X) # F: R(N)-->R(N).
  ! _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO
  !      function, where THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN"
  ! parameterS
  ! ----------
  ! M            --> NUMBER OF ROWS IN A
  ! N            --> NUMBER OF COLUMNS IN A; dimension OF PROBLEM
  ! XPLS(N)      --> NEW ITERATE:  X[K]
  ! D1FCN        --> NAME OF subroutine TO EVALUATE FIRST DERIVATIVES
  ! FPLS(M)      --> _M=1 (OPTIMIZATION) function value AT NEW ITERATE:
  !                       FCN(XPLS)
  !                  _M=N (OPTIMIZATION) value OF FIRST DERIVATIVE
  !                       (GRADIENT) GIVEN BY useR function FCN
  !                  _M=N (SYSTEMS)  function value OF associateD
  !                       MINIMIZATION function
  ! A(NR,N)     <--  FINITE DifFERENCE APPROXIMATION (SEE NOTE).  only
  !                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE returnED
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! RNOISE       --> RELATIVE NOISE IN FCN [F(X)]
  ! Icase        --> =1 OPTIMIZATION (GRADIENT)
  !                  =2 SYSTEMS
  !                  =3 OPTIMIZATION (HESSIAN)
  ! INTERNAL VARIABLES
  ! ------------------
  ! STEPSZ - STEPsize IN THE J-TH VARIABLE DIRECTION
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(inout)  :: xpls(:)
  real(dp),intent(in)      :: fpls(:)
  real(dp),intent(out)     :: a(:,:)
  real(dp),intent(in)      :: sx(:)
  real(dp),intent(in)      :: rnoise
  integer, intent(in)        :: icase
  interface
    subroutine d1fcn(n, x, g)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)    :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: g(:)
    end subroutine d1fcn
  end interface
  ! Workspace
  real(dp) :: fhat(m)
  integer              :: i, j, nm1, jp1
  real(dp)            :: stepsz, xtmpj, sqrtr, rstepsz
  real(dp), parameter :: half = 0.5_dp
  ! FIND J-TH COLUMN OF A
  ! EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J)
  sqrtr = SQRT(rnoise)
  do j = 1,n
    stepsz = sqrtr * MAX(ABS(xpls(j)), one/sx(j))
    xtmpj = xpls(j)
    xpls(j) = xtmpj + stepsz
    call d1fcn(n, xpls, fhat)
    xpls(j) = xtmpj
    rstepsz = one/stepsz
    do i = 1,m
      a(i,j) = (fhat(i) - fpls(i))*rstepsz
    end do
  end do
  if(icase /= 3) return
  ! if COMPUTING HESSIAN, A MUST BE SYMMETRIC
  if(n == 1) return
  nm1 = n-1
  do j = 1,nm1
    jp1 = j+1
    do i = jp1,m
      a(i,j) = (a(i,j) + a(j,i))*half
    end do
  end do
  return
end subroutine estimate_Hessian

subroutine hookdr(n, g, a, udiag, p, sx, stepmx, dlt, &
                    iretcd, amu, dltp, phi, phip0, epsm, itncnt)
  ! PURPOSE
  ! -------
  ! FIND A NEXT NEWTON ITERATE (XPLS) BY THE MORE-HEBdoN METHOD
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
  ! A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN LOWER
  !                  TRIANGULAR PART AND DIAGONAL.
  !                  HESSIAN IN UPPER TRIANGULAR PART AND UDIAG.
  ! UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
  ! P(N)         --> NEWTON STEP
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! DLT         <--> TRUST REGION RADIUS
  ! IRETCD      <--  return CODE
  !                     = 0 SATISFACTORY XPLS FOUND
  !                     = 1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
  !                       DISTINCT FROM X
  ! AMU         <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! DLTP        <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! PHI         <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! PHIP0       <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! EPSM         --> MACHINE EPSILON
  ! ITNCNT       --> ITERATION count
  integer, intent(in)   :: n
  real(dp),intent(in)   :: g(:)
  real(dp),intent(inout):: a(:,:)
  real(dp),intent(in)   :: udiag(:)
  real(dp),intent(in)   :: p(:)
  real(dp),intent(in)   :: sx(:)
  real(dp),intent(in)   :: stepmx
  real(dp),intent(inout):: dlt
  integer, intent(out)  :: iretcd
  real(dp),intent(inout):: amu
  real(dp),intent(inout):: dltp
  real(dp),intent(inout):: phi
  real(dp),intent(inout):: phip0
  real(dp),intent(in)   :: epsm
  integer, intent(in)   :: itncnt
  !
  integer   :: i
  real(dp) :: tmp, rnwtln, alpha, beta
  logical   :: fstime
  iretcd = 4
  fstime = .true.
  tmp = SUM( (sx(1:n)*p(1:n))**2 )
  rnwtln = SQRT(tmp)
  if(itncnt == 1) then
    amu = zero
    !       if FIRST ITERATION AND TRUST REGION NOT PROVIDED BY useR,
    !       COMPUTE INITIAL TRUST REGION.
    if(dlt == -one) then
      alpha = SUM( (g(1:n)/sx(1:n))**2 )
      beta = zero
      do i = 1,n
        tmp = SUM( a(i:n,i)*g(i:n) / (sx(i:n)*sx(i:n)) )
        beta = beta + tmp*tmp
      end do
      dlt = alpha*SQRT(alpha)/beta
      dlt = MIN(dlt, stepmx)
    end if
  end if
  ! FIND NEW STEP BY MORE-HEBdoN ALGORITHM
  call hookst(n, g, a, udiag, p, sx, rnwtln, dlt, amu,  &
              dltp, phi, phip0, fstime, epsm)
  dltp = dlt
  ! CHECK NEW POINT AND UPDATE TRUST REGION
  !     call TREGUP(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,
  !    +         DLT,IRETCD,XPLSP,FPLSP,XPLS,FPLS,MXTAKE,IPR,3,UDIAG)
  ! if(iretcd <= 1) return
  ! GO TO 100
  return
end subroutine hookdr

subroutine hookst( &
! FIND NEW STEP BY MORE-HEBdoN ALGORITHM
& n, g, a, udiag, p, sx, rnwtln, dlt, amu,  &
& dltp, phi, phip0, fstime, epsm)
  ! parameterS
  ! N            --> dimension OF PROBLEM
  ! G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
  ! A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
  !                  LOWER TRIANGULAR PART AND DIAGONAL.
  !                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
  ! UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
  ! P(N)         --> NEWTON STEP
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR N
  ! RNWTLN       --> NEWTON STEP lenGTH
  ! DLT         <--> TRUST REGION RADIUS
  ! AMU         <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! DLTP         --> TRUST REGION RADIUS AT LAST exit FROM THIS ROUTINE
  ! PHI         <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! PHIP0       <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! FSTIME      <--> BOOLEAN. =.true. if FIRST entry TO THIS ROUTINE
  !                  DURING K-TH ITERATION
  ! SC(N)       <--  CURRENT STEP
  ! NWTAKE      <--  BOOLEAN, =.true. if NEWTON STEP TAKEN
  ! EPSM         --> MACHINE EPSILON
  integer, intent(in)   :: n
  real(dp),intent(in)   :: g(:)
  real(dp),intent(inout):: a(:,:)
  real(dp),intent(in)   :: udiag(:)
  real(dp),intent(in)   :: p(:)
  real(dp),intent(in)   :: sx(:)
  real(dp),intent(in)   :: rnwtln
  real(dp),intent(inout):: dlt
  real(dp),intent(inout):: amu
  real(dp),intent(in)   :: dltp
  real(dp),intent(inout):: phi
  real(dp),intent(inout):: phip0
  logical, intent(inout):: fstime
  real(dp),intent(in)   :: epsm
  ! Workspace
  real(dp) :: wrk0(n), sc(n)
  integer   :: i, jp1, j
  real(dp) :: hi, alo, phip
  real(dp) :: amuup, stepln, amulo
  logical   :: done
  ! HI AND ALO ARE CONSTANTS useD IN THIS ROUTINE.
  ! CHANGE HERE if OTHER valueS ARE TO BE SUBSTITUTED.
  hi = 1.5
  alo = .75
  ! -----
  if(rnwtln <= hi*dlt) then
  !       TAKE NEWTON STEP
    sc(1:n) = p(1:n)
    dlt = MIN(dlt,rnwtln)
    amu = zero
    return
  end if
  !     NEWTON STEP NOT TAKEN
  ! SET PHIP TO 1.0 FOR COMPILATION.  THIS subroutine IS NOT CURRENTLY
  ! useD BY TENSOLVE.
  phip = one
  if(amu > zero) then
    amu = amu - (phi+dltp) *((dltp-dlt) + phi) / (dlt*phip)
  end if
  phi = rnwtln - dlt
  if(fstime) then
    wrk0(1:n) = sx(1:n)*sx(1:n)*p(1:n)
  !         SOLVE L*Y = (SX**2)*P
    call forslv(n, a, wrk0, wrk0)
    phip0 = - dnrm2(n, wrk0, 1)**2 / rnwtln
    fstime = .false.
  end if
  phip = phip0
  amulo = -phi/phip
  amuup = SUM( (g(1:n)/sx(1:n))**2 )
  amuup = SQRT(amuup)/dlt
  done = .false.
  !       TEST value OF AMU; GENERATE NEXT AMU if NECESSARY
  100 if(done) return
  if(amu < amulo .or. amu > amuup) then
    amu = MAX(SQRT(amulo*amuup), amuup*1.0E-3)
  end if
  !     COPY (H,UDIAG) TO L
  !     where H <-- H+AMU*(SX**2) [do NOT ACTUALLY CHANGE (H,UDIAG)]
  do j = 1,n
    a(j,j) = udiag(j) + amu*sx(j)*sx(j)
    if(j == n) cycle
    jp1 = j+1
    do i = jp1,n
      a(i,j) = a(j,i)
    end do
  end do
  !     FACTOR H = L(L+)
  call choldc(n, a, zero, SQRT(epsm))
  !     SOLVE H*P = L(L+)*SC = -G
  wrk0(1:n) = -g(1:n)
  call lltslv(n, a, sc, wrk0)
  !     RESET H.  NOTE SINCE UDIAG HAS NOT BEEN DESTROYED WE NEED do
  !     NOTHING HERE.  H IS IN THE UPPER PART AND IN UDIAG, STILL INTACT
  stepln = SUM( sx(1:n)**2 * sc(1:n)**2 )
  stepln = SQRT(stepln)
  phi = stepln - dlt
  wrk0(1:n) = sx(1:n)*sx(1:n)*sc(1:n)
  call forslv(n, a, wrk0, wrk0)
  phip = - dnrm2(n,wrk0,1)**2 /stepln
  if((alo*dlt > stepln .or. stepln > hi*dlt) .and.  &
      (amuup-amulo > zero)) GO TO 170
  !       SC IS ACCEPTABLE HOOKSTEP
  done = .true.
  GO TO 100
  !       SC NOT ACCEPTABLE HOOKSTEP.  select NEW AMU
  170 amulo = MAX(amulo, amu-(phi/phip))
  if(phi < zero) amuup = MIN(amuup,amu)
  amu = amu - (stepln*phi)/(dlt*phip)
  GO TO 100
end subroutine hookst

subroutine hsnint(n, a, sx, method)
  ! PURPOSE
  ! -------
  ! PROVIDE INITIAL HESSIAN WHEN USING SECANT UPDATES
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! A(N,N)      <--  INITIAL HESSIAN (LOWER TRIANGULAR MATRIX)
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! METHOD       --> ALGORITHM TO use TO SOLVE MINIMIZATION PROBLEM
  !                     = 1,2 FACTORED SECANT METHOD useD
  !                     = 3   UNFACTORED SECANT METHOD useD
  integer, intent(in) :: n
  real(dp),intent(out):: a(:,:)
  real(dp),intent(in) :: sx(:)
  integer, intent(in) :: method
  !
  integer :: j
  do j = 1,n
    if(method == 3) a(j,j) = sx(j)*sx(j)
    if(method /= 3) a(j,j) = sx(j)
    if(j == n) cycle
    a(j+1:n,j) = zero
  end do
  return
end subroutine hsnint

subroutine lltslv(n, a, x, b)
  ! PURPOSE
  ! -------
  ! SOLVE AX = B where A HAS THE FORM L(L-TRANSPOSE)
  ! BUT only THE LOWER TRIANGULAR PART, L, IS STORED.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! A(N,N)       --> MATRIX OF FORM L(L-TRANSPOSE).
  !                  ON return A IS UNCHANGED.
  ! X(N)        <--  SOLUTION VECTOR
  ! B(N)         --> RIGHT-HAND SIDE VECTOR
  ! NOTE
  ! ----
  ! if B IS NOT REQUIRED BY callING program, then
  ! B AND X MAY SHARE THE SAME STORAGE.
  integer, intent(in)     :: n
  real(dp),intent(in)   :: a(:,:)
  real(dp),intent(out)  :: x(:)
  real(dp),intent(in)   :: b(:)
  real(dp)  :: wk(n)
  ! FORWARD SOLVE, result IN WK
  call forslv(n, a, wk, b)
  ! BACK SOLVE, result IN X
  call bakslv(n, a, x, wk)
  return
end subroutine lltslv

subroutine optchk( &
& n, x, typsiz, sx, ndigit, epsm,  & !, itnlim, fscale
!& dlt, method, iexp, iagflg, iahflg, stepmx, msg)
& stepmx) !dlt, method, , msg, iexp, iagflg, iahflg
  ! PURPOSE
  ! -------
  ! CHECK INPUT FOR REASONABlenESS
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! X(N)         --> ON entry, ESTIMATE TO ROOT OF FCN
  ! TYPSIZ(N)   <--> TYPICAL size OF EACH COMPONENT OF X
  ! SX(N)       <--  DIAGONAL SCALING MATRIX FOR X
  ! FSCALE      <--> ESTIMATE OF SCALE OF OBJECTIVE function FCN
  ! ITNLIM      <--> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  ! NDIGIT      <--> NUMBER OF GOOD DIGITS IN OPTIMIZATION function FCN
  ! EPSM         --> MACHINE EPSILON
  ! DLT         <--> TRUST REGION RADIUS
  ! METHOD      <--> ALGORITHM INDICATOR
  ! IEXP        <--> EXPENSE FLAG
  ! IAGFLG      <-->  = 1 if ANALYTIC GRADIENT SUPPLIED
  ! IAHFLG      <-->  = 1 if ANALYTIC HESSIAN SUPPLIED
  ! STEPMX      <--> MAXIMUM STEP size
  ! MSG         <--> MESSAGE AND ERROR CODE
  integer, intent(in)   :: n
  real(dp),intent(in)   :: x(:)
  real(dp),intent(inout):: typsiz(:)
  real(dp),intent(out)  :: sx(:)
  !!real(dp),intent(inout):: fscale
  !!integer, intent(inout):: itnlim
  integer, intent(inout):: ndigit
  real(dp),intent(in)   :: epsm
  !!real(dp),intent(inout):: dlt
  !!integer, intent(inout):: method
  !!integer, intent(inout):: iexp
  !!integer, intent(inout):: iagflg
  !!integer, intent(inout):: iahflg
  real(dp),intent(inout)  :: stepmx
  !!integer, intent(inout):: msg
  !
  real(dp) :: stpsiz
  integer   :: i
  ! COMPUTE SCALE MATRIX
  do i = 1,n
    if(typsiz(i) == zero) typsiz(i) = one
    if(typsiz(i) < zero) typsiz(i) = -typsiz(i)
    sx(i) = one/typsiz(i)
  end do
  ! CHECK MAXIMUM STEP size
  stpsiz = zero
  do i = 1, n
    stpsiz = stpsiz + x(i)*x(i)*sx(i)*sx(i)
  end do
  stpsiz =SQRT(stpsiz)
  stepmx = MAX(1.0D3*stpsiz, 1.0D3)
  ! CHECK NUMBER OF DIGITS OF ACCURACY IN function FCN
  ndigit = -LOG10(epsm)
  return
end subroutine optchk

subroutine optdrv( &
& nr, n, x, fcn, d1fcn, d2fcn, typsiz, fscale,  &
& method, iexp, msg, ndigit, itnlim, iagflg, iahflg, &
& dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd)
  ! PURPOSE
  ! -------
  ! DRIVER FOR NON-LINEAR OPTIMIZATION PROBLEM
  ! parameterS
  ! ----------
  ! NR           --> ROW dimension OF MATRIX
  ! N            --> dimension OF PROBLEM
  ! X(N)         --> ON entry: ESTIMATE TO A ROOT OF FCN
  ! FCN          --> NAME OF subroutine TO EVALUATE OPTIMIZATION function
  !                            FCN: R(N) --> R(1)
  ! D1FCN        --> NAME OF subroutine TO EVALUATE GRADIENT OF FCN.
  ! D2FCN        --> NAME OF subroutine TO EVALUATE HESSIAN OF OF FCN.
  ! TYPSIZ(N)    --> TYPICAL size FOR EACH COMPONENT OF X
  ! FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE function
  ! METHOD       --> ALGORITHM TO use TO SOLVE MINIMIZATION PROBLEM
  !                    =1 LINE SEARCH
  !                    =2 doUBLE doGLEG
  !                    =3 MORE-HEBdoN
  ! IEXP         --> =1 if OPTIMIZATION function FCN IS EXPENSIVE TO
  !                  EVALUATE, =0 OTHERWISE.  if SET then HESSIAN WILL
  !                  BE EVALUATED BY SECANT UPDATE INSTEAD OF
  !                  ANALYTIcallY OR BY FINITE DifFERENCES
  ! MSG         <--> ON INPUT:  ( > 0) MESSAGE TO INHIBIT CERTAIN
  !                    AUTOMATIC CHECKS
  !                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR
  ! NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION function FCN
  ! ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  ! IAGFLG       --> =1 if ANALYTIC GRADIENT SUPPLIED
  ! IAHFLG       --> =1 if ANALYTIC HESSIAN SUPPLIED
  ! DLT          --> TRUST REGION RADIUS
  ! GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED close
  !                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! XPLS(N)     <--> ON exit:  XPLS IS LOCAL MINIMUM
  ! FPLS        <--> ON exit:  function value AT SOLUTION, XPLS
  ! GPLS(N)     <--> ON exit:  GRADIENT AT SOLUTION XPLS
  ! ITRMCD      <--  TERMINATION CODE
  ! INTERNAL VARIABLES
  ! ------------------
  ! ANALTL           TOLERANCE FOR COMPARISON OF ESTIMATED AND
  !                  ANALYTICAL GRADIENTS AND HESSIANS
  ! EPSM             MACHINE EPSILON
  ! F                function value: FCN(X)
  ! ITNCNT           CURRENT ITERATION, K
  ! RNF              RELATIVE NOISE IN OPTIMIZATION function FCN.
  !                       NOISE = 10.**(-NDIGIT)
  integer, intent(in)   :: nr
  integer, intent(in)   :: n
  real(dp),intent(inout):: x(:)
  real(dp),intent(inout):: typsiz(:)
  real(dp),intent(inout):: fscale
  integer, intent(inout):: method
  integer, intent(inout):: iexp
  integer, intent(inout):: msg
  integer, intent(inout):: ndigit
  integer, intent(inout):: itnlim
  integer, intent(inout):: iagflg
  integer, intent(inout):: iahflg
  real(dp),intent(inout):: dlt
  real(dp),intent(in)   :: gradtl
  real(dp),intent(inout):: stepmx
  real(dp),intent(in)   :: steptl
  real(dp),intent(inout):: xpls(:)
  real(dp),intent(inout):: fpls
  real(dp),intent(inout):: gpls(:)
  integer, intent(out)  :: itrmcd
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in) :: n
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: f
    end subroutine fcn
    subroutine d1fcn(n, x, g)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in) :: n
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: g(:)
    end subroutine d1fcn
    subroutine d2fcn(nr, n, x, h)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in) :: nr
      integer, intent(in) :: n
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: h(:,:)
    end subroutine d2fcn
  end interface
  ! Workspace
  real(dp) :: a(n,n), udiag(n), g(n), p(n), sx(n), wrk1(n)
  integer :: i, itncnt, iretcd, icscmx
  real(dp):: epsm, f, rnf, dltsav
  real(dp):: amusav, amu, dlpsav, dltp, phisav, phi, phpsav, phip0
  real(dp):: work2(n,1)
  logical :: mxtake
  ! INITIALIZATION
  ! --------------
  p(1:n) = zero
  itncnt = 0
  iretcd = -1
  epsm = EPSILON( 1.D0 )
  
  call optchk( &
  & n, x, typsiz, sx, ndigit, epsm,  & !, fscale, itnlim
  & stepmx) !,dlt, method,  iexp, iagflg, iahflg,, msg
  
  if(msg < 0) return
  
  rnf = MAX(10.0D0**(-ndigit), epsm)
  
  ! EVALUATE FCN(X)
  call fcn(n,x,f)
  
  ! EVALUATE ANALYTIC OR FINITE DifFERENCE GRADIENT AND CHECK ANALYTIC
  ! GRADIENT, if REQUESTED.
  if (iagflg == 1) then
    call d1fcn (n, x, g)
  else
    call fstofd (1, n, x, fcn, g, a, sx, rnf, 1)
  end if
  
  call optstp( &
  & n, x, f, g, wrk1, itncnt, icscmx, itrmcd, gradtl, steptl, sx, &
  & fscale, itnlim, iretcd, mxtake)
  
  if(itrmcd /= 0) GO TO 700
  if(iexp /= 1) GO TO 80
  
  ! if OPTIMIZATION function EXPENSIVE TO EVALUATE (IEXP=1), then
  ! HESSIAN WILL BE OBTAINED BY SECANT UPDATES.  GET INITIAL HESSIAN.
  call hsnint(n, a, sx, method)
  GO TO 100
  
  ! EVALUATE ANALYTIC OR FINITE DifFERENCE HESSIAN AND CHECK ANALYTIC
  ! HESSIAN if REQUESTED (only if useR-SUPPLIED ANALYTIC HESSIAN
  ! ROUTINE D2FCN FILLS only LOWER TRIANGULAR PART AND DIAGONAL OF A).
  80 if (iahflg == 0) then
    if (iagflg == 1) call estimate_Hessian (n, n, x, d1fcn, g, a, sx, rnf, 3)
    if (iagflg /= 1) call sndofd (n, x, fcn, f, a, sx, rnf)
  else
    call d2fcn(nr, n, x, a)              ! Inserted by AJM
  end if
  
  ! ITERATION
  ! ---------
  100 itncnt = itncnt + 1
  ! FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION
  ! (SKIP THIS STEP if LINE SEARCH OR doGSTEP TECHNIQUES BEING useD WITH
  ! SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALreadY OBTAINED FROM SECFAC.)
  if(iexp == 1 .and. method /= 3) GO TO 105
  
  103 call chlhsn(n, a, epsm, sx, udiag)
  ! SOLVE FOR NEWTON STEP:  AP = -G
  105 wrk1(1:n) = -g(1:n)
  call lltslv(n, a, p, wrk1)
  ! DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS = X + P
  ! OR TO CHOOSE XPLS BY A GLOBAL STRATEGY.
  if (iagflg /= 0 .or. method == 1) GO TO 111
  dltsav = dlt
  if (method == 2) GO TO 111
  amusav = amu
  dlpsav = dltp
  phisav = phi
  phpsav = phip0
  
  111 if(method == 2) call dogdrv( &
  & n, x, f, g, a, p, xpls, fpls, fcn, sx, &
  & stepmx, steptl, dlt, iretcd, mxtake)
  
  if(method == 3) call hookdr( &
  & n, g, a, udiag, p, sx, stepmx, dlt, iretcd, &
  & amu, dltp, phi, phip0, epsm, itncnt)
  
  if (method /= 2) xpls(1:n) = x(1:n) + p(1:n) ! Added by AJM
  
  ! if COULD NOT FIND SATISFACTORY STEP AND FORWARD DifFERENCE
  ! GRADIENT WAS useD, RETRY USING CENTRAL DifFERENCE GRADIENT.
  if (iretcd == 1 .and. iagflg == 0) then
  ! SET IAGFLG FOR CENTRAL DifFERENCES
    iagflg = -1
    call fstocd (n, x, fcn, sx, rnf, g)
    if (method == 1) GO TO 105
    dlt = dltsav
    if (method == 2) GO TO 105
    amu = amusav
    dltp = dlpsav
    phi = phisav
    phip0 = phpsav
    GO TO 103
  end if
  
  ! CALCULATE STEP FOR OUTPUT
  p(1:n) = xpls(1:n) - x(1:n)
  ! CALCULATE GRADIENT AT XPLS
  if (iagflg == -1) GO TO 116
  if (iagflg == 0) GO TO 118
  ! ANALYTIC GRADIENT
  call d1fcn (n, xpls, gpls)
  GO TO 120
  ! CENTRAL DifFERENCE GRADIENT
  116 call fstocd (n, xpls, fcn, sx, rnf, gpls)
  GO TO 120
  ! FORWARD DifFERENCE GRADIENT
  118 call fstofd (1, n, xpls, fcn, wrk1, work2, sx, rnf, 1)
  gpls(1:n) = work2(1:n,1)
  ! CHECK WHETHER stopPING CRITERIA SATISFIED
  120 call optstp(n, xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, gradtl, &
                  steptl, sx, fscale, itnlim, iretcd, mxtake)
  if(itrmcd /= 0) GO TO 690
  ! EVALUATE HESSIAN AT XPLS
  if(iexp /= 0) GO TO 150
  if(iahflg == 1) GO TO 140
  if(iagflg == 1) call estimate_Hessian(n, n, xpls, d1fcn, gpls, a, sx, rnf, 3)
  if(iagflg /= 1) call sndofd(n, xpls, fcn, fpls, a, sx, rnf)
  GO TO 150
  140 call d2fcn(nr, n, xpls, a)
  ! X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS
  150 f = fpls
  do i = 1,n
    x(i) = xpls(i)
    g(i) = gpls(i)
  end do
  GO TO 100
  ! TERMINATION
  ! -----------
  ! RESET XPLS,FPLS,GPLS,  if PREVIOUS ITERATE SOLUTION
  690 if(itrmcd /= 3) GO TO 710
  700 fpls = f
  do i = 1,n
    xpls(i) = x(i)
    gpls(i) = g(i)
  end do
  710 msg = 0
  return
end subroutine optdrv

subroutine optif9( &
& nr, n, x, fcn, d1fcn, d2fcn, typsiz, fscale, method,  &
& iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt,  &
& gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd)
  ! PURPOSE
  ! -------
  ! PROVIDE COMPLETE interface TO MINIMIZATION PACKAGE.
  ! useR HAS FULL CONTROL OVER OPTIONS.
  ! parameterS
  ! ----------
  ! NR           --> ROW dimension OF MATRIX
  ! N            --> dimension OF PROBLEM
  ! X(N)         --> ON entry: ESTIMATE TO A ROOT OF FCN
  ! FCN          --> NAME OF subroutine TO EVALUATE OPTIMIZATION function
  !                            FCN: R(N) --> R(1)
  ! D1FCN        --> NAME OF subroutine TO EVALUATE GRADIENT OF FCN.
  ! D2FCN        --> NAME OF subroutine TO EVALUATE HESSIAN OF OF FCN.
  ! TYPSIZ(N)    --> TYPICAL size FOR EACH COMPONENT OF X
  ! FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE function
  ! METHOD       --> ALGORITHM TO use TO SOLVE MINIMIZATION PROBLEM
  !                     = 1 LINE SEARCH
  !                     = 2 doUBLE doGLEG
  !                     = 3 MORE-HEBdoN
  ! IEXP         -->  = 1 if OPTIMIZATION function FCN IS EXPENSIVE TO
  !                  EVALUATE,  = 0 OTHERWISE.  if SET then HESSIAN WILL
  !                  BE EVALUATED BY SECANT UPDATE INSTEAD OF
  !                  ANALYTIcallY OR BY FINITE DifFERENCES
  ! MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN
  !                    AUTOMATIC CHECKS
  !                  ON OUTPUT: (.LT.0) ERROR CODE;  = 0 NO ERROR
  ! NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION function FCN
  ! ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  ! IAGFLG       -->  = 1 if ANALYTIC GRADIENT SUPPLIED
  ! IAHFLG       -->  = 1 if ANALYTIC HESSIAN SUPPLIED
  ! DLT          --> TRUST REGION RADIUS
  ! GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED close
  !                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! XPLS(N)     <--> ON exit:  XPLS IS LOCAL MINIMUM
  ! FPLS        <--> ON exit:  function value AT SOLUTION, XPLS
  ! GPLS(N)     <--> ON exit:  GRADIENT AT SOLUTION XPLS
  ! ITRMCD      <--  TERMINATION CODE
  integer, intent(in)   :: nr
  integer, intent(in)   :: n
  real(dp),intent(inout):: x(:)
  real(dp),intent(inout):: typsiz(:)
  real(dp),intent(inout):: fscale
  integer, intent(inout):: method
  integer, intent(inout):: iexp
  integer, intent(inout):: msg
  integer, intent(inout):: ndigit
  integer, intent(inout):: itnlim
  integer, intent(inout):: iagflg
  integer, intent(inout):: iahflg
  real(dp),intent(inout):: dlt
  real(dp),intent(in)   :: gradtl
  real(dp),intent(inout):: stepmx
  real(dp),intent(in)   :: steptl
  real(dp),intent(inout):: xpls(:)
  real(dp),intent(inout):: fpls
  real(dp),intent(inout):: gpls(:)
  integer, intent(out)  :: itrmcd
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)  :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f
    end subroutine fcn
    subroutine d1fcn(n, x, g)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)  :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: g(:)
    end subroutine d1fcn
    subroutine d2fcn(nr, n, x, h)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)  :: nr
      integer, intent(in)  :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: h(:,:)
    end subroutine d2fcn
  end interface
  call optdrv(nr, n, x, fcn, d1fcn, d2fcn, typsiz, fscale,  &
              method, iexp, msg, ndigit, itnlim, iagflg, iahflg, &
              dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd)
  return
end subroutine optif9

subroutine optstp( &
& n, xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, gradtl, &
& steptl, sx, fscale, itnlim, iretcd, mxtake)
  ! UNCONSTRAINED MINIMIZATION stopPING CRITERIA
  ! --------------------------------------------
  ! FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY OF THE FOLLOWING:
  ! 1) PROBLEM SOLVED WITHIN useR TOLERANCE
  ! 2) CONVERGENCE WITHIN useR TOLERANCE
  ! 3) ITERATION LIMIT REACHED
  ! 4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! XPLS(N)      --> NEW ITERATE X[K]
  ! FPLS         --> function value AT NEW ITERATE F(XPLS)
  ! GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE
  ! X(N)         --> OLD ITERATE X[K-1]
  ! ITNCNT       --> CURRENT ITERATION K
  ! ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX
  !                  [RETAIN value BETWEEN SUCCESSIVE callS]
  ! ITRMCD      <--  TERMINATION CODE
  ! GRADTL       --> TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED close
  !                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE function
  ! ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  ! IRETCD       --> return CODE
  ! MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  integer, intent(in)   :: n
  real(dp),intent(in)   :: xpls(:)
  real(dp),intent(in)   :: fpls
  real(dp),intent(in)   :: gpls(:)
  real(dp),intent(in)   :: x(:)
  integer, intent(in)   :: itncnt
  integer, intent(inout):: icscmx
  integer, intent(out)  :: itrmcd
  real(dp),intent(in)   :: gradtl
  real(dp),intent(in)   :: steptl
  real(dp),intent(in)   :: sx(:)
  real(dp),intent(in)   :: fscale
  integer, intent(in)   :: itnlim
  integer, intent(in)   :: iretcd
  logical, intent(in)   :: mxtake
  integer   :: i, jtrmcd
  real(dp) :: d, rgx
  real(dp) :: relgrd, relstp, rsx
  itrmcd = 0
  ! LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X
  if(iretcd == 1) then
    jtrmcd = 3
    GO TO 600
  end if
  ! FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM.
  ! CHECK WHETHER WITHIN TOLERANCE
  d = MAX(ABS(fpls),fscale)
  rgx = zero
  do i = 1,n
    relgrd = ABS(gpls(i)) * MAX(ABS(xpls(i)), one/sx(i))/d
    rgx = MAX(rgx,relgrd)
  end do
  jtrmcd = 1
  if(rgx <= gradtl) GO TO 600
  if(itncnt == 0) return
  ! FIND DIRECTION IN WHICH RELATIVE STEPsize MAXIMUM
  ! CHECK WHETHER WITHIN TOLERANCE.
  rsx = zero
  do i = 1,n
    relstp = ABS(xpls(i) - x(i)) / MAX(ABS(xpls(i)), one/sx(i))
    rsx = MAX(rsx,relstp)
  end do
  jtrmcd = 2
  if(rsx <= steptl) GO TO 600
  ! CHECK ITERATION LIMIT
  jtrmcd = 4
  if(itncnt >= itnlim) GO TO 600
  ! CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX
  if(.not. mxtake) then
    icscmx = 0
    return
  end if
  icscmx = icscmx + 1
  if(icscmx < 5) return
  jtrmcd = 5
  600 itrmcd = jtrmcd
  return
end subroutine optstp

subroutine sndofd(n, xpls, fcn, fpls, a, sx, rnoise)
  ! PURPOSE
  ! -------
  ! FIND SECOND ORDER FORWARD FINITE DifFERENCE APPROXIMATION "A"
  ! TO THE SECOND DERIVATIVE (HESSIAN) OF THE function DEFINED BY THE SUBP
  ! "FCN" EVALUATED AT THE NEW ITERATE "XPLS"
  ! FOR OPTIMIZATION use THIS ROUTINE TO ESTIMATE
  ! 1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION function
  !    if NO ANALYTICAL useR function HAS BEEN SUPPLIED FOR EITHER
  !    THE GRADIENT OR THE HESSIAN AND if THE OPTIMIZATION function
  !    "FCN" IS INEXPENSIVE TO EVALUATE.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! XPLS(N)      --> NEW ITERATE:   X[K]
  ! FCN          --> NAME OF subroutine TO EVALUATE function
  ! FPLS         --> function value AT NEW ITERATE, F(XPLS)
  ! A(N,N)      <--  FINITE DifFERENCE APPROXIMATION TO HESSIAN
  !                  only LOWER TRIANGULAR MATRIX AND DIAGONAL ARE returnED
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! RNOISE       --> RELATIVE NOISE IN FNAME [F(X)]
  integer, intent(in)      :: n
  real(dp),intent(inout)   :: xpls(:)
  real(dp),intent(in)      :: fpls
  real(dp),intent(out)     :: a(:,:)
  real(dp),intent(in)      :: sx(:)
  real(dp),intent(in)      :: rnoise
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)    :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f
    end subroutine fcn
  end interface
  ! Workspace
  ! STEPSZ(N)    --> WORKSPACE (STEPsize IN I-TH COMPONENT DIRECTION)
  ! ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION)
  real(dp) :: stepsz(n), anbr(n)
  integer   :: i,j,ip1
  real(dp) :: ov3,xtmpi,xtmpj,fhat
  ! FIND I-TH STEPsize AND EVALUATE NEIGHBOR IN DIRECTION OF I-TH UNIT VECTOR.
  ov3 = one/3.0_dp
  do i = 1,n
    stepsz(i) = rnoise**ov3 * MAX(ABS(xpls(i)), one/sx(i))
    xtmpi = xpls(i)
    xpls(i) = xtmpi + stepsz(i)
    call fcn(n,xpls,anbr(i))
    xpls(i) = xtmpi
  end do
  ! CALCULATE COLUMN I OF A
  do i = 1,n
    xtmpi = xpls(i)
    xpls(i) = xtmpi + 2.0*stepsz(i)
    call fcn(n,xpls,fhat)
    a(i,i) = ((fpls - anbr(i)) + (fhat-anbr(i))) / (stepsz(i)*stepsz(i))
  ! CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN
    if(i == n) GO TO 25
    xpls(i) = xtmpi + stepsz(i)
    ip1 = i + 1
    do j = ip1,n
      xtmpj = xpls(j)
      xpls(j) = xtmpj + stepsz(j)
      call fcn(n,xpls,fhat)
      a(j,i) = ((fpls - anbr(i)) + (fhat-anbr(j))) / (stepsz(i)*stepsz(j))
      xpls(j) = xtmpj
    end do
    25 xpls(i) = xtmpi
  end do
  return
end subroutine sndofd

subroutine tregup( &
& n, x, f, g, a, fcn, sc, sx, nwtake, stepmx, steptl, dlt,  &
& iretcd, xplsp, fplsp, xpls, fpls, mxtake, method, udiag)
  ! PURPOSE
  ! -------
  ! DECIDE WHETHER TO ACCEPT XPLS = X+SC AS THE NEXT ITERATE AND UPDATE THE
  ! TRUST REGION DLT.
  ! parameterS
  ! ----------
  ! N            --> dimension OF PROBLEM
  ! X(N)         --> OLD ITERATE X[K-1]
  ! F            --> function value AT OLD ITERATE, F(X)
  ! G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
  ! A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
  !                  LOWER TRIANGULAR PART AND DIAGONAL.
  !                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
  ! FCN          --> NAME OF subroutine TO EVALUATE function
  ! SC(N)        --> CURRENT STEP
  ! SX(N)        --> DIAGONAL SCALING MATRIX FOR X
  ! NWTAKE       --> BOOLEAN, =.true. if NEWTON STEP TAKEN
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! DLT         <--> TRUST REGION RADIUS
  ! IRETCD      <--> return CODE
  !                     = 0 XPLS ACCEPTED AS NEXT ITERATE;
  !                       DLT TRUST REGION FOR NEXT ITERATION.
  !                     = 1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE
  !                       BECAuse XPLS-X .LT. SMALLEST ALLOWABLE STEP lenGTH.
  !                     = 2 F(XPLS) TOO LARGE.  continue CURRENT ITERATION
  !                       WITH NEW REDUCED DLT.
  !                     = 3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL
  !                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO continue
  !                       CURRENT ITERATION WITH NEW doUBLED DLT.
  ! XPLSP(N)    <--> WORKSPACE [value NEEDS TO BE RETAINED BETWEEN
  !                  SUCCESSIVE callS OF K-TH GLOBAL STEP]
  ! FPLSP       <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! XPLS(N)     <--  NEW ITERATE X[K]
  ! FPLS        <--  function value AT NEW ITERATE, F(XPLS)
  ! MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  ! METHOD       --> ALGORITHM TO use TO SOLVE MINIMIZATION PROBLEM
  !                     = 1 LINE SEARCH
  !                     = 2 doUBLE doGLEG
  !                     = 3 MORE-HEBdoN
  ! UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
  integer, intent(in):: n
  real(dp),intent(in):: x(:)
  real(dp),intent(in):: f
  real(dp),intent(in):: g(:)
  real(dp),intent(in):: a(:,:)
  real(dp),intent(in):: sc(:)
  real(dp),intent(in):: sx(:)
  logical, intent(in) :: nwtake
  real(dp),intent(in):: stepmx
  real(dp),intent(in):: steptl
  real(dp),intent(inout):: dlt
  integer, intent(inout):: iretcd
  real(dp),intent(inout):: xplsp(:)
  real(dp),intent(inout):: fplsp
  real(dp),intent(out):: xpls(:)
  real(dp),intent(out):: fpls
  logical, intent(out):: mxtake
  integer, intent(in) :: method
  real(dp),intent(in) :: udiag(:)
  interface
    subroutine fcn(n, x, f)
      use M_Kinds,only:dp
      implicit none
      integer, intent(in)    :: n
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f
    end subroutine fcn
  end interface
  integer   :: i, j, ip1
  real(dp) :: slp, rln
  real(dp) :: dltmp, dltfp, temp, dltf
  mxtake = .false.
  xpls(1:n) = x(1:n) + sc(1:n)
  call fcn(n, xpls, fpls)
  dltf = fpls - f
  slp = dot_product( g(1:n), sc(1:n) )
  ! NEXT STATEMENT ADDED FOR case OF COMPILERS WHICH do NOT OPTIMIZE
  ! EVALUATION OF NEXT "if" STATEMENT (IN WHICH case FPLSP COULD BE UNDEFINED).
  if(iretcd == 4) fplsp = zero
  if(iretcd /= 3 .or. (fpls < fplsp .and. dltf <= 1.e-4*slp)) GO TO 130
  !     if(IRETCD.EQ.3 .and. (FPLS.GE.FPLSP .or. DLTF.GT. 1.E-4*SLP))
  !     then
  !       RESET XPLS TO XPLSP AND TERMINATE GLOBAL STEP
  iretcd = 0
  xpls(1:n) = xplsp(1:n)
  fpls = fplsp
  dlt = .5*dlt
  GO TO 230
  !     else
  !       FPLS TOO LARGE
  130 if(dltf <= 1.e-4*slp) GO TO 170
  !       if(DLTF.GT. 1.E-4*SLP)
  !       then
  rln = zero
  do i = 1,n
    rln = MAX(rln,ABS(sc(i)) / MAX(ABS(xpls(i)),1./sx(i)))
  end do
  if(rln >= steptl) GO TO 150
  !         if(RLN.LT.STEPTL)
  !         then
  !           CANNOT FIND SATISFACTORY XPLS SUFFICIENTLY DISTINCT FROM X
  iretcd = 1
  GO TO 230
  !         else
  !           REDUCE TRUST REGION AND continue GLOBAL STEP
  150 iretcd = 2
  dltmp = -slp*dlt/(2.*(dltf-slp))
  if(dltmp >= .1*dlt) GO TO 155
  !           if(DLTMP.LT. .1*DLT)
  !           then
  dlt = .1*dlt
  GO TO 160
  !           else
  155 dlt = dltmp
  !           end if
  160 GO TO 230
  !         end if
  !       else
  !         FPLS SUFFICIENTLY SMALL
  170 dltfp = zero
  if (method == 3) GO TO 180
  !         if (METHOD .EQ. 2)
  !         then
  do i = 1, n
    temp = dot_product( a(i:n,i), sc(i:n) )
    dltfp = dltfp + temp*temp
  end do
  GO TO 190
  !         else
  180 do i = 1, n
    dltfp = dltfp + udiag(i)*sc(i)*sc(i)
    if (i == n) cycle
    temp = 0
    ip1 = i + 1
    do j = ip1, n
      temp = temp + a(i, j)*sc(i)*sc(j)
    end do
    dltfp = dltfp + 2.0*temp
  end do
  !         end if
  190 dltfp = slp + dltfp/2.0
  if(iretcd == 2 .or. (ABS(dltfp-dltf) > .1*ABS(dltf))  &
      .or. nwtake .or. (dlt > .99*stepmx)) GO TO 210
  !         if(IRETCD.NE.2 .and. (ABS(DLTFP-DLTF) .LE. .1*ABS(DLTF))
  !    +         .and. (.not.NWTAKE) .and. (DLT.LE. .99*STEPMX))
  !         then
  !           doUBLE TRUST REGION AND continue GLOBAL STEP
  iretcd = 3
  xplsp(1:n) = xpls(1:n)
  fplsp = fpls
  dlt = MIN(2.0D0*dlt,stepmx)
  GO TO 230
  !         else
  !           ACCEPT XPLS AS NEXT ITERATE.  CHOOSE NEW TRUST REGION.
  210 iretcd = 0
  if(dlt > .99*stepmx) mxtake = .true.
  if(dltf < .1*dltfp) GO TO 220
  !           if(DLTF.GE. .1*DLTFP)
  !           then
  !             DECREASE TRUST REGION FOR NEXT ITERATION
  dlt = .5*dlt
  GO TO 230
  !           else
  !             CHECK WHETHER TO INCREASE TRUST REGION FOR NEXT ITERATION
  220 if(dltf <= .75*dltfp) dlt = MIN(2.*dlt,stepmx)
  !           end if
  !         end if
  !       end if
  !     end if
  230 return
end subroutine tregup

end module M_TenSolve_UncMin
