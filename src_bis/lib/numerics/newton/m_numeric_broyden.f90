module M_Numeric_Broyden
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_
  implicit none

  private
  !
  public:: Broyden
  !
  logical,public:: TestJacob=.false., ShoResult=.false.
  !
contains

subroutine Broyden( &
& vX,      & !inout= the initial guess, and the root returned
& Residual, Jacobian, &
& TolF,    & !in=    convergence criterion on function values 
& TolX,    & !in=    convergence criterion on dx
!!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
& bFinDif, & !in=    use numeric Jacobian
& MaxIts,  & !in=    maximum number of iterations
& Error_F, & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Gradient, & !out=   
& Nits,    & !out=   number of iterations
& Check,   & !out=   if Check, should check convergence
& iErr)      !out=   error code
!.NumericalRecipes, Broyden’s method embedded in a globally convergent strategy. 
!The length N vector of functions to be zeroed, called fvec in the routine below,
!is returned by a user-supplied routine that must be called funcv and have the declaration function funcv(x).
!use subroutine Jacobian_Numeric and function fmin from newt
!
!check=
!-- false on a normal return,
!-- true if the routine has converged to a local minimum of the function fmin
!-- or if Broyden’s method can make no further progress. 
!In this case try restarting from a different initial guess.
!
  use M_Numeric_Const,      only: Ln10
  use M_Numeric_Tools,only: OuterProd_R,VAbs,Unit_Matrix,Jacobian_Numeric
  use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
  use M_IoTools,    only: OutStrVec
  !use M_NRBroyden,only: &
  !& Get_Diag_Rv,Put_Diag_Rv,Lower_Triangle, &
  !& QRDCMP,QRUPDT,RSOLV
  !!use M_IOTools,only: OutStrVec !for debugging
  !
  real(dp),intent(inout):: vX(:)
  real(dp),intent(in)   :: TolF,TolX !!,TolMin
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X,Gradient
  logical, intent(out)  :: Check
  integer, intent(out)  :: nIts,iErr
  interface
    function Residual(v)
      use M_Kinds
      implicit none
      real(dp),dimension(:),intent(in):: v
      real(dp),dimension(size(v))    :: Residual
    end function Residual
    subroutine Jacobian(v,t)
      use M_Kinds
      implicit none
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
  endinterface
  !
  real(dp),parameter:: STPMX=  100._dp, TolMin= 1.0E-9_dp
  integer  :: I,ITS,K,N
  real(dp) :: F,FOLD,STPMAX,EPS
  !
  real(dp),dimension(size(vX)):: vFunc,vFuncOld,C,D,G,P,S,T,W,XOLD
  real(dp),dimension(size(vX),size(vX)):: QT,R
  logical  :: RESTRT,SING
  !
  EPS= TolX
  N=size(vX)
  vFunc= Residual(vX)
  F=    dot_product(vFunc,vFunc)/2.0D0
  !
  !if (MAXVAL(ABS(vFunc(:))) < 0.01_dp*NewtTolF) then
  !  iErr=1; return !__________________________________________________________return
  !end if
  !
  STPMAX= STPMX*MAX(VABS(vX(:)),real(N,dp))
  RESTRT= .true.
  !
  iErr= -1
  do0: do ITS=1,MAXITS
    nIts=Its
    !
  fNewt_I= fNewt_I +1
  if(fNewtF>0) call OutStrVec(fNewtF,vX(:)/Ln10,Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  if(fNewtR>0) call OutStrVec(fNewtR,vFunc(:),  Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
    !if(DebNewt .and. bOpenNewt) then
    !  call OutStrVec(fNewt1,vX(1:size(vX))/Ln10,Opt_I=Its,Opt_C="F")
    !  call OutStrVec(fNewt2,vFunc(1:size(vX)),  Opt_I=Its,Opt_C="F")
    !end if
    !
    if (RESTRT) then
      !call Jacobian_Numeric(vX,vFunc,R)
      if(bFinDif) then
        call Jacobian_Numeric(Residual,vX,vFunc,R)
      else
        call Jacobian(vX,R)
      end if
      !
      call QRDCMP(R,C,D,SING)
      !
      if (SING) then !call Stop_('SINGULAR JACOBIAN IN Residual')
        iErr=-2; exit do0
      end if
      !
      call UNIT_MATRIX(QT)
      !
      do K=1,N-1
        if(C(K)/=Zero) &
        & QT(K:N,:)= QT(K:N,:) &
        &          - OUTERPROD_R(R(K:N,K),matmul(R(K:N,K),QT(K:N,:)))/C(K)
      end do
      where (LOWER_TRIANGLE(N,N)) R(:,:)=Zero
      call Put_Diag_Rv(D(:),R(:,:))
    else
      S(:)=vX(:)-XOLD(:)
      do I=1,N; T(I)=dot_product(R(I,I:N),S(I:N)); end do
      !
      W(:)=vFunc(:)-vFuncOld(:) -matmul(T(:),QT(:,:))
      !
      where (ABS(W(:)) < EPS*(ABS(vFunc(:))+ABS(vFuncOld(:)))) W(:)=Zero
      if (ANY(W(:) /= Zero)) then
        T(:)=matmul(QT(:,:),W(:))
        S(:)=S(:)/dot_product(S,S)
        call QRUPDT(R,QT,T,S)
        D(:)=GET_DIAG_RV(R(:,:)) !do J=1,size(MAT,1) GET_DIAG_RV(J)=MAT(J,J); end do
        if (ANY(D(:)==Zero)) then !call Stop_('R SINGULAR IN Residual')
          iErr=-2; exit do0
        end if
      end if
    end if
    P(:)=-matmul(QT(:,:),vFunc(:))
    do I=1,N; G(I)=-dot_product(R(1:I,I),P(1:I)); end do
    XOLD(:)=vX(:)
    vFuncOld(:)=vFunc(:)
    FOLD=F
    call RSOLV(R,D,P)
    !
    call LNSRCH(XOLD,FOLD,G,STPMAX,P,vX,F,CHECK) !-> this updates also vFunc !!
    !
    Error_F= MAXVAL(ABS(vFunc(:)))
    Gradient= MAXVAL(ABS(G(:))*MAX(ABS(vX(:)),One)/MAX(F,0.5_dp*N))
    Delta_X= MAXVAL((ABS(vX(:)-XOLD(:)))/MAX(ABS(vX(:)),One))
    !
    if (Error_F < TolF) then
      iErr=0; return !__________________________________________________________return
    end if
    if (CHECK) then
      if(RESTRT .or. Gradient < TolMin) then
        iErr=-3; return !________________________________________________________return
      end if
      RESTRT=.true.
    else
      RESTRT=.false.
      if(Delta_X < TolX) then
        iErr=-4; return !________________________________________________________return
      end if
    end if
  end do do0
  return

contains

subroutine LNSRCH( & !from NR
& vXOld, & !IN= an N-dimensional point
& rFOld, & !IN= the value of the function at point vXOld
& vGrad, & !IN= the gradient of the function at that point
& rStpMax, & !IN=limits the length of the steps
& vDX,   & !INOUT= a direction
& vXnew, & !OUT= a new point vXnew(1:N) along the direction vDX from vXOld
& rFnew, & !OUT= the new function value
& CHECK)   !FALSE= normal exit, TRUE= vXnew is too close to vXOld 
!Given 
!  vXOld(1:N),IN= an N-dimensional point, 
!  rFOld(1:N),IN= the value of the function at vXOld,
!  vGrad(1:N),IN= the gradient of the function at vXOld,
!  vDX(1:N),INOUT= a direction vDX
!finds
!  a new point vXnew(1:N) along the direction vDX from vXOld
!  where the function FUNC has decreased "sufficiently".
!rFnew(1),OUT= the new function value
!rStpMax,IN=
!  an input quantity that limits the length of the steps 
!  so that you do not try to evaluate the function in regions where it is undefined 
!  or subject to overflow. 
!vDX,INOUT= usually the Newton direction.
!Check,OUT=  
!  FALSE on a normal exit, 
!  TRUE when vXnew is too close to vXOld,
!  in a minimization algorithm, this usually signals convergence and can be ignored
!  in a zero-finding algorithm, the calling program should check whether the convergence is spurious.
!Parameters
!  ALF=  ensures sufficient decrease in function value;
!  TolX= the convergence criterion on delta(vXnew)
  !
  real(dp),dimension(:), intent(in)    :: vXOld,vGrad
  real(dp),dimension(:), intent(inout) :: vDX
  real(dp),              intent(in)    :: rFOld,rStpMax
  real(dp),dimension(:), intent(out)   :: vXnew
  real(dp),              intent(out)   :: rFnew
  logical,               intent(out)   :: CHECK
  !
  real(dp),parameter ::&
  & TolX=EPSILON(vXnew), &
  & ALF= EPSILON(vXnew) !1.0E-4_dp  ! !1.0E-6_dp
  real(dp)::&
  & A,Alam,Alam2,AlaMin,B,DISC,F2,PAbs,&
  & RHS1,RHS2,Slope,tmpLam,x
  !
  CHECK=.false.
  PAbs=vAbs(vDX(:))
  if (PAbs>rStpMax) vDX(:)=vDX(:)*rStpMax/PAbs !Scale if attempted step is too big.
  !
  Slope=dot_product(vGrad,vDX) !______________Compute lambda_min
  if (slope>=Zero) call Stop_('roundoff problem in lnsrch') !____________________________stop
  AlaMin=TolX / MAXVAL(ABS(vDX(:))/MAX(ABS(vXOld(:)),One))
  Alam=One !____________________________Always try full Newton step first
  do !__________________________________Start of iteration loop
    vXnew(:)= vXOld(:) +Alam *vDX(:)
    !rFnew=FUNC(vXnew) !rFnew=scalar,vXnew=Vector
    vFunc= Residual(vXnew)
    rFnew= dot_product(vFunc,vFunc)/2.0D0
    if (Alam<AlaMin) then !_____________Convergence on vXnew
      vXnew(:)=vXOld(:) !For zero finding, the calling program should verify the convergence.
      CHECK=.true.
      return !__________________________________________________return!!!
    elseif (rFnew <= rFOld+ALF*Alam*Slope) then
      return !sufficient function decrease______________________return!!!
    else                        !Backtrack.
      if (Alam==One) then; tmpLam=-Slope/(2.0_dp*(rFnew-rFOld-Slope)) !First time.
      else !____________________________Subsequent backtracks
        RHS1=rFnew -rFOld -Alam *Slope
        RHS2=F2-rFOld -Alam2*Slope
        A= (       RHS1/Alam/Alam -      RHS2/Alam2/Alam2) / (Alam-Alam2)
        B= (-Alam2*RHS1/Alam/Alam + Alam*RHS2/Alam2/Alam2) / (Alam-Alam2)
        if (A==Zero) then; tmpLam=-Slope/(2.0_dp*B)
        else
          DISC=B*B-3.0_dp*A*Slope
          if (DISC<Zero)   then; tmpLam=0.5_dp*Alam
          elseif (B<=Zero) then; tmpLam=(-B+SQRT(DISC))/(3.0_dp*A)
          else                 ; tmpLam=-Slope/(B+SQRT(DISC))
          end if
        end if
        if (tmpLam>0.5_dp*Alam) tmpLam=0.5_dp*Alam
      end if
    end if
    Alam2=Alam
    Alam=MAX(tmpLam,0.1_dp*Alam)
    F2=rFnew
  end do
end subroutine LNSRCH

end subroutine Broyden

!!module M_NRBroyden
function Get_Diag_Rv(MAT) !->vector of diagonal elements, used in Broyden method
  use M_Trace,only:Stop_
  real(dp),dimension(:,:),intent(in)::MAT
  real(dp),dimension(size(MAT,1))::Get_Diag_Rv
  !
  integer::J
  !
  if(size(MAT,1)/=size(MAT,2)) call Stop_("size !!! in Get_Diag_Rv") !check that MAT is Square Matrix
  do J=1,size(MAT,1); Get_Diag_Rv(J)=MAT(J,J); end do
end function Get_Diag_Rv

subroutine Put_Diag_Rv(DIAGV,MAT) !used in Broyden method
  use M_Trace,only:Stop_
  real(dp),dimension(:),  intent(in)   :: DIAGV
  real(dp),dimension(:,:),intent(inout):: MAT
  integer::J,N
  !
  if(size(DIAGV) /= MIN(size(MAT,1),size(MAT,2))) call Stop_("size !!! PUT_DIAG_RV")
  N=size(DIAGV)
  do J=1,N; MAT(J,J)=DIAGV(J); end do
end subroutine Put_Diag_Rv

function ARTH_I(First,Increment,N) !used by Broyden
  integer,intent(in)  :: First,Increment,N
  integer,dimension(N):: ARTH_I
  integer             :: K,K2,TEMP
  if(N>0) ARTH_I(1)=FIRST
  if(N<=16) then !(N<=NPAR_ARTH)
    do K=2,N; Arth_I(K)=Arth_I(K-1)+Increment; end do
  else
    do K=2,8; Arth_I(K)=Arth_I(K-1)+Increment; end do !K=2,NPAR2_ARTH
    Temp=Increment*8 !NPAR2_ARTH
    K=8 !NPAR2_ARTH
    do
      if (K >= N) exit
      K2=K+K
      Arth_I(K+1:MIN(K2,N))=Temp+Arth_I(1:MIN(K,N-K))
      Temp=Temp+Temp
      K=K2
    end do
  end if
end function ARTH_I

function OuterDiff_R(A,B)
  real(dp),dimension(:), intent(in)  :: A,B
  real(dp),dimension(size(A),size(B)):: OuterDiff_R
  OuterDiff_R= SPread(A,DIM=2,NCOPIES=size(B)) &
  &          - SPread(B,DIM=1,NCOPIES=size(A))
end function OuterDiff_R

function Outerdiff_I(A,B)
  integer,dimension(:),intent(in) :: A,B
  integer,dimension(size(A),size(B)):: Outerdiff_I
  !
  OUTERDifF_I= SPread(A,DIM=2,NCOPIES=size(B)) &
  &          - SPread(B,DIM=1,NCOPIES=size(A))
end function Outerdiff_I

function Lower_Triangle(J,K,EXTRA) !matricial mask (J,K) of logical
  integer,          intent(in):: J,K
  integer,optional, intent(in):: EXTRA
  logical,dimension(J,K)      :: Lower_Triangle
  integer :: N
  N=0; if (present(EXTRA)) N=EXTRA
  Lower_Triangle=(Outerdiff_I(Arth_I(1,1,J),Arth_I(1,1,K)) > -N)
end function Lower_Triangle

function Assert_Eq4(N1,N2,N3,N4,Str) !"NR"
  use M_Trace,only:Stop_
  character(len=*), intent(in) :: Str
  integer, intent(in)          :: N1,N2,N3,N4
  integer                      :: Assert_Eq4
  if (N1==N2 .and. N2==N3 .and. N3==N4) then
    Assert_Eq4=N1
  else
    call Stop_(Str)
  end if
end function Assert_Eq4

function Assert_EqN(NN,Str) !"NR"
  use M_Trace,only:Stop_
  character(len=*),    intent(in):: Str
  integer,dimension(:),intent(in):: NN
  integer                        :: Assert_Eqn
  if (ALL(NN(2:)==NN(1))) then; Assert_EqN=NN(1)
  else;                         call Stop_(Str)
  end if
end function Assert_EqN

subroutine RSolv(a,d,b) !right triangular backsubstitution
  real(dp),dimension(:,:),intent(in)   :: a
  real(dp),dimension(:),  intent(in)   :: d
  real(dp),dimension(:),  intent(inout):: b
  integer:: i,n
  n=assert_eq4(size(a,1),size(a,2),size(b),size(d),'rsolv')
  b(n)=b(n)/d(n)
  do i=n-1,1,-1; b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i); end do
end subroutine RSolv

subroutine QRUpdt(r,qt,u,v) !from "NR" 
  use M_Numeric_Tools,only:ifirstloc
  real(dp),dimension(:,:),intent(inout):: r,qt
  real(dp),dimension(:),  intent(inout):: u
  real(dp),dimension(:),  intent(in)   :: v
  integer::i,k,n
  n=Assert_EqN((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),size(v)/),'qrupdt')
  k= n+1 -ifirstloc(u(n:1:-1) /= Zero)
  if (k < 1) k=1
  do i=k-1,1,-1
    call Rotate(r,qt,i,u(i),-u(i+1))
    u(i)=pythag(u(i),u(i+1))
  end do
  r(1,:)=r(1,:)+u(1)*v
  do i=1,k-1; call Rotate(r,qt,i,r(i,i),-r(i+1,i)); end do
end subroutine QRUpdt

function Pythag(a,b) !from "NR"
  use M_Kinds
  real(dp),intent(in):: a,b
  real(dp)           :: pythag
  real(dp):: absa,absb
  absa=ABS(a)
  absb=ABS(b)
  if (absa > absb) then
    pythag=absa*SQRT(One+(absb/absa)**2)
  else
    if (absb==Zero) then; pythag=Zero
    else                ; pythag=absb*SQRT(One+(absa/absb)**2)
    end if
  end if
end function Pythag

subroutine Rotate(r,qt,i,a,b) !from "NR"
  real(dp),dimension(:,:),target,intent(inout):: r,qt
  integer,                        intent(in)   :: i
  real(dp),                     intent(in)   :: a,b
  real(dp),dimension(size(r,1)):: temp
  integer:: n
  real:: c,fact,s
  n=assert_eq4(size(r,1),size(r,2),size(qt,1),size(qt,2),'rotate')
  if(a== Zero) then
    c=   Zero
    s=   SIGN(One,b)
  elseif (ABS(a)>ABS(b)) then
    fact=b/a
    c=   sign(One/SQRT(One+fact**2),a)
    s=   fact*c
  else
    fact=a/b
    s=   sign(One/SQRT(One+fact**2),b)
    c=   fact*s
  end if
  temp(i:n)= r(i,i:n)
  r(i,  i:n)= c *temp(i:n) -s *r(i+1,i:n)
  r(i+1,i:n)= s *temp(i:n) +c *r(i+1,i:n)
  temp=       qt(i,:)
  qt(i,  :)= c *temp -s *qt(i+1,:)
  qt(i+1,:)= s *temp +c *qt(i+1,:)
end subroutine Rotate

subroutine QRDcmp(a,c,d,sing) !QR decomposition, from "NR"
  use M_Numeric_Tools,only:OuterProd_R,VAbs
  real(dp),dimension(:,:),intent(inout):: a
  real(dp),dimension(:),  intent(out)  :: c,d
  logical,                 intent(out)  :: sing
  integer :: k,n
  real(dp):: scale,sigma
  n=Assert_Eq4(size(a,1),size(a,2),size(c),size(d),'qrdcmp')
  sing=.false.
  do k=1,n-1
    scale=MAXVAL(ABS(a(k:n,k)))
    if (scale==Zero) then
      sing= .true.
      c(k)= Zero
      d(k)= Zero
    else
      a(k:n,k)=a(k:n,k)/scale
      sigma=   SIGN(vabs(a(k:n,k)),a(k,k))
      a(k,k)=  a(k,k)+sigma
      c(k)=    sigma*a(k,k)
      d(k)=   -scale*sigma
      a(k:n,k+1:n)= a(k:n,k+1:n) &
      &           - OuterProd_R(a(k:n,k),matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
    end if
  end do
  d(n)=a(n,n)
  if (d(n)==Zero) sing=.true.
end subroutine QRDcmp
!!end module M_NRBroyden

end module M_Numeric_Broyden

