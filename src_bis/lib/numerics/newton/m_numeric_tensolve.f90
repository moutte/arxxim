module M_Numeric_Tensolve
  use M_Kinds
  use M_Trace,only: fTrc,Stop_
  implicit none

  private
  !
  public:: TenSolve
  !
contains

subroutine TenSolve( &
& sCod1,    &
& sCod2,    &
& vX,       & 
& Residual, &
& Jacobian, &
& TolF,     & !in=    convergence criterion on function values 
& bFinDif,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
! for TenSolve, iErr is TERMCD,
! i.e. TERMCD<-:  An integer specifying the reason for termination.
! TERMCD= 0 : No termination criterion satisfied
!   occurs if package terminates because of illegal input
! TERMCD= 1 : function tolerance reached.
!   The current iteration is probably a solution.
! TERMCD= 2 : gradient tolerance reached.
!   For nonlinear least squares, the current iteration is probably a solution
!   For nonlinear equations, it could be a solution or a local minimizer
! TERMCD= 3 : Successive iterates within step tolerance.
!   The current iterate may be a solution,
!   or the algorithm is making very slow progress and is not near a solution.
! TERMCD= 4 : Last global step failed to locate a point lower than XP.
!   It is likely that either XP is an approximate solution
!   of the problem or STEPTL is too large.
! TERMCD= 5 : Iteration limit exceeded.
! TERMCD= 6 : Five consecutive steps of length STEPMX have been taken.
  use M_Tensolve
  !
  character(len=*),      intent(in)   :: sCod1,sCod2
  real(dp),dimension(:), intent(inout):: vX
  real(dp),              intent(in)   :: TolF
  logical,               intent(in)   :: bFinDif
  integer,               intent(in)   :: MaxIts
  real(dp),              intent(out)  :: Error_F
  integer,               intent(out)  :: Nits
  integer,               intent(out)  :: iErr
  interface
    function Residual(v)
      use M_Kinds
      implicit none
      real(dp),dimension(:),intent(in):: v
      real(dp),dimension(size(v))     :: Residual
    end function Residual
    subroutine Jacobian(v,t)
      use M_Kinds
      implicit none
      real(dp),dimension(:),              intent(in) :: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
  endinterface
  !
  !MAXM >= M+N+2.
  !MAXN >= N+2.
  !MAXP >= NINT(sqrt(N)).
  integer :: maxm, maxn, maxp
  integer :: m, n, msg, termcd
  integer :: itnlim, jacflg, method, global, ipr
  !integer,parameter:: maxm=120, maxn=50, maxp= 6 !maxm= 100, maxn= 30
  !real(dp):: x0(maxn), xp(maxn), fp(maxm), gp(maxn)
  !real(dp):: typx(maxn),typf(maxm)
  real(dp),allocatable:: x0(:), xp(:), fp(:), gp(:)
  real(dp),allocatable:: typx(:),typf(:)
  real(dp):: gradtl, steptl, ftol, stepmx, dlt
  !
  Error_F = Zero
  Nits = 0
  iErr = 0
  !
  N=    size(vX)
  M=    size(vX)
  MaxN= N+2
  MaxM= M+N+2
  MaxP= NINT(SQRT(real(N))) +1
  allocate( x0(maxn), xp(maxn), fp(maxm), gp(maxn) )
  allocate( typx(maxn),typf(maxm) )
  x0(1:n)= vX(1:n)
  !
  !tsnesi= tensolve with default parameters
  !call tsnesi(maxm, maxn, maxp, x0, m, n, fvec_ts, jac_ts, msg, xp, fp, gp, termcd)
  !
  !set default parameters for module M_Tensolve
  call tsdflt( & 
  & m, n, itnlim, jacflg, gradtl, steptl,  &
  & ftol, method, global, stepmx, dlt, typx, typf, ipr)
  !modify some parameters.
  msg=    16
  ipr=    fTrc
  itnlim= MaxIts
  ftol=   TolF !1.0e-9_dp
  steptl= 1.0e-9_dp
  gradtl= 1.0e-5_dp
  !
  if(bFinDif) then; jacflg= 0
  else            ; jacflg= 1
  end if
  !
  if(trim(sCod1)=="TENSOR")      method= 1  !tensor
  if(trim(sCod1)=="NEWTON")      method= 0  !newton
  if(trim(sCod2)=="TRUSTREGION") global= 1  !trust region
  if(trim(sCod2)=="LINESEARCH")  global= 0  !linesearch
  
  call Stop_("Tensolve Solver Error")
  !--------------------------------------------
  !AM ==> One need to modify tensolve interface
  !--------------------------------------------
  !AM call tsneci( &
  !AM & maxm, maxn, maxp, x0, m, n, typx, typf, itnlim, jacflg,  &
  !AM & gradtl, steptl, ftol, method, global, stepmx, dlt, ipr,  &
  !AM & Residual, Jacobian, msg, xp, fp, gp, Nits, termcd)
  !
  
  vX(1:n)= xp(1:n)
  Error_F= MAXVAL(ABS(fp(1:n)))
  !
  select case(termcd)
    case(0); iErr= -6
    case(1); iErr= 0
    case(2); iErr= -1
    case(3); iErr= -2
    case(4); iErr= -3
    case(5); iErr= -4
    case(6); iErr= -5
  end select
  !
  deallocate( x0, xp, fp, gp )
  deallocate( typx,typf )

end subroutine TenSolve



end module M_Numeric_Tensolve
