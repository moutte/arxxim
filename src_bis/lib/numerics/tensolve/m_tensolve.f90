module M_TenSolve
  use M_Kinds
  use M_TenSolve_BlasPartial
  use M_TenSolve_UncMin
  implicit none
  !
  private
  !
  public:: tsnesi
  public:: tsdflt
  public:: tschki
  public:: tsnesv
  public:: tsneci
  !
  real(dp),allocatable, save:: fc(:), anls(:,:), aja(:,:), s(:,:)
  !FC   : FUCTION valueS AT XC
  !ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
  !AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
  !S    : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  integer,save:: qrank, meqns, nvars
  !Q: NUMERICAL RANK OF JACOBIAN :
  !   Q > P : JACOBIAN IS SINGULAR
  !   Q= P : OTHERWISE
  !
  real(dp),parameter:: ten=10.0_dp, half=0.5_dp, three=3.0_dp

contains
!*******************************************************************************
!ALGORITHM 768, COLLECTED ALGORITHMS FROM ACM.                                 *
!THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,                 *
!VOL. 232, NO. 21, June, 1997, pp. 174--195.                                   *
!*******************************************************************************
!TENSOLVE:  A Software Package for Solving Systems of Nonlinear Equations      *
!           and Nonlinear Least Squares Problems Using Tensor Methods.         *
!AUTHORS:   Ali Bouaricha                                                      *
!           Argonne National Laboratory                                        *
!           MCS Division                                                       *
!           e-mail: bouarich@mcs.anl.gov                                       *
!AND                                                                           *
!           Robert B. Schnabel                                                 *
!           University of colorado at Boulder                                  *
!           Department of computer Science                                     *
!           e-mail: bobby@cs.colorado.edu                                      *
!DATE:      Version of January, 1997                                           *
!*******************************************************************************
!Conversion to Fortran 90 by Alan Miller                                       *
!amiller@bigpond.net.au                                                        *
!*******************************************************************************
!1 Nov 98: Increased dimension of several automatic and allocated arrays
!          to prevent array bound errors.   Several minor changes to reduce
!          the number of fatal errors reported by ELF90.
!21 Jan 99: Fixed a bug in tsqfcn which was causing occasional underflows.
!24 May 01: Changed the intents of several variables in routine tschki.
!           Set termcd= 0 in routine tsnstp.
!           Initialized ierr in routine tsnesv.
!9 Sept 03: Added interface to JAC (renamed from TSDUMJ) in routine TSNESI.
!Latest revision - 9 September 2003
!Purpose of Tensolve:
!     TENSOLVE finds roots of systems of n nonlinear equations in n unknowns,
!     or minimizers of the sum of squares of m > n nonlinear equations in
!     n unknowns.  It allows the user to choose between a tensor method based
!     on a quadratic model and a standard method based on a linear model.
!     Both models calculate an analytic or finite-difference Jacobian matrix
!     at each iteration.  The tensor method augments the linear model with a
!     low-rank, second-order term that is chosen so that the model is hardly
!     more expensive to form, store, or solve than the standard linear model.
!     Either a line search or trust-region step-selection strategy may be
!     selected.  The tensor method is significantly more efficient than the
!     standard method in iterations, function evaluations, and time.  It is
!     especially useful on problems where the Jacobian at the solution is
!     singular.
!     The software can be called with an interface where the user supplies
!     only the function, number of nonlinear equations, number of variables,
!     and starting point; default choices are made for all the other input
!     parameters.  An alternative interface allows the user to specify any
!     input parameters that are different from the defaults.
! List of subroutine and function names called by TENSOLVE:
!    TS2DTR,TSBSLV,TSCHKI,TSCHKJ,TSCPMU,TSCPSS,TSD1SV,TSDFCN,TSDFLT,
!    TSDUMJ,TSFAFA,TSFDFJ,TSFRMT,TSFSCL,TSFSLV,TSJMUV,TSJQTP,TSLMIN,TSLMSP,
!    TSLSCH,TSMAFA,TSMDLS,TSMFDA,TSMFDV,TSMGSA,TSMSDA,TSMSDV,TSMSLV,TSNECI,
!    TSNESI,TSNESV,TSNSTP,TSPRMV,TSRSLT,TSQ1P1,TSQFCN,TSQLFC,TSQMLV,TSQMTS,
!    TSQMUV,TSQRFC,TSRSID,TSSCLF,TSSCLJ,TSSCLX,TSSLCT,TSSMIN,TSSMRD,TSSQP1,
!    TSSSTP,TSSTMX,TSTRUD,TSUDQV,TSUNSF,TSUNSX,TSUPSF,TSUTMD.
!*******************************************************************************
! Packages called by TENSOLVE:
!
!    UNCMIN (R. B. Schnabel, J. E. Koontz, and B. E. Weiss,
!    "A Modular System of Algorithms of Unconstrained Minimization",
!    ACM Trans. Math. Softw., 11 (1985), 419-440).
!
! BLAS called by TENSOLVE:
!    LEVEL 1 BLAS: DNRM2,DSWAP,IDAMAX
!    LEVEL 2 BLAS: DGEMV
!*******************************************************************************
! Parameters and Default Values for the interfaces TSNECI and TSNESI:
!
! Following each variable name in the list below appears a one- or
! two-headed arrow symbol of the form ->, <-, and <-->.
! These symbols signify that the variable is for input, output, and
! input-output, respectively.
! The symbol EPSM in some parts of this section designates the machine
! precision.
!
! X0->: An array of length N that contains an initial estimate
! of the solution x*.
!
! M->: A positive integer specifying the number of nonlinear equations.
!
! N->: A positive integer specifying the number of variables in the problem.
!
! TYPX->:  An array of length N in which the typical size of the
! components of X is specified. The typical component sizes should be
! positive real scalars. If a negative value is specified, its absolute
! value will be used. If 0.0 is specified, 1.0 will be used. This
! vector is used to determine the scaling matrix, Dx. Although the
! package may work reasonably well in a large number of instances without
! scaling, it may fail when the components of x* are of radically
! different magnitude and scaling is not invoked. If the sizes
! of the parameters are known to differ by many orders of magnitude, then
! the scale vector TYPX should definitely be used. For example, if
! it is anticipated that the range of values for the iterates xk would be
!                   x1 in [-10e+10,10e+10]
!                   x2 in [-10e+2,10e+4]
!                   x3 in [-6*10e-6,9*10e-6]
! then an appropriate choice would be TYPX= (1.0e+10,1.0e+3,7.0e-6).
! Module TSDFLT returns TYPX= (1.0,...,1.0).
!
! TYPF->: An array of length M in which the typical size of the components
! of F is specified. The typical component sizes should be positive real
! scalars.  If a negative value is specified, its absolute value will be
! used. If 0.0 is specified, 1.0 will be used. This vector is used to
! determine the scaling matrix DF. TYPF should be chosen so that all
! the components of DF(x) have similar typical magnitudes at points not
! too near a root, and should be chosen in conjunction with FTOL.  It is
! important to supply values of TYPF when the magnitudes of the components
! of F(x) are expected to be very different.  If the magnitudes of the
! components of F(x) are similar, the choice DF= I suffices.  Module
! TSDFLT returns TYPF= (1.0,...,1.0).
! ITNLIM->:  Positive integer specifying the maximum number of iterations to
! be performed before the program is terminated.   Module TSDFLT returns
! ITNLIM= 150. If the user specifies ITNLIM <= 0, the module TSCHKI will
! supply the value 150.
!
! JACFLG->: Integer designating whether or not an analytic Jacobian has
! been supplied by the user.
! JACFLG= 0 : No analytic Jacobian supplied.  The Jacobian is obtained
! by finite differences.
! JACFLG= 1 : Analytic Jacobian supplied.
! The module TSDFLT returns the value 0.  If the user specifies an illegal
! value, the module TSCHKI will supply the value 0.
!
! GRADTL->: Positive scalar giving the tolerance at which the scaled
! gradient of f(x)= 0.5*F(x)-trans*F(x) is considered close enough to
! zero to terminate the algorithm. The scaled gradient is a measure of
! the relative change in F in each direction xj divided by the relative
! change in xj. The module TSDFLT returns the value EPSM**(1/3).  If the
! user specifies a negative value, the module TSCHKI will supply
! the value EPSM**(1/3).
!
! STEPTL->: A positive scalar providing the minimum allowable relative
! step length. STEPTL should be at least as small as 10**(-d), where d
! is the number of accurate digits the user desires in the solution x*.
! The program may terminate prematurely if STEPTL is too large.  Module
! TSDFLT returns the value EPSM**(2/3).  If the user specifies a negative
! value, the module TSCHKI will supply the value EPSM**(2/3).
!
! FTOL->: A positive scalar giving the tolerance at which the scaled
! function DF*F(x) is considered close enough to zero to terminate the
! algorithm. The program is halted if ||DF*F(x)|| (in the infinity norm)
! is <= FTOL. This is the primary stopping condition for nonlinear
! equations; the values of TYPF and FTOL should be chosen so that this
! test reflects the user's idea of what constitutes a solution to the
! problem. The module TSDFLT returns the value EPSM**(2/3). If the
! user specifies a negative value, the module TSCHKI will supply the
! value EPSM**(2/3).
!
! METHOD->: An integer designating which method to use.
!   METHOD= 0 : Newton or Gauss-Newton algorithm is used.
!   METHOD= 1 : Tensor algorithm is used.
!   module TSDFLT returns value 1. If the user specifies an illegal value,
!   module TSCHKI will reset METHOD to 1.
!
! GLOBAL->: An integer designating which global strategy to use.
!   GLOBAL= 0 : Line search is used.
!   GLOBAL= 1 : Two-dimensional trust region is used.
!   Module TSDFLT returns value of 0. If the user specifies an illegal
!   value, module TSCHKI will reset GLOBAL to 0.
!
! STEPMX->: A positive scalar providing the maximum allowable scaled step
! length ||Dx*(x+ - xc)||2, where Dx= diag(1/TYPX_j). STEPMX is used to
! prevent steps that would cause the nonlinear equations problem to
! overflow, and to prevent the algorithm from leaving the area of
! interest in parameter space.  STEPMX should be chosen small enough
! to prevent these occurrences but should be larger than any anticipated
! "reasonable" step. Module TSDFLT returns the value STEPMX= 10e+3.
! If the user specifies a nonpositive value, module TSCHKI sets STEPMX
! to 10e+3.
!
! DLT->: A positive scalar giving the initial trust region radius. When
! the line search strategy is used, this parameter is ignored. For the
! trust region algorithm, if DLT is supplied, its value should reflect
! what the user considers a maximum reasonable scaled step length at
! the first iteration. If DLT= -1.0, the routine uses the length of
! the Cauchy step at the initial iterate instead. The module TSDFLT
! returns the value -1.0. If the user specifies a nonpositive value,
! module TSCHKI sets DLT= -1.0.
!
! IPR->: The unit on which the package outputs information.  TSDFLT returns
! the value 6.
!
! FVEC->: The name of a user-supplied subroutine that evaluates the function F
! at an arbitrary vector X.  The subroutine must conform to the usage
!                      call FVEC(X, F, M, N),
! where X is a vector of length N and F is a vector of length M.  The
! subroutine must not alter the values of X.
!
! JAC->: The name of a user-supplied subroutine that evaluates the first
! derivative (Jacobian) of the function F(x).  The subroutine must conform
! to the usage
!                      call JAC(X, AJA, M, N)
! where X is a vector of length N and the 2-dimensional array AJA of row
! dimension MAXM and column dimension N is the analytic Jacobian of F at
! X.  When using the interface TSNECI, if no analytic Jacobian is supplied
! (JACFLG= 0), the user must use the dummy name TSDUMJ as the value of
! this parameter.
!
! MSG<-->: An integer variable that the user may set on input to inhibit
! certain automatic checks or to override certain default characteristics
! of the package. (For the short call it should be set to 0.) There are
! four "message" features that can be used individually or in combination
! as discussed below.
! MSG= 0 : Values of input parameters, final results, and termination code
! are printed.
! MSG= 2 : Do not check user's analytic Jacobian routine against its
! finite difference estimate.  This may be necessary if the user knows the
! Jacobian is properly coded, but the program aborts because the comparative
! tolerance is too tight.  Do not use MSG= 2 if the analytic Jacobian is
! not supplied.
! MSG= 8 : Suppress printing of the input state, the final results, and
! the stopping condition.
! MSG= 16 : Print the intermediate results; that is, the input state,
! each iteration including the current iterate xk, 0.5*||DF*F(xk)||2**2,
! and grad(f(x))= J(x)-trans*DF**2 F(x), and the final results including
! the stopping conditions.
! The user may specify a combination of features by setting MSG to
! the sum of the individual components. The module TSDFLT returns a value
! of 0. On exit, if the program has terminated because of erroneous
! input, MSG contains an error code indicating the reason.
! MSG= 0   : No error.
! MSG= -1  : Illegal dimension, M <= 0.
! MSG= -2  : Illegal dimension, N <= 0.
! MSG= -3  : Illegal dimension, MAXM < M+N+2.
! MSG= -4  : Illegal dimension, MAXN < N+2.
! MSG= -5  : Illegal dimension, MAXP < NINT(sqrt(N)).
! MSG= -10 : Program asked to override check of analytic Jacobian
! against finite difference estimate, but routine JAC not
! supplied (incompatible input).
! MSG= -11  : Probable coding error in the user's analytic Jacobian
! routine JAC. Analytic and finite difference Jacobian do not agree
! within the assigned tolerance.
!
! XP<-: An array of length N containing the best approximation
! to the solution x*. (If the algorithm has not converged, the final
! iterate is returned).
!
! FP<-: An array of length M containing the function value F(XP).
!
! GP<-: An array of length N containing the gradient of the
! function 0.5*||F(x)||2**2 at XP.
!
! TERMCD<-:  An integer specifying the reason for termination.
! TERMCD= 0 : No termination criterion satisfied (occurs if package
! terminates because of illegal input).
! TERMCD= 1 : function tolerance reached.  The current iteration is
! probably a solution.
! TERMCD= 2 : gradient tolerance reached.  For nonlinear least
! squares, the current iteration is probably a solution; for nonlinear
! equations, it could be a solution or a local minimizer.
! TERMCD= 3 : Successive iterates within step tolerance.  The
! current iterate may be a solution, or the algorithm is making very slow
! progress and is not near a solution.
! TERMCD= 4 : Last global step failed to locate a point lower
! than XP. It is likely that either XP is an approximate solution
! of the problem or STEPTL is too large.
! TERMCD= 5 : Iteration limit exceeded.
! TERMCD= 6 : Five consecutive steps of length STEPMX have been taken.
!///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\
!\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///\\\///

subroutine tsnesi( &
& maxm, maxn, maxp, & !in
& x0,               & !inout
& m, n,             & !in
& fvec, jac,        & !
& msg, xp, fp, gp,  & !inout
& itn,termcd)             !out
!subroutine tsnesi(maxm, maxn, maxp, x0, m, n, fvec, msg, xp, fp, gp, termcd)
!A SIMPLE interface TO THE NONLINEAR EQUATION/NONLINEAR LEAST SQUARES PROBLEMS PACKAGE.
!THE useR HAS NO CONTROL OVER THE PACKAGE OPTIONS.
  integer, intent(in)   :: maxm
  integer, intent(in)   :: maxn
  integer, intent(in)   :: maxp
  real(dp),intent(inout):: x0(:)
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  integer, intent(inout):: msg
  real(dp),intent(inout):: xp(:)
  real(dp),intent(inout):: fp(:)
  real(dp),intent(inout):: gp(:)
  integer, intent(out)  :: itn,termcd
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
    end subroutine fvec
    subroutine jac(x, aja, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: aja(:,:)
      integer, intent(in)    :: m, n
    end subroutine jac
  end interface
  !SUBprogramS callED:
  !  TENSOLVE      ...  TSDFLT,TSCHKI,TSNESV
  real(dp):: typx(n), typf(m), dfn(m), dxn(n)
  integer :: jacflg, itnlim, method
  integer :: global, ipr
  real(dp):: steptl, gradtl, ftol, stepmx, dlt
  integer :: sqrn
  real(dp):: epsm
  ! set default values for each variable to the nonlinear equations/
  ! nonlinear least squares solver
  call tsdflt( &
  & m, n, itnlim, jacflg, gradtl, steptl, ftol, method, global,  &
  & stepmx, dlt, typx, typf, ipr)
  ! check input parameters
  call tschki( &
  & maxm, maxn, maxp, m, n, gradtl, steptl, ftol, itnlim, jacflg,  &
  & method, global, stepmx, dlt, epsm, msg, typx, typf, dxn, dfn,  &
  & sqrn, termcd, ipr)
  if(msg < 0) return
  ! call nonlinear equations/nonlinear least squares solver
  call tsnesv( &
  & maxm, x0, m, n, typx, typf, itnlim, jacflg,  &
  & gradtl, steptl, ftol, method, global, stepmx, dlt, ipr,  &
  & dfn, dxn, epsm, sqrn, fvec, jac, msg, xp, fp, gp, itn, termcd)
  return
end subroutine tsnesi

subroutine tsneci( &
& maxm, maxn, maxp, x0, m, n, typx, typf, itnlim, jacflg, &
& gradtl, steptl, ftol, method, global, stepmx, dlt, ipr, &
& fvec, jac, msg, xp, fp, gp, itn, termcd)
!A COMPLETE interface TO THE NONLINEAR EQUATION/NONLINEAR LEAST SQUARES PACKAGE.
!THE useR HAS FULL CONTROL OVER THE OPTIONS.
!SUBprogramS callED:
!    TENSOLVE      ...  TSCHKI,TSNESV
  integer, intent(in)   :: maxm
  integer, intent(in)   :: maxn
  integer, intent(in)   :: maxp
  real(dp),intent(inout):: x0(:)
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  real(dp),intent(inout):: typx(:)
  real(dp),intent(inout):: typf(:)
  integer, intent(inout):: itnlim
  integer, intent(inout):: jacflg
  real(dp),intent(inout):: gradtl
  real(dp),intent(inout):: steptl
  real(dp),intent(inout):: ftol
  integer, intent(inout):: method
  integer, intent(inout):: global
  real(dp),intent(inout):: stepmx
  real(dp),intent(inout):: dlt
  integer, intent(inout):: ipr
  integer, intent(inout):: msg
  real(dp),intent(out)  :: xp(:)
  real(dp),intent(out)  :: fp(:)
  real(dp),intent(out)  :: gp(:)
  integer, intent(out)  :: itn,termcd
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)  :: m, n
    end subroutine fvec
    subroutine jac(x, aja, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: aja(:,:)
      integer, intent(in)    :: m, n
  end subroutine jac
  end interface
  !
  integer :: sqrn
  real(dp):: epsm, dfn(m), dxn(n)
  !
  ! check input parameters
  call tschki( &
  & maxm, maxn, maxp, m, n, gradtl, steptl, ftol, itnlim, jacflg, &
  & method, global, stepmx, dlt, epsm, msg, typx, typf, dxn, dfn, &
  & sqrn, termcd, ipr)
  if(msg < 0) return
  ! call nonlinear equations/nonlinear least squares solver
  call tsnesv( &
  & maxm, x0, m, n, typx, typf, itnlim, jacflg,  &
  & gradtl, steptl, ftol, method, global, stepmx, dlt, ipr,  &
  & dfn, dxn, epsm, sqrn, fvec, jac, msg, xp, fp, gp, itn, termcd)
  return
end subroutine tsneci

subroutine ts2dtr( &
& aja, shat, anls, dt, g, gbar, xc, method, nwtake, stepmx, &
& steptl, epsm, mxtake, dlt, fq, maxm, m, n, p,    &
& curpos, pivot, pbar, itn, ierr, flag, dxn, dfn, fvec,  &
& fnorm, xpls, fp, fpls, retcd)
!FINDS A NEXT ITERATE BY A 2-dimensionAL TRUST REGION.
  real(dp),intent(in)   :: aja(:,:)
  real(dp),intent(in)   :: shat(:,:)
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(inout):: dt(:)
  real(dp),intent(in)   :: g(:)
  real(dp),intent(inout):: gbar(:)
  real(dp),intent(in)   :: xc(:)
  integer, intent(in)   :: method
  logical, intent(in)   :: nwtake
  real(dp),intent(inout):: stepmx
  real(dp),intent(inout):: steptl
  real(dp),intent(in)   :: epsm
  logical, intent(inout):: mxtake
  real(dp),intent(inout):: dlt
  real(dp),intent(in)   :: fq(:)
  integer, intent(in)   :: maxm
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  integer, intent(in)   :: p
  integer, intent(in)   :: curpos(:)
  integer, intent(in)   :: pivot(:)
  integer, intent(in)   :: pbar(:)
  integer, intent(in)   :: itn
  integer, intent(in)   :: ierr
  integer, intent(in)   :: flag
  real(dp),intent(in)   :: dxn(:)
  real(dp),intent(in)   :: dfn(:)
  real(dp),intent(in)   :: fnorm
  real(dp),intent(out)  :: xpls(:)
  real(dp),intent(out)  :: fp(:)
  real(dp),intent(out)  :: fpls
  integer, intent(out)  :: retcd
  !INPUT parameterS :
  !       AJA    : JACOBIAN AT THE CURRENT ITERATE
  !       SHAT   : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !                AFTER A QL FACTORIZATION
  !       ANLS   : TENSOR TERM MATRIX
  !       DT     : CURRENT STEP
  !       G      : GRADIENT AT CURRENT ITERATE
  !       GBAR   : STEEPEST DESCENT DIRECTION (= -G)
  !       XC     : CURRENT ITERATE
  !       METHOD : METHOD TO use
  !                 =  0  : STANDARD METHOD useD
  !                 =  1  : TENSOR METHOD useD
  !       NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !                NWTAKE =  .true.  : STANDARD STEP TAKEN
  !                NWTAKE =  .false. : TENSOR STEP TAKEN
  !       STEPMX : MAXIMUM STEP ALLOWED
  !       STEPTL : STEP TOLERANCE
  !       EPSM   : MACHINE PRECISION
  !       MXTAKE : BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  !       FQ     : function value AT CURRENT ITERATE MULTIPLIED BY
  !                ORTHOGONAL MATRICES
  !       MAXM   : LEADING dimension OF AJA AND ANLS
  !       M,N    : dimensionS OF PROBLEM
  !       P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !       CURPOS : PIVOT VECTOR (useD DURING THE FACTORIZATION OF THE
  !                JACOBIAN FROM COLUMN 1 TO N-P)
  !       PIVOT  : PIVOT VECTOR ( useD DURING THE FACTORIZATION OF THE
  !                JACOBIAN FROM COLUMN N-P+1 TO N)
  !       PBAR   : PIVOT VECTOR (useD DURING THE FACTORIZATION OF THE
  !                JACOBIAN if IT IS SINGULAR
  !       FNORM  :  0.5 * || FC ||**2
  !       ITN    : ITERATION NUMBER
  !       IERR   : return CODE FROM THE QRP FACTORIZATION ROUTINE:
  !                IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !                IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !       FLAG   : return CODE WITH THE FOLLOWING MEANINGS :
  !                FLAG =  0 : NO SINGULARITY DETECTED DURING FACTORIZATION OF
  !                             THE JACOBIAN FROM COLUMN 1 TO N
  !                FLAG =  1 : SINGULARITY DETECTED DURING FACTORIZATION
  !                             OF THE JACOBIAN FROM COLUMN 1 TO N-P
  !                FLAG =  2 : SINGULARITY DETECTED DURING FACTORIZATION
  !                             OF THE JACOBIAN FROM COLUMN N-P+1 TO N
  !        DXN   : DIAGONAL SCALING MATRIX FOR X
  !        DFN   : DIAGONAL SCALING MATRIX FOR F
  !        FVEC  : subroutine TO EVALUATE THE useR'S function
  !INPUT-OUTPUT parameterS :
  !       DLT    : INITIAL TRUST RADIUS (= -1.0D0) if IT IS NOT SUPPLIED
  !                BY THE useR ON entry AND CURRENT TRUST RADIUS ON exit
  !OUTPUT parameterS :
  !       XPLS   : NEXT ITERATE
  !       FP     : function value AT NEXT ITERATE
  !       FPLS   : 0.5 * || FP ||**2
  !       RETCD  : return CODE FROM subroutine (SEE subroutine TSTRUD FOR MEANING)
  !       SUBprogramS callED:
  !       LEVEL 1 BLAS  ...  DNRM2
  !       TENSOLVE      ...  TSPRMV,TSUTMD,TSJMUV,TSUDQV,TSSMIN,TSRSID,
  !       TENSOLVE      ...  TSTRUD
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds,only: dp
      implicit none
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: f(:)
      integer, intent(in) :: m, n
    end subroutine fvec
  end interface
  ! Workspace
  real(dp):: &
  & d(m), xplsp(n), adt(n), ag(n), temp(m), vn(m), vnp(m), vns(m), &
  & const1(p), const2(p)
  real(dp):: res, alph, sum, resg, optim
  real(dp):: scres, fplsp, rresg
  logical :: dtaken
  !
  dtaken= .false.
  retcd= 4
  if(dlt== -one) then
    ! set DLT to length of Cauchy step
    alph= dnrm2(n, g, 1)
    alph= alph**2
    call tsprmv(vn, g, pivot, n, 1)
    if(ierr== 0) then
      call tsutmd(aja, vn, m, n, vnp)
    else
      call tsprmv(vns, vn, pbar, n, 1)
      call tsutmd(aja, vns, m+n, n, vnp)
    end if
    dlt= alph*SQRT(alph) / dnrm2(n, vnp, 1)**2
    if(dlt > stepmx) then
      dlt= stepmx
    end if
  end if
  ! form an orthonormal basis for the two-dimensional subspace
  gbar(1:n)= -g(1:n)
  res= dnrm2(n, dt, 1)
  sum= -dot_product( gbar(1:n), dt(1:n) ) / res**2
  gbar(1:n)= gbar(1:n) + sum * dt(1:n)
  resg= dnrm2(n, gbar, 1)
  if(resg > zero) then
    rresg= one/resg
    gbar(1:n)= rresg * gbar(1:n)
  end if
  dt(1:n)= dt(1:n) / res
  ! compute Jacobian times DT
  call tsjmuv( &
  & itn, method, dt, curpos, pivot, pbar, aja, shat,  &
  & flag, ierr, m, n, p, vn, adt)
  ! compute Jacobian times GBAR
  call tsjmuv( &
  & itn, method, gbar, curpos, pivot, pbar, aja, shat,  &
  & flag, ierr, m, n, p, vnp, ag)
  if(.not. nwtake) then
    ! compute SHAT times VN
    call tsudqv(shat, vn, n, p, const1)
    ! compute SHAT times VNP
    call tsudqv(shat, vnp, n, p, const2)
  end if
  ! normalize DT
  70 if(res <= dlt) then
    dtaken= .true.
    d(1:n)= dt(1:n)*res
    dlt= res
    else
      ! find the global minimizer of one-variable function in the
      ! interval (-dlt, dlt)
      call tssmin( &
      & anls, fq, adt, ag, const1, const2, dlt, m, n,  &
      & p, nwtake, ierr, epsm, optim)
      ! compute the global step
      d(1:n)= optim*dt(1:n) + SQRT(dlt**2 - optim**2) * gbar(1:n)
    end if
  ! compute the tensor model residual
  call tsrsid( &
  & itn, method, fq, d, curpos, pivot, pbar, aja, anls,  &
  & shat, flag, nwtake, ierr, maxm, m, n, p, scres)
  ! check whether the global step is acceptable
  call tstrud( &
  & m, n, xc, fnorm, g, d, dtaken, stepmx, steptl, dlt, mxtake, &
  & dxn, dfn, fvec, scres, retcd, xplsp, fplsp, temp, xpls, fp, fpls)
  if(retcd >= 2) GO TO 70
  return
end subroutine ts2dtr

subroutine tsbslv(r, m, n, b, y)
!doES A BACKWARD SOLVE.
  real(dp),intent(in)   :: r(:,:)
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  real(dp),intent(in)   :: b(:)
  real(dp),intent(out)  :: y(:)
  !INPUT parameterS :
  !    R  : UPPER TRIANGULAR MATRIX OBTAINED FROM A QR FACTORIZATION
  !         OF AN M BY N MATRIX A. DIAG(R) IS STORED IN ROW M+2. THIS
  !         IS THE STORAGE SCHEME useD IN STEWART, G. W., III(1973)
  !         "INTRODUCTION TO MATRIX COMPUTATION", ACADEMIC PRESS, NEW YORK
  !    M  : ROW dimension OF MATRIX A
  !    N  : COLUMN dimension OF MATRIX A
  !    B  : RIGHT HAND SIDE
  !OUTPUT parameterS :
  !    Y :  VECTOR SOLUTION ON exit
  integer :: j
  ! solve R Y= B
  y(n)= b(n) / r(m+2,n)
  if(n > 2) then
    y(1:n-1)= zero
    do j= n-1,2,-1
      y(1:j)= y(1:j) + y(j+1) * r(1:j,j+1)
      y(j)=   (b(j) - y(j))/r(m+2,j)
    enddo
    y(1)= y(1) + r(1,2)*y(2)
    y(1)= (b(1) - y(1)) /r(m+2,1)
  else if(n== 2) then
    y(1)= (b(1) - (r(1,2) * y(2))) / r(m+2,1)
  end if
  return
end subroutine tsbslv

subroutine tschki( &
!CHECKS INPUT FOR REASONABlenESS.
& maxm, maxn, maxp, m, n, gradtl, steptl, ftol, itnlim, &
& jacflg, method, global, stepmx, dlt, epsm, msg, typx, &
& typf, dxn, dfn, sqrn, termcd, ipr)
  integer, intent(in)   :: maxm
  integer, intent(in)   :: maxn
  integer, intent(in)   :: maxp
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  real(dp),intent(inout):: gradtl
  real(dp),intent(inout):: steptl
  real(dp),intent(inout):: ftol
  integer, intent(inout):: itnlim
  integer, intent(inout):: jacflg
  integer, intent(inout):: method
  integer, intent(inout):: global
  real(dp),intent(inout):: stepmx
  real(dp),intent(inout):: dlt
  real(dp),intent(out)  :: epsm
  integer, intent(inout):: msg
  real(dp),intent(inout):: typx(:)
  real(dp),intent(inout):: typf(:)
  real(dp),intent(out)  :: dxn(:)
  real(dp),intent(out)  :: dfn(:)
  integer, intent(out)  :: sqrn
  integer, intent(out)  :: termcd
  integer, intent(in)   :: ipr
  ! N.B. intent of gradtl, steptl, ftol, stepmx, dlt, msg changed to IN/OUT
  !      24 May 2001
  !INPUT parameterS :
  !    MAXM : LEADING dimension OF WORKSPACE WRKNEM
  !           (SEE TOP OF THIS FILE )
  !    MAXN : LEADING dimension OF WORKSPACE WRKNEN
  !           (SEE TOP OF THIS FILE )
  !    MAXP : LEADING dimension OF WORKSPACE WRKUNC
  !           (SEE TOP OF THIS FILE )
  !    M,N  : dimensionS OF PROBLEM
  !    IPR  : DEVICE TO WHICH TO Send OUTPUT
  !INPUT-OUTPUT parameterS :
  !    GRADTL : TOLERANCE AT WHICH GRADIENT CONSIDERED close
  !             ENOUGH TO ZERO TO TERMINATE ALGORITHM
  !    STEPTL : TOLERANCE AT WHICH SUCCESSIVE ITERATES
  !             CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  !    FTOL   : TOLERANCE AT WHICH function value CONSIDERED
  !             close ENOUGH TO ZERO
  !    ITNLIM : MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  !    STEPMX : MAXIMUM STEP ALLOWED IN TRUST REGION
  !    DLT    : TRUST RADIUS
  !    JACFLG : JACOBIAN FLAG WITH THE FOLLOWING MEANINGS :
  !             JACFLG= 1 : ANALYTIC JACOBIAN SUPPLIED
  !             JACFLG= 0 : ANALYTIC JACOBIAN NOT SUPPLIED
  !    METHOD : METHOD TO use
  !             METHOD= 0 : STANDARD METHOD IS useD
  !             METHOD= 1 : TENSOR METHOD IS useD
  !    GLOBAL : GLOBAL STRATEGY useD
  !             GLOBAL= 0 : LINE SEARCH useD
  !             GLOBAL= 1 : 2-dimensionAL TRUST REGION useD
  !    TYPX   : TYPICAL size FOR EACH COMPONENT OF X
  !    TYPF   : TYPICAL size FOR EACH COMPONENT OF F
  !    MSG    : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  !OUTPUT parameterS :
  !    TERMCD: TERMINATION CODE
  !    DXN   : DIAGONAL SCALING MATRIX FOR X
  !    DFN   : DIAGONAL SCALING MATRIX FOR F
  !    SQRN  : MAXIMUM COLUMN dimension OF S AND FV
  real(dp),parameter:: thous= 1000._dp
  integer :: i, leng
  real(dp):: temp
  ! check that parameters only take on acceptable values
  ! if not, set them to default values
  ! set TERMCD to zero in case we abort prematuraly
  termcd= 0
  ! compute machine precision
  epsm= EPSILON( 1.0D0 )
  ! check dimensions of the problem
  if(m <= 0) then
    write(ipr,'(A,I5)') '  TSCHKI     ILLEGAL dimension M=',m
    msg= -1
    return
  end if
  if(n <= 0) then
    write(ipr,'(A,I5)') '  TSCHKI     ILLEGAL dimension N=',n
    msg= -2
    return
  end if
  ! check leading dimensions of the problem
  leng= m+n+2
  if(maxm < leng) then
    write(ipr,603) maxm,leng
    msg= -3
    return
  end if
  leng= n+2
  if(maxn < leng) then
    write(ipr,604) maxn,leng
    msg= -4
    return
  end if
  temp= SQRT(DBLE(n))
  sqrn= nint(temp)
  if(maxp < sqrn) then
    write(ipr,605) maxp,sqrn
    msg= -5
    return
  end if
  !
  ! check JACFLG, METHOD, and GLOBAL
  if(jacflg /= 1) jacflg= 0
  if(method /= 0 .and. method /= 1) method= 1
  if(global /= 0 .and. global /= 1) global= 0
  if(MOD(msg/2,2)== 1 .and. jacflg== 0) then
    write(ipr,610) msg,jacflg
    msg= -10
    return
  end if
  !
  ! check scale matrices
  do i= 1,n
    if(typx(i)== zero) typx(i)= one
    if(typx(i) < zero) typx(i)= -typx(i)
    dxn(i)= one/typx(i)
  enddo
  do i= 1,m
    if(typf(i)== zero) typf(i)= one
    if(typf(i) < zero) typf(i)= -typf(i)
    dfn(i)= one/typf(i)
  enddo
  !
  ! check gradient, step, and function tolerances
  temp= one/three
  if(gradtl < zero) gradtl= epsm**temp
  if(steptl < zero) steptl= epsm**(two*temp)
  if(ftol   < zero) ftol=   epsm**(two*temp)
  !
  ! check iteration limit
  if(itnlim <= 0)   itnlim= 150
  ! check STEPMX and DLT
  if(stepmx < zero) stepmx= thous
  if(dlt <= zero) then
    dlt= -one
    if(dlt > stepmx) dlt= stepmx
  end if
  return
  !601 format('  TSCHKI     ILLEGAL dimension M=',i5)
  !602 format('  TSCHKI     ILLEGAL dimension N=',i5)
  603 format('  TSCHKI     ILLEGAL dimension MAXM=',i5,' < M+N+2=',i5)
  604 format('  TSCHKI     ILLEGAL dimension MAXN=',i5,' < N+2=',i5)
  605 format('  TSCHKI     ILLEGAL dimension MAXP=',i5,' <',  &
             '  NINT(SQRT (N))=',i5)
  610 format('  TSCHKI     useR REQUESTS THAT ANALYTIC JACOBIAN BE',  &
             ' ACCEPTED AS PROPERLY CODED (MSG=',i5,')'/  &
             '  TSCHKI     BUT ANALYTIC JACOBIAN NOT SUPPLIED',  &
             ' (JACFLG=',i5,')')
end subroutine tschki

!-------------------------------------------------------------------------------
!-- subroutine tschkj --
!-- CHECKS THE ANALYTIC JACOBIAN AGAINST ITS FINITE DifFERENCE APPROXIMATION.
subroutine tschkj( &
& ajanal, xc, fc, m, n, epsm, dfn, dxn, typx, ipr, fvec, msg)
  real(dp),intent(in)   :: ajanal(:,:)
  real(dp),intent(inout):: xc(:)
  real(dp),intent(inout):: fc(:)
  integer,  intent(in)  :: m
  integer,  intent(in)  :: n
  real(dp),intent(in)   :: epsm
  real(dp),intent(in)   :: dfn(:)
  real(dp),intent(in)   :: dxn(:)
  real(dp),intent(in)   :: typx(:)
  integer, intent(in)   :: ipr
  integer, intent(inout):: msg
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds,only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
    end subroutine fvec
  end interface
  !
  !INPUT parameterS :
  !     AJANAL : ANALYTIC JACOBIAN AT XC
  !     XC   : CURRENT ITERATE
  !     FC   : function value AT XC
  !     M,N  : dimensionS OF PROBLEM
  !     EPSM : MACHINE PRECISION
  !     DFN  : DIAGONAL SCALING MATRIX FOR F
  !     DXN  : DIAGONAL SCALING MATRIX FOR X
  !     TYPX : TYPICAL size FOR EACH COMPONENT OF X
  !     IPR  : DEVICE TO WHICH TO Send OUTPUT
  !     FVEC : subroutine TO EVALUATE THE useR'S function
  !INPUT-OUTPUT parameterS :
  !     MSG : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  !     SUBprogramS callED:
  !     LEVEL 1 BLAS  ...  IDAMAX
  !     TENSOLVE      ...  TSUNSX,TSUNSF,TSSCLX,TSSCLF
  !     useR          ...  FVEC
  !
  ! Workspace
  real(dp):: fhat(m), wrk1(m)
  integer :: i, j
  real(dp):: ndigit, rnoise, sqrns, stepsz, xtmpj, dinf, rstpsz
  real(dp):: tol
  real(dp),parameter:: quart= 0.25_dp
  !
  ! unscale XC and FC
  call tsunsx(xc, dxn, n)
  call tsunsf(fc, dfn, m)
  !
  ! compute the finite difference Jacobian and check it against the analytic one
  ndigit= -LOG10(epsm)
  rnoise= MAX(ten**(-ndigit),epsm)
  sqrns = SQRT(rnoise)
  tol= epsm**quart
  do j= 1,n
    stepsz= sqrns*MAX(ABS(xc(j)),one)
    xtmpj= xc(j)
    xc(j)= xtmpj + stepsz
    call fvec(xc, fhat, m, n)
    xc(j)= xtmpj
    rstpsz= one/stepsz
    do i= 1,m
      wrk1(i)= (fhat(i) - fc(i))*rstpsz
    enddo
    do i= 1,m
      wrk1(i)= wrk1(i)*dfn(i)*typx(j)
    enddo
    dinf= ABS(wrk1(idamax(m, wrk1, 1)))
    do i= 1,m
      if(ABS(ajanal(i,j) - wrk1(i)) > tol*dinf) then
        write(ipr,50)
        msg= -11
        return
      end if
    enddo
  enddo
  ! scale back XC and FC
  call tssclx(xc,dxn,n)
  call tssclf(fc,dfn,m)
  50 format(/'  TSCHKJ      PROBABLE ERROR IN CODING OF ANALYTIC JACOBIAN')
  return
end subroutine tschkj

subroutine tscpmu(r, n, epsm, mu)
!COMPUTES A SMALL PERTURBATION MU.
!MU IS useD IN THE COMPUTATION OF THE LEVENBERG-MARQUARDT STEP.
  real(dp),intent(in)   :: r(:,:)
  integer, intent(in)     :: n
  real(dp),intent(in)   :: epsm
  real(dp),intent(out)  :: mu
  !INPUT parameterS :
  !    R  : UPPER TRIANGULAR MATRIX
  !    N  : COLUMN dimension OF R
  !    EPSM :  MACHINE PRECISION
  !OUTPUT parameterS :
  !    MU= SQRT(L1 NORM OF R * INFINITY NORM OF R * N * EPSM * 100)
  integer              :: i,j
  real(dp)            :: aifnrm, total, al1nrm
  real(dp), parameter :: hund= 100.0_dp
  ! compute the infinity norm of R
  aifnrm= zero
  do i= 1,n
    total= SUM( ABS(r(i,i:n)) )
    aifnrm= MAX(aifnrm,total)
  enddo
  ! compute the l1 norm of R
  al1nrm= zero
  do j= 1,n
    total= SUM( r(1:j,j) )
    al1nrm= MAX(al1nrm,total)
  enddo
  ! compute MU
  mu= SQRT(aifnrm*al1nrm*n*epsm*hund)
  return
end subroutine tscpmu

subroutine tscpss( &
& s, m, n, p, method, global, epsm, fcq, qhat, anls,  &
& dn, fqq, ptilda, curpos, pbar, zero1, ierr, resnew, flag)
!
!THIS ROUTINE COMPUTES THE STANDARD STEP.  NOTE THAT AT THIS STAGE
!THE JACOBIAN MATRIX (QHAT) HAS ALreadY BEEN FACTORED FROM COLUMNS 1
!TO N-P DURING THE TENSOR STEP COMPUTATION.  THIS ROUTINE FACTORS
!THE MATRIX QHAT FROM COLUMN N-P+1 TO N, THEREBY OBTAINING A QR
!FACTORIZATION OF THE FULL MATRIX QHAT, then COMPUTES THE STANDARD
!STEP BY PREMULTIPLYING THE RIGH-HAND SIDE FCQ BY AN ORTHOGONAL
!MATRIX AND BY PERFORMING A BACKWARD SOLVE.
!
  real(dp),intent(in)   :: s(:,:)
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  integer, intent(in)   :: p
  integer, intent(in)   :: method
  integer, intent(in)   :: global
  real(dp),intent(in)   :: epsm
  real(dp),intent(in)   :: fcq(:)
  real(dp),intent(inout):: qhat(:,:)
  real(dp),intent(inout):: anls(:,:)
  real(dp),intent(out)  :: dn(:)
  real(dp),intent(out)  :: fqq(:)
  integer, intent(out)  :: ptilda(:)
  integer, intent(in)   :: curpos(:)
  integer, intent(out)  :: pbar(:)
  integer, intent(out)  :: zero1
  integer, intent(inout):: ierr
  real(dp),intent(out)  :: resnew
  integer, intent(out)  :: flag
  !INPUT parameterS :
  !     S    : FACTORED MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !            (OBTAINED FROM TSQLFC subroutine)
  !     M,N  : dimensionS OF PROBLEM
  !     P    : COLUMN dimension OF MATRIX S
  !   METHOD : METHOD useD :
  !            METHOD= 0 : STANDARD METHOD IS useD
  !            METHOD= 1 : TENSOR METHOD IS useD
  !   GLOBAL : GLOBAL STRATEGY useD
  !            GLOBAL= 0 : LINE SEARCH IS useD
  !            GLOBAL= 1 : 2-dimensionAL TRUST REGION IS useD
  !   EPSM   : MACHINE PRECISION
  !   FCQ    : function value AT CURRENT ITERATE MULTIPLIED BY AN
  !            ORTHOGONAL MATRIX
  !   CURPOS : PIVOT VECTOR FOR THE FACTORIZATION OF QHAT FROM COLUMN 1 TO N-P
  !     INPUT-OUTPUT parameterS :
  !     QHAT  : FACTORED MATRIX FROM COLUMN 1 TO N-P
  !             ON entry AND FACTORED MATRIX FROM 1 TO N ON exit
  !     ANLS  : TENSOR TERM MATRIX ON entry AND ANLS MULTIPLIED BY
  !             ORTHOGONAL MATRICES ON exit (THIS IS PERFORMED IN THE
  !             case where THE GLOBAL STRATEGY useD IS THE 2-dimensionAL
  !             TRUST REGION)
  !OUTPUT parameterS :
  !     DN    : STANDARD STEP
  !     FQQ   : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES (THIS IS useD IN THE case where
  !             THE GLOBAL STRATEGY useD IS THE 2-dimensionAL
  !             TRUST REGION)
  !     PTILDA: PIVOT VECTOR FOR THE FACTORIZATION OF THE
  !             MATRIX QHAT FROM N-P+1 TO N
  !     PBAR  : PIVOT VECTOR FOR THE FACTORIZATION OF THE
  !             TRANSFORMED MATRIX QHAT FROM 1 TO N
  !             IN case OF SINGULARITY
  !     ZERO1 : FIRST ZERO COLUMN OF MATRIX QHAT IN case OF SINGULARITY
  !     IERR  : returnED CODE WITH THE FOLLOWING MEANING :
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 0 : OTHERWISE
  !     RESNEW: RESIDUAL OF THE STANDARD MODEL
  !     FLAG  : returnED CODE WITH THE FOLLOWING MEANINGS :
  !             FLAG= 0 : NO SINGULARITY DETECTED
  !             FLAG= 1 : SINGULARITY DETECTED DURING QR FACTORIZATION
  !                        OF QHAT FROM COLUMN 1 TO N-P
  !             FLAG= 2 : SINGULARITY DETECTED DURING QR FACTORIZATION
  !                        OF QHAT FROM COLUMN N-P+1 TO N
  !SUBprogramS callED:
  !     TENSOLVE      ...  TSQRFC,TSQMUV,TSQMTS,TSSMRD,TSBSLV,TSPRMV
  !     TENSOLVE      ...  TSQMLV,TSCPMU
  real(dp):: y(n), w(m+n), fqt(m+n)
  integer :: zerotm, i, j
  real(dp):: mu
  flag= 0
  ! initialization
  fqq(1:m+n)= zero
  w(1:m)= -fcq(1:m)
  ! if the Jacobian is singular then compute the Levenberg-Marquardt
  ! step (label 20)
  if(ierr== 1) then
    flag= 1
    GO TO 20
  end if
  ! factor the matrix QHAT from column n-p+1 to n
  call tsqrfc(qhat, n, m, n-p+1, n, ierr, epsm, ptilda, zero1)
  if(m== n .and. ierr== 0) then
    zerotm= zero1 - 1
  else
    zerotm= zero1
  end if
  ! premultiply W by the orthogonal matrix resulting from the QR
  ! factorization of QHAT
  call tsqmuv(qhat, w, fqq, m, n-p+1, zerotm, .false.)
  if(method== 1 .and. global== 1) &
  ! premultiply ANLS by the orthogonal matrix resulting from the QR
  ! factorization of QHAT
  & call tsqmts(anls, qhat, m, m, p, n-p+1, zerotm)
  if(ierr== 1) then
    flag= 2
    GO TO 20
  end if
  ! compute the residual of the standard model
  call tssmrd(fqq, resnew, dn, mu, ierr, m, n)
  ! if QHAT is nonsingular perform a backward solve to obtain Y
  call tsbslv(qhat, m, n, fqq, y)
  ! pivot Y
  call tsprmv(dn, y, ptilda, n, 0)
  if(n /= 1) then
    call tsprmv(y, dn, curpos, n, 0)
    ! premultiply Y by the orthogonal matrix resulting from the QL
    ! factorization of S
    call tsqmlv(n, p, s, y, dn, .true.)
  end if
  if(global== 1) then
    ierr= 0
    fqq(1:m)= -fqq(1:m)
  end if
  return
  !                    @   SINGULAR case   @
  ! solve ( QHAT-trans QHAT + MU I ) DN= -QHAT-trans W
  ! put the diagonal elements stored in row m+2 of QHAT into their
  ! propre positions and zero out the unwanted portions of QHAT
  20 do j= 1, zero1-1
    qhat(j,j)= qhat(m+2,j)
    qhat(j+1:m+n,j)= zero
  enddo
  do j= zero1, n
    qhat(zero1:m+n,j)= zero
  enddo
  ! compute a small perturbation MU
  call tscpmu(qhat, n, epsm, mu)
  ! form the augmented matrix QHAT by adding an (n*n) diag(MU) in the bottom
  do i= m+1,m+n
    qhat(i,i-m)= mu
  enddo
  ! factor the transformed matrix QHAT from 1 to n
  call tsqrfc(qhat, n, m+n, 1, n, ierr, epsm, pbar, zero1)
  if(method== 1 .and. global== 1) then
    ! premultiply ANLS by the orthogonal matrix resulting from the QR
  ! factorization of QHAT
      call tsqmts(anls, qhat, m+n, m, p, 1, zero1)
  end if
  ! compute the Levenberg-Marquardt step and the residual of the standard model
  if(flag== 1) then
    call tsqmuv(qhat, w, fqq, m+n, 1, n+1, .false.)
    call tsbslv(qhat, m+n, n, fqq, y)
    call tsprmv(dn, y, pbar, n, 0)
    call tsprmv(y, dn, curpos, n, 0)
    call tsqmlv(n, p, s, y, dn, .true.)
    call tssmrd(fqq, resnew, dn, mu, ierr, m, n)
    if(global== 1) then
      ierr= 1
      fqq(1:m+n)= -fqq(1:m+n)
    end if
    return
  else
    call tsqmuv(qhat, fqq, fqt, m+n, 1, n+1, .false.)
    call tsbslv(qhat, m+n, n, fqt, dn)
    call tsprmv(y, dn, pbar, n, 0)
    call tsprmv(dn, y, ptilda, n, 0)
    call tsprmv(y, dn, curpos, n, 0)
    call tsqmlv(n, p, s, y, dn, .true.)
    call tssmrd(fqt, resnew, dn, mu, ierr, m, n)
    if(global== 1) then
      ierr= 1
      fqq(1:m+n)= -fqt(1:m+n)
    end if
  end if
  return
end subroutine tscpss
subroutine tsd1sv(aja, s, anls, fn, x, maxm, m, n, p, q, epsm, pivot, d1)
!SOLVES THE FIRST N-Q LINEAR EQUATIONS IN N-P UNKNOWNS OF THE TENSOR MODEL.
  real(dp),intent(out)     :: aja(:,:)
  real(dp),intent(inout)  :: s(:,:)
  real(dp),intent(inout)  :: anls(:,:)
  real(dp),intent(in)      :: fn(:)
  real(dp),intent(inout)  :: x(:)
  integer, intent(in)        :: maxm
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  integer, intent(in)        :: p
  integer, intent(in)        :: q
  real(dp),intent(in)      :: epsm
  integer, intent(inout)    :: pivot(:)
  real(dp),intent(out)     :: d1(:)
  !INPUT parameterS :
  !    AJA : JACOBIAN MATRIX AT CURRENT ITERATE
  !    S   : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !    ANLS: TENSOR TERM MATRIX AT CURRENT ITERATE
  !    FN  : function value AT CURRENT ITERATE
  !     X  : SOLUTION OF THE LOWER M-N+Q QUADRATIC EQUATIONS IN P
  !          UNKNOWNS OF THE TENSOR MODEL
  !    MAXM: LEADING dimension OF AJA AND ANLS
  !    M,N : dimensionS OF PROBLEM
  !    P   : COLUMN dimension OF S AND ANLS
  !    Q   : NUMERICAL RANK OF JACOBIAN :
  !          Q > P : JACOBIAN IS SINGULAR
  !          Q= P : OTHERWISE
  !    EPSM: MACHINE PRECISION
  !OUTPUT parameterS :
  !     PIVOT : PIVOT VECTOR
  !     D1 : SOLUTION OF THE N-Q LINEAR EQUATIONS IN N-P UNKNOWNS OF
  !          THE TENSOR MODEL
  !SUBprogramS callED:
  !    LEVEL 2 BLAS  ...  DGEMV
  !    TENSOLVE      ...  TSSTMX,TSBSLV,TSQRFC,TSPRMV
  !    TENSOLVE      ...  TSFSLV,TSQMUV
  integer   :: zero1, i, j, ierr, icol
  real(dp) :: wrk1(n), wrk2(n)
  real(dp) :: epsm1
  real(dp), parameter :: alpha= 1.e-04_dp
  ! compute the top right (n-q) x p submatrix of AJA times X
  call dgemv('N', n-q, p, one, aja(:,n-p+1:), maxm, x, 1, zero, d1, 1)
  ! compute S-trans times X
  call tsstmx(s, x, n, p, wrk2)
  ! compute 0.5 * (S-trans times X)**2
  wrk1(1:p)= half * wrk2(1:p)**2
  ! compute 0.5 * (top (n-q) x p submatrix of ANLS) * (S-trans times X)**2
  call dgemv('N', n-q, p, one, anls(:,:), maxm, wrk1, 1, zero, wrk2, 1)
  wrk1(1:n-q)= -fn(1:n-q) - d1(1:n-q) - wrk2(1:n-q)
  ! if the Jacobian is nonsingular then solve for the first
  ! n-p components of the tensor step and return
  if(p== q) then
    call tsbslv(aja, m, n-p, wrk1, d1)
    return
  end if
  wrk2(n-q+1:n-p)= zero
  ! copy top left (n-q) x (n-p) submatrix of AJA into bottom of AJA
  aja(m+3:m+n-q+2,1:n-p)= aja(1:n-q,1:n-p)
  ! copy the transpose of the top left (n-q) x (n-p) submatrix of AJA
  ! into top of AJA
  do j= 1,n-q
    aja(j,j)= aja(m+2,j)
    do i= j+1,n-p
      aja(i,j)= aja(j,i)
    enddo
  enddo
  ! zero out the upper triangular (n-q) x (n-q) triangular part of
  ! the transpose of the top left (n-q) x (n-p) submatrix of AJA
  do j= 1,n-q
    aja(1:j-1,j)= zero
  enddo
  ! factor the transpose of the top left (n-q) x (n-p) submatrix of AJA
  epsm1= epsm*alpha
  call tsqrfc(aja, n-q, n-p, 1, n-q, ierr, epsm1, pivot, zero1)
  if(ierr== 0) then
    icol= n-q
  else
    icol= zero1-1
  end if
  call tsprmv(d1, wrk1, pivot, n-q, 0)
  ! solve for the first n-p components of the tensor step
  call tsfslv(aja, d1, n-p, icol, wrk2)
  call tsqmuv(aja, wrk2, d1, n-p,1, zero1, .true.)
  ! copy the (n-q) x (n-p) submatrix of AJA from bottom back to top of AJA
  aja(1:n-q,1:n-p)= aja(m+3:m+n-q+2,1:n-p)
  return
end subroutine tsd1sv

subroutine tsdfcn(p, x, g)
!COMPUTES THE ANALYTIC GRADIENT OF THE function GIVEN BY subroutine TSQFCN.
  integer, intent(in)     :: p
  real(dp),intent(in)   :: x(:)
  real(dp),intent(out)  :: g(:)
  !INPUT parameterS :
  !    P    : COLUMN dimension OF ANLS AND S
  !    X    : POINT AT WHICH GRADIENT IS EVALUATED
  !OUTPUT parameterS :
  !    G : GRADIENT AT X
  !SUBprogramS callED:
  !    LEVEL 2 BLAS  ...  DGEMV
  !    TENSOLVE      ...  TSSTMX
  integer :: i, j, l
  real(dp):: wrk1(meqns), wrk2(p), wrk3(p), wrk4(meqns), wrk5(meqns)
  !
  ! compute the lower right (m-n+q) x p submatrix of AJA times X
  call dgemv('N', meqns-nvars+qrank, p, one, aja(nvars-qrank+1:,nvars-p+1:), &
             meqns, x, 1, zero, wrk1, 1)
  ! compute S-trans times X
  call tsstmx(s, x, nvars, p, wrk3)
  ! compute 0.5 * (S-trans times X)**2
  wrk2(1:p)= half * wrk3(1:p)**2
  ! compute 0.5 * (lower (m-n+q) x p submatrix of ANLS) *
  ! (S-trans times X)**2
  call dgemv('N', meqns-nvars+qrank, p, one, anls(nvars-qrank+1:,:), meqns,  &
             wrk2, 1, zero, wrk4, 1)
  do i= 1,meqns-nvars+qrank
    wrk4(i)= wrk4(i) + fc(nvars-qrank+i) + wrk1(i)
  enddo
  ! compute AJA-trans * WRK4
  call dgemv('T', meqns-nvars+qrank, p, one, aja(nvars-qrank+1:,nvars-p+1:), &
             meqns, wrk4, 1, zero, wrk1, 1)
  ! compute ANLS-trans * WRK4
  call dgemv('T', meqns-nvars+qrank, p, one, anls(nvars-qrank+1:,:), meqns,  &
             wrk4, 1, zero, wrk5, 1)
  ! compute S * diag(S-trans * WRK3) * WRK5
  wrk2(1:p)= zero
  l= p+1
  do j= 1,p
    l= l-1
    wrk2(l)= s(nvars+2,l)
    do i= l+1,p
      wrk2(i)= s(nvars-p+j,i)
    enddo
    wrk2(1:p)= wrk2(1:p)*wrk3(1:p)
    g(j)=      dot_product( wrk2(1:p), wrk5(1:p) )
  enddo
  g(1:p)= g(1:p) + wrk1(1:p)
  return
end subroutine tsdfcn

subroutine tsdflt(&
& m, n, itnlim, jacflg, gradtl, steptl, ftol, method,  &
& global, stepmx, dlt, typx, typf, ipr)
! SETS default valueS FOR EACH INPUT VARIABLE TO THE NONLINEAR EQUATION ALGORITHM.
! N.B. Argument MSG has been removed.
  integer, intent(in) :: m
  integer, intent(in) :: n
  integer, intent(out):: itnlim
  integer, intent(out):: jacflg
  real(dp),intent(out):: gradtl
  real(dp),intent(out):: steptl
  real(dp),intent(out):: ftol
  integer, intent(out):: method
  integer, intent(out):: global
  real(dp),intent(out):: stepmx
  real(dp),intent(out):: dlt
  real(dp),intent(inout):: typx(:)
  real(dp),intent(inout):: typf(:)
  integer, intent(out):: ipr
  !
  real(dp):: eps
  real(dp),parameter :: thous= 1000._dp
  !
  jacflg= 0
  eps=    EPSILON(one)
  gradtl= eps**(one/three)
  steptl= eps**(two/three)
  ftol=   eps**(two/three)
  itnlim= 150
  method= 1
  global= 0
  stepmx= thous
  dlt= -one
  !msg= 0
  ipr= 6
  typx(1:n)= one
  typf(1:m)= one
  return
end subroutine tsdflt

subroutine tsdumj(x, aja, m, n)
!DUMMY ROUTINE TO PREVENT UNSATISFIED external DIAGNOSTIC
!WHEN SPECifIC ANALYTIC JACOBIAN IS NOT SUPPLIED.
  real(dp),intent(in)   :: x(:)
  real(dp),intent(out)  :: aja(:,:)
  integer, intent(in)     :: m, n
  !INPUT parameterS:
  !   X   : POINT AT WHICH JACOBIAN IS EVALUATED
  !   AJA : JACOBIAN MATRIX
  aja(m,n)= x(1)
  return
end subroutine tsdumj

subroutine d2fcn(nr, n, x, h)
!THIS IS A DUMMY ROUTINE TO PREVENT UNSATISFIED external DIAGNOSTIC WHEN
!A REFERENCE TO THE MATRIX OF SECOND DERIVATIVES IS passED TO UNCMIN.
  integer, intent(in)    :: nr
  integer, intent(in)    :: n
  real(dp),intent(in)  :: x(:)
  real(dp),intent(out) :: h(:,:)
  h(nr,n)= x(1)
  return
end subroutine d2fcn

subroutine tsfafa( &
& anls, fq, adt, ag, const1, const2, alpha, dlt, m, n, p,  &
& nwtake, ierr, vn, fn_val)
!COMPUTES || F + J*D + 0.5*A*D**2 ||**2 IN THE L2 NORM SENS,
!where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
!N.B. Changed to a subroutine by AJM
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: fq(:)
  real(dp),intent(in)   :: adt(:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: alpha
  real(dp),intent(in)   :: dlt
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  logical, intent(in)     :: nwtake
  integer, intent(in)     :: ierr
  real(dp),intent(out)  :: vn(:)
  real(dp),intent(out)  :: fn_val
  !INPUT parameterS
  !    ANLS   : TENSOR TERM MATRIX
  !    FQ     : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES
  !    ADT    : JACOBIAN MATRIX TIMES DT
  !     AG    : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1 : SHAT-TRANS TIMES DT
  !    CONST2 : SHAT-TRANS TIMES GBAR
  !    ALPHA  : POINT AT WHICH TO EVALUATE THE subroutine TSFAFA
  !    DLT    : CURRENT TRUST RADIUS
  !    M,N    : dimensionS OF THE PROBLEM
  !    P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !    NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !             NWTAKE= .true.  : STANDARD STEP TAKEN
  !             NWTAKE= .false. : TENSOR STEP TAKEN
  !    IERR   : return CODE FROM QRP FACTORIZATION ROUTINE:
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS
  !    VN     : F + J*D + 0.5*A*D**2
  !    TSFAFA :  || F + J*D + 0.5*A*D**2 ||**2
  !             where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2)
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSMAFA
  integer              :: len
  call tsmafa( &
  & anls, fq, adt, ag, const1, const2, alpha, dlt, m, n, p,  &
  & nwtake, ierr, vn)
  len= m
  if(ierr > 0) len= m + n
  fn_val= half * SUM( vn(1:len)**2 )
  return
end subroutine tsfafa

subroutine tsfdfj(xc, fc, m, n, epsm, fvec, aja)
!COMPUTES THE FINITE DifFERENCE JACOBIAN AT THE CURRENT ITERATE XC.
  real(dp),intent(inout)  :: xc(:)
  real(dp),intent(in)      :: fc(:)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(in)      :: epsm
  real(dp),intent(out)     :: aja(:,:)
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
  end subroutine fvec
  end interface
  !INPUT parameterS :
  !    XC   : CURRENT ITERATE
  !    FC   : function value AT XC
  !    M,N  : dimensionS OF PROBLEM
  !    EPSM : MACHINE PRECISION
  !    FVEC : subroutine TO EVALUATE THE useR'S function
  !OUTPUT parameterS :
  !    AJA : FINITE DifFERENCE JACOBIAN AT XC
  !    SUBprogramS callED:
  !    useR   ...  FVEC
  real(dp):: fhat(m)
  integer :: j
  real(dp):: ndigit, rnoise, stepsz, xtmpj, sqrtr, rstpsz
  ndigit= -LOG10(epsm)
  rnoise= MAX(ten**(-ndigit),epsm)
  sqrtr= SQRT(rnoise)
  do j= 1,n
    stepsz= sqrtr*MAX(ABS(xc(j)),one)
    xtmpj= xc(j)
    xc(j)= xtmpj+stepsz
    call fvec(xc, fhat, m, n)
    xc(j)= xtmpj
    rstpsz= one/stepsz
    aja(1:m,j)= (fhat(1:m) - fc(1:m))*rstpsz
  enddo
  return
end subroutine tsfdfj

subroutine tsfrmt(shat, s, aja, fv, fn, m, n, p, idp, a)
!FORMS THE TENSOR TERM MATRIX OF THE TENSOR MODEL.
  real(dp),intent(inout)  :: shat(:,:)
  real(dp),intent(in)      :: s(:,:)
  real(dp),intent(in)      :: aja(:,:)
  real(dp),intent(in)      :: fv(:,:)
  real(dp),intent(in)      :: fn(:)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  integer, intent(in)        :: p
  integer, intent(in)        :: idp(:)
  real(dp),intent(out)     :: a(:,:)
  !INPUT parameterS :
  !    SHAT: MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !    S   : MATRIX OF PREVIOUS DIRECTIONS
  !    AJA : JACOBIAN MATRIX AT CURRENT ITERATE
  !    FV  : MATRIX OF PAST function valueS
  !    FN  : function value AT CURRENT ITERATE
  !    M   : ROW dimension OF MATRICES A,FV,AND AJA
  !    N   : COLUMN dimension OF JACOBIAN MATRIX
  !    P   : COLUMN dimension OF MATRIX SHAT
  !    IDP : VECTOR WHICH KEEPS TRACK OF LINEARLY INDEPendENT
  !          DIRECTION POSITIONS WITHIN THE MATRIX S
  !OUTPUT parameterS :
  !    A   : TENSOR TERM MATRIX
  !    SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  !    UNCMIN        ...  CHOLDC,LLTSLV
  real(dp) :: am(p,p), x(p), b(p), scale(p)
  integer   :: i, j, jj
  real(dp) :: sum, sc, tol, diagmx
  ! scale the matrix SHAT and save scaling in SCALE
  do j= 1,p
    sc= one/dnrm2(n, shat(:,j), 1)
    shat(1:n,j)= sc * shat(1:n,j)
    scale(j)= sc**2
  enddo
  ! form the matrix AM= (Si Sj)**2
  do j= 1,p
    jj= idp(j)
    do i= 1,p
      am(i,j)= dot_product( s(1:n,idp(i)), s(1:n,jj) )**2
    enddo
  enddo
  ! scale the matrix AM
  do i= 1,p
    do j= 1,p
      am(i,j)= scale(i)*scale(j)*am(i,j)
    enddo
  enddo
  ! perform a Cholesky decomposition of AM
  tol= zero
  diagmx= zero
  call choldc(p, am, diagmx, tol)
  ! form the tensor term matrix A
  do i= 1,m
    do j= 1,p
      jj= idp(j)
      sum= dot_product( aja(i,1:n), s(1:n,jj) )
      b(j)= two*(fv(i,jj) - fn(i) - sum)
      b(j)= scale(j)*b(j)
    enddo
    ! solve AM*X= B
      call lltslv(p, am, x, b)
    ! copy X into row i of A
      a(i,1:p)= x(1:p)
    enddo
  return
end subroutine tsfrmt

subroutine tsfscl(x, dx, df, m, n, fvec, f)
!EVALUATES THE function AT THE CURRENT ITERATE X then SCALES ITS value.
  real(dp),intent(inout)  :: x(:)
  real(dp),intent(in)      :: dx(:)
  real(dp),intent(in)      :: df(:)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(out)     :: f(:)
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
  end subroutine fvec
  end interface
  !INPUT parameterS :
  !     X  : CURRENT ITERATE
  !     DX : DIAGONAL SCALING MATRIX FOR X
  !     DF : DIAGONAL SCALING MATRIX FOR F
  !    M,N : dimensionS OF PROBLEM
  !   FVEC : subroutine TO EVALUATE function
  !OUTPUT parameterS :
  !     F  : SCALED function value AT CURRENT ITERATE X
  !SUBprogramS callED:
  !     TENSOLVE      ...  TSUNSX,TSSCLF,TSSCLX
  !     useR          ...  FVEC
  call tsunsx(x, dx, n)
  call fvec(x,f, m, n)
  call tssclf(f, df, m)
  call tssclx(x, dx, n)
  return
end subroutine tsfscl

subroutine tsfslv(l, b, m, n, y)
!doES A FORWARD SOLVE.
  real(dp),intent(in) :: l(:,:)
  real(dp),intent(in) :: b(:)
  integer, intent(in) :: m
  integer, intent(in) :: n
  real(dp),intent(out):: y(:)
  !INPUT parameterS :
  !    L   : THE TRANSPOSE OF THE UPPER TRIANGULAR MATRIX OBTAINED
  !          FROM A QR FACTORIZATION OF AN M BY N MATRIX A. DIAG(L)
  !          IS STORED IN ROW M+2. THIS IS THE STORAGE SCHEME useD
  !          IN STEWART, G. W., III(1973) "INTRODUCTION TO MATRIX
  !          COMPUTATION", ACADEMIC PRESS,NEW YORK
  !    B   : RIGHT HAND SIDE
  !     M  : ROW dimension OF MATRIX A
  !     N  : COLUMN dimension OF MATRIX A
  !OUTPUT parameterS :
  !     Y  : VECTOR SOLUTION ON exit
  integer   :: j
  real(dp) :: s
  ! solve L Y= B
  y(1)= b(1) / l(m+2,1)
  if(n > 1) then
    s= l(1,2) * y(1)
    y(2)= (b(2) - s) / l(m+2,2)
    do j= 3,n
      s= dot_product( l(1:j-1,j), y(1:j-1) )
      y(j)= (b(j) - s) / l(m+2,j)
    enddo
  end if
  return
end subroutine tsfslv

subroutine tsjmuv( &
& itn, method, v, curpos, pivot, pbar, aja, shat,  &
& flag, ierr, m, n, p, vn, av)
!CALCULATES THE PRODUCT JACOBIAN TIMES A VECTOR.
  integer, intent(in)     :: itn
  integer, intent(in)     :: method
  real(dp),intent(in)   :: v(:)
  integer, intent(in)     :: curpos(:)
  integer, intent(in)     :: pivot(:)
  integer, intent(in)     :: pbar(:)
  real(dp),intent(in)   :: aja(:,:)
  real(dp),intent(in)   :: shat(:,:)
  integer, intent(in)     :: flag
  integer, intent(in)     :: ierr
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  real(dp),intent(out)  :: vn(:)
  real(dp),intent(out)  :: av(:)
  !INPUT parameterS
  !     ITN    : CURRENT ITERATION NUMBER
  !     METHOD : METHOD TO BE useD
  !     V      : VECTOR TO BE MULTIPLIED BY AJA
  !     CURPOS : PIVOT VECTOR (useD DURING THE FACTORIZATION OF AJA
  !              FROM COLUMN 1 TO N-P)
  !     PIVOT  : PIVOT VECTOR (useD DURING THE FACTORIZATION OF AJA
  !              FROM COLUMN N-P+1 TO N)
  !     PBAR   : PIVOT VECTOR (useD DURING THE FACTORIZATION OF AJA
  !              if IT IS SINGULAR
  !     AJA    : JACOBIAN MATRIX AT CURRENT ITERATE
  !     SHAT   : MATRIX OF LINEARLY INDEPendENT DIRECTIONS AFTER
  !              A QL FACTORIZATION
  !     FLAG   : return CODE WITH THE FOLLOWING MEANINGS:
  !             FLAG= 0 : NO SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN 1 TO N
  !             FLAG= 1 : SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN 1 TO N-P
  !             FLAG= 2 : SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN N-P+1 TO N
  !     IERR   : return CODE FROM QRP FACTORIZATION ROUTINE:
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !     M,N    : dimensionS OF THE PROBLEM
  !     P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !OUTPUT parameterS
  !     VN     : ? previously erroneously described as workspace
  !     AV     : JACOBIAN TIMES V
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSPRMV,TSQMLV,TSUTMD
  real(dp) :: wrk1(n), wrk2(n)
  integer :: len
  if(itn== 1 .or. method== 0) then
    call tsprmv(wrk1, v, pivot, n, 1)
    if(ierr== 1) then
      call tsprmv(wrk2, wrk1, pbar, n, 1)
      wrk1= wrk2
    end if
  else if(n== 1) then
    vn(1)= v(1)
  else
    call tsqmlv(n, p, shat, v, vn, .false.)
    call tsprmv(wrk2, vn, curpos, n, 1)
    if(flag== 0) then
      call tsprmv(wrk1, wrk2, pivot, n, 1)
    else if(flag== 1) then
      call tsprmv(wrk1, wrk2, pbar, n, 1)
    else if(flag== 2 ) then
      call tsprmv(wrk1, wrk2, pivot, n, 1)
      call tsprmv(wrk2, wrk1, pbar, n, 1)
      wrk1= wrk2
    end if
  end if
  len= m
  if(ierr > 0) len= m + n
  call tsutmd(aja, wrk1, len, n, av)
  return
end subroutine tsjmuv

subroutine tsjqtp(q, n, m, p, aja)
!GETS J*(Q-TRANS) BY COMPUTING EACH ROW OF THE resultING MATRIX AS FOLLOWS:
!(J*Q-TRANS)I-TH ROW<--Q*(J)I-TH ROW.
  real(dp),intent(in)      :: q(:,:)
  integer, intent(in)        :: n
  integer, intent(in)        :: m
  integer, intent(in)        :: p
  real(dp),intent(inout)  :: aja(:,:)
  !INPUT parameterS :
  !    Q    : resultING MATRIX FROM A QL FACTORIZATION
  !    M,N  : dimensionS OF PROBLEM
  !    P    : COLUMN dimension OF MATRIX Q
  !INPUT-OUTPUT parameterS :
  !    AJA : JACOBIAN MATRIX ON entry AND JACOBIAN MULTIPLIED BY THE
  !          ORTHOGONAL MATRIX Q ON exit
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSQMLV
  real(dp) :: wrk1(n), wrk2(n)
  integer   :: i
  do i= 1,m
    ! copy the i-th row of AJA into WRK1
      wrk1= aja(i,1:n)
      call tsqmlv(n, p, q, wrk1, wrk2, .false.)
    ! form the i-th row of AJA*(Q-trans)
      aja(i,1:n)= wrk2
    enddo
  return
end subroutine tsjqtp

subroutine tslmin( &
& xc, xp, p1, q, anls, fq, adt, ag, const1, const2,  &
& dlt, m, n, p, nwtake, ierr, tol, xplus)
!FINDS A LOCAL MINIMIZER OF A ONE-VARIABLE function IN AN INTERVAL [XC XP].
  real(dp),intent(inout)  :: xc
  real(dp),intent(in)      :: xp
  real(dp),intent(inout)  :: p1
  real(dp),intent(out)     :: q
  real(dp),intent(in)      :: anls(:,:)
  real(dp),intent(in)      :: fq(:)
  real(dp),intent(in)      :: adt(:)
  real(dp),intent(in)      :: ag(:)
  real(dp),intent(in)      :: const1(:)
  real(dp),intent(in)      :: const2(:)
  real(dp),intent(in)      :: dlt
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  integer, intent(in)        :: p
  logical, intent(in)        :: nwtake
  integer, intent(in)        :: ierr
  real(dp),intent(in)      :: tol
  real(dp),intent(out)     :: xplus
  !INPUT parameterS :
  !    XC,XP  : LOWER AND UPPER BOUND OF INTERVAL IN WHICH THE SEARCH IS PERFORMED
  !    P1,Q   : FIRST DERIVATIVES OF THE ONE-VARIABLE function
  !    ANLS   : TENSOR TERM MATRIX
  !    FQ     : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES
  !    ADT    : JACOBIAN TIMES THE STEP DT (SEE subroutine TS2DTR)
  !    AG     : JACOBIAN TIMES THE GRADIENT G (SEE subroutine TS2DTR)
  !    CONST1 : SHAT-TRANS * DT  (SEE subroutine TS2DTR)
  !    CONST2 : SHAT-TRANS * GBAR (SEE subroutine TS2DTR)
  !    DLT    : TRUST RADIUS
  !    M,N    : dimensionS OF PROBLEM
  !    P      : COLUMN dimension OF MATRIX ANLS
  !    NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !             NWTAKE= .true.  : STANDARD STEP TAKEN
  !             NWTAKE= .false. : TENSOR STEP TAKEN
  !    IERR   : return CODE FROM QRP FACTORIZATION ROUTINE:
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : OTHERWISE
  !    TOL    : SMALL TOLERANCE
  !OUTPUT parameterS :
  !    XPLUS  :  LOCAL MINIMIZER OF THE ONE-VARIABLE function
  !SUBprogramS callED :
  !    TENSOLVE      ...  TSMSDA,TSFAFA,TSLMSP,TSMFDA
  real(dp) :: vn(n+m)
  integer              :: itercd, retcd, itncnt
  real(dp)            :: aleft, aright, t, e, s, sinit, tmp
  real(dp), parameter :: ott= 1.0e-04_dp, small= 2.0e-20_dp
  logical              :: skip
  retcd= 0
  aleft= MIN(xc,xp)
  aright= MAX(xc,xp)
  itncnt= 0
  t= ABS(xc-xp)
  skip= .false.
  ! compute the second derivative value at the current point
  call tsmsda( &
  & anls, fq, adt, ag, const1, const2, xc, dlt, m, n, p,  &
  & nwtake, ierr, skip, vn, e)
  10 if(e > zero) then
    s= -p1/e
    if(ABS(s) > two*t) then
      if (s < zero) then
        s= -two*t
      else
        s= two*t
      end if
    end if
  else
    if (p1 > zero) then
      s= -t
    else
      s= t
    end if
  end if
  if(xc+s > aright) s= aright - xc
  if(xc+s < aleft)  s= aleft - xc
  sinit= ABS(s)
  ! compute a next iterate XPLUS
  20 call tsfafa( &
  & anls, fq, adt, ag, const1, const2, xc+s, dlt, m, n, p,  &
  & nwtake, ierr, vn, tmp)
  if (tmp > q + ott*s*p1) then
    s= s/2
    if(ABS(s) < small*sinit .or. s== zero) then
      retcd= 1
    else
      GO TO 20
    end if
  end if
  xplus= xc + s
  itncnt= itncnt + 1
  ! check stopping criteria
  call tslmsp( &
  & xc, xplus, itncnt, retcd, itercd, anls, adt, ag,  &
  & const1, const2, dlt, m, n, p, nwtake, ierr, tol, vn)
  if(itercd > 0) return
  ! update XC
  xc= xplus
  ! compute function and derivative values at the new point
  call tsfafa(anls, fq, adt, ag, const1, const2, xc, dlt, m, n, p, nwtake, &
              ierr, vn, q)
  p1= tsmfda(anls, adt, ag, const1, const2, xc, dlt, m, n, p, nwtake, &
              ierr, vn)
  skip= .true.
  call tsmsda(anls, fq, adt, ag, const1, const2, xc, dlt, m, n, p,  &
              nwtake, ierr, skip, vn, e)
  GO TO 10
end subroutine tslmin

subroutine tslmsp( &
!CHECKS THE stopPING CRITERIA FOR A LOCAL MINIMIZER.
& xc, xp, itncnt, retcd, itercd, anls, adt, ag, const1,  &
& const2, dlt, m, n, p, nwtake, ierr, tol, vn)
  real(dp),intent(in)  :: xc
  real(dp),intent(in)  :: xp
  integer, intent(in)    :: itncnt
  integer, intent(in)    :: retcd
  integer, intent(out)   :: itercd
  real(dp),intent(in)  :: anls(:,:)
  real(dp),intent(in)  :: adt(:)
  real(dp),intent(in)  :: ag(:)
  real(dp),intent(in)  :: const1(:)
  real(dp),intent(in)  :: const2(:)
  real(dp),intent(in)  :: dlt
  integer, intent(in)    :: m
  integer, intent(in)    :: n
  integer, intent(in)    :: p
  logical, intent(in)    :: nwtake
  integer, intent(in)    :: ierr
  real(dp),intent(in)  :: tol
  real(dp),intent(in)  :: vn(:)
  !INPUT parameterS :
  !    XC       : CURRENT ITERATE (FROM SEARCH subroutine)
  !    XP       : NEXT ITERATE (FROM SEARCH subroutine)
  !    ITNCNT   : ITERATION LIMIT
  !    RETCD    : return CODE FROM LINE SEARCH
  !    DLT      : TRUST RADIUS
  !    AJA      : JACOBIAN AT THE CURRENT ITERATE
  !    NR       : LEADING dimension OF THE JACOBIAN MATRIX
  !    M,N      : dimensionS OF THE PROBLEM
  !    P        : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !    NWTAKE   : logical VARIABLE WITH THE FOLLOWING MEANINGS :
  !               NWTAKE= .true.  : STANDARD STEP TAKEN
  !               NWTAKE= .false. : TENSOR STEP TAKEN
  !    IERR     : return CODE FROM THE QRP FACTORIZATION ROUTINE :
  !               IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !               IERR= 1 : OTHERWISE
  !    VN       : WORKING VECTOR
  !    TOL      : SMALL TOLERANCE
  !    METHOD   : METHOD TO use
  !            = 0   : STANDARD METHOD useD
  !            = 1   : TENSOR METHOD useD
  !OUTPUT parameterS :
  !    ITERCD  : return CODE WITH FOLLOWING MEANINGS :
  !              ITERCD= 1 : FIRST DERIVATIVE AT THE CURRENT POINT
  !                           close TO 0
  !              ITERCD= 2 : SUCCESSIVE ITERATES WITHIN TOLERANCE
  !              ITERCD= 3 : LINE SEARCH FAILED TO LOCATE A POINT
  !                           LOWER THAT THE CURRENT POINT
  !              ITERCD= 4 : ITERATION LIMIT EXCEEDED
  real(dp) :: grdt
  grdt= SQRT(tol)
  itercd= 0
  if(retcd== 1) then
    itercd= 3
  else if(ABS(tsmfda( &
  & anls, adt, ag, const1, const2, xp, dlt,  &
  & m, n, p, nwtake, ierr, vn)) < grdt) then
    itercd= 1
  else if(xp /= zero .and. ABS(xp-xc)/ABS(xp) <= tol) then
    itercd= 2
  else if(itncnt >= 150) then
    itercd= 4
  end if
  return
end subroutine tslmsp

subroutine tslsch( &
& m, n, xc, d, g, steptl, dx, df, fvec,  &
& mxtake, stepmx, xp, fp, fcnorm, fpnorm, retcd)
!FINDS A NEXT ITERATE USING A STANDARD LINE SEARCH METHOD.
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(in)      :: xc(:)
  real(dp),intent(inout)  :: d(:)
  real(dp),intent(in)      :: g(:)
  real(dp),intent(in)      :: steptl
  real(dp),intent(in)      :: dx(:)
  real(dp),intent(in)      :: df(:)
  logical, intent(out)       :: mxtake
  real(dp),intent(in)      :: stepmx
  real(dp),intent(out)     :: xp(:)
  real(dp),intent(out)     :: fp(:)
  real(dp),intent(in)      :: fcnorm
  real(dp),intent(out)     :: fpnorm
  integer, intent(out)       :: retcd
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
  end subroutine fvec
  end interface
  !INPUT parameterS :
  !       M,N : dimensionS OF PROBLEM
  !       XC  : CURRENT ITERATE
  !       D   : SEARCH DIRECTION
  !       G   : GRADIENT AT CURRENT ITERATE
  !    STEPTL : RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                ARE CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  !       DX  : DIAGONAL SCALING MATRIX FOR X
  !       DF  : DIAGONAL SCALING MATRIX FOR F
  !       FVEC: subroutine TO EVALUATE THE function
  !     STEPMX: MAXIMUM ALLOWABLE STEP size
  !OUTPUT parameterS :
  !    MXTAKE: BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  !       XP : NEXT ITARATE
  !       FP : function value AT NEXT ITERATE
  !   FCNORM : 0.5 * || F(XC) ||**2
  !   FPNORM : 0.5 * || F(XP) ||**2
  !    RETCD : return CODE WITH THE FOLLOWING MEANING :
  !                RETCD= 0 : SATISFACTORY LOCATION OF A NEW ITERATE
  !                RETCD= 1 : NO SATISFACTORY POINT FOUND SUFFICIENTLY
  !                            DISTINCT FROM X
  !SUBprogramS callED:
  !       LEVEL 1 BLAS  ...  DNRM2
  !       TENSOLVE      ...  TSFSCL
  !       useR          ...  FVEC
  integer   :: i
  real(dp) :: slope, releng, temp1, temp2, almda, temp, almdat, almdam, sln, scl
  real(dp), parameter :: alpha= 1.0D-4, tenth= 0.1_dp, z99= 0.99_dp
  mxtake= .false.
  sln= dnrm2(n, d, 1)
  if(sln > stepmx) then
    ! step longer than maximum allowed
      scl= stepmx/sln
    d(1:n)= scl * d(1:n)
    sln= stepmx
  end if
  ! compute SLOPE =  G-trans * D
  slope= dot_product( g(1:n), d(1:n) )
  ! initialization of RETCD
  retcd= 0
  ! compute the smallest value allowable for the damping
  ! parameter ALMDA, i.e, ALMDAM
  releng= zero
  do i= 1,n
    temp1= MAX(ABS(xc(i)), one)
    temp2= ABS(d(i))/temp1
    releng= MAX(releng,temp2)
  enddo
  almdam= steptl/releng
  almda= one
  ! compute the next iterate XP
  40 xp(1:n)= xc(1:n) + almda*d(1:n)
  ! evaluate the objective function at XP and its residual
  call tsfscl(xp, dx, df, m, n, fvec, fp)
  fpnorm= half*dnrm2(m, fp, 1)**2
  ! test whether the full step produces enough decrease in the l2 norm of
  ! the objective function.  If not update ALMDA and compute a new step
  if (fpnorm > (fcnorm + (alpha* almda * slope))) then
    almdat= ((-almda**2)*slope) / (two*(fpnorm - fcnorm - almda*slope))
    temp= almda/ten
    almda= MAX(temp,almdat)
    if(almda < almdam) then
      retcd= 1
      return
    end if
    GO TO 40
  else
    if(almda== tenth .and. sln > z99*stepmx) mxtake=.true.
  end if
  return
end subroutine tslsch

subroutine tsmafa( &
& anls, f, adt, ag, const1, const2, alpha, dlt,  &
& m, n, p, nwtake, ierr, vn)
!COMPUTES THE VECTOR VN= F(XC) + J(XC)*D + 0.5*A*D**2,
!where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: f(:)
  real(dp),intent(in)   :: adt(:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: alpha
  real(dp),intent(in)   :: dlt
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  logical, intent(in)     :: nwtake
  integer, intent(in)     :: ierr
  real(dp),intent(out)  :: vn(:)
  !INPUT parameterS :
  !    ANLS  : TENSOR TERM MATRIX
  !     ADT  : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !      AG  : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1: SHAT-TRANS * DT (SEE subroutine TS2DTR)
  !    CONST2: SHAT-TRABS * GBAR (SEE subroutine TS2DTR)
  !    ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
  !      DLT : CURRENT TRUST RADIUS
  !      M,N : dimensionS OF THE PROBLEM
  !      P   : COLUMN dimension OF THE MATRIX ANLS
  !   NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS
  !               NWTAKE= .true.  : STANDARD STEP TAKEN
  !               NWTAKE= .false. : TENSOR STEP TAKEN
  !     IERR : return CODE FROM THE QRP FACTORIZATION ROUTINE :
  !            IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !            IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS :
  !     VN  : F + J*D + 0.5*A*D**2, where
  !           D= ALPHA*DT + SQRT(DLT**2-ALPHA**2)
  integer   :: i, j, len
  real(dp) :: expr, const
  expr= SQRT(dlt**2 - alpha**2)
  vn(1:n)= alpha*adt(1:n) + expr*ag(1:n)
  vn(n+1:n+m)= zero
  len= m
  if(ierr > 0) len= m + n
  do i= 1, len
    vn(i)= vn(i) + f(i)
  enddo
  if(nwtake) return
  do j= 1,p
    const= half*(alpha*const1(j) + expr*const2(j))**2
    vn(1:len)= vn(1:len) + const * anls(1:len,j)
  enddo
  return
end subroutine tsmafa
subroutine tsmdls( &
& aja, shat, anls, xc, m, n, p, dt, g, dx, df, &
& fvec, method, steptl, global, stepmx, epsm, fq, dn, fqq, &
& pivot, curpos, pbar, mxtake, xp, fp, fcnorm, fpnorm,   &
& zero1, retcd, ierr)
!FINDS A NEXT ITERATE USING A LINE SEARCH METHOD.
!IT TRIES THE FULL TENSOR STEP FIRST. if THIS IS NOT SUCCESSFUL then
!IT COMPUTES THE STANDARD DIRECTION AND COMPUTES A STEP IN THAT
!DIRECTION. NEXT, if THE TENSOR DIRECTION IS DESCENT, IT COMPUTES
!A STEP IN THE TENSOR DIRECTION.  THE ITERATE THAT PRODUCES
!THE LOWER RESIDUAL IS THE NEXT ITERATE FOR THE NONLINEAR ALGORITHM.
  real(dp),intent(inout)  :: aja(:,:)
  real(dp),intent(in)      :: shat(:,:)
  real(dp),intent(inout)  :: anls(:,:)
  real(dp),intent(in)      :: xc(:)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  integer, intent(in)        :: p
  real(dp),intent(inout)  :: dt(:)
  real(dp),intent(in)      :: g(:)
  real(dp),intent(in)      :: dx(:)
  real(dp),intent(in)      :: df(:)
  integer, intent(in)        :: method
  real(dp),intent(in)      :: steptl
  integer, intent(in)        :: global
  real(dp),intent(in)      :: stepmx
  real(dp),intent(in)      :: epsm
  real(dp),intent(in)      :: fq(:)
  real(dp),intent(out)     :: dn(:)
  real(dp),intent(out)     :: fqq(:)
  integer, intent(out)       :: pivot(:)
  integer, intent(in)        :: curpos(:)
  integer, intent(out)       :: pbar(:)
  logical, intent(out)       :: mxtake
  real(dp),intent(out)     :: xp(:)
  real(dp),intent(out)     :: fp(:)
  real(dp),intent(in)      :: fcnorm
  real(dp),intent(out)     :: fpnorm
  integer, intent(out)       :: zero1
  integer, intent(out)       :: retcd
  integer, intent(inout)    :: ierr
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
    end subroutine fvec
  end interface
  !INPUT parameterS
  !    AJA    : JACOBIAN AT CURRENT ITERATE
  !    SHAT   : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !             AFTER A QL FACORIZATION
  !    ANLS   : TENSOR TERM MATRIX
  !    XC     : CURRENT ITERATE
  !    M,N    : dimensionS OF THE PROBLEM
  !    P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !    DT     : TENSOR STEP
  !    G      : GRADIENT AT CURRENT ITERATE
  !    DX     : DIAGONAL SCALING MATRIX FOR X
  !    DF     : DIAGONAL SCALING MATRIX FOR F
  !    GBAR   : STEEPEST DESCENT DIRECTION (= -G)
  !    METHOD : METHOD TO use
  !            = 0  : STANDARD METHOD useD
  !            = 1  : TENSOR METHOD useD
  !    STEPTL : STEP TOLERANCE
  !    GLOBAL : GLOBAL STRATEGY useD
  !               =  0 : LINE SEARCH IS useD
  !               =  1 : 2-dimensionAL TRUST REGION IS useD
  !    STEPMX : MAXIMUM ALLOWABLE STEP size
  !    EPSM   : MACHINE PRECISION
  !    FQ     : function value AT CURRENT ITERATE MULTIPLIED BY AN
  !             ORTHOGOL MATRIX
  !OUTPUT parameterS
  !    DN     : NEWTON STEP
  !    FQQ    : FQ MULTIPLIED BY AN ORTHOGONAL MATRIX
  !    CURPOS : PIVOT VECTOR (useD DURING THE FACTORIZATION OF THE
  !             JACOBIAN FROM COLUMN 1 TO N-P)
  !    PIVOT  : PIVOT VECTOR (useD DURING THE FACTORIZATION OF THE
  !             JACOBIAN FROM COLUMN N-P+1 TO N)
  !    PBAR   : PIVOT VECTOR (useD DURING THE FACTORIZATION OF THE
  !             JACOBIAN if IT IS SINGULAR
  !    MXTAKE : BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  !    XP     : NEXT ITERATE
  !    FP     : function value AT NEXT ITERATE
  !    FCNORM :  0.5 * || F(XC) ||**2
  !    FPNORM :  0.5 * || F(XP) ||**2
  !    ZERO1  : FIRST ZERO COLUMN OF THE JACOBIAN IN case OF SINGULARITY
  !    RETCD  : return CODE WITH THE FOLLOWING MEANING :
  !             RETCD =  0 : SATISFACTORY LOCATION OF A NEW ITERATE
  !             RETCD =  1 : NO SATISFACTORY POINT FOUND SUFFICIENTLY
  !                           DISTINCT FROM X
  !    IERR   : return CODE FROM THE QRP FACTORIZATION ROUTINE
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  !    TENSOLVE      ...  TSFSCL,TSCPSS,TSLSCH
  real(dp) :: wrk1(n), wrk2(m)
  integer  :: i, flag= 0, retcd1
  real(dp) :: slope, releng, temp1, temp2, almda, resnew, f1n, dtnorm, gnorm
  real(dp) :: sln, scl, beta, temp, almdat, almdam
  real(dp), parameter :: alpha= 1.0e-4_dp, tenth= 0.1_dp, z99= 0.99_dp
  mxtake= .false.
  sln= dnrm2(n, dt, 1)
  if(sln > stepmx) then
    ! step longer than maximum allowed
    scl= stepmx/sln
    dt(1:n)= scl * dt(1:n)
    sln= stepmx
  end if
  ! compute SLOPE= G-Trans * DT
  slope= dot_product( g(1:n), dt(1:n) )
  ! initialization of RETCD
  retcd= 0
  ! compute the smallest value allowable for the damping
  ! parameter ALMDA, i.e, ALMDAM
  releng= zero
  do i= 1,n
    temp1= MAX(ABS(xc(i)), one)
    temp2= ABS(dt(i))/temp1
    releng= MAX(releng, temp2)
  enddo
  almdam= steptl/releng
  almda= one
  ! compute the next iterate XP
  xp(1:n)= xc(1:n) + almda*dt(1:n)
  ! evaluate the objective function at XP and its residual
  call tsfscl(xp, dx, df, m, n, fvec, fp)
  fpnorm= half*dnrm2(m, fp, 1)**2
  ! test whether the full tensor step produces enough decrease in the
  ! l2 norm of of the objective function
  if (fpnorm < fcnorm + alpha* almda * slope) return
  ! compute the standard direction
  call tscpss(shat, m, n, p, method, global, epsm, fq, aja,  &
              anls, dn, fqq, pivot, curpos, pbar, zero1, ierr, resnew, flag)
  ! compute a step in the standard direction
  call tslsch(m, n, xc, dn, g, steptl, dx, df, fvec,  &
              mxtake, stepmx, wrk1, wrk2, fcnorm, f1n, retcd1)
  ! test whether the tensor direction is descent
  dtnorm= dnrm2(n, dt, 1)
  gnorm= dnrm2(n, g, 1)
  if(m > n) then
    beta= tenth
  else
    beta= alpha
  end if
  temp1= -beta*dtnorm*gnorm
  ! compute a step in the tensor direction
  if(slope <= temp1) then
    50 almdat= ((-almda**2)*slope)/(two*(fpnorm - fcnorm - almda*slope))
    temp= almda/ten
    almda= MAX(temp, almdat)
    if(almda < almdam) then
      if(retcd1== 1) then
        retcd= 1
        GO TO 70
      end if
    end if
    xp(1:n)= xc(1:n) + almda*dt(1:n)
    call tsfscl(xp, dx, df, m, n, fvec, fp)
    fpnorm= half*dnrm2(m, fp, 1)**2
    if (fpnorm > (fcnorm + (alpha* almda * slope))) GO TO 50
    if(almda== tenth .and. sln > z99*stepmx) mxtake=.true.
    ! select the next iterate that produces the lower function value
      70 if(f1n < fpnorm) then
      xp(1:n)= wrk1(1:n)
      fp(1:m)= wrk2(1:m)
      fpnorm =  f1n
    end if
  else
    xp(1:n)= wrk1(1:n)
    fp(1:m)= wrk2(1:m)
    fpnorm= f1n
  end if
  return
end subroutine tsmdls

function tsmfda( &
& anls, adt, ag, const1, const2, alpha, dlt,  &
& m, n, p, nwtake, ierr, vn) &
& result(fn_val)
!COMPUTES THE DERIVATIVE OF || F + J*D + 0.5*A*D**2 ||**2
!IN THE L2 NORM SENS, where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
  real(dp),intent(in)  :: anls(:,:)
  real(dp),intent(in)  :: adt(:)
  real(dp),intent(in)  :: ag(:)
  real(dp),intent(in)  :: const1(:)
  real(dp),intent(in)  :: const2(:)
  real(dp),intent(in)  :: alpha
  real(dp),intent(in)  :: dlt
  integer, intent(in)    :: m
  integer, intent(in)    :: n
  integer, intent(in)    :: p
  logical, intent(in)    :: nwtake
  integer, intent(in)    :: ierr
  real(dp),intent(in)  :: vn(:)
  real(dp)              :: fn_val
  !INPUT parameterS
  !    ANLS   : TENSOR MATRIX
  !    FQ     : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES
  !    ADT    : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !    AG     : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1 : SHAT-TRANS TIMES DT (SEE subroutine TS2DTR)
  !    CONST2 : SHAT-TRANS TIMES GBAR (SEE subroutine TS2DTR)
  !    ALPHA  : POINT AT WHICH TO EVALUATE THE DERIVATIVE OF function
  !    DLT    : CURRENT TRUST RADIUS
  !    M,N    : dimensionS OF THE PROBLEM
  !    P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !    NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !             NWTAKE= .true.  : STANDARD STEP TAKEN
  !             NWTAKE= .false. : TENSOR STEP TAKEN
  !    IERR   : return CODE FROM QRP FACTORIZATION ROUTINE:
  !             IERR=0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR=1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS
  !    VN     : F + J*D + 0.5*A*D**2
  !    TSMFDA : DERIVATIVE IN ALPHA OF || F + J*D + 0.5*A*D**2 ||**2
  !             where D= ALPHA*DT + SQRT(DLT**2 - ALPHA**2)
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSMFDV
  ! Workspace
  !       VNP    : DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
  real(dp) :: vnp(n+m)
  integer   :: len
  call tsmfdv( &
  & anls, adt, ag, const1, const2, alpha, dlt, m, n, p, nwtake, &
  & ierr, vnp)
  len= m
  if(ierr > 0) len= m + n
  fn_val= dot_product( vnp(1:len), vn(1:len) )
  return
end function tsmfda

subroutine tsmfdv( &
& anls, adt, ag, const1, const2, alpha, dlt,  &
& m, n, p, nwtake, ierr, vnp)
  !COMPUTES THE DERIVATIVE IN ALPHA OF THE VECTOR
  !VN= F(XC) + J(XC)*D + 0.5*A*D**2,
  !where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: adt(:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: alpha
  real(dp),intent(in)   :: dlt
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  logical, intent(in)     :: nwtake
  integer, intent(in)     :: ierr
  real(dp),intent(out)  :: vnp(:)
  !INPUT parameterS :
  !    ANLS  : TENSOR TERM MATRIX
  !     ADT  : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !      AG  : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1: SHAT-TRANS TIMES DT (SEE subroutine TS2DTR)
  !    CONST2: SHAT-TRANS TIMES GBAR (SEE subroutine TS2DTR)
  !    ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
  !      DLT : CURRENT TRUST RADIUS
  !      M,N : dimensionS OF THE PROBLEM
  !      P   : COLUMN dimension OF THE MATRIX ANLS
  !   NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS :
  !              NWTAKE= .true.  : STANDARD STEP TAKEN
  !              NWTAKE= .false. : TENSOR STEP TAKEN
  !     IERR : return CODE FROM THE QRP FACTORIZATION ROUTINE
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS :
  !     VNP  : THE DERIVATIVE IN ALPHA OF VN= F(XC) + J(XC)*D +
  !            0.5*A*D**2, where D= ALPHA*DT +  SQRT(DLT**2-ALPHA**2)
  integer   :: j, len
  real(dp) :: quant1, quant2, expr, const
  quant1= SQRT(dlt**2 - alpha**2)
  expr= - alpha/quant1
  vnp(1:n)= adt(1:n) + expr*ag(1:n)
  vnp(n+1:n+m)= zero
  if(nwtake) return
  quant2= quant1 - alpha**2/quant1
  len= m
  if(ierr > 0) len= m + n
  do j= 1,p
    const= half &
    &    *(two*alpha*(const1(j)**2 - const2(j)**2) + two*const1(j)*const2(j)*quant2)
    vnp(1:len)= vnp(1:len) + const * anls(1:len,j)
  enddo
  return
end subroutine tsmfdv

subroutine tsmgsa(s, n, sqrn, itn, shat, p, idp)
!FINDS A SET OF LINEARLY INDEPendENT DIRECTIONS USING
!THE MODIFIED GRAM-SCHMIDT ALGORITHM.
  real(dp),intent(in)   :: s(:,:)
  integer, intent(in)     :: n
  integer, intent(in)     :: sqrn
  integer, intent(in)     :: itn
  real(dp),intent(out)  :: shat(:,:)
  integer, intent(out)    :: p
  integer, intent(out)    :: idp(:)
  !INPUT parameterS :
  !    S   : MATRIX OF PAST DIRECTIONS
  !    N   : ROW dimension OF MATRIX S AND SHAT
  !    SQRN: MAXIMUM COLUMN dimension OF SHAT
  !    ITN : CURRENT ITERATION NUMBER
  !OUTPUT parameterS :
  !    SHAT: MATRIX OF LINEARLY INDEPendENT DIRECTIONS
  !    P   : COLUMN dimension OF THE MATRIX SHAT
  !    IDP : VECTOR THAT KEEPS TRACK OF THE INDICES CORRESPONDING TO
  !          THE LINEARLY INDEPendENT DIRECTIONS IN THE MATRIX S
  !    SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  integer :: j, k, l
  real(dp):: tol, tj, sj, sum, rtjs
  if(sqrn < itn) then
    k= sqrn
  else
    k= itn-1
  end if
  tol= SQRT(two)/two
  shat(1:n,1:k)= s(1:n,1:k)
  p= 0
  do j= 1,k
    tj= dnrm2(n, shat(:,j), 1)
    sj= dnrm2(n, s(:,j), 1)
    if(tj/sj > tol) then
      p= p + 1
      idp(p)= j
      rtjs= one/tj**2
      do l= j+1,k
        sum= -rtjs * dot_product( shat(1:n,l), shat(1:n,j) )
        shat(1:n,l)= shat(1:n,l) + sum * shat(1:n,j)
      enddo
    end if
  enddo
  do j= 1,p
    shat(1:n,j)= s(1:n,idp(j))
  enddo
  return
end subroutine tsmgsa

subroutine tsmsda( &
& anls, fq, adt, ag, const1, const2, alpha, dlt, m, n, p, &
& nwtake, ierr, skip, vn, fn_val)
!COMPUTES THE SECOND DERIVATIVE OF
!  || F + J*D + 0.5*A*D**2 ||**2 IN THE L2 NORM SENS,
!where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
!N.B. Changed from a function to a subroutine by AJM
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: fq(:)
  real(dp),intent(in)   :: adt(:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: alpha
  real(dp),intent(in)   :: dlt
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  integer, intent(in)   :: p
  logical, intent(in)   :: nwtake
  integer, intent(in)   :: ierr
  logical, intent(in)   :: skip
  real(dp),intent(out)  :: vn(:)
  real(dp),intent(out)  :: fn_val
  !INPUT parameterS
  !    ANLS   : TENSOR TERM MATRIX AT CURRENT ITERATE
  !    FQ     : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES
  !    ADT    : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !     AG    : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1 : SHAT-TRANS TIMES DT (SEE subroutine TS2DTR)
  !    CONST2 : SHAT-TRANS TIMES GBAR (SEE subroutine TS2DTR)
  !    ALPHA  : POINT AT WHICH TO EVALUATE THE SECOND DERIVATIVE OF function
  !    DLT    : CURRENT TRUST RADIUS
  !    M,N    : dimensionS OF THE PROBLEM
  !    P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !    NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !             NWTAKE= .true.  : STANDARD STEP TAKEN
  !             NWTAKE= .false. : TENSOR STEP TAKEN
  !    IERR   : return CODE FROM QRP FACTORIZATION ROUTINE
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS
  !    VN     : F + J*D + 0.5*A*D**2
  !    TSMSDA : SECOND DERIVATIVE IN ALPHA OF || F + J*D + 0.5*A*D**2 ||**2
  !             where D=ALPHA*DT + SQRT(DLT**2-ALPHA**2)
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSMAFA,TSMFDV,TSMSDV
  ! Workspace
  !       VNP    : DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
  !       VNS    : SECOND DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
  real(dp) :: vnp(n+m), vns(n+m)
  integer   :: len
  if(.not. skip) then
    call tsmafa( &
    & anls, fq, adt, ag, const1, const2, alpha, dlt, m, n, p,  &
    & nwtake, ierr, vn)
    call tsmfdv( &
    & anls, adt, ag, const1, const2, alpha, dlt, m, n, p,   &
    & nwtake, ierr, vnp)
  end if
  call tsmsdv(anls, ag, const1, const2, alpha, dlt, m, n, p, nwtake, ierr, vns)
  len= m
  if(ierr > 0) len= m + n
  fn_val= SUM( vnp(1:len)**2 ) + dot_product( vns(1:m), vn(1:m) )
  return
end subroutine tsmsda

subroutine tsmsdv( &
& anls, ag, const1, const2, alpha, dlt, m, n, p, nwtake,  &
& ierr, vns)
!THIS ROUTINE COMPUTES THE SECOND DERIVATIVE IN ALPHA OF THE VECTOR
!VN= F(XC) + J(XC)*D + 0.5*A*D**2,
!where D= ALPHA*DT + SQRT(DLT**2-ALPHA**2).
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: alpha
  real(dp),intent(in)   :: dlt
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  logical, intent(in)     :: nwtake
  integer, intent(in)     :: ierr
  real(dp),intent(out)  :: vns(:)
  !INPUT parameterS :
  !    ANLS  : TENSOR TERM MATRIX
  !     ADT  : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !      AG  : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1: SHAT-TRANS * DT (SEE subroutine TS2DTR)
  !    CONST2: SHAT-TRABS * GBAR (SEE subroutine TS2DTR)
  !    ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
  !      DLT : CURRENT TRUST RADIUS
  !      NR  : LEADING dimension OF ANLS
  !      M,N : dimensionS OF THE PROBLEM
  !      P   : COLUMN dimension OF THE MATRIX ANLS
  !   NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS :
  !               NWTAKE= .true.  : STANDARD STEP TAKEN
  !               NWTAKE= .false. : TENSOR STEP TAKEN
  !     IERR : return CODE FROM THE QRP FACTORIZATION ROUTINE :
  !             IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !             IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !OUTPUT parameterS :
  !     VNS  : THE SECOND DERIVATIVE IN ALPHA OF VN= F(XC) + J(XC)*D
  !            + 0.5*A*D**2, where D= ALPHA*DT +  SQRT(DLT**2-ALPHA**2)
  integer :: j, len
  real(dp):: quant1, expr, const, quant2
  real(dp), parameter :: onepf= 1.5_dp
  quant1= dlt**2 - alpha**2
  expr=  -dlt**2 * SQRT(quant1) / quant1**2
  vns(1:n)=  expr*ag(1:n)
  vns(n+1:n+m)= zero
  if(nwtake) return
  quant2= -three*alpha/SQRT(quant1) - alpha**3/quant1**onepf
  len= m
  if(ierr > 0) len= m + n
  do j= 1,p
    const= half &
    &    *(two*(const1(j)**2 - const2(j)**2) +two*const1(j)*const2(j)*quant2)
    vns(1:len)= vns(1:len) + const *anls(1:len,j)
  enddo
  return
end subroutine tsmsdv

subroutine tsmslv( &
!FINDS THE TENSOR AND STANDARD STEPS.
& p, shat, maxm, sqrn, m, n, epsm, method, &
& global, x, typxu, xpls, gpls, curpos, pbar, pivot, fq, fqq, &
& dn, dt, restns, resnew, itrmcd, flag, zero1, ierr)
  integer, intent(in)        :: p
  real(dp),intent(inout)  :: shat(:,:)
  integer, intent(in)        :: maxm
  integer, intent(in)        :: sqrn
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(in)      :: epsm
  integer, intent(in)        :: method
  integer, intent(in)        :: global
  real(dp),intent(out)     :: x(:)
  real(dp),intent(inout)  :: typxu(:)
  real(dp),intent(out)     :: xpls(:)
  real(dp),intent(inout)  :: gpls(:)
  integer, intent(out)       :: curpos(:)
  integer, intent(out)       :: pbar(:)
  integer, intent(out)       :: pivot(:)
  real(dp),intent(out)     :: fq(:)
  real(dp),intent(out)     :: fqq(:)
  real(dp),intent(out)     :: dn(:)
  real(dp),intent(out)     :: dt(:)
  real(dp),intent(out)     :: restns
  real(dp),intent(out)     :: resnew
  integer, intent(out)       :: itrmcd
  integer, intent(out)       :: flag
  integer, intent(out)       :: zero1
  integer, intent(out)       :: ierr
  !INPUT parameterS :
  !    P      : COLUMN dimension OF MATRICES ANLS AND S
  !    MAXM   : LEADING dimension OF AJA AND ANLS
  !    SQRN   : LEADING dimension OF MATRICES A AND WRK
  !    M,N    : dimensionS OF PROBLEM
  !    EPSM   : MACHINE PRECISION
  !    X      : ESTIMATE TO A ROOT OF FCN (useD BY UNCMIN)
  !    TYPXU  : TYPICAL size FOR EACH COMPONENT OF X (useD BY UNCMIN)
  !    METHOD : METHOD TO use
  !             METHOD= 0 : STANDARD METHOD IS useD
  !             METHOD= 1 : TENSOR METHOD IS useD
  !    GLOBAL : GLOBAL STRATEGY useD
  !INPUT/OUTPUT parameterS :
  !    SHAT   : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !OUTPUT parameterS :
  !    DN     : STANDARD STEP
  !    DT     : TENSOR STEP
  !    FLAG   : returnED CODE WITH THE FOLLOWING MEANING :
  !             FLAG= 0 : NO SINGULARITY DETECTED WHEN FACTORIZING AJA
  !             FLAG= 1 : SINGULARITY DETECTED WHEN FACTORIZING AJA
  !                        FROM 1 TO N-P
  !             FLAG= 2 : SINGULARITY DETECTED WHEN FACTORIZING AJA
  !                        FROM N-P TO N
  !    IERR   : returnED CODE WITH THE FOLLOWING MEANING :
  !               IERR= 0 : NO SINGULARITY DETECTED WHEN FACTORIZING AJA
  !               IERR= 1 : SINGULARITY DETECTED WHEN FACTORIZING AJA
  !    XPLS   : LOCAL MINIMUM OF OPTIMIZATION function FCN (useD BY UNCMIN)
  !    FPLS   : function value AT SOLUTION OF OPTIMIZATION function FCN
  !             (useD IN UNCMIN)
  !    GPLS   : GRADIENT AT SOLUTION XPLS (useD BY UNCMIN)
  !    CURPOS,PIVOT,PBAR : PIVOT VECTORS
  !    RESTNS : TENSOR RESIDUAL
  !    RESNEW : STANDARD RESIDUAL
  !    ITRMCD : TERMINATION CODE (FOR UNCMIN)
  !    SUBprogramS callED:
  !    TENSOLVE      ...  TSQLFC,QTRNS,TSQRFC,TSQMTS,TSQMUV,TSSQP1
  !    TENSOLVE      ...  TSQ1P1,TSD1SV,TSPRMV,TSQMLV,TSCPSS
  !    UNCMIN        ...  DFAUT,OPTif9
  real(dp) :: wrk3(m), wrk4(m)
  integer   :: msg, itnlim, ipr
  integer   :: q, meth, iexp, ndigit, iagflg, iahflg
  real(dp) :: root, typfu, dlt, gradlt, stepmx, steptl, fpls
  itrmcd= 0
  if(n== 1) then
    shat(2,1)= one
    shat(3,1)= one
    curpos(1)= 1
    fq= fc(1:m)
  else
    ! perform a QL decomposition of S
    call tsqlfc(shat, n, p, ierr)
    ! compute AJA times Q-trans
    call tsjqtp(shat, n, m, p, aja)
    ! perform a QR factorization of AJA
    call tsqrfc(aja, n, m, 1, n-p, ierr, epsm, curpos, zero1)
    if(ierr== 1) then; q= n - zero1 + 1
    else             ; q= p
    end if
    qrank= q
    call tsqmts(anls, aja, m, m, p, 1, zero1)
    call tsqmuv(aja, fc, fq, m, 1, zero1, .false.)
  end if
  ! Minimize the lower m-n+q quadratic equations in p unknowns of the tensor
  ! model.  The minimization is performed analytically if p=1,q>1, or
  ! p=1,q=1,m>n, or n=1,m>n.  Otherwise an unconstrained minimization package,
  ! UNCMIN, is used.
  if((p== 1 .and. q > 1) .or. (p== 1 .and. q== 1 .and. m > n)  &
        .or. (n== 1 .and. m > n)) then
    call tssqp1(aja, anls, shat, fq, m, n, q, root, restns)
    xpls(1)= root
  else if(m== n .and. p== 1 .and. q== 1 .or. (m== 1 .and. n== 1)) then
    call tsq1p1(aja, anls, shat, fq, n, root, restns)
    xpls(1)= root
  else
    call dfaut( &
    & p, typxu, typfu, meth, iexp, msg, ndigit, itnlim,  &
    & iagflg, iahflg, ipr, dlt, gradlt, stepmx, steptl)
    iagflg= 1
    iahflg= 0
    iexp=   0
    meth=   2
    msg=    9
    x(1:p)= zero
    call optif9( &
    & sqrn, p, x, tsqfcn, tsdfcn, d2fcn, typxu, typfu, meth, iexp,  &
    & msg, ndigit, itnlim, iagflg, iahflg, dlt, gradlt,  &
    & stepmx, steptl, xpls, fpls, gpls, itrmcd)
    ! compute the tensor residual
    restns= SQRT(two*fpls)
  end if
  wrk4(n-p+1:n)= xpls(1:p)
  if(n== 1) then
    dt(1)= wrk4(1)
  else
    ! compute the first n-p components of the tensor step
    call tsd1sv(aja, shat, anls, fq, xpls, maxm, m, n, p, q, epsm, pivot, wrk3)
    call tsprmv(wrk4, wrk3, curpos, n-p, 0)
    ! premultiply the tensor step by the orthogonal matrix resulting
    ! from the QL factorization of S
    call tsqmlv(n, p, shat, wrk4, dt, .true.)
  end if
  ! compute the standard step if needed
  if(global== 1 .or. (m > n .and. global== 0)) then
    call tscpss( &
    & shat, m, n, p, method, global, epsm, fq, aja, anls, &
    & dn, fqq, pivot, curpos, pbar, zero1, ierr, resnew, flag)
  end if
  return
end subroutine tsmslv

subroutine tsnesv( &
!THE DRIVER FOR NONLINEAR EQUATIONS/NONLINEAR LEAST SQUARES PROBLEMS.
& maxm, xc, m, n, typx, typf, itnlim, jacflg,  &
& gradtl, steptl, ftol, method, global, stepmx, dlt, ipr,  &
& dfn, dxn, epsm, sqrn, fvec, jac, msg, xp, fp, gp, itn, termcd)
  integer, intent(in)     :: maxm
  real(dp),intent(inout)  :: xc(:)
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  real(dp),intent(in)     :: typx(:)
  real(dp),intent(in)     :: typf(:)
  integer, intent(in)     :: itnlim
  integer, intent(in)     :: jacflg
  real(dp),intent(in)     :: gradtl
  real(dp),intent(inout)  :: steptl
  real(dp),intent(in)     :: ftol
  integer, intent(in)     :: method
  integer, intent(in)     :: global
  real(dp),intent(inout)  :: stepmx
  real(dp),intent(inout)  :: dlt
  integer, intent(in)     :: ipr
  real(dp),intent(in)     :: dfn(:)
  real(dp),intent(in)     :: dxn(:)
  real(dp),intent(in)     :: epsm
  integer, intent(in)     :: sqrn
  integer, intent(inout)  :: msg
  real(dp),intent(out)    :: xp(:)
  real(dp),intent(out)    :: fp(:)
  real(dp),intent(out)    :: gp(:)
  integer, intent(out)    :: itn,termcd
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds,only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)  :: m, n
    end subroutine fvec
    subroutine jac(x, aja, m, n)
      use M_Kinds,only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: aja(:,:)
      integer, intent(in)  :: m, n
    end subroutine jac
  endinterface
  !INPUT parameterS :
  !      MAXM   : LEADING dimension OF AJA, ANLS, AND FV
  !      XC     : INITIAL ESTIMATE OF SOLUTION
  !      M,N    : dimensionS OF PROBLEM
  !      TYPX   : TYPICAL size FOR EACH COMPONENT OF X
  !      TYPF   : TYPICAL size FOR EACH COMPONENT OF F
  !      ITNLIM : MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  !      JACFLG : JACOBIAN FLAG WITH THE FOLLOWING MEANINGS:
  !               JACFLG= 1 if ANALYTIC JACOBIAN SUPPLIED
  !               JACFLG= 0 if ANALYTIC JACOBIAN NOT SUPPLIED
  !      GRADTL : TOLERANCE AT WHICH GRADIENT IS CONSIDERED close ENOUGH
  !               TO ZERO TO TERMINATE ALGORITHM
  !      STEPTL : TOLERANCE AT WHICH SUCCESSIVE ITERATES ARE CONSIDERED
  !               close ENOUGH TO TERMINATE ALGORITHM
  !      FTOL   : TOLERANCE AT WHICH function value IS CONSIDERED close
  !               ENOUGH TO ZERO
  !      METHOD : METHOD TO use
  !               METHOD= 0 : STANDARD METHOD IS useD
  !               METHOD= 1 : TENSOR METHOD IS useD
  !      GLOBAL : GLOBAL STRATEGY TO use
  !               GLOBAL= 0 : LINE SEARCH
  !               GLOBAL= 1 : 2-dimensionAL TRUST REGION
  !      STEPMX : MAXIMUM ALLOWABLE STEP size
  !      DLT    : TRUST REGION RADIUS
  !      IPR    : DEVICE TO WHICH TO Send OUTPUT
  !      DFN    : DIAGONAL SCALING MATRIX FOR F
  !      DXN    : DIAGONAL SCALING MATRIX FOR X
  !      EPSM   : MACHINE PRECISION
  !      SQRN   : MAXIMUM COLUMN dimension OF ANLS, S, AND SHAT
  !      FVEC   : NAME OF subroutine TO EVALUATE function
  !      JAC    : (optional) NAME OF subroutine TO EVALUATE JACOBIAN.
  !               MUST BE DECLARED external IN callING ROUTINE
  !INPUT-OUTPUT parameterS :
  !      MSG : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  !OUTPUT parameterS :
  !      XP : SOLUTION TO THE SYSTEM OF NONLINEAR EQUATIONS
  !      FP : function value AT THE SOLUTION
  !      GP : GRADIENT AT THE SOLUTION
  !      TERMCD : TERMINATION CODE
  !SUBprogramS callED:
  !      LEVEL 1 BLAS  ...  DNRM2
  !      LEVEL 2 BLAS  ...  DGEMV
  !      TENSOLVE      ...  TSSCLX,TSFSCL,TSSCLJ,TSCHKJ,TSNSTP,TSSSTP,
  !      TENSOLVE      ...  TSLSCH,TS2DTR,TSRSLT,TSMGSA,TSFRMT,TSMSLV,
  !      TENSOLVE      ...  TSSLCT,TSMDLS,TSUPSF
  !Local variables (were dummy arguments in F77 version)
  !      X      : ESTIMATE TO A ROOT OF FCN ( useD BY UNCMIN)
  !      TYPXU  : TYPICAL size FOR EACH COMPONENT OF X (useD BY UNCMIN)
  !      XPLS   : LOCAL MINIMUM OF OPTIMIZATION function FCN useD BY UNCMIN
  !      GPLS   : GRADIENT AT SOLUTION XPLS (useD BY UNCMIN)
  !      DN     : STANDARD STEP
  !      DT     : TENSOR STEP
  !      SHAT   : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !      CURPOS,PIVOT,PBAR : PIVOT VECTORS
  ! Workspace
  real(dp) :: x(n), typxu(n), xpls(n), gpls(n), dn(n), dt(n), shat(m+2,sqrn)
  integer  :: curpos(n), pivot(n), pbar(n)
  real(dp) :: df(n), gbar(n), fq(m), fqq(m+n), fhat(m), fv(m,n)
  integer  :: p, flag, retcd, zero1, ierr, itrmcd, icscmx
  real(dp) :: fpls, fnorm, restns, resnew
  logical  :: nwtake, mxtake
  ! initialization
  itn= 0
  ierr= 0
  nwtake= .true.
  !
  allocate( fc(m), anls(m+2,sqrn), aja(m+n+2,n), s(m+2,sqrn) )
  !
  meqns= m
  nvars= n
  !
  call tssclx(xc, dxn, n)
  !
  if(MOD(msg/8,2) /= 1) then
    write(ipr,896)            !896 format('  TSNESV      TYPICAL X')
    write(ipr,900) typx(1:n)  !900 format(100('  TSNESV     ', 3(g20.13, "   ")/))
    write(ipr,897)            !897 format('  TSNESV      DIAGONAL SCALING MATRIX FOR X')
    write(ipr,900) dxn(1:n)   !900 format(100('  TSNESV     ', 3(g20.13, "   ")/))
    write(ipr,898)            !898 format('  TSNESV      TYPICAL F')
    write(ipr,900) typf(1:m)  !900 format(100('  TSNESV     ', 3(g20.13, "   ")/))
    write(ipr,899)            !899 format('  TSNESV      DIAGONAL SCALING MATRIX FOR F')
    write(ipr,900) dfn(1:m)   !900 format(100('  TSNESV     ', 3(g20.13, "   ")/))
    write(ipr,901) jacflg     !901 format('  TSNESV      JACOBIAN FLAG     = ', i1)
    write(ipr,902) method     !902 format('  TSNESV      METHOD useD       = ', i1)
    write(ipr,903) global     !903 format('  TSNESV      GLOBAL STRATEGY   = ', i1)
    write(ipr,904) itnlim     !904 format('  TSNESV      ITERATION LIMIT   = ', i5)
    write(ipr,905) epsm       !905 format('  TSNESV      MACHINE EPSILON   = ', g20.13)
    write(ipr,906) steptl     !906 format('  TSNESV      STEP TOLERANCE    = ', g20.13)
    write(ipr,907) gradtl     !907 format('  TSNESV      GRADIENT TOLERANCE= ', g20.13)
    write(ipr,908) ftol       !908 format('  TSNESV      function TOLERANCE= ', g20.13)
    write(ipr,909) stepmx     !909 format('  TSNESV      MAXIMUM STEP size = ', g20.13)
    write(ipr,910) dlt        !910 format('  TSNESV      TRUST REG RADIUS  = ', g20.13)
  end if
  ! Evaluate analytic or finite difference Jacobian and check analytic
  ! Jacobian, if requested
  call tsfscl(xc, dxn, dfn, m, n, fvec, fc)
  call tssclj(xc, dxn, typx, fc, dfn, m, n, epsm, jacflg, fvec, jac, aja)
  if(jacflg== 1) then
    if(MOD(msg/2,2)== 0) then
      call tschkj(aja, xc, fc, m, n, epsm, dfn, dxn, typx, ipr, fvec, msg)
      if(msg < 0) then
        deallocate( fc, anls, aja, s )
        return
      end if
    end if
  end if
  ! compute the gradient at the current iterate XC
  call dgemv('T', m, n, one, aja, maxm, fc, 1, zero, gp, 1)
  ! compute the residual of FC
  fnorm= half*dnrm2(m, fc, 1)**2
  ! check stopping criteria for input XC
  call tsnstp( &
  & gp, xc, fc, xc, steptl, gradtl, retcd, ftol, itn,  &
  & itnlim, icscmx, mxtake, m, n, msg, ipr, fnorm, termcd)
  if(termcd > 0) then
    fpls= fnorm
    GO TO 120
  end if
  ! iteration 1
  itn= 1
  ! compute the standard step
  fhat(1:m)= fc(1:m)
  call tssstp(aja, fhat, m, n, epsm, global, dn, fqq, pivot, pbar, ierr)
  ! choose next iterate XP by a global strategy
  if(global== 0) then
    call tslsch( &
    & m, n, xc, dn, gp, steptl, dxn, dfn, fvec,  &
    & mxtake, stepmx, xp, fp, fnorm, fpls, retcd)
  else
    shat(1:n,1:sqrn)= zero
    call ts2dtr( &
    & aja, shat, anls, dn, gp, gbar, xc, method, nwtake, stepmx,   &
    & steptl, epsm, mxtake, dlt, fqq, maxm, m, n, sqrn,      &
    & curpos, pivot, pbar, itn, ierr, flag, dxn, dfn, fvec, fnorm, &
    & xp, fp, fpls, retcd)
  end if
  if(MOD(msg/8,2)== 0) call tsrslt(n, xc, fnorm, gp, 0, ipr)
  ! evaluate the Jacobian at the new iterate XP
  call tssclj(xp, dxn, typx, fp, dfn, m, n, epsm, jacflg, fvec, jac, aja)
  ! compute the gradient at the new iterate XP
  call dgemv('T', m, n, one, aja, maxm, fp, 1, zero, gp, 1)
  ! check stopping criteria for the new iterate XP
  call tsnstp( &
  & gp, xp, fp, xc, steptl, gradtl, retcd, ftol, itn,  &
  & itnlim, icscmx, mxtake, m, n, msg, ipr, fpls, termcd)
  if(termcd > 0) GO TO 120
  if(MOD(msg/16,2)== 1) call tsrslt(n, xp, fpls, gp, itn, ipr)
  ! update S and FV
  s(1:n,1)= xc(1:n) - xp(1:n)
  fv(1:m,1)= fc(1:m)
  ! update XC and FC
  xc(1:n)= xp(1:n)
  fc(1:m)= fp(1:m)
  fnorm= fpls
  ! iteration > 1
  80 itn= itn + 1
  ! if the standard method is selected then compute the standard step
  if(method== 0) then
    fhat(1:m)= fc(1:m)
    call tssstp(aja, fhat, m, n, epsm, global, df, fqq, pivot, pbar, ierr)
  end if
  ! if the tensor method is selected then form the tensor model
  if(method== 1) then
    ! select the past linearly independent directions
    call tsmgsa(s, n, sqrn, itn, shat, p, curpos)
    !
    ! form the tensor term
    call tsfrmt(shat, s, aja, fv, fc, m, n, p, curpos, anls)
    !
    ! solve the tensor model for the tensor step DT and compute DN
    ! as a by-product if the global strategy selected is the
    ! two-dimensional trust region or M > N
    call tsmslv( &
    & p, shat, maxm, sqrn, m, n, epsm, method, &
    & global, x, typxu, xpls, gpls, curpos, pbar, pivot,  &
    & fq, fqq, dn, dt, restns, resnew, itrmcd, flag, zero1, ierr)
    !
    ! decide which step to use (DN or DT)
    if(global== 1 .or. (m > n .and. global== 0)) &
    & call tsslct(restns, resnew, itrmcd, fc, m, n, dn, dt, gp, df, nwtake)
    !
  end if
  ! choose the next iterate XP by a global strategy
  if(global== 0) then
    if(method== 0) then
      call tslsch( &
      & m, n, xc, df, gp, steptl, dxn, dfn, fvec,  &
      & mxtake, stepmx, xp, fp, fnorm, fpls, retcd)
    else if(m== n) then
      call tsmdls( &
      & aja, shat, anls, xc, m, n, p, dt, gp, dxn, dfn, &
      & fvec, method, steptl, global, stepmx, epsm, fq, dn, fqq,   &
      & pivot, curpos, pbar, mxtake, xp, fp, fnorm, fpls, zero1,   &
      & retcd, ierr)
    else
      call tslsch( &
      & m, n, xc, df, gp, steptl, dxn, dfn, fvec,  &
      & mxtake, stepmx, xp, fp, fnorm, fpls, retcd)
    end if
  else
    call ts2dtr( &
    & aja, shat, anls, df, gp, gbar, xc, method, nwtake, stepmx,   &
    & steptl, epsm, mxtake, dlt, fqq, maxm, m, n, p, curpos, &
    & pivot, pbar, itn, ierr, flag, dxn, dfn, fvec, fnorm, xp, fp, &
    & fpls, retcd)
  end if
  ! evaluate the Jacobian at the new iterate XP
  call tssclj(xp, dxn, typx, fp, dfn, m, n, epsm, jacflg, fvec, jac, aja)
  ! evaluate the gradient at the new iterate XP
  call dgemv('T', m, n, one, aja, maxm, fp, 1, zero, gp, 1)
  ! check stopping criteria for the new iterate XP
  call tsnstp( &
  & gp, xp, fp, xc, steptl, gradtl, retcd, ftol, itn,  &
  & itnlim, icscmx, mxtake, m, n, msg, ipr, fpls, termcd)
  !
  if(termcd > 0) GO TO 120
  !
  if(MOD(msg/16,2)== 1) call tsrslt(n, xp, fpls, gp, itn, ipr)
  !
  ! if tensor method is selected then update the matrices S and FV
  if(method== 1) call tsupsf(fc, xc, xp, sqrn, itn, m, n, s, fv)
  !
  ! update XC, FC, and FNORM
  xc(1:n)= xp(1:n)
  fc(1:m)= fp(1:m)
  fnorm= fpls
  !
  GO TO 80
  !
  ! termination
  120 if(MOD(msg/8,2)== 0) then
    if(itn /= 0) then
      call tsrslt(n, xp, fpls, gp, itn, ipr)
    else
      fpls= half*dnrm2(m, fc, 1)**2
      call tsrslt(n, xc, fpls, gp, itn, ipr)
    end if
  end if
  !
  deallocate( fc, anls, aja, s )
  !
  return
  !
  896 format('  TSNESV      TYPICAL X')
  897 format('  TSNESV      DIAGONAL SCALING MATRIX FOR X')
  898 format('  TSNESV      TYPICAL F')
  899 format('  TSNESV      DIAGONAL SCALING MATRIX FOR F')
  900 format(100('  TSNESV     ', 3(g20.13, "   ")/))
  901 format('  TSNESV      JACOBIAN FLAG     = ', i1)
  902 format('  TSNESV      METHOD useD       = ', i1)
  903 format('  TSNESV      GLOBAL STRATEGY   = ', i1)
  904 format('  TSNESV      ITERATION LIMIT   = ', i5)
  905 format('  TSNESV      MACHINE EPSILON   = ', g20.13)
  906 format('  TSNESV      STEP TOLERANCE    = ', g20.13)
  907 format('  TSNESV      GRADIENT TOLERANCE= ', g20.13)
  908 format('  TSNESV      function TOLERANCE= ', g20.13)
  909 format('  TSNESV      MAXIMUM STEP size = ', g20.13)
  910 format('  TSNESV      TRUST REG RADIUS  = ', g20.13)
end subroutine tsnesv

subroutine tsnstp( &
!DECIDES WHETHER TO TERMINATE THE NONLINEAR ALGORITHM.
& g, xplus, fplus, xc, steptl, gradtl, retcd, ftol, itn,  &
& itnlim, icscmx, mxtake, m, n, msg, ipr, fnorm, termcd)
  ! N.B. Added TERMCD= 0  24 May 2001
  real(dp),intent(in)    :: g(:)
  real(dp),intent(in)    :: xplus(:)
  real(dp),intent(in)    :: fplus(:)
  real(dp),intent(in)    :: xc(:)
  real(dp),intent(in)    :: steptl
  real(dp),intent(in)    :: gradtl
  integer, intent(in)      :: retcd
  real(dp),intent(in)    :: ftol
  integer, intent(in)      :: itn
  integer, intent(in)      :: itnlim
  integer, intent(inout)  :: icscmx
  logical, intent(in)      :: mxtake
  integer, intent(in)      :: m
  integer, intent(in)      :: n
  integer, intent(in)      :: msg
  integer, intent(in)      :: ipr
  real(dp),intent(in)    :: fnorm
  integer, intent(out)     :: termcd
  !INPUT parameterS :
  !    G     : GRADIENT AT XC
  !    XPLUS : NEW ITERATE
  !    FPLUS : function value AT XPLUS
  !    XC    : CURRENT ITERATE
  !    STEPTL: STEP TOLERANCE
  !    GRADTL: GRADIENT TOLERANCE
  !    RETCD : return CODE WITH THE FOLLOWING MEANINGS :
  !            RETCD= 0 : SUCCESSFUL GLOBAL STRATEGY
  !            RETCD= 1 : UNSUCCESSFUL GLOBAL STRATEGY
  !    FTOL  : function TOLERANCE
  !    ITN   : ITERATION NUMBER
  !    ITNLIM: ITERATION NUMBER LIMIT
  !    ICSCMX: NUMBER CONSECUTIVE STEPS >= STEPMX
  !    MXTAKE: BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH
  !    M     : dimension OF FPLUS
  !    N     : dimension OF G, XC, AND XPLUS
  !    MSG   : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  !    IPR   : DEVICE TO WHICH TO Send OUTPUT
  !OUTPUT parameterS :
  !    TERMCD: return CODE WITH THE FOLLOWING MEANINGS :
  !           TERMCD= 0 NO TERMINATION CRITERION SATISFIED
  !           TERMCD > 0 : SOME TERMINATION CRITERION SATISFIED
  !           TERMCD= 1 : NORM OF SCALED function value IS LESS THAN FTOL
  !           TERMCD= 2 : GRADIENT TOLERANCE REACHED
  !           TERMCD= 3 : SCALED DISTANCE BETWEEN LAST TWO STEPS < STEPTL
  !           TERMCD= 4 : UNSUCCESSFUL GLOBAL STRATEGY
  !           TERMCD= 5 : ITERATION LIMIT EXCEEDED
  !           TERMCD= 6 : 5 CONSECUTIVE STEPS OF lenGTH STEPMX HAVE BEEN TAKEN
  !SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  IDAMAX
  integer   :: i
  real(dp) :: res, d, rgx, relgrd, rsx, relstp
  termcd= 0
  ! check whether scaled function is within tolerance
  res= ABS(fplus(idamax(m, fplus, 1)))
  if(res <= ftol) then
    termcd= 1
    if(MOD(msg/8,2)== 0) write(ipr,701)
    return
  end if
  ! check whether scaled gradient is within tolerance
  d= one/MAX(fnorm, DBLE(n/2))
  rgx= zero
  do i= 1,n
    relgrd= ABS(g(i)) * MAX(ABS(xplus(i)), one)*d
    rgx= MAX(rgx,relgrd)
  enddo
  if(rgx <= gradtl) then
    termcd= 2
    if(MOD(msg/8,2)== 0) write(ipr,702)
    return
  end if
  if(itn== 0) return
  if(retcd== 1) then
    termcd= 4
    if(MOD(msg/8,2)== 0)  write(ipr,703)
    return
  end if
  ! check whether relative step length is within tolerance
  rsx= zero
  do i= 1,n
    relstp= ABS(xplus(i) - xc(i))/MAX(xplus(i), one)
    rsx= MAX(rsx, relstp)
  enddo
  if(rsx <= steptl) then
    termcd= 3
    if(MOD(msg/8,2)== 0) write(ipr,704)
    return
  end if
  ! check iteration limit
  if(itn >= itnlim) then
    termcd= 5
    if(MOD(msg/8,2)== 0) write(ipr,705)
  end if
  ! check number of consecutive steps .ge. stepmx
  if(mxtake) then
    icscmx= icscmx + 1
    if(icscmx >= 5) then
      termcd= 6
      if(MOD(msg/8,2)== 0) write(ipr,706)
    end if
  else
    icscmx= 0
  end if
  return
  701 format(/,'  TSNSTP      function value close TO ZERO')
  702 format(/,'  TSNSTP      RELATIVE GRADIENT close TO ZERO')
  703 format(/,'  TSNSTP      LAST GLOBAL STEP FAILED TO LOCATE A',/  &
      '  TSNSTP      POINT LOWER THAN THE CURRENT ITERATE')
  704 format(/,'  TSNSTP      SUCCESSIVE ITERATES WITHIN TOLERANCE',/  &
      '  TSNSTP      CURRENT ITERATE IS PROBABLY SOLUTION')
  705 format(/,'  TSNSTP      ITERATION LIMIT EXCEEDED',/  &
      '  TSNSTP      ALGORITHM FAILED')
  706 format(/,'  TSNSTP      MAXIMUM STEP size EXCEEDED 5',  &
      ' CONSECUTIVE TIMES',/  &
      '  TSNSTP      EITHER THE function IS UNBOUNDED BELOW',/  &
      '  TSNSTP      BECOMES ASYMPTOTIC TO A FINITE value',/  &
      '  TSNSTP      FROM ABOVE IN SOME DIRECTION',/  &
      '  TSNSTP      OR STEPMX IS TOO SMALL')
end subroutine tsnstp

subroutine tsprmv(x, y, pivot, n, job)
!THIS subroutine PERFORMS A VECTOR PERMUTATION.
  real(dp),intent(out)  :: x(:)
  real(dp),intent(in)   :: y(:)
  integer, intent(in)     :: pivot(:)
  integer, intent(in)     :: n
  integer, intent(in)     :: job
  !INPUT parameterS :
  !      Y :  VECTOR TO TSPRMV
  !  PIVOT :  PIVOT VECTOR
  !      N :  dimension OF THE VECTORS Y AND PIVOT
  !OUTPUT parameterS :
  !      X : PIVOTED VECTOR
  if(job== 0) then; x(pivot(1:n))= y(1:n) ! permute Y
  else            ; x(1:n)= y(pivot(1:n)) ! reverse permute of Y
  end if
  return
end subroutine tsprmv

subroutine tsrslt(n, xp, fval, gp, itn, ipr)
!printS INformatION.
  use M_Numeric_Tools,only: fNewt_I,fNewtF,fNewtR
  use M_IOtools
  !
  integer, intent(in)    :: n
  real(dp),intent(in)  :: xp(:)
  real(dp),intent(in)  :: fval
  real(dp),intent(in)  :: gp(:)
  integer, intent(in)    :: itn
  integer, intent(in)    :: ipr
  !INPUT parameterS :
  !    M,N  : dimensionS OF PROBLEM
  !    XP   : NEXT ITERATE
  !    FVAL : SUM OF SQUARES OF F(XP)
  !    GP   : GRADIENT AT XP
  !    ITN  : ITERATION NUMBER
  !    IPR  : DEVICE TO WHICH TO Send OUTPUT
  !
!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!
!!ADDED J.MOUTTE 19/02/2008 10:50
!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!
  fNewt_I= fNewt_I +1
  !if(fNewtF>0) call OutStrVec(fNewtF,vX(:)/Ln10,Opt_I=fNewt_I,Opt_J=Its,Opt_C="F")
  if(fNewtF>0) call OutStrVec(fNewtF,xp(1:N),     Opt_I=fNewt_I,Opt_J=Itn,Opt_C="G")
  if(fNewtR>0) call OutStrVec(fNewtR,ABS(gp(1:N)),Opt_I=fNewt_I,Opt_J=Itn,Opt_C="G")
!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!
!!ADDED J.MOUTTE end
!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!
  !
  write(ipr,801) itn
  write(ipr,802)
  write(ipr,803) xp(1:n)
  write(ipr,804)
  write(ipr,805) fval
  write(ipr,806)
  write(ipr,807) gp(1:n)
  !
  801 format(/,'  TSRSLT    ITERATION K  = ',i5)
  802 format('  TSRSLT    X(K)')
  803 format(100('  TSRSLT    ', 3(g20.13, "   "),/))
  804 format('  TSRSLT    function AT X(K)')
  805 format('  TSRSLT       ', g20.13)
  806 format('  TSRSLT    GRADIENT AT X(K)')
  807 format(100('  TSRSLT    ', 3(g20.13, "   "),/))
  return
end subroutine tsrslt

subroutine tsq1p1(aja, anls, s, f, n, root, restns)
  !SOLVES THE LOWER M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
  !OF THE TENSOR MODEL WHEN Q= 1 AND P= 1.
  real(dp),intent(in)   :: aja(:,:)
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: s(:,:)
  real(dp),intent(in)   :: f(:)
  integer, intent(in)     :: n
  real(dp),intent(out)  :: root
  real(dp),intent(out)  :: restns
  !INPUT parameterS :
  !    AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
  !    ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
  !    S    : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !    F    : function value AT CURRENT ITERATE MULTIPIED BY AN
  !           ORTHOGONAL MATRIX
  !    N    : COLUMN dimension OF AJA
  !OUTPUT parameterS :
  !    ROOT   : SOLUTION TO THE SYSTEM
  !    RESTNS : TENSOR RESIDUAL
  real(dp):: delta, t1, t2
  ! find the roots of the equation:
  ! F(N) + AJA(N,N)*D + 0.5*ANLS(N,1)*(S(N+2,1)*D)**2
  t1= aja(n,n)
  t2= anls(n,1) * s(n+2,1)**2
  if(anls(n,1)== zero) then
    root= -f(n)/t1
  else
    delta= t1**2 - two*f(n)*t2
    if(delta >= zero) then
      root= (-t1 + SIGN(one,t1) * SQRT(delta))/t2
    else
      root= -t1/t2
    end if
  end if
  ! compute tensor residual
  restns= ABS(f(n) + aja(n,n)*root + half*anls(n,1)*(s(n+2,1)**2)* (root**2))
  return
end subroutine tsq1p1

subroutine tsqfcn(p, x, sum)
  ! THIS ROUTINE IS useD TO EVALUATE THE RESIDUAL OF THE LAST M-N+P
  ! QUADRATIC EQUATIONS IN P UNKNOWNS OF THE TENSOR MODEL. NOTE THAT
  ! THIS ROUTINE IS callED BY UNCMIN TO SOLVE THE NONLINEAR LEAST SQUARES
  ! PART OF THE TENSOR MODEL.
  integer, intent(in)     :: p
  real(dp),intent(in)   :: x(:)
  real(dp),intent(out)  :: sum
  !INPUT parameterS :
  !    P : dimension OF THE PROBLEM SOLVED BY UNCMIN
  !INPUT-OUTPUT parameterS :
  !    X : NULL VECTOR ON entry AND APPROXIMATION OF THE SOLUTION
  !        TO THE SYSTEM OF M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
  !        OF THE TENSOR MODEL ON exit
  !    OUTPUT parameterS :
  !    SUM : RESIDUAL OF THE LAST M-N+P QUADRATIC EQUATIONS IN P
  !          UNKNOWNS OF THE TENSOR MODEL
  !SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  !    LEVEL 2 BLAS  ...  DGEMV
  !    TENSOLVE      ...  TSSTMX
  real(dp):: wrk1(meqns), wrk2(p), wrk3(p), wrk4(meqns), wrk5(meqns)
  integer :: i
  real(dp):: small
  small= 4.0_dp * SQRT( TINY(1.0_dp) )
  ! compute the lower right (m-n+q) x p submatrix of AJA times X
  call dgemv('N', meqns-nvars+qrank, p, one, aja(nvars-qrank+1:,nvars-p+1:), &
             meqns, x, 1, zero, wrk1, 1)
  ! compute S-trans times X
  call tsstmx(s, x, nvars, p, wrk3)
  ! compute 0.5 * (S-trans times X)**2
  do i= 1, p
    if (ABS(wrk3(i)) > small) then
      wrk2(i)= half * wrk3(i)**2
    else
      wrk2(i)= zero
    end if
  enddo
  ! compute 0.5 * (down (m-n+q) x p submatrix of ANLS) * (S-trans times X)**2
  call dgemv('N', meqns-nvars+qrank, p, one, anls(nvars-qrank+1:,:), meqns,  &
             wrk2, 1, zero, wrk4, 1)
  do i= 1,meqns-nvars+qrank
    wrk5(i)= wrk4(i) + fc(nvars-qrank+i) + wrk1(i)
  enddo
  sum= half*dnrm2(meqns-nvars+qrank, wrk5, 1)**2
  return
end subroutine tsqfcn

subroutine tsqlfc(ql, m, n, ierr)
!PERFORMS A QL DECOMPOSITION.
  real(dp),intent(inout)  :: ql(:,:)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  integer, intent(out)       :: ierr
  !INPUT parameterS :
  !     M   : ROW dimension OF QL
  !     N   : COLUMN dimension OF QL
  !INPUT-OUTPUT parameterS :
  !     QL : INPUT MATRIX ON entry AND FACTORED MATRIX ON exit
  !OUTPUT parameterS :
  !     IERR : return CODE WITH THE FOLLOWING MEANINGS :
  !            IERR= 1 : SINGULARITY DETECTED
  !            IERR= 0 : OTHERWISE
  integer   :: i, j, k
  real(dp) :: nu, sigma, tau, rnu
  ! zero out rows m+1 and m+2 of matrix QL
  do j= 1,n
    ql(m+1,j)= zero
    ql(m+2,j)= zero
  enddo
  ierr= 0
  k= 1
  20 if(k < m .and. k <= n) then
    ! find NU= max element of col K on or above l-diagonal
      nu= zero
    do i= 1, m+1-k
      nu= MAX(nu, ABS(ql(i,k)))
    enddo
    if(nu /= zero) then
      ! normalize col K on or above l-diagonal
      rnu= one/nu
      ql(1:m-k+1,k)= rnu * ql(1:m-k+1,k)
      ! code to find SIGMA= SGN(QL(M+1-K,K))*l2-norm of col K on or
      ! above l-diagonal
      sigma= dnrm2(m-k+1, ql(:,k), 1)
      sigma= SIGN(sigma, ql(m+1-k,k))
      ! store last element(1st in normal QR) of U-vector in QL(M+1-K,K)
      ql(m+1-k,k)= ql(m+1-k,k)+sigma
      ! store value of <U,U>/2 in QL(M+1,K)
      ql(m+1,k)= sigma*ql(m+1-k,k)
      if(ql(m+1,k)== zero) then
        ierr= 1
        return
      end if
      ! store L(M+1-K,K) in QL(M+2,K)
      ql(m+2,k)= -nu*sigma
      ! code to get (I-2U*UT/<U,U>)*QL for kth iteration
      if(k < n) then
        do j= k+1,n
          ! loop to get TAU= <U,J-TH COL OF QL>
          tau= dot_product( ql(1:m-k+1,k), ql(1:m-k+1,j) ) - tau/ql(m+1,k)
          ! loop to get (I-2U*UT/<U,U>)*j-th col of QL
          ql(1:m-k+1,j)= ql(1:m-k+1,j) + tau * ql(1:m-k+1,k)
        enddo
      end if
      k= k + 1
    else
      ierr= 1
      return
    end if
      GO TO 20
    end if
  if(m== n) ql(m+2,n)= ql(1,n)
  if(ql(m+2,n)== zero) then
    ierr= 1
  end if
  return
end subroutine tsqlfc

subroutine tsqmlv(n, p, q, v, w, trans)
!MULTIPLIES AN ORTHOGONAL MATRTIX Q OR ITS TRANSPOSE BY A VECTOR.
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  real(dp),intent(in)   :: q(:,:)
  real(dp),intent(in)   :: v(:)
  real(dp),intent(out)  :: w(:)
  logical, intent(in)     :: trans
  !INPUT parameterS :
  !    N  : dimension OF VECTORS V AND W
  !    P  : COLUMN dimension OF MATRIX Q
  !    Q  : ORTHOGONAL MATRIX (OBTAINED FROM TSQLFC subroutine)
  !    V  : VECTOR TO BE MULTIPLIED BY Q
  !    TRANS : BOOLEAN parameter:
  !           =.true. : PERFORM Q-TRANS*V
  !           =.false.: PERFORM Q*V
  !OUTPUT parameterS :
  !    W  : VECTOR Q*V OR Q-TRANS*V ON exit
  integer   :: j, k
  real(dp) :: tau, const
  w(1:n)= v(1:n)
  do j= 1,p
    if(trans) then
      k= p + 1 - j
    else
      k= j
    end if
    tau= dot_product( q(1:n-k+1,k), w(1:n-k+1) )
    const= -tau/q(n+1,k)
    w(1:n-k+1)= w(1:n-k+1) + const * q(1:n-k+1,k)
  enddo
  return
end subroutine tsqmlv

subroutine tsqmts(anls, qhat, mj, m, p, lb, zero1)
!MULTIPLIES AN ORTHOGONAL MATRIX QHAT BY THE TENSOR MATRIX ANLS.
  real(dp),intent(inout):: anls(:,:)
  real(dp),intent(in):: qhat(:,:)
  integer, intent(in):: mj
  integer, intent(in):: m
  integer, intent(in):: p
  integer, intent(in):: lb
  integer, intent(in):: zero1
  !INPUT parameterS :
  !   QHAT : ORTHOGONAL MATRIX (OBTAINED FROM TSQRFC subroutine)
  !     MJ : ROW dimension OF QHAT
  !     M  : ROW dimension OF MATRIX ANLS
  !     P  : COLUMN dimension OF MATRIX ANLS
  !     LB : STARTING COLUMN FROM WHICH QR DECOMPOSITION WAS PERFORMED
  !   ZERO1: FIRST ZERO COLUMN OF MATRIX QHAT IN case OF SINGULARITY
  !INPUT-OUTPUT parameterS :
  !     ANLS : MATRIX TO BE MULTIPLIED BY AN ORTHOGONAL MATRIX
  !     ON entry AND THE MATRIX QHAT*ANLS ON exit
  !     SUBprogramS callED:
  !     TENSOLVE      ...  TSQMUV
  !Workspace
  real(dp) :: wrk1(m)
  integer   :: j
  do j= 1,p
    call tsqmuv(qhat, anls(:,j), wrk1, mj, lb, zero1, .false.)
    anls(1:m,j)= wrk1(1:m)
  enddo
  return
end subroutine tsqmts

subroutine tsqmuv(q, v, w, m, lb, zero1, transp)
!MULTIPLIES AN ORTHOGONAL MATRIX BY A VECTOR.
  real(dp),intent(in) :: q(:,:)
  real(dp),intent(in) :: v(:)
  real(dp),intent(out):: w(:)
  integer, intent(in) :: m
  integer, intent(in) :: lb
  integer, intent(in) :: zero1
  logical, intent(in) :: transp
  !INPUT parameterS :
  !      Q  : FACTORED MATRIX (OBTAINED FROM subroutine TSQRFC)
  !      V  : VECTOR TO BE MULTIPLIED BY THE ORTHOGONAL MATRIX Q
  !      M  : ROW dimension OF MATRIX Q
  !      LB : STARTING COLUMN FROM WHICH QR DECOMPOSITION WAS PERFORMED
  !    ZERO1: FIRST ZERO COLUMN OF THE MATRIX Q
  !  TRANSP : BOOLEAN parameter :
  !              =.true. : PERFORM Q-TRANS*V
  !              =.false.: PERFORM Q*V
  !OUTPUT parameterS :
  !      W : Q*V OR Q-TRANS*V ON exit
  integer   :: ub, a, b, c, k
  real(dp) :: tau, const
  ! copy the vector V to W
  !! print '(A,3I3)',"v,w,m",size(W),size(v),m ; pause
  w(1:m)= v(1:m)
  ub= zero1-1
  if(transp) then
    a= ub
    b= lb
    c= -1
  else
    a= lb
    b= ub
    c= 1
  end if
  do k= a,b,c
    tau= dot_product( q(k:m,k), w(k:m) )
    const= -tau/q(m+1,k)
    w(k:m)= w(k:m) + const * q(k:m,k)
  enddo
  return
end subroutine tsqmuv

subroutine tsqrfc(qr, n, m, lb, ub, ierr, epsm, curpos, zero1)
!PERFORMS COLUMN-PIVOTED QR DECOMPOSITION ON AN M*N MATRIX.
!THE DECOMPOSITION IS doNE BETWEEN COLS LB AND UB.
  real(dp),intent(inout)  :: qr(:,:)
  integer, intent(in)        :: n
  integer, intent(in)        :: m
  integer, intent(in)        :: lb
  integer, intent(in)        :: ub
  integer, intent(out)       :: ierr
  real(dp),intent(in)      :: epsm
  integer, intent(out)       :: curpos(:)
  integer, intent(out)       :: zero1
  !INPUT parameterS :
  !      N  : COLUMN dimension OF MATRIX QR
  !      M  : ROW dimension OF MATRIX QR
  !   LB,UB : SUBSPACE OF QR DECOMPOSITION
  !   EPSM  : MACHINE PRECISION
  !INPUT-OUTPUT parameterS :
  !      QR  : INPUT MATRIX ON entry AND FACTORED MATRIX ON exit
  !OUTPUT parameterS :
  !      IERR : return CODE WITH TH FOLLOWING MEANINGS:
  !             IERR =  1 : SINGULARITY DETECTED
  !             IERR =  0 : OTHERWISE
  !      CURPOS :  PIVOT VECTOR
  !      ZERO1  :  FIRST ZERO COLUMN OF MATRIX QR IN case OF SINGULARITY
  !      SUBprogramS callED:
  !      LEVEL 1 BLAS  ...  DNRM2,DSWAP,IDAMAX
  real(dp) :: al2nrm(m)
  integer   :: colpiv, i, j, k, l
  real(dp) :: colmax, sigma, tau, amax
  real(dp) :: nu, rnu
  ! zero rows m+1 and m+2 of QR matrix
  do j= 1,n
    curpos(j)= j
  enddo
  do j= lb,ub
    qr(m+1,j)= zero
    qr(m+2,j)= zero
  enddo
  ierr= 0
  zero1= ub+1
  k= lb
  !  get L2NORM**2 of columns (LB to UB)
  do j= k,ub
    al2nrm(j)= dnrm2(m-k+1, qr(k:,j), 1)**2
  enddo
  40 if(k < m .and. k <= ub) then
      amax= zero
    do j= k,ub
      if(al2nrm(j) >= amax) then
        amax= al2nrm(j)
        colpiv= j
      end if
    enddo
      if(amax== zero) then
      ierr= 1
      zero1= k
      return
    end if
      if(k== lb) then
      colmax= amax
    end if
      if(al2nrm(colpiv) <= epsm*colmax) then
      ierr= 1
      zero1= k
      return
    end if
      if(colpiv /= k) then
      call dswap(m+2, qr(:,colpiv), 1, qr(:,k), 1)
      l= curpos(k)
      curpos(k)= curpos(colpiv)
      curpos(colpiv)= l
      call dswap(1, al2nrm(colpiv:), 1, al2nrm(k:), 1)
    end if
    ! find NU= max element of col K on or below diagonal
      l= idamax(m-k+1, qr(k:,k), 1) + k - 1
    nu= ABS(qr(l,k))
      if(nu== zero) then
      ierr= 1
      zero1= k
      return
    end if
    ! normalize col K on or below diagonal
    rnu= one/nu
    qr(k:m,k)= rnu * qr(k:m,k)
    ! code to find SIGMA= SGN(QR(K,K))*l2-norm of col K on or
    ! below diagonal
    sigma= dnrm2(m-k+1, qr(k:,k), 1)
    sigma= SIGN(sigma, qr(k,k))
    ! store 1st element of U-vector in QR(K,K)
    qr(k,k)= qr(k,k)+sigma
    ! store value of <U,U>/2 in QR(M+1,K)
    qr(m+1,k)= sigma*qr(k,k)
    if(qr(m+1,k)== zero) then
      ierr= 1
      zero1= k
      return
    end if
    ! store R(K,K) in QR(M+2,K)
    qr(m+2,k)= -nu*sigma
    if(qr(m+2,k)== zero) then
      ierr= 1
      zero1= k
      return
    end if
    ! code to get (I-2U*UT/<U,U>)*QR for kth iteration
    if(k < n) then
      do j= k+1,n
        ! loop to get UT*J-TH col of QR
              tau= -dot_product( qr(k:m,k), qr(k:m,j) ) / qr(m+1,k)
        ! loop to get (I-2U*UT/<U,U>)*j-th col of QR
              qr(k:m,j)= qr(k:m,j) + tau * qr(k:m,k)
      enddo
    end if
    ! update l2norm**2 (K+1 to UB)
    do i= k+1,ub
      al2nrm(i)= al2nrm(i) - qr(k,i)**2
    enddo
      k= k+1
    GO TO 40
    else
      if(lb== ub) colmax= al2nrm(lb)
    end if
  if(m== ub) qr(m+2,ub)= qr(m,m)
  if(ABS(qr(m+2,ub)) <= epsm*colmax) then
    ierr= 1
    zero1= ub
  end if
  return
end subroutine tsqrfc

subroutine tsrsid( &
!COMPUTES || F + J*D + 0.5*A*D**2 ||**2 IN THE L2 NORM SENS AT A GIVEN STEP D.
& itn, method, fq, d, curpos, pivot, pbar, aja, anls,  &
& shat, flag, nwtake, ierr, maxm, m, n, p, scres)
  integer, intent(in):: itn
  integer, intent(in):: method
  real(dp),intent(in):: fq(:)
  real(dp),intent(in):: d(:)
  integer, intent(in):: curpos(:)
  integer, intent(in):: pivot(:)
  integer, intent(in):: pbar(:)
  real(dp),intent(in):: aja(:,:)
  real(dp),intent(in):: anls(:,:)
  real(dp),intent(in):: shat(:,:)
  integer, intent(in):: flag
  logical, intent(in):: nwtake
  integer, intent(in):: ierr
  integer, intent(in):: maxm
  integer, intent(in):: m
  integer, intent(in):: n
  integer, intent(in):: p
  real(dp),intent(out):: scres
  !INPUT parameterS
  !     ITN   : CURRENT ITERATION NUMBER
  !     METHOD: METHOD TO BE useD
  !     FQ    : function value AT CURRENT ITERATE MULTIPLIED BY
  !             ORTHOGONAL MATRICES
  !     D     : STEP AT WHICH TO EVALUATE THE TENSOR MODEL
  !     CURPOS: PIVOT VECTOR (useD DURING THE FACTORIZATION OF AJA
  !             FROM COLUMN 1 TO N-P)
  !     PIVOT : PIVOT VECTOR ( useD DURING THE FACTORIZATION OF AJA
  !             FROM COLUMN N-P+1 TO N)
  !     PBAR  : PIVOT VECTOR (useD DURING THE FACTORIZATION OF AJA
  !             if IT IS SINGULAR
  !     AJA   : JACOBIAN MATRIX AT CURRENT ITERATE
  !     ANLS  : TENSOR TERM MATRIX AT CURRENT ITERATE
  !     SHAT  : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS AFTER
  !             A QL FACTORIZATION
  !     FLAG  : return CODE WITH THE FOLLOWING MEANINGS:
  !             FLAG= 0 : NO SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN 1 TO N
  !             FLAG= 1 : SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN 1 TO N-P
  !             FLAG= 2 : SINGULARITY DETECTED DURING FACTORIZATION
  !                        OF THE JACOBIAN FROM COLUMN N-P+1 TO N
  !     NWTAKE: logical VARIABLE WITH THE FOLLOWING MEANINGS:
  !             NWTAKE= .true.  : NEWTON STEP TAKEN
  !             NWTAKE= .false. : TENSOR STEP TAKEN
  !     IERR  : return CODE FROM QRP FACTORIZATION ROUTINE:
  !             IERR= 0 : NO SINGULARITY DETECTED
  !             IERR= 1 : SINGULARITY DETECTED
  !    MAXM   : LEADING dimension OF AJA AND ANLS
  !    M,N    : dimensionS OF THE PROBLEM
  !    P      : COLUMN dimension OF THE MATRICES SHAT AND ANLS
  !OUTPUT parameterS
  !     SCRES :  || F + J*D + 0.5*A*D**2 ||**2
  !    SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  !    LEVEL 2 BLAS  ...  DGEMV
  !    TENSOLVE      ...  TSJMUV,TSUDQV
  real(dp):: wrk1(n+m), vn(n+m), vnp(m), vns(m)
  integer :: len
  call tsjmuv( &
  & itn, method, d, curpos, pivot, pbar, aja, shat, flag,  &
  & ierr, m, n, p, vns, wrk1)
  wrk1(n+1:n+m)= zero
  len= m
  if(ierr > 0) len= m + n
  vn(1:len)= wrk1(1:len) + fq(1:len)
  if( .not. nwtake) then
    call tsudqv(shat, vns, n, p, vnp)
    vnp(1:p)= vnp(1:p) ** 2
    call dgemv('N', len, p, half, anls, maxm, vnp, 1, one, vn, 1)
  end if
  scres= dnrm2(len, vn, 1)
  return
end subroutine tsrsid

subroutine tssclf(f,df,m)
!SCALES A function value F.
  real(dp),intent(inout) :: f(:)
  real(dp),intent(in):: df(:)
  integer, intent(in):: m
  !INPUT parameterS :
  !    DF : DIAGONAL SCALING MATRIX FOR F
  !    M  : dimension OF F
  !INPUT-OUTPUT parameterS :
  !    F  : UNSCALED function value ON entry AND SCALED function
  !         value ON exit
  f(1:m)= df(1:m)*f(1:m)
  return
end subroutine tssclf

subroutine tssclj(x, dx, typx, f, df, m, n, epsm, jacflg, fvec, jac, aja)
!COMPUTES THE JACOBIAN MATRIX AT THE CURRENT ITERATE X AND SCALES ITS value.
  real(dp),intent(inout):: x(:)
  real(dp),intent(in):: dx(:)
  real(dp),intent(in):: typx(:)
  real(dp),intent(inout):: f(:)
  real(dp),intent(in):: df(:)
  integer, intent(in):: m
  integer, intent(in):: n
  real(dp),intent(in):: epsm
  integer, intent(in):: jacflg
  real(dp),intent(inout):: aja(:,:)
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
    end subroutine fvec
    subroutine jac(x, aja, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in) :: x(:)
      real(dp),intent(out):: aja(:,:)
      integer, intent(in) :: m, n
    end subroutine jac
  end interface
  !INPUT parameterS :
  !     X    : SCALED CURRENT ITERATE
  !     DX   : DIAGONAL SCALING MATRIX FOR X
  !     F    : SCALED function value AT X
  !     DF   : DIAGONAL SCALING MATRIX FOR F
  !     M,N  : dimensionS OF PROBLEM
  !     EPSM : MACHINE PRECISION
  !   JACFLG : JACOBIAN FLAG
  !     FVEC : subroutine TO EVALUATE function
  !     JAC  : subroutine TO EVALUATE ANALYTIC JACOBIAN
  !INPUT-OUTPUT parameterS :
  !     AJA  : SCALED JACOBIAN AT CURRENT ITERATE
  !SUBprogramS callED:
  !     TENSOLVE      ...  TSUNSX,TSUNSF,TSFDFJ,TSSCLF,TSSCLX
  !     useR          ...  FVEC,JAC
  integer :: i,j
  ! unscale X AND F
  call tsunsx(x, dx, n)
  call tsunsf(f, df, m)
  ! compute the finite difference or analytic Jacobian at X
  if(jacflg== 0) then
    call tsfdfj(x, f, m, n, epsm, fvec, aja)
  else
    call jac(x, aja, m, n)
  end if
  ! scale the Jacobian matrix
  do j= 1,n
    do i= 1,m
      aja(i,j)= aja(i,j)*df(i)*typx(j)
    enddo
  enddo
  ! scale back X AND F
  call tssclf(f, df, m)
  call tssclx(x, dx, n)
  return
end subroutine tssclj

subroutine tssclx(x, dx, n)
!SCALES A VECTOR X.
  real(dp),intent(inout)  :: x(:)
  real(dp),intent(in)      :: dx(:)
  integer, intent(in)        :: n
  !INPUT parameterS :
  !    DX : DIAGONAL SCALING MATRIX FOR X
  !    N  : dimension OF X
  !OUTPUT parameterS :
  !    X  : SCALED VECTOR X
  x(1:n)= dx(1:n)*x(1:n)
  return
end subroutine tssclx

subroutine tsslct(residt, residn, itrmcd, fc, m, n, dn, dt, g, df, nwtake)
! THIS ROUTINE DECIDES WHICH DIRECTION TO CHOOSE: THE TENSOR OR THE
! STANDARD DIRECTION. THE STANDARD DIRECTION IS CHOSEN WHENEVER
! THE TENSOR DIRECTION IS NOT DESCENT OR THE TENSOR DIRECTION IS TO A
! MINIMIZER OF THE TENSOR MODEL AND doESN'T PROVIDE ENOUGH DECREASE
! IN THE TENSOR MODEL, OR THE QUADRATIC SYSTEM OF Q EQUATIONS IN P
! UNKNOWNS CANNOT BE SOLVED BY UNCMIN WITHIN 150 ITERATIONS.
  real(dp),intent(in)   :: residt
  real(dp),intent(in)   :: residn
  integer, intent(in)     :: itrmcd
  real(dp),intent(in)   :: fc(:)
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  real(dp),intent(in)   :: dn(:)
  real(dp),intent(in)   :: dt(:)
  real(dp),intent(in)   :: g(:)
  real(dp),intent(out)  :: df(:)
  logical, intent(out)    :: nwtake
  !INPUT parameterS :
  !    RESIDT : TENSOR RESIDUAL
  !    RESIDN : NEWTON RESIDUAL
  !    ITRMCD : UNCMIN TERMINATION CODE
  !    FC : function value AT XC
  !    M,N: dimensionS OF PROBLEM
  !    DN : STANDARD STEP
  !    DT : TENSOR STEP
  !    G  : GRADIENT value AT XC
  !OUTPUT parameterS :
  !    DF     : EITHER THE STANDARD OR TENSOR STEP ON exit
  !    NWTAKE : BOOLEAN value WITH THE FOLLOWING MEANINGS:
  !           =.true.  : STANDARD STEP IS TAKEN
  !           =.false. : TENSOR STEP IS TAKEN
  !    SUBprogramS callED:
  !    LEVEL 1 BLAS  ....  DNRM2
  real(dp) :: anrmfc, dtnorm, gnorm
  real(dp) :: temp, temp1, beta, gama
  real(dp), parameter:: tenth= 0.1_dp, onett= 1.0e-04_dp
  nwtake= .false.
  anrmfc= dnrm2(m, fc, 1)
  dtnorm= dnrm2(n, dt, 1)
  gnorm= dnrm2(n, g, 1)
  temp= dot_product( dt(1:n), g(1:n) )
  gama= half
  if(m > n) then
    beta= tenth
  else
    beta= onett
  end if
  temp1= -beta*dtnorm*gnorm
  if(residt >= gama*(anrmfc + residn) .or. temp > temp1 .or. itrmcd== 4) then
    df(1:n)= dn(1:n)
    nwtake= .true.
  else
    df(1:n)= dt(1:n)
  end if
  return
end subroutine tsslct

subroutine tssmin( &
& anls, fq, adt, ag, const1, const2, dlt, m, n,  &
& p, nwtake, ierr, epsm, sol)
!THIS ROUTINE MINIMIZES THE TENSOR MODEL OVER THE SUBSPACE SPANNED BY
!THE TENSOR STEP AND THE STEEPEST DECENT DIRECTION. if THE NEWTON STEP
!WERE CHOSEN, IT WILL MINIMIZE THE NEWTON MODEL OVER THE SUBSPACE
!SPANNED BY THE NEWTON STEP AND THE STEEPEST DESCENT DIRECTION.
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: fq(:)
  real(dp),intent(in)   :: adt(:)
  real(dp),intent(in)   :: ag(:)
  real(dp),intent(in)   :: const1(:)
  real(dp),intent(in)   :: const2(:)
  real(dp),intent(in)   :: dlt
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  logical, intent(in)     :: nwtake
  integer, intent(in)     :: ierr
  real(dp),intent(in)   :: epsm
  real(dp),intent(out)  :: sol
  !INPUT parameterS :
  !     ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
  !     FQ   : function value AT CURRENT ITERATE MULTIPLIED BY
  !            ORTHOGONAL MATRICES
  !     ADT  : JACOBIAN MATRIX TIMES DT (SEE subroutine TS2DTR)
  !      AG  : JACOBIAN MATRIX TIMES GBAR (SEE subroutine TS2DTR)
  !    CONST1: SHAT-TRANS * DT  (SEE subroutine TS2DTR)
  !    CONST2: SHAT-TRANS * GBAR (SEE subroutine TS2DTR)
  !    ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
  !      DLT : CURRENT TRUST RADIUS
  !       M,N: dimensionS OF THE PROBLEM
  !         P: COLUMN dimension OF MATRIX ANLS
  !   NWTAKE : logical VARIABLE WITH THE FOLLOWING MEANINGS :
  !            NWTAKE= .true.  : STANDARD STEP TAKEN
  !            NWTAKE= .false. : TENSOR STEP TAKEN
  !   IERR   : return CODE FROM QRP FACTORIZATION ROUTINE
  !            IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !            IERR= 1 : SINGULARITY OF JACOBIAN DETECTED
  !   EPSM   : MACHINE PRECISION
  !OUTPUT parameterS :
  !     SOL   : GLOBAL MINIMIZER OF THE ONE VARIABLE function IN ALPHA
  !             ||F + J*(ALPHA*DT + BETA*GBAR) + 0.5*A*(ALPHA*DT +
  !             BETA*GBAR)**2||**2 where BETA= SQRT(DLT**2 - ALPHA**2)
  !SUBprogramS callED:
  !     TENSOLVE      ...  TSFAFA,TSMFDA,TSLMIN
  real(dp):: vn(n+m)
  integer :: INT
  real(dp):: tol, dl, s, sp, c, a
  real(dp):: d, s1, b, q, bc, optim, ac, glopt, bloop, eloop, incr
  real(dp),parameter:: ohund= 0.01_dp, tenth= 0.1_dp
  logical   :: first
  !
  first= .true.
  tol= epsm**(two/three)
  INT= 40
  dl= tol
  if(dlt < tol) then
    dl= tol*tenth
  else if(dlt > tol*ten) then
    dl= tol*ten
  end if
  if(dlt <= ohund) then
    INT= 10
  end if
  ! find global minimizer of FALPHA
  bloop= -dlt+dl
  eloop= dlt*(INT-2)/INT
  incr= two*(dlt-dl)/INT
  s= bloop
  10 sp= s
  s1= s+incr
  ! evaluate FALPHA(SP) and the derivative of FALPHA at SP
  if(first) then
    call tsfafa( &
    & anls, fq, adt, ag, const1, const2, sp, dlt, m, n, p,  &
    & nwtake, ierr, vn, c)
    a= tsmfda( &
    & anls, adt, ag, const1, const2, sp, dlt, m, n, p, nwtake,  &
    & ierr, vn)
  else
    c= d
    a= b
  end if
  ! evaluate FALPHA(S1) and the derivative of FALPHA at S1
  call tsfafa(anls, fq, adt, ag, const1, const2, s1, dlt, m, n, p, nwtake, &
             ierr, vn, d)
  b= tsmfda(anls, adt, ag, const1, const2, s1, dlt, m, n, p, nwtake, ierr, vn)
  ! minimize FALPHA in the subinterval [SP,S1]
  if(a <= zero .and. b >= zero) then
    if(c > d) then
      q= d
      bc= b
      call tslmin( &
      & s1, sp, bc, q, anls, fq, adt, ag, const1, const2,  &
      & dlt, m, n, p, nwtake, ierr, tol, optim)
    else
      q= c
      ac= a
      call tslmin( &
      & sp, s1, ac, q, anls, fq, adt, ag, const1, const2,  &
      & dlt, m, n, p, nwtake, ierr, tol, optim)
    end if
  else if(a <= zero .and. b <= zero) then
    if(c <= d) then
      q= c
      ac= a
      call tslmin( &
      & sp, s1, ac, q, anls, fq, adt, ag, const1, const2,  &
      & dlt, m, n, p, nwtake, ierr, tol, optim)
    else
      optim= s1
      q= d
    end if
  else if(a >= zero .and. b >= zero) then
    if(c >= d) then
      q= d
      bc= b
      call tslmin( &
      & s1, sp, bc, q, anls, fq, adt, ag, const1, const2,  &
      & dlt, m, n, p, nwtake, ierr, tol, optim)
    else
      optim= sp
      q= c
    end if
  end if
  ! update the global minimizer of FALPHA so far
  if(first) then
    if(a > zero .and. b < zero) then
      glopt= c
      sol= sp
      if(c > d) then
        glopt= d
        sol= s1
      end if
    else
      glopt= q
      sol= optim
    end if
    first= .false.
  else if(glopt >= q) then
    glopt= q
    sol= optim
  end if
  s= s + incr
  if(s <= eloop) GO TO 10
  return
end subroutine tssmin

subroutine tssmrd(vect, resnew, x, mu, ierr, m, n)
!COMPUTES THE RESIDUAL OF THE STANDARD MODEL.
  real(dp),intent(in)   :: vect(:)
  real(dp),intent(out)  :: resnew
  real(dp),intent(in)   :: x(:)
  real(dp),intent(in)   :: mu
  integer, intent(in)     :: ierr
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  !INPUT parameterS :
  !    VECT : RIGHT HAND SIDE VECTOR OF THE NEWTON/GAUSS-NEWTON
  !           EQUATIONS AFTER BEING MULTIPLIED BY ORTHOGONAL MATRICES
  !           (SEE subroutine TSCPSS)
  !    X    : STANDARD STEP COMPUTED BY THE subroutine TSCPSS
  !    MU   : A SMALL PERTURBATION useD IN COMPUTING THE STANDARD
  !           STEP WHEN THE JACOBIAN IS SINGULAR
  !    IERR : return CODE WITH THE FOLLOWING MEANINGS :
  !           IERR= 0 : NO SINGULARITY OF JACOBIAN DETECTED
  !           IERR= 1 : OTHERWISE
  !    M,N  : dimension OF PROBLEM
  !OUTPUT parameterS :
  !    ------------------
  !    RESNEW : RESIDUAL OF THE STANDARD MODEL
  !    SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  real(dp) :: temp, prod
  if(ierr== 0) then
    if(m== n) then
      resnew= zero
    else
      resnew= dnrm2(m-n, vect(n+1:), 1)
    end if
  else
    temp= dnrm2(m, vect(n+1:), 1)**2
    prod= mu * dnrm2(n, x, 1)**2
    resnew= SQRT(temp-prod)
  end if
  return
end subroutine tssmrd

subroutine tssqp1(aja, anls, s, f, m, n, q, root, restns)
!SOLVES THE LOWER M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
!OF THE TENSOR MODEL WHEN Q > 1 AND P= 1.
  real(dp),intent(in)   :: aja(:,:)
  real(dp),intent(in)   :: anls(:,:)
  real(dp),intent(in)   :: s(:,:)
  real(dp),intent(in)   :: f(:)
  integer, intent(in)     :: m
  integer, intent(in)     :: n
  integer, intent(in)     :: q
  real(dp),intent(out)  :: root
  real(dp),intent(out)  :: restns
  !INPUT parameterS :
  !    AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
  !    ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
  !    S    : MATRIX OF PAST LINEARLY INDEPendENT DIRECTIONS
  !    F    : function value AT CURRENT ITERATE MULTIPLIED BY AN
  !           ORTHOGONAL MATRIX
  !    M,N  : ROW AND COLUMN dimensionS OF AJA
  !    Q    : NUMERICAL RANK OF JACOBIAN :
  !           Q > P : JACOBIAN IS SINGULAR
  !           Q= P : OTHERWISE
  !OUTPUT parameterS :
  !    ROOT   : SOLUTION TO THE SYSTEM
  !    RESTNS : TENSOR RESIDUAL
  !SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  integer   :: i
  real(dp) :: temp, a, b, c, d, res1, res2, res3, res, s1, s2, s3
  real(dp) :: t, t0, t1, t2, t3, a1, a2, a3, theta, onetrd
  real(dp), parameter:: quart= 0.25_dp, four= 4.0_dp, nine= 9.0_dp
  real(dp), parameter:: tseven= 27.0_dp, ffour= 54.0_dp
  real(dp), parameter:: pi= 3.141592653589793_dp
  ! compute the coefficients of a third degree polynomial
  onetrd= one/three
  a= zero
  b= zero
  c= zero
  temp= dnrm2(m-n+q, f(n-q+1:), 1)**2
  d= two * dot_product( aja(n-q+1:m,n), f(n-q+1:m) )
  t0= s(n+2,1)**2
  t1= t0**2
  do i= n-q+1,m
    t2= aja(i,n)
    t3= anls(i,1) * t0
    c= c + two * (t2**2 + f(i) * t3)
    b= b + three * t2 * t3
    a= a + anls(i,1)**2 * t1
  enddo
  ! compute the roots of the third degree polynomial
  if(a== zero) then
    if(b /= zero) then
      t0= SQRT(c**2 - four*b*d)
      t1= two*b
      s1= (-c + t0)/t1
      s2= (-c - t0)/t1
      res1= ABS(temp + d*s1 + half*c*(s1**2) + onetrd*b*(s1**3))
      res2= ABS(temp + d*s2 + half*c*(s2**2) + onetrd*b*(s2**3))
      if(res1 > res2) then
        root=  s2
        res =  res2
      else
        root=  s1
        res =  res1
      end if
      restns =  SQRT(res)
      return
    else if(c /= zero) then
      root= -d/c
      res= ABS(temp + d*root + half*c*(root**2))
      restns= SQRT(res)
      return
    else
      root= zero
      restns= SQRT(temp)
      return
    end if
  else if(d== zero) then
    root= zero
    restns= SQRT(temp)
    return
  end if
  t3= d
  a1= b/a
  a2= c/a
  a3= d/a
  t0= (three*a2 - a1**2)/nine
  t1= (nine*a1*a2 - tseven*a3 - two*a1**3)/ffour
  d= t0**3 + t1**2
  if(d > zero) then
    t2= t1 + SQRT(d)
    t= t1 - SQRT(d)
    if(t < zero) then
      t= -(-t)**onetrd
    else
      t= t**onetrd
    end if
    if(t2 < zero)then
      t2= -(-t2)**onetrd
    else
      t2= t2**onetrd
    end if
    s1= t2 + t - a1/three
    s3= s1
    s2= s1
  else
    t= t1/SQRT(-t0**3)
    theta= ACOS(t)
    theta= theta/three
    t= two*SQRT(-t0)
    s1= t*COS(theta) - a1/three
    s2= t*COS(theta + pi*two/three) - a1/three
    s3= t*COS(theta + pi*four/three) - a1/three
  end if
  ! compute the tensor residual for each root
  res1= ABS(temp + t3*s1 + half*c*(s1**2) + onetrd*b*(s1**3)+ quart*a*(s1**4))
  res2= ABS(temp + t3*s2 + half*c*(s2**2) + onetrd*b*(s2**3)+ quart*a*(s2**4))
  res3= ABS(temp + t3*s3 + half*c*(s3**2) + onetrd*b*(s3**3)+ quart*a*(s3**4))
  ! select the root that produces the smallest tensor residual
  res= res1
  root= s1
  if(res > res2) then
    res= res2
    root= s2
  end if
  if(res > res3) then
    res= res3
    root= s3
  end if
  restns= SQRT(res)
  return
end subroutine tssqp1

subroutine tssstp(aja, fn, m, n, epsm, iglobl, dn, fq, pivot, pbar, ierr)
!FINDS THE STANDARD STEP WHEN THE ITERATION NUMBER IS EQUAL TO 1
!OR THE INPUT parameter "METHOD" IS SET TO 0.
  real(dp),intent(inout):: aja(:,:)
  real(dp),intent(inout):: fn(:)
  integer, intent(in)   :: m
  integer, intent(in)   :: n
  real(dp),intent(in)   :: epsm
  integer, intent(in)   :: iglobl
  real(dp),intent(out)  :: dn(:)
  real(dp),intent(out)  :: fq(:)
  integer, intent(out)  :: pivot(:)
  integer, intent(out)  :: pbar(:)
  integer, intent(out)  :: ierr
  !INPUT parameterS :
  !    AJA   : JACOBIAN MATRIX AT CURRENT ITERATE
  !    FN    : function value AT CURRENT ITERATE
  !    M,N   : dimensionS OF PROBLEM
  !    EPSM  : MACHINE EPSILON
  !    IGLOBL: GLOBAL STRATEGY useD :
  !            = 0 : LINE SEARCH useD
  !            = 1 : 2-dimensionAL TRUST REGION useD
  !OUTPUT parameterS :
  !    DN : STANDARD STEP
  !    FQ : function value AT CURRENT ITERATE MULTIPLIED BY
  !           ORTHOGONAL MATRICES (THIS IS useD IN THE case where
  !           THE GLOBAL STRATEGY IS THE 2-dimensionAL TRUST REGION)
  !    PIVOT,PBAR : PIVOT VECTORS
  !    IERR : returnED CODE WITH THE FOLLOWING MEANING :
  !           IERR =  1 : SINGULARITY OF JACOBIAN DETECTED (ZERO1 IS useD TO
  !                        KEEP TRACK OF THE FIRST COLUMN where SINGULARITY IS
  !                        DETECTED)
  !           IERR =  0 : OTHERWISE
  !SUBprogramS callED:
  !    TENSOLVE      ...  TSQRFC,TSQMUV,TSBSLV,TSPRMV,TSCPMU
  real(dp) :: y(n), w(m+n)
  integer   :: zero1, zerotm, i, j
  real(dp) :: mu
  ! perform a QR factorization of AJA
  call tsqrfc(aja, n, m, 1, n, ierr, epsm, pivot, zero1)
  fn(1:m)= -fn(1:m)
  if(ierr== 0) then
    if(m== n) then
      zerotm= zero1-1
    else
      zerotm= zero1
    end if
    ! multiply (-FN) by the orthogonal matrix resulting from the QR
    ! decomposition of AJA
    call tsqmuv(aja, fn, w, m, 1, zerotm, .false.)
    ! solve AJA*DN =  W
    call tsbslv(aja, m, n, w, y)
    call tsprmv(dn, y, pivot, n, 0)
    if(iglobl== 1) then
      ierr= 0
      fq(1:m)= -w(1:m)
    end if
    return
  else
    ! AJA is singular
    call tsqmuv(aja, fn, w, m, 1, zero1, .false.)
    ! solve ( AJA-trans AJA + MU I ) DN= -AJA-trans FN
    ! put the diagonal elements stored in row m+2 of AJA into their
    ! proper positions and zero out the unwanted portions of AJA
    do j= 1, zero1-1
      aja(j,j)= aja(m+2,j)
      aja(j+1:m+n,j)= zero
    enddo
    do j= zero1, n
      aja(zero1:m+n,j)= zero
    enddo
    ! compute a perturbation MU
    call tscpmu(aja, n, epsm, mu)
    ! form the augmented Jacobian matrix by adding an nxn diag(mu) in
    ! the bottom of AJA
    do i= m+1, m+n
      aja(i,i-m)= mu
    enddo
    ! factorize the transformed matrix AJA from 1 to n and compute
    ! the standard step DN
    call tsqrfc(aja, n, m+n, 1, n, ierr, epsm, pbar, zero1)
    call tsqmuv(aja, w, fq, m+n, 1, n+1, .false.)
    call tsbslv(aja, m+n, n, fq, dn)
    call tsprmv(y, dn, pbar, n, 0)
    call tsprmv(dn, y, pivot, n, 0)
  end if
  if(iglobl== 1) then
    ierr= 1
    fq(1:m+n)= - fq(1:m+n)
  end if
  return
end subroutine tssstp

subroutine tsstmx(s, x, n, p, wrk2)
!COMPUTES SHAT-TRANS * X, where SHAT IS AN UPSIDE doWN
!TRIANGULAR MATRIX resultING FROM A QL FACTORIZATION OF A MATRIX
!A AND X IS A VECTOR.
  real(dp),intent(in)   :: s(:,:)
  real(dp),intent(in)   :: x(:)
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  real(dp),intent(out)  :: wrk2(:)
  !INPUT parameterS :
  !   SHAT  : UPSIDE doWN TRIANGULAR MATRIX resultING FROM A QL FACTORIZATION
  !   X     : VECTOR TO BE MULTIPLIED BY SHAT
  !   N     : ROW dimension OF THE MATRIX A
  !   P     : COLUMN dimension OF SHAT
  !OUTPUT parameterS :
  !   WRK2  :  SHAT-TRANS * X
  real(dp) :: wrk1(p)
  integer   :: col
  wrk1(1:p)= zero
  wrk2(1)= s(n+2,1) * x(p)
  if(p > 1) then
    wrk1(p)= s(n,2)
    wrk1(p-1)= s(n+2,2)
    wrk2(2)= dot_product( wrk1, x(1:p) )
    do col= 3, p
      wrk1(p-col+2:p)= s(n-col+2:n,col)
      wrk1(p-col+1)= s(n+2,col)
      wrk2(col)= dot_product( wrk1, x(1:p) )
    enddo
  end if
  return
end subroutine tsstmx

subroutine tstrud( &
!DECIDES WHETHER TO ACCEPT XPLS=X+SC AS THE NEXT ITERATE
!AND UPDATES THE TRUST REGION DLT.
& m, n, x, f, g, sc, nwtake, stepmx, steptl, dlt, mxtake, dxn, &
& dfn, fvec, scres, iretcd, xplsp, fplsp, fprev, xpls, fp, fpls)
  integer, intent(in)        :: m
  integer, intent(in)        :: n
  real(dp),intent(in)      :: x(:)
  real(dp),intent(in)      :: f
  real(dp),intent(in)      :: g(:)
  real(dp),intent(in)      :: sc(:)
  logical, intent(in)        :: nwtake
  real(dp),intent(in)      :: stepmx
  real(dp),intent(in)      :: steptl
  real(dp),intent(inout)  :: dlt
  logical, intent(out)       :: mxtake
  real(dp),intent(in)      :: dxn(:)
  real(dp),intent(in)      :: dfn(:)
  real(dp),intent(in)      :: scres
  integer, intent(inout)    :: iretcd
  real(dp),intent(inout)  :: xplsp(:)
  real(dp),intent(inout)  :: fplsp
  real(dp),intent(inout)  :: fprev(:)
  real(dp),intent(out)     :: xpls(:)
  real(dp),intent(out)     :: fp(:)
  real(dp),intent(out)     :: fpls
  interface
    subroutine fvec(x, f, m, n)
      use M_Kinds, only: dp
      implicit none
      real(dp),intent(in)  :: x(:)
      real(dp),intent(out) :: f(:)
      integer, intent(in)    :: m, n
    end subroutine fvec
  end interface
  ! parameterS
  ! ----------
  ! M,N          --> dimensionS OF PROBLEM
  ! X(N)         --> OLD ITERATE X[K-1]
  ! F            --> 0.50D0 * || FC ||**2
  ! G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
  ! SC(N)        --> CURRENT STEP
  ! NWTAKE       --> BOOLEAN,=.true. if INPUT STEP TAKEN
  ! STEPMX       --> MAXIMUM ALLOWABLE STEP size
  ! STEPTL       --> RELATIVE STEP size AT WHICH SUCCESSIVE ITERATES
  !                  CONSIDERED close ENOUGH TO TERMINATE ALGORITHM
  ! DLT         <--> TRUST REGION RADIUS
  ! MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM lenGTH useD
  ! DXN         --->DIAGONAL SCALING MATRIX FOR X
  ! DFN         --->DIAGONAL SCALING MATRIX FOR F
  ! FVEC        --->subroutine TO EVALUATE function
  ! IRETCD      <--> return CODE
  !                   =0 XPLS ACCEPTED AS NEXT ITERATE;
  !                       DLT TRUST REGION FOR NEXT ITERATION.
  !                   =1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE
  !                       BECAuse XPLS-X < SMALLEST ALLOWABLE STEP lenGTH.
  !                   =2 F(XPLS) TOO LARGE.  continue CURRENT ITERATION
  !                       WITH NEW REDUCED DLT.
  !                   =3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL
  !                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO continue
  !                       CURRENT ITERATION WITH NEW doUBLED DLT.
  ! XPLSP(N)    <--> WORKSPACE [value NEEDS TO BE RETAINED BETWEEN
  !                  SUCCESSIVE callS OF K-TH GLOBAL STEP]
  ! FPLSP       <--> [RETAIN value BETWEEN SUCCESSIVE callS]
  ! FPREV       ---> WORKING VECTOR
  ! XPLS(N)     <--  NEW ITERATE X[K]
  ! FP          <--  function value AT NEXT ITERATE
  ! FPLS        <--  function value AT NEW ITERATE, F(XPLS)
  !
  !SUBprogramS callED:
  !    LEVEL 1 BLAS  ...  DNRM2
  !    TENSOLVE      ...  TSFSCL
  integer :: i
  real(dp):: stepln, dltfp, slope, dltf, slp, pq, rln, dltmp
  real(dp),parameter :: tenth= 0.1_dp, znn= 0.99_dp
  !
  mxtake= .false.
  xpls(1:n)= x(1:n) + sc(1:n)
  stepln= dnrm2(n, sc, 1)
  call tsfscl(xpls, dxn, dfn, m, n, fvec, fp)
  fpls= half*dnrm2(m, fp, 1)**2
  dltf= fpls - f
  slope= dot_product( g(1:n), sc(1:n) )
  slp= half*scres**2 - f
  ! next statement added for case of compilers which do not optimize
  ! evaluation of next "if" statement (in which case fplsp could be undefined)
  if(iretcd== 4) fplsp= zero
  if(iretcd /= 3 .or. (fpls < fplsp .and. dltf <= 1.e-4_dp*slp)) GO TO 130
  !       reset XPLS to XPLSP and terminate global step
  iretcd= 0
  xpls(1:n)= xplsp(1:n)
  fpls= fplsp
  dlt= half*dlt
  fp(1:m)= fprev(1:m)
  GO TO 230
  !       FPLS too large
  130 if(dltf <= 1.e-4_dp*slp) GO TO 170
  pq= one
  rln= zero
  do i= 1,n
    rln= MAX(rln,ABS(sc(i)) / MAX(ABS(xpls(i)),one/pq))
  enddo
  if(rln >= steptl) GO TO 150
  !           cannot find satisfactory XPLS sufficiently distinct from X
  iretcd= 1
  GO TO 230
  !           reduce trust region and continue global step
  150 iretcd= 2
  dltmp= -slope*stepln/(two*(dltf-slope))
  if(dltmp >= tenth*dlt) GO TO 155
  dlt= tenth*dlt
  GO TO 230
  155 if(dltmp > half*dlt) then
    dlt= half*dlt
  else
    dlt= dltmp
  end if
  !         FPLS sufficiently small
  170 dltfp= half*scres**2-f
  if(iretcd== 2 .or. (ABS(dltfp-dltf) > tenth*ABS(dltf)  &
      .and. dltfp > slope) .or. nwtake .or. dlt > znn*stepmx) GO TO 210
  !           double trust region and continue global step
  iretcd= 3
  xplsp(1:n)= xpls(1:n)
  fplsp= fpls
  dlt= MIN(two*dlt,stepmx)
  fprev(1:m)= fp(1:m)
  GO TO 230
  !           accept XPLS as next iterate.  Choose new trust region.
  210 iretcd= 0
  if(dlt > znn*stepmx) mxtake= .true.
  if(dltf < tenth*dltfp) GO TO 220
  !             Decrease trust region for next iteration
  dlt= half*dlt
  GO TO 230
  !             Check whether to increase trust region for next iteration
  220 if(dltf <= .75_dp*dltfp) dlt= MIN(two*dlt, stepmx)
  230 return
end subroutine tstrud

subroutine tsudqv(shat, wrk1, n, p, const1)
! THIS ROUTINE COMPUTES SHAT-TRANS * WRK1, where SHAT IS AN UPSIDE
! doWN TRIANGULAR MATRIX resultING FROM A QL FACTORIZATION OF A
! MATRIX A AND WRK1 IS A VECTOR OF lenGTH N.
  real(dp),intent(in)   :: shat(:,:)
  real(dp),intent(in)   :: wrk1(:)
  integer, intent(in)     :: n
  integer, intent(in)     :: p
  real(dp),intent(out)  :: const1(:)
  ! INPUT parameterS
  !      SHAT : UPSIDE doWN TRIANGULAR MATRIX resultING FROM A QL FACTORIZATION
  !      WRK1 : VECTOR TO BE MULTIPLIED BY SHAT
  !      N    : dimension OF MATRIX A
  !      P    : COLUMN dimension OF SHAT
  ! OUTPUT parameterS
  !      CONST1 : SHAT * WRK1
  integer   :: j
  const1(1)= shat(n+2,1) * wrk1(n)
  if(p > 1) then
    const1(2)= shat(n,2) * wrk1(n) + shat(n+2,2) * wrk1(n-1)
  end if
  do j= 3,p
    const1(j)= dot_product( shat(n-j+2:n,j), wrk1(n-j+2:n) ) +  &
                shat(n+2,j)*wrk1(n-j+1)
  enddo
  return
end subroutine tsudqv

subroutine tsunsf(f,df,m)
!UNSCALES A function value F.
  real(dp),intent(inout):: f(:)
  real(dp),intent(in)   :: df(:)
  integer, intent(in)   :: m
  !INPUT parameterS :
  !    DF : DIAGONAL SCALING MATRIX FOR F
  !    M  : dimension OF F
  !INPUT-OUTPUT parameterS :
  !    F  : SCALED function value ON entry AND UNSCALED function
  !         value ON exit
  f(1:m)= f(1:m) / df(1:m)
  return
end subroutine tsunsf

subroutine tsunsx(x, dx, n)
!UNSCALES A VECTOR X.
  real(dp),intent(inout):: x(:)
  real(dp),intent(in)   :: dx(:)
  integer, intent(in)   :: n
  !INPUT parameterS :
  !    DX : DIAGONAL SCALING MATRIX FOR X
  !    N  : dimension OF X
  !OUTPUT parameterS :
  !    X  : UNSCALED VECTOR X
  x(1:n)= x(1:n) / dx(1:n)
  return
end subroutine tsunsx

subroutine tsupsf(fc, xc, xp, sqrn, itn, m, n, s, fv)
!UPDATE THE MATRIX S OF PAST DIRECTIONS AND THE MATRIX FV OF function valueS.
  real(dp),intent(inout):: fc(:)
  real(dp),intent(in):: xc(:)
  real(dp),intent(in):: xp(:)
  integer, intent(in):: sqrn
  integer, intent(in):: itn
  integer, intent(in):: m
  integer, intent(in):: n
  real(dp),intent(out):: s(:,:)
  real(dp),intent(inout):: fv(:,:)
  !INPUT parameterS :
  !    FC  : function value AT CURRENT ITERATE
  !    XC  : CURRENT ITERATE X[K-1]
  !    XP  : NEW ITERATE X[K]
  !    SQRN: MAXIMUM COLUMN dimension OF S AND FV
  !    ITN : ITERATION NUMBER
  !    M   : ROW dimension OF MATRIX FV
  !    N   : ROW dimension OF MATRIX S
  !    STEP: WORKING VECTOR
  !INPUT-OUTPUT parameterS :
  !    S   :  MATRIX OF PAST DIRECTIONS (XK - XC)
  !    FV  :  MATRIX OF PAST functionS valueS
  !
  integer   :: j, l
  real(dp) :: step(n)
  ! update FV
  if(sqrn < itn) then; l= sqrn
  else               ; l= itn
  end if
  do j= l-1, 1, -1
    fv(1:m,j+1)= fv(1:m,j)
  enddo
  fv(1:m,1)= fc(1:m)
  ! update S
  step(1:n)= xc(1:n) - xp(1:n)
  do j= l-1,1,-1
    s(1:n,j+1)= s(1:n,j) + step(1:n)
  enddo
  s(1:n,1)= step(1:n)
  return
end subroutine tsupsf

subroutine tsutmd(aja, d, m, n, v)
!MULTIPLIES AN UPPER TRIANGULAR MATRIX (AS STORED IN STEWART) BY A VECTOR D.
  real(dp),intent(in):: aja(:,:)
  real(dp),intent(in):: d(:)
  integer, intent(in):: m
  integer, intent(in):: n
  real(dp),intent(out):: v(:)
  !INPUT
  !  AJA  : JACOBIAN AT CURRENT ITERATE
  !  D    : VECTOR TO BE MULTIPLIED BY AJA
  !  M,N  : dimensionS OF PROBLEM
  !OUTPUT
  !  V  : VECTOR AJA * D ON exit
  integer :: j
  v(1)= aja(m+2,1) * d(1) + aja(1,2) * d(2)
  v(2)= aja(m+2,2) * d(2)
  do j= 3, n
    v(j)= aja(m+2,j) * d(j)
    v(1:j-1)= v(1:j-1) + d(j) * aja(1:j-1,j)
  enddo
  return
end subroutine tsutmd

end module M_TenSolve
