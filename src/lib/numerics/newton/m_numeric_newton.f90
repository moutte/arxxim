module M_Numeric_Newton
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_
  implicit none

  private
  !
  public:: NewtLnsrch
  public:: Newton_Press
  public:: Newton_Walker
  public:: Newton_Kelley
  public:: Newton
  public:: NewtonChess
  !
  logical,public:: TestJacob=.false.
  logical,public:: ShoResult=.false.
  !
  logical:: DeltaX_Ok= .true.
  !
contains

subroutine NewtLnsrch( & !from Press et al, Numerical Recipes
& vX,       & !inout= the initial guess, and the root returned
!!& vLPos,  & !in=    enforce vX(vLPos=TRUE)>0
& Residual, & !
& Jacobian, & !
& Converge, & !
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
!!& TolMin, & !in=    whether spurious convergence to a minimum of fmin has occurred
& bFinDif,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   maxval(abs(vFunc(:)))
& Delta_X,  & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Gradient, & !out=
& Nits,     & !out=   number of iterations
& Check,    & !out=   if Check, should check convergence
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
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  use M_Numeric_Const,only: Ln10 !for debugging
  use M_IoTools,      only: GetUnit, OutStrVec !for debugging
  !
  !!logical, dimension(:), intent(in)   :: vLPos
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
      real(dp),dimension(size(v))     :: Residual
    end function Residual
    subroutine Jacobian(v,t)
      use M_Kinds
      implicit none
      real(dp),dimension(:),              intent(in) :: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    logical function Converge(vF,vTolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: vF(:),vTolF(:)
    end function Converge
  endinterface
  !
  !integer, parameter :: MAXITS=NewtMaxIts
  !NewtTolF=   convergence criterion on function values
  !NewtTolMin= criterion for spurious convergence
  !NewtTolX=   convergence criterion on dx
  !
  real(dp),parameter:: STPMX= 100.0
  !scaled maximum step length allowed in line searches
  real(dp),parameter:: Alf= Epsilon(vX) !1.0D-12
  !
  logical :: bSingul
  integer :: ITS
  real(dp):: D,F,rFOld,rStpMax,rSlope,Norm_vDX !,FMin
  !
  real(dp),dimension(size(vX)):: vFunc,vGrad,vDX,vXOld
  integer, dimension(size(vX)):: vIndex
  real(dp),dimension(size(vX),size(vX)):: tJacob,tJacobNum
  !
  real(dp):: A,Alam,Alam2,AlaMin,B,DISC,F2
  real(dp):: RHS1,RHS2,tmpLam
  !
  vFunc= Residual(vX)
  F= 0.5_dp*dot_product(vFunc,vFunc)
  !
  iErr=-1
  !Error_F= maxval(abs(vFunc(:))) !-> error on function
  !if (Error_F < 0.01_dp*TolF) then !Test for initial guess being a root
  !  nIts=1
  !  iErr=1 !-----------Use more stringent test than simply NewtTolF ???
  !  return
  !end if
  !
  if(ShoResult) print '(/,A,/)',"iTs,Error_F,Gradient,Delta_X"
  !
  rStpMax=STPMX*MAX(vAbs(vX(:)),real(size(vX),dp)) !Calculate StpMax for line searches
  !
  do its=1,MaxIts
    !--------------------------------------------start of iteration loop
    !
    nIts=its
    !
    call Newton_Sho(vX,vFunc,Its)
    !
    if(bFinDif) then
      call Jacobian_Numeric(Residual,vX,vFunc,tJacob)
    else
      call Jacobian(vX,tJacob)
      if(TestJacob) then
      !for test, calc. numeric Jacobian and compare with analytical one
        call Jacobian_Numeric(Residual,vX,vFunc,tJacobNum)
        call TestJacob_Sho(tJacob,tJacobNum)
        TestJacob=.false. !calculate only on first call to NewtLnsrch
      end if
    end if
    !
    vGrad(:)= matmul(vFunc(:),tJacob(:,:)) ! compute grad(f) for the line search.
    vXOld(:)= vX(:)                        ! store vX
    rFOld=    F                            ! store F
    vDX(:)=  -vFunc(:)                     ! right-hand side for linear equations.
    !
    !--- solve linear equations by LU decomposition
    call LU_Decomp(tJacob,vIndex,D,bSingul)
    !
    if(bSingul) then
      iErr= -2 !-> "SINGULAR JACOBIAN, cf Log File"
      call Jacobian_Sho(tJacob)
      return !System not solved --------------------------------- return
    end if
    !
    call LU_BakSub(tJacob,vIndex,vDX)
    !---/ solve linear equations by LU decomposition
    !
    Norm_vDX= vAbs(vDX(:))
    !-- Scale if attempted step is too big --
    if (Norm_vDX>rStpMax) vDX(:)= rStpMax *vDX(:) /Norm_vDX
    !
    rSlope= dot_product(vGrad,vDX)
    !
    if (rSlope>=Zero) then
      if(idebug>1) then
        write(fTrc,'(/,A,/)') "ROUNdoFF PROBLEM -> vGrad, vDX ="
        call OutStrVec(fTrc,vGrad(1:size(vX)),OPt_C="G")
        call OutStrVec(fTrc,vDX(1:size(vX)),  OPt_C="G")
      end if
      iErr= -3 !-> 'roundoff problem in linesearch, cf log file'
      return ! System not solved -------------------------------- return
    end if
    !
    !----------------------------------------linesearch and backtracking
    !
    CHECK= .false.
    AlaMin= TolX /maxval(abs(vDX(:))/MAX(abs(vXOld(:)),One))
    Alam=   One !always try full Newton step first
    !
    DoLineSearch: do !Start of iteration loop
      !
      !!write(51,'(I3,G15.6)') Its, Alam
      !
      vX(:)= vXOld(:) +Alam*vDX(:)
      !
      vFunc= Residual(vX)
      F= dot_product(vFunc,vFunc)/2.0D0
      !
      if (Alam<AlaMin) then !Convergence on vX
        !
        vX(:)=vXOld(:)
        !-- For zero finding, the calling program should verify the convergence
        CHECK=.true.
        exit DoLineSearch !-----------------------------------------exit
        !
      elseif (F <= rFOld + ALF*Alam*rSlope) then
        !
        exit DoLineSearch !sufficient function decrease ----------- exit
        !
      else !--< Backtrack.
        !
        if (Alam==One) then !---------------------------------First time
          tmpLam= -rSlope /(2.0_dp*(F-rFOld-rSlope))
          !
        else !-------------------------------------Subsequent backtracks
          !print *,Alam2   ; pause_
          RHS1= F -rFOld -Alam *rSlope
          RHS2= F2-rFOld -Alam2*rSlope
          A= (       RHS1/Alam**2 -      RHS2/Alam2**2) /(Alam-Alam2)
          B= (-Alam2*RHS1/Alam**2 + Alam*RHS2/Alam2**2) /(Alam-Alam2)
          !
          if (A==Zero) then
            tmpLam= -rSlope/(2.0_dp*B)
            !
          else
            if(B<-1.0D100) B=-1.0D100
            if(B> 1.0D100) B= 1.0D100
            !print *,A     ; pause_
            !print *,B     ; pause_
            !print *,"B*B",B*B   ; pause_
            DISC= B*B -3.0_dp*A*rSlope
            !print *,DISC  ; pause_
            if     (DISC<Zero) then  ; tmpLam= 0.5_dp*Alam
            elseif (B<=Zero)   then  ; tmpLam= (-B+SQRT(DISC))/3.0_dp/A
            else                     ; tmpLam= -rSlope/(B+SQRT(DISC))
            end if
          end if
          !
          if (tmpLam>0.5_dp*Alam) tmpLam=0.5_dp*Alam
          !
        end if
        !
      end if
      !
      Alam2= Alam
      Alam=  max(tmpLam,0.1_dp*Alam)
      !
      F2=F
      !
    end do DoLineSearch
    !
    !---------------------------------------/linesearch and backtracking
    !
    Error_F=  maxval( abs(vFunc(:)) ) !-> error on function value
    !
    !-- Press-> max variation on vX
    !Delta_X= maxval( abs(vX(:)-vXOld(:))/MAX(abs(vX(:)),One) )
    !
    !-- Archim-> max relative variation on vX
    Delta_X=  maxval( abs( (vX(:)-vXOld(:)) /vX(:)) )
    !
    Gradient= maxval( abs(vGrad(:)) *MAX(abs(vX(:)),One) /MAX(F,0.5_dp*size(vX)) )
    !
    if(ShoResult) print '(I3,3(G15.6,A))',iTs,Error_F,"=F ",Gradient,"=G ",Delta_X,"=X"
    !
    !-------------------------------- Test for convergence on function values --
    !! if(Converge(vFunc)) then
    if(Error_F < TolF) then
      iErr= 0
      call Newton_Sho(vX,vFunc,Its) !mod add 19/06/2008 08:22
      return !----------------------------------------------------return
    end if
    !! if(Check) then !Check for gradient of f being zero
    !!   Check= (Gradient < TolMin ) !convergence on dx
    !!   iErr= -4
    !!   return !-------------------------------------------------return
    !! end if
    if(DeltaX_Ok .and. Delta_X < TolX) then
    !-> while Error_F still high, current vX very close to previous vX
    !-> arrived to a stationary point ??
      iErr= -5 !
      return !----------------------------------------------------return
    end if
    !
    !----------------------------------------------end of iteration loop
  end do
end subroutine NewtLnsrch

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine Newton_Press( & !from "NR"
& vX,       & !inout= the initial guess, and the root returned
& vIsPlus,  & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
& Residual, & !interface
& Jacobian, & !interface
& Converge, & !interface
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& bFinDif,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   maxval(abs(vFunc(:)))
& Delta_X,  & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
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
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Tools,only: fNewtJ
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  !
  real(dp),intent(inout):: vX(:)
  logical, intent(in)   :: vIsPlus(:)
  real(dp),intent(in)   :: TolF,TolX !!,TolMin
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X !,Gradient
  !! logical, intent(out)  :: Check
  integer, intent(out)  :: nIts,iErr
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
    logical function Converge(vF,vTolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: vF(:),vTolF(:)
    end function Converge
  endinterface
  !
  !integer, parameter :: MAXITS=NewtMaxIts
  ! NewtTolF=   convergence criterion on function values
  ! NewtTolMin= criterion for spurious convergence
  ! NewtTolX=   convergence criterion on dx
  !
  ! scaled maximum step length allowed in line searches
  real(dp),parameter:: STPMX= 100.0
  !
  real(dp),parameter:: Alf= Epsilon(vX) !1.0D-12
  !
  logical :: Check
  logical :: bSingul
  logical :: ForcePositive
  integer :: ITS
  real(dp):: D,F,rFOld,rStpMax,rSlope,Norm_vDX,Gradient !,FMin
  !
  real(dp),dimension(size(vX)):: vFunc,vDX,vX0,vGrad
  integer, dimension(size(vX)):: vIndex
  real(dp),dimension(size(vX),size(vX)):: tJacob
  real(dp),dimension(size(vX)):: vTolF
  !
  real(dp):: A,Alam,Alam2,AlaMin,B,DISC,F2,X
  real(dp):: RHS1,RHS2,tmpLam
  !
  ForcePositive= count(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  vFunc= Residual(vX)
  F= dot_product(vFunc,vFunc) *0.5_dp
  !
  iErr=-1
  !
  if(ShoResult) print '(/,A,/)',"iTs,Error_F,Gradient,Delta_X"
  !
  ! Calculate StpMax for line searches
  rStpMax= STPMX*MAX(vAbs(vX(:)),real(size(vX),dp))
  !
  do its=1,MaxIts
    !
    nIts=its
    !
    call Newton_Sho(vX,vFunc,Its)
    !
    if(bFinDif) then
      call Jacobian_Numeric(Residual,vX,vFunc,tJacob)
    else
      call Jacobian(vX,tJacob)
    end if
    !
    vGrad(:)= matmul(vFunc(:),tJacob(:,:)) ! compute grad(f) for the line search.
    vX0(:)= vX(:)                          ! store vX
    rFOld=    F                            ! store F
    vDX(:)=  -vFunc(:)                     ! right-hand side for linear equations.
    !
    !-------------------------solve linear equations by LU decomposition
    call LU_Decomp(tJacob,vIndex,D,bSingul)
    !
    if(bSingul) then
      iErr= -2 !-> "SINGULAR JACOBIAN, cf Log File"
      !call Jacobian_Sho(fNewtJ,tJacob)
      return !-------------------------------System not solved -> return
    end if
    !
    call LU_BakSub(tJacob,vIndex,vDX)
    !------------------------/solve linear equations by LU decomposition
    !
    !! if(fNewtJ>0) then
      !! call Jacobian_Sho(fNewtJ,tJacob)
      !! fNewtJ= 0  ;  close(fNewtJ)
    !! end if
    !
    Norm_vDX= vAbs(vDX(:))
    !-- Scale if attempted step is too big --
    if (Norm_vDX>rStpMax) vDX(:)= rStpMax *vDX(:) /Norm_vDX
    !
    rSlope= dot_product(vGrad,vDX)
    !
    if (rSlope>=Zero) then
      iErr= -3 !-> 'roundoff problem in linesearch, cf log file'
      return !-------------------------------System not solved -> return
    end if
    !
    !-------------------------------------------------enforce positivity
    if(ForcePositive) then
      !--- under relaxation technique from Bethke
      !--- if (Xi/2 + dXi)-0
      X= maxval(-vDX(:)/vX0(:)*2.0D0,mask=vIsPlus(:))
      if(iDebug>2 .and. X>1.0D0) write(69,*) " Press-UR, ", X
      if(X>1.0D0) vDX(:)= vDX(:) /X
      !do
      !  X= MINVAL(vX(:),mask=vIsPlus(:))
      !  if(X>Zero) exit
      !  vDX(:)= vDX(:)*0.5D0
      !  vX(:)= vX0(:) + vDX(:)
      !end do
    end if
    !------------------------------------------------/enforce positivity
    !
    vX(:)= vX0(:) +vDX(:)
    !
    !----------------------------------------linesearch and backtracking
    !
    CHECK= .false.
    AlaMin= TolX /maxval(abs(vDX(:))/MAX(abs(vX0(:)),One))
    Alam=   One !always try full Newton step first
    !
    DoLineSearch: do !Start of iteration loop
      !
      !!write(51,'(I3,G15.6)') Its, Alam
      !
      vX(:)= vX0(:) +Alam*vDX(:)
      !
      vFunc= Residual(vX)
      F= dot_product(vFunc,vFunc)/2.0D0
      !
      if (Alam<AlaMin) then !Convergence on vX
        !
        vX(:)=vX0(:)
        !-- For zero finding, the calling program should verify the convergence
        CHECK=.true.
        exit DoLineSearch !-----------------------------------------exit
        !
      elseif (F <= rFOld + ALF*Alam*rSlope) then
        !
        exit DoLineSearch !sufficient function decrease ------------exit
        !
      else !--Backtrack.
        !
        if (Alam==One) then !---------------------------------First time
          tmpLam= -rSlope /(2.0_dp*(F-rFOld-rSlope))
          !
        else !-------------------------------------Subsequent backtracks
          RHS1= F -rFOld -Alam *rSlope
          RHS2= F2-rFOld -Alam2*rSlope
          A= (       RHS1/Alam**2 -      RHS2/Alam2**2) /(Alam-Alam2)
          B= (-Alam2*RHS1/Alam**2 + Alam*RHS2/Alam2**2) /(Alam-Alam2)
          !
          if (A==Zero) then
            tmpLam= -rSlope/(2.0_dp*B)
            !
          else
            B= min(B, 1.0D100)
            B= max(B,-1.0D100)
            DISC= B*B -3.0_dp*A*rSlope
            if     (DISC<Zero) then  ; tmpLam= 0.5_dp*Alam
            elseif (B<=Zero)   then  ; tmpLam= (-B+SQRT(DISC))/3.0_dp/A
            else                     ; tmpLam= -rSlope/(B+SQRT(DISC))
            end if
          end if
          !
          tmpLam= min(tmpLam,0.5_dp*Alam)
          !
        end if
        !
      end if
      !
      Alam2= Alam
      Alam=  max(tmpLam,0.1_dp*Alam)
      !
      F2=F
      !
    end do DoLineSearch
    !
    !---------------------------------------/linesearch and backtracking
    !
    Error_F=  maxval( abs(vFunc(:)) ) !-> error on function value
    !
    !-- Press-> max relative variation on vX
    Delta_X= maxval( abs(vX(:)-vX0(:))/MAX(abs(vX(:)),One) )
    !
    !-- Archim-> max relative variation on vX
    ! Delta_X=  maxval( abs( (vX(:)-vX0(:)) /vX(:)) )
    !
    Gradient= maxval( abs(vGrad(:)) *MAX(abs(vX(:)),One) /MAX(F,0.5_dp*size(vX)) )
    !
    if(ShoResult) print '(I3,3(G15.6,A))',iTs,Error_F,"=F ",Gradient,"=G ",Delta_X,"=X"
    !
    !----------------------------test for convergence on function values
    if(Converge(vFunc,vTolF)) then
    !if(Error_F < TolF) then
      iErr= 0
      call Newton_Sho(vX,vFunc,Its) !mod add 19/06/2008 08:22
      return !----------------------------------------------------return
    end if
    !---/
    !! if(Check) then !Check for gradient of f being zero
    !!   Check= (Gradient < TolMin ) !convergence on dx
    !!   iErr= -4
    !!   return !-------------------------------------------------return
    !! end if
    !------------------------------------------test for stationary point
    if(DeltaX_Ok .and. Delta_X < TolX) then
    !-> while Error_F still high, current vX very close to previous vX
    !-> arrived to a stationary point ??
      iErr= -5 !
      return !----------------------------------------------------return
    end if
    !---/
    !
    !----------------------------------------------end of iteration loop
  end do
end subroutine Newton_Press

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine Newton_Walker( & !
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
& bFinDif,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   maxval(abs(fVec(:)))
& Delta_X,  & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!--
!-- translated from matlab code by Homer Walekr
!-- http://users.wpi.edu/~walker/MA590
!-- cf http://users.wpi.edu/~walker/MA590/HANdoUTS/newton-backtracking.pdf
!--
!-- Given an initial guess x for a root in n dimensions,
!-- take Nits Newton-Raphson steps to improve the root
!-- Stop if the root converges,
!-- in either summed absolute variable increments NewtTolX
!-- or summed absolute function values NewtTolF
!--
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Tools,only: fNewtJ
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  !
  !use M_Numeric_Const,  only: Ln10
  !use M_IOTools,only: OutStrVec
  !
  real(dp),intent(inout):: vX(:)
  logical, intent(in)   :: vIsPlus(:)
  real(dp),intent(in)   :: TolF
  real(dp),intent(in)   :: TolX
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X
  integer, intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    logical function Converge(vF,vTolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: vF(:),vTolF(:)
    end function Converge
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vFunc0,vDX,vX0
  real(dp),dimension(size(vX))         :: vTolF
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  !
  integer :: Its
  real(dp):: D
  logical :: bSingul
  logical :: ForcePositive
  !
  real(dp),parameter:: &
  & SigmaMax= 0.5_dp,  &
  & SigmaMin= 0.1_dp,  &
  & Tau=      1.0D-4
  integer, parameter:: iArmMax= 12
  integer :: iArm
  real(dp):: Norm_vF0,Norm_vF,delta,lambda,Sigma
  real(dp):: X
  !
  ForcePositive= count(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  iErr=-1
  !
  do Its=1,MaxIts
    !
    vX0(:)= vX(:)
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(sum(vFunc0(:)*vFunc0(:)))
    !
    call Newton_Sho(vX0,vFunc0,Its)
    !
    if(bFinDif) then
      call Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
    else
      call Jacobian(vX0,tJac)
    end if
    !
    !! if(fNewtJ>0) then
      !! call Jacobian_Sho(fNewtJ,tJac)
      !! fNewtJ= 0  ;  close(fNewtJ)
    !! end if
    !
    !----------------------solve linear equations using LU decomposition
    call LU_Decomp(tJac, vIndex, D, bSingul)
    if(bSingul) then
      iErr=-2  ;  exit !--------------------------------------------exit
    end if
    vDX(:)= -vFunc0(:)
    call LU_BakSub(tJac,vIndex,vDX)
    !---------------------/solve linear equations using LU decomposition
    !
    if(ForcePositive) then
      X= maxval(-vDX(:)/vX0(:)*2.0D0,mask=vIsPlus(:))
      if(iDebug>2 .and. X>1.0D0) write(69,*) "Walker-UR, ", X
      if(X>1.0D0) vDX(:)= vDX(:) /X
      !do
      !  X= MINVAL(vX(:),mask=vIsPlus(:))
      !  if(X>Zero) exit
      !  vDX(:)= vDX(:)*0.5D0
      !  vX(:)= vX0(:) + vDX(:)
      !end do
    end if
    !
    vX(:)= vX0(:) + vDX(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(sum(vFunc(:)*vFunc(:)))
    !
    lambda= 1.0d0
    iArm=   0
    !
    !---------------------------Test the step and backtrack as necessary
    do while ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      iArm= iArm +1
      !
      if(iArm > iArmMax) then !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  exit
      end if
      !
      delta= (Norm_vF /Norm_vF0)**2 -1.0d0 +2.0d0*lambda
      !
      if (delta > Zero) then
        Sigma=  lambda/delta
        Sigma=  min(Sigma,SigmaMax)
        Sigma=  max(Sigma,SigmaMin)
      else
        Sigma=  SigmaMax
      end if
      !
      vDX(:)=  Sigma*vDX(:)
      lambda=  Sigma*lambda
      vX(:)=   vX0(:) + vDX(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(sum(vFunc(:)*vFunc(:)))
      !
    end do
    !--------------------------/Test the step and backtrack as necessary
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vDX(:)))
    !
    !-- convergence --
    !if (Error_F < TolF) then
    if(Converge(vFunc,vTolF)) then  ;  iErr= 0  ;  exit !-----------exit
    end if
    !
    !-- stationary --
    if (Delta_X < TolX) then        ;  iErr=-5  ;  exit !-----------exit
    end if
    !
  end do
  !
  nIts=Its
  !
end subroutine Newton_Walker

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine Newton_Kelley( & !
& vX,        & !inout= the initial guess, and the root returned
& vIsPlus,   & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
& Residual,  & !
& Jacobian,  & !
& Converge,  & !
& TolF,      & !in=    convergence criterion on function values
& TolX,      & !in=    convergence criterion on dx
& bFinDif,   & !in=    use numeric Jacobian
& MaxIts,    & !in=    maximum number of iterations
& Error_F,   & !out=   maxval(abs(fVec(:)))
& Delta_X,   & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Nits,      & !out=   number of iterations
& iErr)        !out=   error code
!--
!-- translated from MATLAB code nsold.m,
!-- http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsold.m
!-- cf CT Kelley, Solving nonlinear equations with Newton's method, SIAM, 2003
!--
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Tools,only: fNewtJ
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  !
  real(dp),intent(inout):: vX(:)
  logical, intent(in)   :: vIsPlus(:)
  real(dp),intent(in)   :: TolF,TolX
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X
  integer, intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    logical function Converge(F,tolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: F(:),tolF(:)
    end function Converge
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vFunc0,vDX,vX0,vStep,vTolF
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  !
  integer :: Its
  real(dp):: D,X
  logical :: bSingul
  logical :: ForcePositive
  !
  !-- Armijo parameters --
  real(dp),parameter:: SigmaMax= 0.5_dp
  real(dp),parameter:: SigmaMin= 0.1_dp
  real(dp),parameter:: Tau=      1.0D-4
  integer, parameter:: iArmMax= 10
  !--
  !
  integer :: iArm
  real(dp):: Norm_vF0,Norm_vF,Norm_Ratio
  real(dp):: lambda,lam_m,lam_c
  real(dp):: ff_0, ff_m, ff_c
  !
  ForcePositive= count(vIsPlus(:))>0
  !
  vTolF(:)= TolF
  !
  iErr=-1
  Norm_Ratio= One
  !
  do Its=1,MaxIts
    !
    vX0(:)= vX(:)
    !
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(dot_product(vFunc0(:),vFunc0(:)))
    !
    call Newton_Sho(vX0,vFunc0,Its)
    !
    !----------------------solve linear equations using LU decomposition
    !!if(Its==1 .or. Norm_Ratio > 0.5_dp) then
      !-- evaluate and factorize Jacobian
      !-- only when norm(Residual) strongly decreases
      !
      if(bFinDif) then
        call Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
      else
        call Jacobian(vX0,tJac)
      end if
      !
      call LU_Decomp(tJac, vIndex, D, bSingul)
      if(bSingul) then
        iErr=-2  ;  exit
      end if
      !
    !!end if
    !
    vDX(:)= -vFunc0(:)
    call LU_BakSub(tJac,vIndex,vDX)
    !-------------------------------------------------------------/solve
    !
    !! if(fNewtJ>0) then
      !! call Jacobian_Sho(fNewtJ,tJac)
      !! fNewtJ= 0  ;  close(fNewtJ)
    !! end if
    !
    if(ForcePositive) then
      X= maxval(-vDX(:)/vX0(:)*2.0D0,mask=vIsPlus(:))
      if(iDebug>2 .and. X>1.0D0) write(69,*) "Kelley-UR, ", X
      if(X>1.0D0) vDX(:)= vDX(:) /X
      !do
      !  X= MINVAL(vX(:),mask=vIsPlus(:))
      !  if(X>Zero) exit
      !  vStep(:)= vStep(:)*0.5D0
      !  vX(:)= vX0(:) + vStep(:)
      !end do
    end if
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
    Norm_vF= SQRT(dot_product(vFunc(:),vFunc(:)))
    !
    ff_0= Norm_vF0**2  ! initial value
    ff_c= Norm_vF**2   ! current value
    ff_m= ff_c         ! previous value in linesearch
    !
    !---------------------------test the step and backtrack as necessary
    do while ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      if(iArm==0) then
        lambda=  SigmaMax*lambda
      else
        call Parab3p(SigmaMin,SigmaMax, lam_c,lam_m, ff_0,ff_c,ff_m, lambda)
      end if
      !
      iArm= iArm + 1
      if(iArm > iArmmax) then !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  exit
      end if
      !
      lam_m= lam_c    ! lam_c=  current steplength
      lam_c= lambda   ! lam_m=  previous steplength
      !
      vStep(:)= lambda *vDX(:)
      vX(:)= vX0(:) + vStep(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(dot_product(vFunc(:),vFunc(:)))
      !
      ff_m= ff_c
      ff_c= Norm_vF**2
      !
    end do
    !--------------------------/test the step and backtrack as necessary
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vStep(:)))
    !
    Norm_Ratio= Norm_vF /Norm_vF0
    !! write(71,'(A,G16.6)') "Kelley, Norm_Ratio ", Norm_Ratio
    !
    !if (Error_F < TolF) then
    if(Converge(vFunc,vTolF)) then  ;  iErr= 0   ;  exit
    end if
    !
    if (Delta_X<TolX) then          ;  iErr=-5   ;  exit
    end if
    !
  end do
  !
  !! write(71,*)
  nIts=Its
  !
end subroutine Newton_Kelley

subroutine parab3p(sigma0,sigma1, lambdac,lambdam, ff0,ffc,ffm, lambdap)
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
  real(dp),intent(in) :: sigma0,sigma1
  real(dp),intent(in) :: lambdac, lambdam, ff0, ffc, ffm
  real(dp),intent(out):: lambdap
  !
  real(dp):: c1, c2
  !
  !-------------------- Compute coefficients of interpolation polynomial
  ! p(lambda)=  ff0 + (-c1 lambda + c2 lambda^2)/d1
  ! d1=  (lambdac - lambdam)*lambdac*lambdam < 0
  ! so, if c2 > 0 we have negative curvature and default to
  ! lambdap=  sigma1 * lambda.
  !
  c2= lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
  !
  if(c2 >=0.0d0) then
    lambdap=  sigma1*lambdac
  else
    c1= lambdam *lambdam *(ffc-ff0) - lambdac *lambdac *(ffm-ff0)
    lambdap= c1 *0.5d0 /c2
    lambdap= max(lambdap, sigma0*lambdac)
    lambdap= min(lambdap, sigma1*lambdac)
  end if
  !
  return
end subroutine parab3p

subroutine Newton_Walker_old( & !
& vX,       & !inout=  the initial guess, and the root returned
& Residual, & !
& Jacobian, & !
& Converge, & !
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& bFinDif,  & !in=    use numeric Jacobian
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   maxval(abs(fVec(:)))
& Delta_X,  & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!--
!-- translated from matlab code by Homer Walekr
!-- http://users.wpi.edu/~walker/MA590
!-- cf http://users.wpi.edu/~walker/MA590/HANdoUTS/newton-backtracking.pdf
!--
!-- Given an initial guess x for a root in n dimensions,
!-- take Nits Newton-Raphson steps to improve the root
!-- Stop if the root converges,
!-- in either summed absolute variable increments NewtTolX
!-- or summed absolute function values NewtTolF
!--
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  !
  use M_Numeric_Const,  only: Ln10
  use M_IOTools,only: OutStrVec
  !
  real(dp),intent(inout):: vX(:)
  real(dp),intent(in)   :: TolF,TolX
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X
  integer, intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    logical function Converge(vF,vTolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: vF(:),vTolF(:)
    end function Converge
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vFunc0,vDX,vX0,vTolF
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  !
  integer :: Its
  real(dp):: D
  logical :: bSingul
  !
  real(dp),parameter:: &
  & SigmaMax= 0.5_dp,  &
  & SigmaMin= 0.1_dp,  &
  & Tau=      1.0D-4
  integer, parameter:: iArmMax= 12
  integer :: iArm
  real(dp):: Norm_vF0,Norm_vF,delta,lambda,Sigma
  !
  vTolF(:)= TolF
  !
  iErr=-1
  !
  do Its=1,MaxIts
    !
    vX0(:)= vX(:)
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(sum(vFunc0(:)*vFunc0(:)))
    !
    call Newton_Sho(vX0,vFunc0,Its)
    !
    if(bFinDif) then
      call Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
    else
      call Jacobian(vX0,tJac)
    end if
    !
    !---solve linear equations using LU decomposition --
    call LU_Decomp(tJac, vIndex, D, bSingul)
    if(bSingul) then
      iErr=-2  ;  exit
    end if
    vDX(:)= -vFunc0(:)
    call LU_BakSub(tJac,vIndex,vDX)
    !---/solve
    !
    vX(:)= vX0(:) + vDX(:)
    !
    vFunc(:)= Residual(vX)
    Norm_vF= SQRT(sum(vFunc(:)*vFunc(:)))
    !
    lambda= 1.0d0
    iArm=   0
    !
    !---------------------------Test the step and backtrack as necessary
    do while ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      iArm= iArm +1
      !
      if(iArm > iArmMax) then !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  exit
      end if
      !
      delta= (Norm_vF /Norm_vF0)**2 -1.0d0 +2.0d0*lambda
      !
      if (delta > Zero) then
        Sigma=  lambda/delta
        Sigma=  min(Sigma,SigmaMax)
        Sigma=  max(Sigma,SigmaMin)
      else
        Sigma=  SigmaMax
      end if
      !
      vDX(:)=  Sigma*vDX(:)
      lambda=  Sigma*lambda
      vX(:)=   vX0(:) + vDX(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(sum(vFunc(:)*vFunc(:)))
      !
    end do
    !--------------------------/Test the step and backtrack as necessary
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vDX(:)))
    !
    !if (Error_F < TolF) then
    if(Converge(vFunc,vTolF)) then
      iErr= 0  ;  exit
    end if
    if (Delta_X < TolX) then !-- stationary ==
      iErr=-5  ;  exit
    end if
    !
  end do
  !
  nIts=Its
  !
end subroutine Newton_Walker_old

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine Newton_Kelley_old( & !
& vX,        & !inout= the initial guess, and the root returned
& Residual,  & !
& Jacobian,  & !
& Converge,  & !
& TolF,      & !in=    convergence criterion on function values
& TolX,      & !in=    convergence criterion on dx
& bFinDif,   & !in=    use numeric Jacobian
& MaxIts,    & !in=    maximum number of iterations
& Error_F,   & !out=   maxval(abs(fVec(:)))
& Delta_X,   & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Nits,      & !out=   number of iterations
& iErr)        !out=   error code
!--
!-- translated from MATLAB code nsold.m,
!-- http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsold.m
!-- cf CT Kelley, Solving nonlinear equations with Newton's method, SIAM, 2003
!--
  use M_Numeric_Tools,only: vAbs,Jacobian_Numeric
  use M_Numeric_Mat,  only: LU_BakSub,LU_Decomp
  !
  use M_Numeric_Const,only: Ln10
  use M_IOTools,      only: OutStrVec
  !
  real(dp),intent(inout):: vX(:)
  real(dp),intent(in)   :: TolF,TolX
  logical, intent(in)   :: bFinDif
  integer, intent(in)   :: MaxIts
  real(dp),intent(out)  :: Error_F,Delta_X
  integer, intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    logical function Converge(F,tolF)
      use M_Kinds
      implicit none
      real(dp),intent(in):: F(:),tolF(:)
    end function Converge
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vFunc0,vDX,vX0,vStep,vTolF
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  !
  integer :: Its
  real(dp):: D
  logical :: bSingul
  !
  !--- Armijo parameters --
  real(dp),parameter:: SigmaMax= 0.5_dp
  real(dp),parameter:: SigmaMin= 0.1_dp
  real(dp),parameter:: Tau=      1.0D-4
  integer, parameter:: iArmMax= 10
  !---/
  !
  integer :: iArm
  real(dp):: Norm_vF0,Norm_vF,Norm_Ratio
  real(dp):: lambda,lam_m,lam_c
  real(dp):: ff_0, ff_m, ff_c
  !
  vTolF(:)= TolF
  !
  iErr=-1
  Norm_Ratio= One
  !
  do Its=1,MaxIts
    !
    vX0(:)= vX(:)
    !
    vFunc0(:)= Residual(vX0)
    Norm_vF0= SQRT(dot_product(vFunc0(:),vFunc0(:)))
    !
    call Newton_Sho(vX0,vFunc0,Its)
    !
    !if(DebNewt .and. bOpenNewt) then
    !  call OutStrVec(fNewt1,X(1:size(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  call OutStrVec(fNewt2,vFunc0(1:size(X)),  Opt_I=Its,Opt_C="F")
    !end if
    !
    !-------------------------- solve linear equations using LU decomposition --
    !!if(Its==1 .or. Norm_Ratio > 0.5_dp) then
      !-- evaluate and factorize Jacobian
      !-- only when norm(Residual) strongly decreases
      !
      if(bFinDif) then
        call Jacobian_Numeric(Residual,vX0,vFunc0,tJac)
      else
        call Jacobian(vX0,tJac)
      end if
      !
      call LU_Decomp(tJac, vIndex, D, bSingul)
      if(bSingul) then
        iErr=-2  ;  exit
      end if
      !
    !!end if
    !
    vDX(:)= -vFunc0(:)
    call LU_BakSub(tJac,vIndex,vDX)
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
    Norm_vF= SQRT(dot_product(vFunc(:),vFunc(:)))
    !
    ff_0= Norm_vF0**2
    ff_c= Norm_vF**2
    ff_m= Norm_vF**2
    !
    !---------------------------test the step and backtrack as necessary
    do while ( Norm_vF > (1.0d0-tau*lambda)*Norm_vF0 )
      !
      !!write(52,'(I3,G15.6)') Its, lambda
      !
      if(iArm==0) then
        lambda=  SigmaMax*lambda
      else
        call Parab3p(SigmaMin,SigmaMax, lam_c,lam_m, ff_0,ff_c,ff_m, lambda)
      end if
      !
      iArm=  iArm + 1
      if(iArm > iArmmax) then !! error('Maximum number of backtracks reached.')
        iErr= -4  ;  exit
      end if
      !
      ! lam_c=  current steplength
      ! lam_m=  previous steplength
      lam_m= lam_c
      lam_c= lambda
      !
      vStep(:)= lambda *vDX(:)
      vX(:)= vX0(:) + vStep(:)
      vFunc(:)= Residual(vX)
      Norm_vF= SQRT(dot_product(vFunc(:),vFunc(:)))
      !
      ff_m= ff_c
      ff_c= Norm_vF**2
      !
    end do
    !------------------------------/ test the step and backtrack as necessary --
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vStep(:)))
    !
    Norm_Ratio= Norm_vF /Norm_vF0
    !! write(71,'(G16.6)') Norm_Ratio
    !
    !if (Error_F < TolF) then
    if(Converge(vFunc,vTolF)) then
      iErr= 0; exit
    end if
    if (Delta_X<TolX) then
      iErr=-5; exit
    end if
    !
  end do
  !
  nIts=Its
  !
end subroutine Newton_Kelley_old

subroutine parab3p_old(sigma0,sigma1, lambdac,lambdam, ff0,ffc,ffm, lambdap)
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
  real(dp),intent(in) :: sigma0,sigma1
  real(dp),intent(in) :: lambdac, lambdam, ff0, ffc, ffm
  real(dp),intent(out):: lambdap
  !
  real(dp):: c1, c2
  !
  !---------------------Compute coefficients of interpolation polynomial
  ! p(lambda)=  ff0 + (-c1 lambda + c2 lambda^2)/d1
  ! d1=  (lambdac - lambdam)*lambdac*lambdam < 0
  ! so, if c2 > 0 we have negative curvature and default to
  ! lambdap=  sigma1 * lambda.
  !
  c2= lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
  !
  if(c2 >=0.0d0) then
    lambdap=  sigma1*lambdac
  else
    c1= lambdam *lambdam *(ffc-ff0) - lambdac *lambdac *(ffm-ff0)
    lambdap= c1 *0.5d0 /c2
    lambdap= max(lambdap, sigma0*lambdac)
    lambdap= min(lambdap, sigma1*lambdac)
  end if
  !
  return
end subroutine parab3p_old

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine NewtonChess( &
& vX, &      !initial guess x for a root in N dimensions
& Residual, Jacobian, &
& TolF,    & !in=    convergence criterion on function values
& TolX,    & !in=    convergence criterion on dx
!!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
& MaxIts,  & !in=    maximum number of iterations
& Error_F, & !out=   maxval(abs(fVec(:)))
& Delta_X, & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
!& Gradient, & !out=
& Nits,    & !out=   number of iterations
!& Check,   & !out=   if Check, should check convergence
& iErr)      !out=   error code
!Given an initial guess x for a root in N dimensions,
!take Nits Newton-Raphson steps to improve the root
!Stop if the root converges,
!in either summed absolute variable increments NewtTolX
!or summed absolute function values NewtTolF.
  use M_Numeric_Tools,   only: vAbs, iMaxLoc_R
  use M_Numeric_Mat,only: LU_BakSub,LU_Decomp
  !
  use M_Numeric_Const,  only: Ln10
  use M_IOTools,only: OutStrVec
  !
  real(dp),dimension(:), intent(inout):: vX
  real(dp),              intent(in)   :: TolF,TolX !!,TolMin
  integer,               intent(in)   :: MaxIts
  real(dp),              intent(out)  :: Error_F,Delta_X !,Gradient
  !logical,               intent(out)  :: Check
  integer,               intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vDX
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  integer :: Its
  real(dp):: D
  logical ::bSingul
  !
  !Chess/Concepts, p71, Newton with "polishing factor"
  integer :: I
  real(dp):: R,Alfa
  real(dp),parameter::  A=0.5D0, B=3.0D0, C=-0.9D0
  !iErr=-1: MaxIts reached
  !iErr=-2: singular Jacobian
  !iErr=-5:
  !
  iErr=-1
  !
  do Its=1,MaxIts
    !
    vFunc(:)= Residual(vX)
    call Jacobian(vX,tJac)
    !
    !if(DebNewt .and. bOpenNewt) then
    !  call OutStrVec(fNewt1,X(1:size(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  call OutStrVec(fNewt2,vFunc(1:size(X)),  Opt_I=Its,Opt_C="F")
    !end if
    !
    !--- Solve linear equations using LU decomposition --
    call LU_Decomp(tJac,vIndex,D,bSingul)
    if(bSingul) then
      iErr=-2; exit
    end if
    vDX(:)= -vFunc(:)
    call LU_BakSub(tJac,vIndex,vDX) !solve vDX= -vFunc.inv(tJac)
    !---/ Solve
    !
    R= One
    !Chess/Concepts, p71, Newton with "polishing factor"
    Alfa= maxval(abs(vDX(:)/vX(:)))  !alfa_i= |dX_i|/X_i
    if(Alfa > A) then !-> far from root
      I= iMaxLoc_R(abs(vDX(:)/vX(:)))
      if(vDX(I)>Zero) R=   (Alfa*B - A*A)/(B + Alfa -2.0D0*A)*vX(I)/vDX(I)
      if(vDX(I)<Zero) R= C*(Alfa   - A*A)/(B + Alfa -2.0D0*A)*vX(I)/vDX(I)
    end if
    !
    vX= vX + R *vDX
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vDX(:)))
    !
    if (Error_F<TolF) then
      iErr= 0; exit
    end if
    !
    if (Delta_X<TolX) then
      iErr=-5; exit
    end if
  end do
  !
  nIts=Its
end subroutine NewtonChess

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine Newton( &
& vX, &      ! inout= the initial guess, and the root returned
& Residual, Jacobian, &
& TolF,     & !in=    convergence criterion on function values
& TolX,     & !in=    convergence criterion on dx
& MaxIts,   & !in=    maximum number of iterations
& Error_F,  & !out=   maxval(abs(fVec(:)))
& Delta_X,  & !out=   maxval( abs(vX(:)-vXOld(:)) / max(abs(vX(:)),One) )
& Nits,     & !out=   number of iterations
& iErr)       !out=   error code
!-----------------------------------------------------------------------
! Given an initial guess x for a root in n dimensions,
! take Nits Newton-Raphson steps to improve the root
! Stop if the root converges,
! in either summed absolute variable increments NewtTolX
! or summed absolute function values NewtTolF
!-----------------------------------------------------------------------
  use M_Numeric_Tools,   only: vAbs
  use M_Numeric_Mat,only: LU_BakSub,LU_Decomp
  !
  use M_Numeric_Const,  only: Ln10
  use M_IOTools,only: OutStrVec
  !
  real(dp),dimension(:), intent(inout):: vX
  real(dp),              intent(in)   :: TolF,TolX !!,TolMin
  integer,               intent(in)   :: MaxIts
  real(dp),              intent(out)  :: Error_F,Delta_X !,Gradient
  !logical,               intent(out)  :: Check
  integer,               intent(out)  :: nIts,iErr
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
      real(dp),dimension(:), intent(in):: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
  endinterface
  !
  real(dp),dimension(size(vX))         :: vFunc,vDX
  real(dp),dimension(size(vX),size(vX)):: tJac
  integer, dimension(size(vX))         :: vIndex
  integer :: Its
  real(dp):: D
  logical :: bSingul
  real(dp),parameter:: Alfa= 0.5_dp
  !
  iErr=-1
  do Its=1,MaxIts
    !
    vFunc(:)= Residual(vX)
    call Jacobian(vX,tJac)
    !
    !if(DebNewt .and. bOpenNewt) then
    !  call OutStrVec(fNewt1,X(1:size(X))/Ln10,Opt_I=Its,Opt_C="F")
    !  call OutStrVec(fNewt2,vFunc(1:size(X)),  Opt_I=Its,Opt_C="F")
    !end if
    !
    !--- solve linear equations using LU decomposition --
    call LU_Decomp(tJac,vIndex,D,bSingul)
    if(bSingul) then
      iErr=-2  ;  exit
    end if
    vDX(:)= -vFunc(:)
    call LU_BakSub(tJac,vIndex,vDX)
    !---/ solve
    !
    vX= vX + vDX *Alfa
    !
    Error_F= maxval(abs(vFunc(:)))
    Delta_X= maxval(abs(vDX(:)))
    !
    if (Error_F < TolF) then
      iErr= 0; exit
    end if
    if (Delta_X<TolX) then
      iErr=-5; exit
    end if
  end do
  !
  nIts=Its
  !
  return
end subroutine Newton

!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!!///\\\!
!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!!\\\///!

subroutine TestJacob_Sho(tJacob,tJacobNum)
  use M_IoTools,only: GetUnit, OutStrVec
  !
  real(dp),dimension(:,:),intent(in):: tJacob,tJacobNum
  !
  real(dp),dimension(size(tJacob,1))  :: xx,xr
  integer:: ff,I
  !
  write(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN"
  write(fTrc,'(A)') "results in zzDynam_TestJacob.log"
  call GetUnit(ff)
  open(ff,file="debug_dynam_jacob.log")
  do I=1,size(tJacob,1)
    call OutStrVec(ff,tJacob(I,:),   Opt_I=I,Opt_S="JacobAna",Opt_C="G")
    call OutStrVec(ff,tJacobNum(I,:),Opt_I=I,Opt_S="JacobDif",Opt_C="G")
    !
    xx(:)= tJacob(I,:) -tJacobNum(I,:)
    call OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltJac",   Opt_C="G")
    !
    xr(:)= tJacob(I,:) +tJacobNum(I,:)
    where(xr(:)>1.D-16) ; xx(:)= xx(:) / xr(:)
    elsewhere           ; xx(:)= Zero
    endwhere
    call OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltRel",   Opt_C="G")
  end do
  write(ff,*)
  close(ff)
  write(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN end"
end subroutine TestJacob_Sho

subroutine Newton_Sho(vX,vFunc,Its) !,nAq)
  use M_Numeric_Const,        only: Ln10
  use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
  use M_IoTools,      only: OutStrVec
  !
  real(dp),intent(in):: vX(:),vFunc(:)
  integer, intent(in):: Its
  !
  fNewt_I= fNewt_I +1
  if(fNewtF>0) call OutStrVec(fNewtF,vX(:)/Ln10,   Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  if(fNewtR>0) call OutStrVec(fNewtR,abs(vFunc(:)),Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
  !
end subroutine Newton_Sho

subroutine Jacobian_Sho(t)
  real(dp),intent(in):: t(:,:)
  integer:: I,J
  do I=1,size(t,1)
    do J=1,size(t,2)
      if(t(I,J)/=0) then  ;  write(fTrc,'(F7.1,A1)',advance="no") t(I,J), T_
      else                ;  write(fTrc,'(A1,  A1)',advance="no") "0",    T_
      end if
    end do
    write(fTrc,*)
  end do
  return
end subroutine Jacobian_Sho

end module M_Numeric_Newton
