module M_Equil_Solve
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Stop_,Pause_
  implicit none
  !
  private
  !
  public:: Equil_Solve
  !
contains

subroutine Equil_Solve(NewtIts,Newt_iErr)
!--
!-- basic procedure for speciation, including loop on gammas
!--
  !---------------------------------------------------------------------
  use M_Numeric_Const,only: Ln10,TinyDP
  use M_Safe_Functions
  use M_IoTools,      only: OutStrVec
  use M_Dtb_Const,    only: T_CK
  use M_SolModel_Calc, only: Solmodel_CalcGamma
  use M_Basis
  !---------------------------------------------------------------------
  use M_Global_Vars, only: vSpc,SolModel
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_Basis_Vars,  only: isW,nAx,vOrdAq,vLAx,MWSv,nCi,nAs
  !
  use M_Equil_Vars,  only: vMolF,vFasMole
  use M_Equil_Vars,  only: vLnGam,vLnAct,vDGapp_As
  use M_Equil_Vars,  only: LogForAqu,DirectSub,vMolSec,cMethod,cEquMode
  use M_Equil_Vars,  only: vTooLow,nEvalFunc
  use M_Equil_Vars,  only: fActiz,GamMaxIts,TolGam,TolAct
  use M_Equil_Vars,  only: vEquFas,nEquFas
  !---------------------------------------------------------------------
  
  integer,intent(out):: NewtIts,Newt_iErr
  !
  real(dp),dimension(:),allocatable:: vX,vLnGamOld,vLnActOld !,vXsec
  logical, dimension(:),allocatable:: vXisPlus
  !real(dp),dimension(:),allocatable:: vTolF,vTolX,
  integer :: nAq,iGam,iCount,I,J,Newt_nIts
  real(dp):: r1,r2
  real(dp):: OsmoSv
  logical :: OkConverged

  integer:: nF !nr of aqu'species in array vX of residual
  integer:: nM !nr of non-aqu'species in array vX of residual

  nAq= count(vSpc(:)%Typ(1:3)=="AQU")
  !~ OsmoSv= One

  !allocate(vTolF(1:nF+nM))   ;  vTolF= NewtTolF !; vTolF(nCi+1:nAs)= 1.0E-06_dp
  !allocate(vTolX(1:nF+nM))   ;  vTolX= NewtTolX

  allocate(vLnGamOld(1:nAq)) ;  vLnGamOld(1:nAq)= vLnGam(1:nAq)
  allocate(vLnActOld(1:nAq)) ;  vLnActOld(1:nAq)= vLnAct(1:nAq)

  nM= nEquFas        ! new version, using vEquFas, etc

  if(DirectSub) then  ;  nF= nCi
  else                ;  nF= nCi +nAs
  end if
  !
  allocate(vMolSec(nAs))
  allocate(vX(1:nF+nM))
  allocate(vXisPlus(1:nF+nM))
  !
  !~ if(DirectSub) allocate(vXsec(1:nAs))
  !
  vX(1:nF)= vMolF(vOrdAq(1:nF))
  vXisPlus(:)= .false.
  if(.not. LogForAqu) vXisPlus(1:nF)= .true.
  !
  if(LogForAqu) vX(1:nF)= FSafe_vLog(vX(1:nF))
  !
  if(nM>0) then

    do I=1,nEquFas
      vX(nF+I)= vEquFas(I)%MolFs
      !vX(nF+I)= vFasMole(I)
    end do

  end if
  !
  iCount=  0
  iGam=    0
  NewtIts= 0
  !
  OkConverged= .false.
  !
  if(iDebug>2) &
  & write(72,'(A)') "=====================================Equil_Solve=="
  !-------------------------------------------------- loop on Gamma's --
  DoGamma: do
    !
    iCount=iCount+1
    !
    call Update_DGapp_As(vDGapp_As,vLnGam)
    !
    !--------------------------------------------------------- solver --
    call Equil_Newton(cMethod,vX,vXisPlus,Newt_nIts,Newt_iErr)
    !--------------------------------------------------------/ solver --
    !
    NewtIts= NewtIts +Newt_Nits
    !
    if(Newt_iErr<0) exit DoGamma !--======Newton failed, exit DoGamma ==
    !
    if(LogForAqu) then
      !--- Back to Mole Quantities
      if(DirectSub) then
        !~ call Compute_SecondSpecies(vX,vXsec)
        vMolF(vOrdAq(1:nCi))= exp(vX(1:nF))
        !~ vMolF(vOrdAq(nCi+1:nCi+nAs))= exp(vXsec(1:nAs))
        vMolF(vOrdAq(nCi+1:nCi+nAs))= vMolSec(1:nAs)
      else
        vMolF(vOrdAq(1:nF))= exp(vX(1:nF))
      end if
    else
      if(DirectSub) then
        !~ call Compute_SecondSpecies(vX,vXsec)
        vMolF(vOrdAq(1:nCi))= vX(1:nF)
        !~ vMolF(vOrdAq(nCi+1:nCi+nAs))= vXsec(1:nAs)
        vMolF(vOrdAq(nCi+1:nCi+nAs))= vMolSec(1:nAs)
      else
        vMolF(vOrdAq(1:nF))= vX(1:nF)
      end if
    end if
    !
    !--- mobile aqu'species located at end of vOrdAq !!
    do J=1,nAx
      vMolF(vOrdAq(nCi+nAs+J))= MWSv*vMolF(isW) & !
      *exp(vLnAct(vOrdAq(nCi+nAs+J)) -vLnGam(vOrdAq(nCi+nAs+J)))
    end do
    !
    if(OkConverged) exit DoGamma
    !
    vLnGamOld= vLnGam
    vLnActOld= vLnAct
    !
    !------------------------------------------- update activ'coeff's --
    call Solmodel_CalcGamma( & !
    & TdgK,Pbar,     & !
    & SolModel,      & !IN
    & isW,vSpc,      & !IN
    & vLAx,          & !IN
    & vMolF(1:nAq),  & !IN
    & vLnAct(1:nAq), & !OUT
    & vLnGam(1:nAq), & !OUT
    & vTooLow(1:nAq),OsmoSv) !OUT
    !
    if(fActiZ>0) call OutStrVec(fActiz,vLnAct(1:nAq)/Ln10,Opt_I=iGam)
    !------------------------------------------/ update activ'coeff's --
    !
    !if(all(ABS(vLnGamOld-vLnGam)<TolGam)) exit !exit when convergence on gammas
    ! r1= MAXVAL(ABS((vLnGamOld-vLnGam)/vLnGam))
    ! r2= MAXVAL(ABS((vLnActOld-vLnAct)/vLnAct))
    r1= MAXVAL(ABS(vLnGamOld-vLnGam))
    r2= MAXVAL(ABS(vLnActOld-vLnAct))
    !
    !if(iDebug>0) write(fTrc,'(2(A,I3),2(A,G15.6))') &
    !& "iGam=",iGam,"/ Newt_nIts=",Newt_nIts,"/ deltaGamma=",r1,"/ OsmoSv=",OsmoSv
    !if(iDebug>3) print '(2(A,I3),2(A,G15.6))', &
    !& "iGam=",iGam,"/ Newt_nIts=",Newt_nIts,"/ deltaGamma=",r1,"/ OsmoSv=",OsmoSv
    !
    if(iDebug>2) write(72,'(A,3I4,2(A,G15.6))') &
    & "iGam /Newt_nIts /nEvalFunc=",iGam,Newt_nIts,nEvalFunc, &
    & "/ deltaGamma=",r1,"/ deltaLnAct=",r2
    if(iDebug==4) print '(A,3I4,2(A,G15.6))', &
    & "iGam /Newt_nIts /nEvalFunc=",iGam,Newt_nIts,nEvalFunc, &
    & "/ deltaGamma=",r1,"/ deltaLnAct=",r2
    !
    !------------------------------------------------- if Convergence --
    ! test convergence on activ'coeffs and on activities
    ! -> update mole nrs according to current activ'coeffs
    !OkConverged= (r1< TolGam *MAX(One,MAXVAL(ABS(vLnGam))) &
    !&       .and. r2< TolAct *MAX(One,MAXVAL(ABS(vLnAct))) )
    OkConverged= (r1<TolGam .and. r2<TolAct)
    !------------------------------------------------/ if Convergence --
    !
    iGam=iGam+1
    if(iGam>=GamMaxIts) exit DoGamma !=< TOO MANY LOOPS, exit DoGamma ==
    !
  end do DoGamma
  !-------------------------------------------------/ loop on Gamma's --
  !
  if(iGam>=GamMaxIts) Newt_iErr=-4
  !
  !--- retrieve mole numbers of mineral or gas --
  if(nEquFas>0) then

    do I=1,nEquFas
      vEquFas(I)%MolFs= vX(nF+I)
      vFasMole(I)= vX(nF+I)
    end do

  end if
  !---/
  !
  deallocate(vMolSec)
  deallocate(vX)
  deallocate(vXisPlus)
  !~ if(DirectSub) deallocate(vXsec)
  deallocate(vLnGamOld)
  deallocate(vLnActOld)
  !
end subroutine Equil_Solve

subroutine Update_DGapp_As(vDGapp_As,vLnGam)
!--
!-- update the "apparent" deltaG's of second'aqu'species
!-- (include the gamma's on solute species)
!--
  use M_Global_Vars,only: vSpc
  use M_Basis_Vars, only: vOrdAq,vOrdPr,vOrdAs,nCi,nAs,tNuAs
  real(dp),intent(in) :: vLnGam(:)
  real(dp),intent(out):: vDGapp_As(:)
  !
  integer:: I
  do I=1,nAs
    vDGapp_As(I)= &
    &   vSpc(vOrdAs(I))%G0rt  - dot_product(tNuAs(I,:),    vSpc(vOrdPr(:))%G0rt ) &
    & + vLnGam(vOrdAq(nCi+I)) - dot_product(tNuAs(I,2:nCi),vLnGam(vOrdAq(2:nCi)))
  end do
end subroutine Update_DGapp_As

subroutine Compute_SecondSpecies(vX,vXsec)
  use M_Basis_Vars,only: nAs,nAx,nMx,tNuAs,nCi,isW,MWSv,vOrdPr
  use M_Equil_Vars,only: vDGapp_As,vLnGam,vLnAct,LogForAqu
  !
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vXsec(:)
  !
  integer :: iAs,iCi,i
  real(dp):: X
  !
  if(LogForAqu) then

    do iAs=1,nAs
      !
      X=   vLnAct(isW) *tNuAs(iAs,isW) &
      &  + (One - SUM(tNuAs(iAs,2:nCi)))*(vX(isW) +log(MWSv)) &
      &  - vDGapp_As(iAs)
      do iCi=1,nCi
        if(iCi/=isW) X= X + vX(iCi) *tNuAs(iAs,iCi)
      end do
      do I=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+I)) *tNuAs(iAs,nCi+I)
      end do
      !
      vXsec(iAs)= X
      !
    end do
    !
    if(iDebug>2) then
      write(71,'(A)') "====================SOLVE=="
      do i=1,size(vX)
        write(71,'(A,G15.6)') "resPrim=",exp(vX(i))
      end do
      do i=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",exp(vXsec(i))
      end do
    end if

  else

    do iAs=1,nAs
      !
      X= vLnAct(isW)*tNuAs(iAs,isW) -vDGapp_As(iAs)
      ! note: vDGapp_As includes the delta on activity coeff's
      do i=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+I))*tNuAs(iAs,nCi+I)
      end do
      X= exp(X)
      !
      do iCi=1,nCi
        if(iCi/=isW .and. tNuAs(iAs,iCi) /= Zero) &
        !& X= X *(vMol(vOrdAqu(iCi))/vMol(isW)/MWSv)**tNuAs(iAs,iCi)
        & X= X *(vX(iCi)/vX(isW)/MWSv)**tNuAs(iAs,iCi)
      end do
      !vMol(nCi+iAs)= X *vMol(isW) *MWSv
      vXsec(iAs)= X *vX(isW) *MWSv
      !
    end do
    !
    if(iDebug>2) then
      write(71,'(A)') "====================SOLVE=="
      do i=1,size(vX)
        write(71,'(A,G15.6)') "resPrim=",vX(i)
      end do
      do i=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",vXsec(i)
      end do
    end if

  end if
  !pause
  !
end subroutine Compute_SecondSpecies

subroutine Equil_Newton(cMethod,vX,vXisPlus,Newt_nIts,Newt_iErr)
  use M_Numeric_Newton
  use M_Equil_Residual
  use M_Equil_Jacobian
  !
  use M_Equil_Vars,only: NewtTolF,NewtTolX,NewtMaxIts
  use M_Equil_Vars,only: LogForAqu,DirectSub,bFinDif,nEvalFunc
  !
  !integer, intent(in)   :: iMethod
  character(len=*),intent(in)   :: cMethod
  real(dp),        intent(inout):: vX(:)
  logical,         intent(in)   :: vXisPlus(:)
  integer,         intent(out)  :: Newt_nIts
  integer,         intent(out)  :: Newt_iErr
  !
  real(dp):: Error_F,Gradient,Delta_X
  logical :: Check
  logical :: JacobNumeric
  ! logical :: ForcePositive
  !
  Error_F= Zero
  Delta_X= Zero
  Gradient= Zero
  Check= .false.
  !
  allocate(vTolCoef(size(vX)))
  vTolCoef(:)= One
  !
  JacobNumeric= bFinDif
  !
  if(iDebug>2) nEvalFunc= 0
  !
  !select case(iMethod)
  select case(trim(cMethod))
  
  case("NEWTONKELLEY") !(0)
    call Newton_Kelley( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,       & !in=    convergence criterion on function values
    & NewtTolX,       & !in=    convergence criterion on dx
    & JacobNumeric,   & !in=    use numeric Jacobian
    & NewtMaxIts,     & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr)    !out=   error code
  !
  case("NEWTONWALKER") !(1)
    call Newton_Walker( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,      & !in=    convergence criterion on function values
    & NewtTolX,      & !in=    convergence criterion on dx
    & JacobNumeric,  & !in= use numeric Jacobian
    & NewtMaxIts,    & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr   & !out=   error code
    & )
  !
  case("NEWTONPRESS") !(2)
    call Newton_Press( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,      & !in=    convergence criterion on function values
    & NewtTolX,      & !in=    convergence criterion on dx
    & JacobNumeric,  & !in= use numeric Jacobian
    & NewtMaxIts,    & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr   & !out=   error code
    & )
  !
  case("NEWTLNSRCH") !(3)
    call NewtLnsrch( & !from "NR"
    & vX,          & !inout= the initial guess, and the root returned

    & Equil_Residual, Equil_Jacobian, Equil_Converge, &
    & NewtTolF,    & !in=    convergence criterion on function values
    & NewtTolX,    & !in=    convergence criterion on dx

    & JacobNumeric, & !in=   use numeric Jacobian
    & NewtMaxIts,  & !in=    maximum number of iterations
    & Error_F,     & !out=   MAXVAL(ABS(fVec(:)))

    & Delta_X,     & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient,    & !out=
    & Newt_Nits,   & !out=   number of iterations
    & Check,       & !out=   if Check, should check convergence
    & Newt_iErr)     !out=   error code
  !
  end select
  !
  deallocate(vTolCoef)
  !
  return
end subroutine Equil_Newton

end module M_Equil_Solve

