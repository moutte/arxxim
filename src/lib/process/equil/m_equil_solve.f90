MODULE M_Equil_Solve
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Solve
  !
CONTAINS

SUBROUTINE Equil_Solve(NewtIts,Newt_iErr)
!--
!-- basic procedure for speciation, including loop on gammas
!--
  !---------------------------------------------------------------------
  USE M_Numeric_Const,ONLY: Ln10,TinyDP
  USE M_Safe_Functions
  USE M_IoTools,      ONLY: OutStrVec
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_SolModel_Calc, ONLY: Solmodel_CalcGamma
  USE M_Basis
  !---------------------------------------------------------------------
  USE M_Global_Vars, ONLY: vSpc,SolModel
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_Basis_Vars,  ONLY: isW,nAx,vOrdAq,vLAx,MWSv,nCi,nAs
  !
  USE M_Equil_Vars,  ONLY: vMolF,vFasMole
  USE M_Equil_Vars,  ONLY: vLnGam,vLnAct,vDGapp_As
  USE M_Equil_Vars,  ONLY: LogForAqu,DirectSub,vMolSec,cMethod,cEquMode
  USE M_Equil_Vars,  ONLY: vTooLow,nEvalFunc
  USE M_Equil_Vars,  ONLY: fActiz,GamMaxIts,TolGam,TolAct
  USE M_Equil_Vars,  ONLY: vEquFas,nEquFas
  !---------------------------------------------------------------------
  
  INTEGER,INTENT(OUT):: NewtIts,Newt_iErr
  !
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vX,vLnGamOld,vLnActOld !,vXsec
  LOGICAL, DIMENSION(:),ALLOCATABLE:: vXisPlus
  !REAL(dp),DIMENSION(:),ALLOCATABLE:: vTolF,vTolX,
  INTEGER :: nAq,iGam,iCount,I,J,Newt_nIts
  REAL(dp):: r1,r2
  REAL(dp):: OsmoSv
  LOGICAL :: OkConverged

  INTEGER:: nF !nr of aqu'species in array vX of residual
  INTEGER:: nM !nr of non-aqu'species in array vX of residual

  nAq= COUNT(vSpc(:)%Typ(1:3)=="AQU")
  !~ OsmoSv= One

  !ALLOCATE(vTolF(1:nF+nM))   ;  vTolF= NewtTolF !; vTolF(nCi+1:nAs)= 1.0E-06_dp
  !ALLOCATE(vTolX(1:nF+nM))   ;  vTolX= NewtTolX

  ALLOCATE(vLnGamOld(1:nAq)) ;  vLnGamOld(1:nAq)= vLnGam(1:nAq)
  ALLOCATE(vLnActOld(1:nAq)) ;  vLnActOld(1:nAq)= vLnAct(1:nAq)

  nM= nEquFas        ! new version, using vEquFas, etc

  IF(DirectSub) THEN  ;  nF= nCi
  ELSE                ;  nF= nCi +nAs
  ENDIF
  !
  ALLOCATE(vMolSec(nAs))
  ALLOCATE(vX(1:nF+nM))
  ALLOCATE(vXisPlus(1:nF+nM))
  !
  !~ IF(DirectSub) ALLOCATE(vXsec(1:nAs))
  !
  vX(1:nF)= vMolF(vOrdAq(1:nF))
  vXisPlus(:)= .FALSE.
  IF(.NOT. LogForAqu) vXisPlus(1:nF)= .TRUE.
  !
  IF(LogForAqu) vX(1:nF)= FSafe_vLog(vX(1:nF))
  !
  IF(nM>0) THEN

    DO I=1,nEquFas
      vX(nF+I)= vEquFas(I)%Mole
      !vX(nF+I)= vFasMole(I)
    ENDDO

  ENDIF
  !
  iCount=  0
  iGam=    0
  NewtIts= 0
  !
  OkConverged= .FALSE.
  !
  IF(iDebug>2) &
  & WRITE(72,'(A)') "=====================================Equil_Solve=="
  !-------------------------------------------------- loop on Gamma's --
  DoGamma: DO
    !
    iCount=iCount+1
    !
    CALL Update_DGapp_As(vDGapp_As,vLnGam)
    !
    !--------------------------------------------------------- solver --
    CALL Equil_Newton(cMethod,vX,vXisPlus,Newt_nIts,Newt_iErr)
    !--------------------------------------------------------/ solver --
    !
    NewtIts= NewtIts +Newt_Nits
    !
    IF(Newt_iErr<0) EXIT DoGamma !========Newton failed, exit DoGamma ==
    !
    IF(LogForAqu) THEN
      !--- Back to Mole Quantities
      IF(DirectSub) THEN
        !~ CALL Compute_SecondSpecies(vX,vXsec)
        vMolF(vOrdAq(1:nCi))= EXP(vX(1:nF))
        !~ vMolF(vOrdAq(nCi+1:nCi+nAs))= EXP(vXsec(1:nAs))
        vMolF(vOrdAq(nCi+1:nCi+nAs))= vMolSec(1:nAs)
      ELSE
        vMolF(vOrdAq(1:nF))= EXP(vX(1:nF))
      ENDIF
    ELSE
      IF(DirectSub) THEN
        !~ CALL Compute_SecondSpecies(vX,vXsec)
        vMolF(vOrdAq(1:nCi))= vX(1:nF)
        !~ vMolF(vOrdAq(nCi+1:nCi+nAs))= vXsec(1:nAs)
        vMolF(vOrdAq(nCi+1:nCi+nAs))= vMolSec(1:nAs)
      ELSE
        vMolF(vOrdAq(1:nF))= vX(1:nF)
      ENDIF
    ENDIF
    !
    !--- mobile aqu'species located at end of vOrdAq !!
    DO J=1,nAx
      vMolF(vOrdAq(nCi+nAs+J))= MWSv*vMolF(isW) & !
      *EXP(vLnAct(vOrdAq(nCi+nAs+J)) -vLnGam(vOrdAq(nCi+nAs+J)))
    ENDDO
    !
    IF(OkConverged) EXIT DoGamma
    !
    vLnGamOld= vLnGam
    vLnActOld= vLnAct
    !
    !------------------------------------------- update activ'coeff's --
    CALL Solmodel_CalcGamma( & !
    & TdgK,Pbar,     & !
    & SolModel,      & !IN
    & isW,vSpc,      & !IN
    & vLAx,          & !IN
    & vMolF(1:nAq),  & !IN
    & vLnAct(1:nAq), & !OUT
    & vLnGam(1:nAq), & !OUT
    & vTooLow(1:nAq),OsmoSv) !OUT
    !
    IF(fActiZ>0) CALL OutStrVec(fActiz,vLnAct(1:nAq)/Ln10,Opt_I=iGam)
    !------------------------------------------/ update activ'coeff's --
    !
    !IF(ALL(ABS(vLnGamOld-vLnGam)<TolGam)) EXIT !EXIT when convergence on gammas
    ! r1= MAXVAL(ABS((vLnGamOld-vLnGam)/vLnGam))
    ! r2= MAXVAL(ABS((vLnActOld-vLnAct)/vLnAct))
    r1= MAXVAL(ABS(vLnGamOld-vLnGam))
    r2= MAXVAL(ABS(vLnActOld-vLnAct))
    !
    !IF(iDebug>0) WRITE(fTrc,'(2(A,I3),2(A,G15.6))') &
    !& "iGam=",iGam,"/ Newt_nIts=",Newt_nIts,"/ deltaGamma=",r1,"/ OsmoSv=",OsmoSv
    !IF(iDebug>3) PRINT '(2(A,I3),2(A,G15.6))', &
    !& "iGam=",iGam,"/ Newt_nIts=",Newt_nIts,"/ deltaGamma=",r1,"/ OsmoSv=",OsmoSv
    !
    IF(iDebug>2) WRITE(72,'(A,3I4,2(A,G15.6))') &
    & "iGam /Newt_nIts /nEvalFunc=",iGam,Newt_nIts,nEvalFunc, &
    & "/ deltaGamma=",r1,"/ deltaLnAct=",r2
    IF(iDebug==4) PRINT '(A,3I4,2(A,G15.6))', &
    & "iGam /Newt_nIts /nEvalFunc=",iGam,Newt_nIts,nEvalFunc, &
    & "/ deltaGamma=",r1,"/ deltaLnAct=",r2
    !
    !------------------------------------------------- if Convergence --
    ! test convergence on activ'coeffs and on activities
    ! -> update mole nrs according to current activ'coeffs
    !OkConverged= (r1< TolGam *MAX(One,MAXVAL(ABS(vLnGam))) &
    !&       .AND. r2< TolAct *MAX(One,MAXVAL(ABS(vLnAct))) )
    OkConverged= (r1<TolGam .AND. r2<TolAct)
    !------------------------------------------------/ if Convergence --
    !
    iGam=iGam+1
    IF(iGam>=GamMaxIts) EXIT DoGamma !=< TOO MANY LOOPS, EXIT DoGamma ==
    !
  ENDDO DoGamma
  !-------------------------------------------------/ loop on Gamma's --
  !
  IF(iGam>=GamMaxIts) Newt_iErr=-4
  !
  !--- retrieve mole numbers of mineral or gas --
  IF(nEquFas>0) THEN

    DO I=1,nEquFas
      vEquFas(I)%Mole= vX(nF+I)
      vFasMole(I)= vX(nF+I)
    ENDDO

  ENDIF
  !---/
  !
  DEALLOCATE(vMolSec)
  DEALLOCATE(vX)
  DEALLOCATE(vXisPlus)
  !~ IF(DirectSub) DEALLOCATE(vXsec)
  DEALLOCATE(vLnGamOld)
  DEALLOCATE(vLnActOld)
  !
ENDSUBROUTINE Equil_Solve

SUBROUTINE Update_DGapp_As(vDGapp_As,vLnGam)
!--
!-- update the "apparent" deltaG's of second'aqu'species
!-- (include the gamma's on solute species)
!--
  USE M_Global_Vars,ONLY: vSpc
  USE M_Basis_Vars, ONLY: vOrdAq,vOrdPr,vOrdAs,nCi,nAs,tNuAs
  REAL(dp),INTENT(IN) :: vLnGam(:)
  REAL(dp),INTENT(OUT):: vDGapp_As(:)
  !
  INTEGER:: I
  DO I=1,nAs
    vDGapp_As(I)= &
    &   vSpc(vOrdAs(I))%G0rt  - DOT_PRODUCT(tNuAs(I,:),    vSpc(vOrdPr(:))%G0rt ) &
    & + vLnGam(vOrdAq(nCi+I)) - DOT_PRODUCT(tNuAs(I,2:nCi),vLnGam(vOrdAq(2:nCi)))
  ENDDO
END SUBROUTINE Update_DGapp_As

SUBROUTINE Compute_SecondSpecies(vX,vXsec)
  USE M_Basis_Vars,ONLY: nAs,nAx,nMx,tNuAs,nCi,isW,MWSv,vOrdPr
  USE M_Equil_Vars,ONLY: vDGapp_As,vLnGam,vLnAct,LogForAqu
  !
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vXsec(:)
  !
  INTEGER :: iAs,iCi,i
  REAL(dp):: X
  !
  IF(LogForAqu) THEN

    DO iAs=1,nAs
      !
      X=   vLnAct(isW) *tNuAs(iAs,isW) &
      &  + (One - SUM(tNuAs(iAs,2:nCi)))*(vX(isW) +LOG(MWSv)) &
      &  - vDGapp_As(iAs)
      DO iCi=1,nCi
        IF(iCi/=isW) X= X + vX(iCi) *tNuAs(iAs,iCi)
      ENDDO
      DO I=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+I)) *tNuAs(iAs,nCi+I)
      ENDDO
      !
      vXsec(iAs)= X
      !
    END DO
    !
    if(iDebug>2) then
      write(71,'(A)') "====================SOLVE=="
      do i=1,size(vX)
        write(71,'(A,G15.6)') "resPrim=",EXP(vX(i))
      enddo
      do i=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",EXP(vXsec(i))
      enddo
    endif

  ELSE

    DO iAs=1,nAs
      !
      X= vLnAct(isW)*tNuAs(iAs,isW) -vDGapp_As(iAs)
      ! note: vDGapp_As includes the delta on activity coeff's
      DO i=1,nAx+nMx
        X= X + vLnAct(vOrdPr(nCi+I))*tNuAs(iAs,nCi+I)
      ENDDO
      X= EXP(X)
      !
      DO iCi=1,nCi
        IF(iCi/=isW .AND. tNuAs(iAs,iCi) /= Zero) &
        !& X= X *(vMol(vOrdAqu(iCi))/vMol(isW)/MWSv)**tNuAs(iAs,iCi)
        & X= X *(vX(iCi)/vX(isW)/MWSv)**tNuAs(iAs,iCi)
      END DO
      !vMol(nCi+iAs)= X *vMol(isW) *MWSv
      vXsec(iAs)= X *vX(isW) *MWSv
      !
    END DO
    !
    if(iDebug>2) then
      write(71,'(A)') "====================SOLVE=="
      do i=1,size(vX)
        write(71,'(A,G15.6)') "resPrim=",vX(i)
      enddo
      do i=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",vXsec(i)
      enddo
    endif

  ENDIF
  !pause
  !
END SUBROUTINE Compute_SecondSpecies

SUBROUTINE Equil_Newton(cMethod,vX,vXisPlus,Newt_nIts,Newt_iErr)
  USE M_Numeric_Newton
  USE M_Equil_Residual
  USE M_Equil_Jacobian
  !
  USE M_Equil_Vars,ONLY: NewtTolF,NewtTolX,NewtMaxIts
  USE M_Equil_Vars,ONLY: LogForAqu,DirectSub,bFinDIF,nEvalFunc
  !
  !INTEGER, INTENT(IN)   :: iMethod
  CHARACTER(LEN=*),INTENT(IN)   :: cMethod
  REAL(dp),        INTENT(INOUT):: vX(:)
  LOGICAL,         INTENT(IN)   :: vXisPlus(:)
  INTEGER,         INTENT(OUT)  :: Newt_nIts
  INTEGER,         INTENT(OUT)  :: Newt_iErr
  !
  REAL(dp):: Error_F,Gradient,Delta_X
  LOGICAL :: Check
  LOGICAL :: JacobNumeric
  ! LOGICAL :: ForcePositive
  !
  Error_F= Zero
  Delta_X= Zero
  Gradient= Zero
  Check= .FALSE.
  !
  ALLOCATE(vTolCoef(SIZE(vX)))
  vTolCoef(:)= One
  !
  JacobNumeric= bFinDIF
  !
  IF(iDebug>2) nEvalFunc= 0
  !
  !SELECT CASE(iMethod)
  SELECT CASE(TRIM(cMethod))
  
  CASE("NEWTONKELLEY") !(0)
    CALL Newton_Kelley( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,       & !in=    convergence criterion on FUNCTION values
    & NewtTolX,       & !in=    convergence criterion on dx
    & JacobNumeric,   & !in=    use numeric Jacobian
    & NewtMaxIts,     & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr)    !out=   error code
  !
  CASE("NEWTONWALKER") !(1)
    CALL Newton_Walker( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,      & !in=    convergence criterion on FUNCTION values
    & NewtTolX,      & !in=    convergence criterion on dx
    & JacobNumeric,  & !in= use numeric Jacobian
    & NewtMaxIts,    & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr   & !out=   error code
    & )
  !
  CASE("NEWTONPRESS") !(2)
    CALL Newton_Press( & !
    & vX,             & !inout= the initial guess, and the root returned
    & vXisPlus,       & !

    & Equil_Residual, & !
    & Equil_Jacobian, & !
    & Equil_Converge, & !

    & NewtTolF,      & !in=    convergence criterion on FUNCTION values
    & NewtTolX,      & !in=    convergence criterion on dx
    & JacobNumeric,  & !in= use numeric Jacobian
    & NewtMaxIts,    & !in=    maximum number of iterations

    & Error_F,    & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Newt_Nits,  & !out=   number of iterations
    & Newt_iErr   & !out=   error code
    & )
  !
  CASE("NEWTLNSRCH") !(3)
    CALL NewtLnsrch( & !from "NR"
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
  END SELECT
  !
  DEALLOCATE(vTolCoef)
  !
  RETURN
END SUBROUTINE Equil_Newton

ENDMODULE M_Equil_Solve

