module M_Equil_Jacobian
  use M_Trace,        only: iDebug,fTrc,T_,Stop_,Pause_
  use M_Numeric_Const,only: MinExpDP,MaxExpDP
  use M_Kinds
  implicit none
  !
  private
  !
  public:: Equil_Jacobian
  !
contains

subroutine Equil_Jacobian(vX,tJac)
!--
!-- EquJacobian of the system EquResidual; used by Newton
!-- NB: in vX and in tJac, solvent has index 1 !!!!
!--
  use M_Safe_Functions
  !
  use M_Global_Vars,only: SolModel
  use M_Basis_Vars, only: tAlfPr,tAlfAs,tNuAs,tAlfFs,tNuFas
  use M_Basis_Vars, only: vOrdPr,vOrdAq,nAx,nAs,nCi
  !
  use M_Equil_Vars,only: vMolal,vLnAct,vLnGam,vMolSec
  !! use M_Equil_Vars,only: vFasYes
  use M_Equil_Vars,only: cEquMode
  use M_Equil_Vars,only: vAffScale,vFasAff,dXi
  !! use M_Equil_Vars,only: OsmoSv,UseOsmoSv
  use M_Equil_Vars,only: DebNewt,DebJacob
  use M_Equil_Vars,only: LogForAqu,DirectSub
  use M_Equil_Vars,only: tAlfEq,tNuEq,nEquFas
  !---------------------------------------------------------------------
  real(dp),dimension(:),                intent(in) :: vX
  real(dp),dimension(size(vX),size(vX)),intent(out):: tJac
  !---------------------------------------------------------------------
  real(dp),dimension(1:nCi):: dRatdX
  integer :: iCi,jCi,iAs,iAx,I,isW
  real(dp):: MWsv
  real(dp):: X
  real(dp):: AffScale
  real(dp),allocatable:: vLnXf(:),vXf(:)
  real(dp),allocatable:: vSumNuAs(:)
  !
  integer:: nF !nr of aqu'species in array vX of residual
  integer:: nM !nr of non-aqu'species in array vX of residual
  !---------------------------------------------------------------------
  tJac=  Zero
  !
  !if(DirectSub) call Stop_("No Analytic Jacobian ...")
  !
  if(DirectSub) then  ;   nF= nCi
  else                ;   nF= nCi +nAs
  end if
  !
  isW= SolModel%iSolvent
  MWSv = SolModel%MolWeitSv
  nM= nEquFas
  !
  allocate(vXf(nF))  ;  vXf= Zero
  !
  allocate(vLnXf(nF))
  if(LogForAqu) then
    vLnXf(1:nF)= vX(1:nF)
    vXf(:)= FSafe_vExp(vLnXf(:)) ! avoid overflow with EXP
  else
    vXf(1:nF)= vX(1:nF)
    vLnXf(1:nF)= FSafe_vLog(vX(1:nF))
  end if
  !
  allocate(vSumNuAs(nAs))
  do iAs=1,nAs
    X= Zero
    do I=1,nCi
      if(I/=isW) X= X+tNuAs(iAs,I)
    end do
    vSumNuAs(iAs)= X
  end do
  !
  !do I=1,nAx
  !  vMolal(vOrdPr(nCi+I))= exp(vLnAct(vOrdPr(nCi+I))-vLnGam(vOrdPr(nCi+I)))
  !end do
  !if(nAx>0) vMolal(vOrdPr(nCi+1:nCi+nAx))= &
  !& exp(vLnAct(vOrdPr(nCi+1:nCi+nAx))-vLnGam(vOrdPr(nCi+1:nCi+nAx))) !-> i.e. molality
  !
  if(LogForAqu) then
    !------------------------------------------ material conservation --
    !--------------------------- rows 1:nCi - "internal" prim'species --
    do iCi=1,nCi
      !--- solvent
      tJac(iCi,1)= tAlfPr(iCi,1)*vXf(1)
      !-------- prim'mobile'aqu'species -> add corresponding moles of Solvent --
      !! if(nAx>0) &
      !! & tJac(iCi,1)= tJac(iCi,1) &
      !! & + vXf(isW) *MWSv &
      !! &   *dot_product(tAlfPr(iCi,nCi+1:nCi+nAx),vMolal(vOrdPr(nCi+1:nCi+nAx)))
      do iAx=1,nAx
        tJac(iCi,1)= tJac(iCi,1) &
        &          + vXf(isW) *MWSv *tAlfPr(iCi,nCi+iAx) *vMolal(vOrdPr(nCi+iAx))
      end do
      !--- prim'intern'species
      do jCi=2,nCi
        tJac(iCi,jCi)= tAlfPr(iCi,jCi)*vXf(jCi)
      end do
      !
      !--- second'species
      if(DirectSub) then
        do iAs=1,nAs
          tJac(iCi,1)= tJac(iCi,1) &
          &          + vMolSec(iAs) *tAlfAs(iCi,iAs) *(One - vSumNuAs(iAs))
        end do
        do jCi=2,nCi
          do iAs=1,nAs
            tJac(iCi,jCi)= tJac(iCi,jCi) &
            &            + vMolSec(iAs) *tAlfAs(iCi,iAs) *tNuAs(iAs,jCi)
          end do
        end do
      else
        do iAs=1,nAs
          tJac(iCi,nCi+iAs)= vXf(nCi+iAs) *tAlfAs(iCi,iAs)
        end do
      end if
      !
    end do
    !---/"internal" prim'species --
    !
    !! if(bMolal) then
      !! tJac(isW,:)= Zero
      !! tJac(isW,1)= vXf(isW)
    !! end if
    !
    if(DirectSub) then
      !
    else
      !--- rows nCi+1:nAq- sec'species --
      do iAs=1,nAs
        tJac(nCi+iAs,1)= One-sum(tNuAs(iAs,2:nCi)) !for Solvent
        tJac(nCi+iAs,2:nCi)=     tNuAs(iAs,2:nCi)  !for other aqu.species
      end do
      do iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= -One
      end do
      !---/
    end if
    !
    !! if(UseOsmoSv) then
      !! !
      !! SumXi= sum(vXf(1:nAs))-vXf(isW)
      !! !
      !! do iAs=1,nAs !add terms for osmotic coeff'
        !! tJac(nCi+iAs,1  )= tJac(nCi+iAs,1) &
        !! &                + tNuAs(iAs,isW) *OsmoSv *SumXi /vXf(isW) != d/dXw(LnActW)
        !! do I=2,nAs !for other aqu.species
          !! tJac(nCi+iAs,I)= tJac(nCi+iAs,I) &
          !! &              - tNuAs(iAs,isW) *OsmoSv *vXf(I)/vXf(isW) != d/dXi(LnActW)
        !! end do
      !! end do
      !! !
    !! end if

  else
    !-------------------------------------------------- NOT LogForAqu --
    !! if(DirectSub) call Stop_("No Analytic Jacobian ...")
    !
    !------------------------------------------ material conservation --
    do iCi=1,nCi
      !
      !--- solvent --
      tJac(iCi,1)= tAlfPr(iCi,1)
      do iAx=1,nAx
        tJac(iCi,1)= tJac(iCi,1) &
        &          + tAlfPr(iCi,nCi+iAx)*vMolal(vOrdPr(nCi+iAx))*MWsv
      end do
      !---/
      !--- prim'intern'species
      tJac(iCi,2:nCi)= tAlfPr(iCi,2:nCi)
      !
      if(DirectSub) then
        !--- solvent --
        do iAs=1,nAs
          if(tNuAs(iAs,iCi)/= 0.0D0) then
            ! RESIDUAL ! X= X + tAlfAs(iCi,iAs) *vMolal(iAs)*(vX(isW)*MWSv)
            tJac(iCi,1)= tJac(iCi,1) &
            &          + tAlfAs(iCi,iAs) *vMolal(iAs)*MWSv
          end if
        end do
      else
        !--- second'species
        tJac(iCi,nCi+1:nCi+nAs)= tAlfAs(iCi,1:nAs)
      end if
      !
    end do
    !-----------------------------------------/ material conservation --
    !
    if(.not. DirectSub) then
      !-------------------------------------- sec'species equilibrium --
      !------------------------ deriv' second'species vs prim'species --
      do iAs=1,nAs
        !! tJac(nCi+iAs,1)= vMolal(iAs)*(vX(isW)*MWSv) *(One -vSumNuAs(iAs)) /vX(isW)
        tJac(nCi+iAs,1)= vMolSec(iAs) *(One -vSumNuAs(iAs)) /vX(isW)
        do iCi=1,nCi
          if(iCi/=isW .and. tNuAs(iAs,iCi)/= 0.0D0) &
          !! & tJac(nCi+iAs,iCi)= vMolal(iAs)*(vX(isW)*MWSv) *tNuAs(iAs,iCi) /vX(iCi)
          & tJac(nCi+iAs,iCi)= vMolSec(iAs) *tNuAs(iAs,iCi) /vX(iCi)
        end do
      end do
      !---------------------- deriv' second'species vs second'species --
      do iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= -One
      end do
      !-------------------------------------/ sec'species equilibrium --
    end if
    !
    !-------------------------------------------------/ NOT LogForAqu --
  end if
  !
  if(nEquFas>0) then

    select case(cEquMode)

    case("EQ2")
      do I=1,nEquFas
        do iCi=1,nCi
          tJac(iCi,nF+I)= tAlfEq(iCi,I)
        end do
        if(LogForAqu) then
          tJac(nF+I,isW)=   sum(tNuEq(I,2:nCi))
          tJac(nF+I,2:nCi)= - tNuEq(I,2:nCi)
        else
          tJac(nF+I,isW)=   sum(tNuEq(I,2:nCi)) / vX(isW)
          tJac(nF+I,2:nCi)= - tNuEq(I,2:nCi) / vX(2:nCi)
        end if
      end do

    case("EQ1")
      do I=1,nEquFas
        AffScale= sum(abs(tNuEq(I,:)))
        dRatdX(isW)= - sum(tNuEq(I,2:nCi)) &
        &            *FSafe_Exp(- vFasAff(I)/AffScale) &
        &            /AffScale
        dRatdX(2:nCi)= tNuEq(I,2:nCi) &
        &            *FSafe_Exp(- vFasAff(I)/AffScale) &
        &            /AffScale
        do iCi=1,nCi
          tJac(iCi,1:nCi)= tJac(iCi,1:nCi) &
          &              + tAlfEq(iCi,I) *dXI *dRatdX(1:nCi)
        end do
        !
        tJac(nF+I,1:nCi)= - dXi *dRatdX(1:nCi)
        tJac(nF+I,nF+I)=    One
      end do

    end select

  end if
  !
  if(DebNewt .and. DebJacob) then
    call ShoJacob(nCi,nF,nAs,tJac)
    DebJacob=.false.
  end if
  !
  deallocate(vSumNuAs)
  deallocate(vXf)
  deallocate(vLnXf)
  !
end subroutine Equil_Jacobian

subroutine ShoJacob(nCi,nF,nAs,tJac)
  use M_IoTools,only: GetUnit
  use M_Files,  only: DirLog
  use M_Global_Vars, only: vSpc
  !
  integer, intent(in):: nCi,nAs,nF
  real(dp),intent(in):: tJac(:,:)
  !
  integer::I,J,f
  !
  write(fTrc,'(/,A)') "< ShoJacob"
  write(fTrc,'(A)') "Isaac Newton_&_Carl_Gustav_Jacob_Jacobi -> file=JACOB.LOG"
  call GetUnit(f)
  open(f,file=trim(DirLog)//"debug_jacob.log")
  do I=1,nCi; write(f,'(A7,A1)',advance="no") vSpc(I)%NamSp,    T_; end do
  do I=1,nAs; write(f,'(A7,A1)',advance="no") vSpc(nCi+I)%NamSp,T_; end do
  write(f,*)
  do I=1,nF
    do J=1,nF
      if(tJac(I,J)/=Zero) then; write(f,'(E8.1,A1)',advance="no") tJac(I,J), T_
      else;                     write(f,'(A1,A1)',  advance="no") ".",       T_
      end if
    end do
    write(f,'(I3)') I
  end do
  write(f,'(A7)') "_______"
  do I=1,nF
    do J=1,nF
      if(tJac(I,J)/=Zero) then; write(f,'(G15.4,A1)',advance="no") tJac(I,J), T_
      else                    ; write(f,'(A1,A1)',   advance="no") "0",       T_
      end if
    end do
    write(f,'(I3)') I
  end do
  close(f)
  write(fTrc,'(A,/)') "</ ShoJacob"
end subroutine ShoJacob

end module M_Equil_Jacobian


