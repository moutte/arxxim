module M_Dynam_Jacobian
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_, Stop_
  implicit none
  !
  private
  !
  public:: Dynam_Jacobian
  !
contains

subroutine BuildVecs( & !
& LogForAqu,          & !IN
& LogForMin,          & !IN
& nF,nEqA,nMkA,       & !IN
& vX,                 & !IN
& vXf,vLnXf,vXeq,vXmk)
  !
  use M_Safe_Functions
  use M_Numeric_Const,only: MinExpDP,MaxExpDP,Ln10
  !
  logical, intent(in) :: LogForAqu
  logical, intent(in) :: LogForMin
  integer, intent(in) :: nF,nEqA,nMkA
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vXf(:),vLnXf(:)
  real(dp),intent(out):: vXeq(:),vXmk(:)
  !
  integer:: N
  !
  if(LogForAqu) then
    vXf(1:nF)= FSafe_vExp(vX(1:nF)) !-> avoid overflow
    vLnXf(:)= vX(1:nF)
  else
    vXf(1:nF)= vX(1:nF)
    vLnXf(:)= FSafe_vLog(vX(1:nF))
  end if
  !
  if(nEqA>0) then
    vXeq=Zero
    N= nF
    if(LogForMin) then
      != works with log(mineral amount)
      vXeq(1:nEqA)= FSafe_vExp(vX(N+1:N+nEqA))
    else
      != works directly with mineral amounts
      vXeq(1:nEqA)=vX(N+1:N+nEqA)
    end if
  end if
  !
  if(nMkA>0) then
    vXmk= Zero
    N= nF+nEqA
    if(LogForMin) then
      != works with log(mineral amount)
      vXmk(1:nMkA)= FSafe_vExp(vX(N+1:N+nMkA))
    else
      != works directly with mineral amounts
      vXmk(1:nMkA)=vX(N+1:N+nMkA)
    end if
  end if
  !
end subroutine BuildVecs

subroutine Dynam_Jacobian(vX,tJac)
!--
!-- tJac(1:n,1:n)- Jacobian of Dynam_Residual(1:n)
!--
  use M_Global_Vars,only: vKinFas
  use M_Basis_Vars, only: nCi,nAs,MWsv
  use M_Basis_Vars, only: tAlfAs,tNuAs,tAlfPr
  !
  use M_Dynam_Solve_Vars
  use M_Dynam_Vars, only: &
  & tAlfKin,tNu_Kin, &
  & DirectSub,  &
  & vMolSec,      &
  & LogForMin,LogForAqu, &
  & vMolarVol,  &
  & vLKinActiv,vLEquActiv, &
  & vSurfK,vVmQsK,vVmAct, &
  & UDarcy,dTime,VBox,dX,PhiF,FOut, &
  & DebNewt,DebJacob, &
  & CoupledCoores
  !
  real(dp),dimension(:),                intent(in) :: vX
  real(dp),dimension(size(vX),size(vX)),intent(out):: tJac
  !
  real(dp),dimension(:),allocatable:: vXf,vLnXf,vXmk,vXeq,vSumNuAs
  !
  real(dp):: dum(size(vX))
  real(dp):: Coeff_PhiF,TotFF
  integer :: iPr,jPr,iAs,jAs,iMk,jMk,iMkA,jMkA,I,J
  integer :: nF,nMk,nMkA,nEqA,N
  !
  !nAi= nCi +nAs
  !nF= nAi
  !
  nMk=  size(vKinFas)
  nMkA= count(vLKinActiv)
  nEqA= count(vLEquActiv)
  !
  if(nMkA +nEqA >0) Coeff_PhiF= UDarcy*dTime/dX /PhiF/PhiF /VBox
  !
  if(DirectSub) then  ;  nF= nCi
  else                ;  nF= nCi +nAs
  end if
  !
  allocate(vXf(nF))
  allocate(vLnXf(nF))
  allocate(vXmk(nMkA))
  allocate(vXeq(nEqA))
  allocate(vSumNuAs(nAs))
  !
  do iAs=1,nAs
    vSumNuAs(iAs)= SUM(tNuAs(iAs,2:nCi))
  enddo
  !
  if(nMkA>0) then
    !iMkA=0
    do iMk=1,nMk
      if(vLKinActiv(iMk)) then
        !iMkA= iMkA + 1
        !-------------------------------------- derivative of VM-Surf*VmA*VmQ --
        dVMdLnXf(iMk,1:nCi)= &
        & vSurfK(iMk) *   vVmAct(iMk)*dVmQdLnXf(iMk,1:nCi) ! &
        ! &    + vVmQsK(iMk)*dVmAdLnXf(iMk,1:nF) )
        !
        dVMdLnXm(iMk,1:nMk)= &
        & vVmQsK(iMk)*vVmAct(iMk) * dSRdLnXm(iMk,1:nMk)
        !
        dVMdXm(iMk,1:nMk)= &
        & vVmQsK(iMk)*vVmAct(iMk) * dSRdXm(iMk,1:nMk)
        !-------------------------------------/ derivative of VM-Surf*VmA*VmQ --
      else !not(vLKinActiv(iMk))
        dVMdLnXf(iMk,1:nCi)= Zero
        dVMdLnXm(iMk,1:nMk)= Zero
        dVMdXm  (iMk,1:nMk)= Zero
      end if
    enddo
  end if
  !
  !-----------------------------------------------------------------------------
  call BuildVecs( & !
  & LogForAqu,LogForMin,nF,nEqA,nMkA,vX, & !
  & vXf,vLnXf,vXeq,vXmk)
  !
  ! call BuildVecs(LogForMin,nF,nEqA,nMkA,vX,vXf,vXeq,vXmk)
  !-----------------------------------------------------------------------------
  !
  tJac=Zero
  !
  !-----------------------------------------------------------------------------
  !derive Dynam_Residual(iPr)= & !variable terms only
  !& (One + UDarcy*dTime/dX/ PhiF) * &
  !& ( dot_product(tAlfPr(iPr,1:nCp),vY(1:nCp))   & !# moles in fluid, from prim'species
  !& + dot_product(tAlfAs(iPr,1:nAs),vY(nCp+1:nCp+nAs))) & !# moles in fluid
  !- vTotF(iPr) & !mole nrs in fluid at previous step
  !- vTotInj(iPr)*UDarcy*dTime/dX &
  !+ dTime * dot_product(tAlfKin(iPr,1:nMk),vVm(1:nMk))
  !
  !! PhiF= One - dot_product(vY(nF+1:nF+nMk),vMolarVol(1:nMk)) /VBox  !porosity at time t
  !
  !----------------------------------------------------------------- Jacobian --
  !
  !------------------------------------------------- rows 1:nCi: Prim'Species --
  do iPr=1,nCi
    !
    !----------------------------------------------- 1:nCi-- Aqu'Prim'Species --
    do jPr=1,nCi
      !
      !! NOTE: Fout=  UDarcy *Vbox /dX /PhiF 
      !! m3.s-1 flow related to porous volume
      !
      if(LogForAqu) then
        tJac(iPr,jPr)= &
        & (One + dTime *FOut /VBox) *tAlfPr(iPr,jPr) *vXf(jPr)
        !! + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdLnXf(1:nMk,jPr))
      else
        tJac(iPr,jPr)= &
        & (One + dTime *FOut /VBox) *tAlfPr(iPr,jPr)
      end if
      !
      if(nMkA>0) then
        do iMk=1,nMk
          if(vLKinActiv(iMk)) then
            tJac(iPr,jPr)= tJac(iPr,jPr) &
            &            + dTime *tAlfKin(iPr,iMk) *dVMdLnXf(iMk,jPr)
          end if
        enddo
      end if
      !
    enddo
    !----------------------------------------------/ 1:nCi-- Aqu'Prim'Species --
    !
    if(DirectSub) then
      do jAs=1,nAs
        tJac(iPr,1)= tJac(iPr,1) &
        &          + (One + dTime *FOut /VBox) &
        &            *vMolSec(jAs) *tAlfAs(iPr,jAs) *(One - vSumNuAs(jAs))
      enddo
      do jPr=2,nCi
        do jAs=1,nAs
          tJac(iPr,jPr)= tJac(iPr,jPr) &
          &            + (One + dTime *FOut /VBox) &
          &              *vMolSec(jAs) *tAlfAs(iPr,jAs) *tNuAs(jAs,jPr)
        enddo
      enddo
    else
    !--------------------------------------------- 1:nCi-- Aqu'Second'species --
      do jAs=1,nAs
        ! same as prim'sp,
        ! but we write again because tAlfa is split into tAlfPr U tAlfAs
        if(LogForAqu) then
          tJac(iPr,nCi+jAs)= &
          & (One + dTime *FOut /VBox) *vXf(nCi+jAs) *tAlfAs(iPr,jAs)
        else
          tJac(iPr,nCi+jAs)= &
          & (One + dTime *FOut /VBox) *tAlfAs(iPr,jAs)
        end if
        !! + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdLnXf(1:nMk,jAs))
        !<removed: VM does not depend on second'species>
        !if(nMkA>0) then
        !  do iMk=1,nMk
        !    if(vLKinActiv(iMk)) then
        !      tJac(iPr,jAs)= tJac(iPr,jAs) &
        !      &            + dTime *tAlfKin(iPr,iMk) *dVMdLnXf(iMk,jAs)
        !    end if
        !  enddo
        !end if
        !</removed>
      enddo
    end if
    !----------------------------------------------------/ Aqu'Second'species --
    !
    if(nMkA +nEqA >0) then
      if(DirectSub) then
        TotFF= dot_product(tAlfPr(iPr,1:nCi),vXf(1    :nCi)) &
        &    + dot_product(tAlfAs(iPr,1:nAs),vMolSec(1:nAs))
      else
        TotFF= dot_product(tAlfPr(iPr,1:nCi),vXf(1    :nCi)) &
        &    + dot_product(tAlfAs(iPr,1:nAs),vXf(nCi+1:nCi+nAs))
      end if
    end if
    !
    !-------------------------------------------------- 1:nCi-- equil' phases --
    if(nEqA>0) then
      !
      N= nF
      !
      J= 0
      do iMk=1,nMk
        if(vLEquActiv(iMk)) then
          J=J+1
          !! vFunc(I)= vFunc(I) + tAlfKin(iPr,iMk)*vXeq(J)
          !! CAVEAT: option LogForMin not implemented !!
          if(CoupledCoores) then
            tJac(iPr,N+J)= tAlfKin(iPr,iMk)
          else
            !dPhi/dXeq= - vMolarVol(iMk) /VBox
            !dPsi/dXeq= TotFF*(UDarcy*dTime/dX) *(-1/PhiF/PhiF) *dPhi/dXeq
            tJac(iPr,N+J)= Coeff_PhiF *TotFF *vMolarVol(iMk) + tAlfKin(iPr,iMk)
          end if
        end if
      enddo
    end if
    !---------------------------------------------------------/ equil' phases --
    !
    !------------------------------------------------ 1:nCi -- kinetic phases --
    if(nMkA>0) then
      !
      N= nF+nEqA
      !
      J=0
      do iMk=1,nMk
        !
        if(vLKinActiv(iMk)) then
          !
          J=J+1
          if(LogForMin) then
            !
            if(CoupledCoores) then
              tJac(iPr,N+J)= &
              & dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
            else !-> Phi implicit
              tJac(iPr,N+J)= &
              & Coeff_PhiF *TotFF * vXmk(J) *vMolarVol(iMk) &
              + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
              !**! + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdLnXm(1:nMk,iMk))
            end if
            !
          else
            !
            if(CoupledCoores) then
              tJac(iPr,N+J)= &
              & dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
              !**! & dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
            else !-> Phi implicit
              tJac(iPr,N+J)= &
              & Coeff_PhiF *TotFF *vMolarVol(iMk) &
              + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
              !**! + dTime * dot_product(tAlfKin(iPr,1:nMk),dVMdXm(1:nMk,iMk))
            end if
            !
          end if
          !
        end if !if(vLKinActiv(iMk))
        !
      enddo
    end if
    !--------------------------------------------------------/ kinetic phases --
  enddo
  !-----------------------------------------------/ rows 1:nCi (Prim'Species) --
  !
  !----------------------------------- rows nCi+1:nCi+nAs: Aqu'Second'Species --
  !
  !derive Psi(nCp+iAs)=
  ! dot_product(tNuAs(iAs,1:nCp),vX(1:nCp)) - vX(nCp+iAs) &
  !- (SUM(tNuAs(iAs,1:nCp)) - One) *(vX(isW) + log(MWSv)) + Constant Terms
  !
  if(.not. DirectSub) then
    !
    if(LogForAqu) then
      !------------------------------------------------------------ LogForAqu --
      !-------------------------------- deriv' second'species vs prim'species --
      do iAs=1,nAs
        tJac(nCi+iAs,1)=     One - vSumNuAs(iAs) !Solvent
        do iPr=2,nCi
          tJac(nCi+iAs,iPr)= tNuAs(iAs,iPr)
        enddo
      enddo
      !---/
      !------------------------------ deriv' second'species vs second'species --
      do iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= - One !-> matrix block = -(Identity)
      enddo
      !---/
      !------------------------------------------------------------/LogForAqu --
      !
    else
      !
      !-------------------------------------------------------- NOT LogForAqu --
      !-------------------------------- deriv' second'species vs prim'species --
      do iAs=1,nAs
        tJac(nCi+iAs,1)= vMolSec(iAs) *(One -vSumNuAs(iAs)) /vX(1)
        do iPr=2,nCi
          if(tNuAs(iAs,iPr)/= 0.0D0) &
          & tJac(nCi+iAs,iPr)= vMolSec(iAs) *tNuAs(iAs,iPr) /vX(iPr)
        end do
      end do
      !---/
      !------------------------------ deriv' second'species vs second'species --
      do iAs=1,nAs
        tJac(nCi+iAs,nCi+iAs)= - One !-> matrix block = -(Identity)
      enddo
      !---/
      !--------------------------------------------------------/NOT LogForAqu --
    end if
  end if
  !
  !------------------------------------------ rows nF+1:nF+nEqA: equil'phases --
  if(nEqA>0) then
    !
    N= nF
    !
    J= 0
    do iMk=1,nMk
      if(vLEquActiv(iMk)) then
        J= J +1
        tJac(N +J,1)= - SUM(tNu_Kin(iMk,2:nCi))
        tJac(N +J,2:nCi)=   tNu_Kin(iMk,2:nCi)
      end if
    enddo
  end if
  !-----------------------------------------/ rows nF+1:nF+nEqA: equil'phases --
  !
  !----------------------------------- rows nF+nEqA+1:nF+nEqA+nMk: kin'phases --
  !------ Dynam_Residual(nF+i)- X(nF+i)) -Xi(nF+i) - dTime*vSurfK(i)*vVmF_(i) --
  !
  if(nMkA>0) then
    !
    N= nF +nEqA
    !
    iMkA=0
    do iMk=1,nMk !!!CHECK!!!
      !
      if(vLKinActiv(iMk)) then
        !
        iMkA=iMkA+1
        !
        tJac(N+iMkA,1:nCi)= -dTime * dVMdLnXf(iMk,1:nCi)
        !
        jMkA=0
        do jMk=1,nMk
          if(vLKinActiv(jMk)) then
            !
            jMkA=jMkA+1
            !
            if(jMkA==iMkA) then
              if(LogForMin) then
                tJac(N+iMkA,N+jMkA)= vXmk(iMkA) - dTime * dVMdLnXm(iMk,jMk)
              else
                tJac(N+iMkA,N+jMkA)= One - dTime * dVMdXm(iMk,jMk)
              end if
            else
              if(LogForMin) then
                tJac(N+iMkA,N+jMkA)= - dTime * dVMdLnXm(iMk,jMk)
              else
                tJac(N+iMkA,N+jMkA)= - dTime * dVMdXm(iMk,jMk)
              end if
            end if
            !
          end if
        enddo
      end if
      !
    enddo
  end if !if(nMkA>0)
  !----------------------------------/ rows nF+nEqA+1:nF+nEqA+nMk: kin'phases --
  !
  deallocate(vXf)
  deallocate(vLnXf)
  deallocate(vXmk)
  deallocate(vXeq)
  deallocate(vSumNuAs)
  !
  !-----------------------------------------------------------/ DynamJacobian --
  !
  !DebNewt= (iDebug>2)
  if(DebNewt .and. DebJacob) then !
    !
    !
    write(31,'(A)') "Dynamic Jacobian"
    !
    !--- print max value of each column
    dum=MAXVAL(ABS(tJac),DIM=2)
    do I=1,nF+nEqA+nMkA; write(31,'(E8.1,A1)',advance="no") dum(I), T_; enddo
    write(31,'(A)') "MAXVAL(COL)"
    !---/
    !
    do I=1,nF+nEqA+nMkA
      do J=1,nF+nEqA+nMkA
        if(tJac(I,J)/=Zero) then; write(31,'(F7.1,A1)',advance="no") tJac(I,J), T_
        else;                     write(31,'(A1,  A1)',advance="no") ".",       T_
        end if
      enddo
      write(31,'(I3)') I
    enddo
    !
    DebJacob=.false.
    !
  end if
  !
end subroutine Dynam_Jacobian

end module M_Dynam_Jacobian
