module M_Box_Residual_Dynam

  !==============================================
  ! Residual Computation Similar to Dynam Model
  !==============================================

  use M_Kinds
  use M_Trace

  use M_Box_Residual_Vars
 
  use M_Box_Vars  
  use M_Box_RefState_Vars, only : vG0rt_Aq, vLnK_As
  use M_Box_Kinetics_Vars, only : vKinRateK, Reaction_Limiter
  use M_Box_Solver_Vars,   only : ThetaRate
  use M_Box_System_Vars
  use M_Box_Thermo_Vars, only : vLnAct, vLnGam

  implicit none
  private

  !-- Public Functions
  public:: Box_Residual_Dynam

contains

    !-------

  function Box_Residual_Dynam(vX)
    !==================================================
    ! Purpose : compute the residual for the box 
    ! problem used in reaction transport coupled mode
    !==================================================
    use M_Box_Utils
    use M_Box_Prepare
    !--
    
    implicit none

    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(size(vX))     :: Box_Residual_Dynam
    
    !-- local variables
    integer :: iPr,iAs, iMk, iAq
    real(dp):: dni
    real(dp):: Lnai
    real(dp):: LnQ, LnK
   
    !---------------------------------------------------
    Box_Residual_Dynam = Zero

    !---- Prepare Terms --------------------------------
    
    call Box_Prepare_Residual(vX)
    
    !---- Compute total of component in fluid 
    vTotf(1:nCp) = Zero
    
    do iPr = 1,nPr
       iAq = iPr
       vTotf(1:nCp) = vTotf(1:nCp) + tAlfPr(1:nCp,iPr)*vMolf(iAq)
    end do
    
    do iAs = 1,nAs
       iAq = nPr + iAs
       vTotf(1:nCp) = vTotf(1:nCp) + tAlfAs(1:nCp,iAs)*vMolf(iAq)
    end do    
    
    !---- Mass balance equations for equilibrium components
    vEqBal(1:nCp) = Zero
    
    do iPr = 1,nPr
       iAq = iPr       
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfPr(1:nCp,iPr)* &
            &  ( One + dTime*Vout ) * vXf(iAq)
    end do
    
    do iAs = 1,nAs
       iAq = nPr + iAs       
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfAs(1:nCp,iAs)* &
            &  ( One + dTime*Vout ) * vXf(iAq)
    end do
    
    do iMk = 1,nMk
       dni =  dTime * ( ThetaRate*vVm(iMk) + (One-ThetaRate)*vKinRateK(iMk) )
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfMk(1:nCp,iMk)* dni
    end do
    
    vEqBal(1:nCp) = vEqBal(1:nCp) - vTotf(1:nCp)
    vEqBal(1:nCp) = vEqBal(1:nCp) + dTime* vQInj(1:nCp)

    !---- Mass action law equations for secondary species      
    
    vEqEquil(1:nAs) = Zero

    do iAs = 1,nAs
      
       !-- activity contribution 
       LnQ = Zero
       do iPr= 1,nPr
          iAq = iPr
          LnQ = LnQ + tNuAs(iAs,iPr) * vLnai(iAq)
       end do
       
       iAq = nPr + iAs
       LnQ = LnQ - One*vLnai(iAq)
       
       !-- reference state potential contribution 
       LnK = vLnK_As(iAs)
       
       !-- add the contributions
       vEqEquil(iAs) = LnQ - LnK
       
    end do
    
    !---- Mass balance equations for kinetics components
    
    vEqKin(1:nMk) = Zero
        
    do iMk=1,nMk
       ! iSpc = nPr + nAs + iMk
       vEqKin(iMk) = (vXm(iMk) - vMolK(iMk)) &
            - dTime*( ThetaRate*vVm(iMk) + (One- ThetaRate)*vKinRateK(iMk) )
    end do

    !---- Store the result in vEqSystem and form Residual 
    vEqSystem = Zero
    vEqSystem(1:nCp)         = vEqBal(1:nCp)
    vEqSystem(nCp+1:nCp+nAs) = vEqEquil(1:nAs)
    vEqSystem(nAq+1:nAq+nMk) = vEqKin(1:nMk)
    
    BOX_Residual_Dynam = vEqSystem
    
    call Box_Info_Residual (vX)   
    
  end function Box_Residual_Dynam

  !---
  
  subroutine Box_Info_Residual(vX)
    !==================================================
    ! Purpose   : Info Debug Residual
    !==================================================

    use M_Box_Info
    use M_EchoMat
    use M_Box_Debug_Vars, only : LDebug_Residual
    implicit none
    real(dp),dimension(:),intent(in) :: vX
    integer :: fres = 44
    
    if (LDebug_Residual) then
       !// Debug  
       call Box_Info_Variable(vX,'Variable vX')
       call Box_Info_Equation(vEqSystem,'Residual vEqSystem')
       call Box_Info_TotF(vTotf,'Total Fluid vTotF')
    end if

    if (LDebug_Residual) then
    Write(fres,*) " >> Info Residual ..."
    Write(fres,*) "================================="
    Write(fres,*) "Vout =", Vout 
    Write(fres,*) "dTime=", dTime
    Write(fres,*) "================================="

    call EchoVec_Real(vTotf,       size(vTotf),    TRTAG='T', Unit=fres)
    call EchoVec_Real(vQInj,       size(vQInj),    TRTAG='T', Unit=fres)
    Write(fres,*) "---------------------------------"
    call EchoVec_Real(vLnAi,       size(vLnAi),       TRTAG='T', Unit=fres)
    call EchoVec_Real(vLnAct,      size(vLnAct),     TRTAG='T', Unit=fres)
    call EchoVec_Real(vLnGam,      size(vLnGam),     TRTAG='T', Unit=fres)
    call EchoVec_Real(vLnK_As,     size(vLnK_As),     TRTAG='T', Unit=fres)
    Write(fres,*) "---------------------------------"
    call EchoVec_Real(vXf,         size(vXf),      TRTAG='T', Unit=fres)
    call EchoVec_Real(vMolF ,      size(vXf),      TRTAG='T', Unit=fres)
    call EchoVec_Real(vMolF - vXf, size(vXf),      TRTAG='T', Unit=fres)
    Write(fres,*) "---------------------------------"
    call EchoVec_Real(vXm ,        size(vXm),         TRTAG='T', Unit=fres)
    call EchoVec_Real(vVm ,        size(vVm),         TRTAG='T', Unit=fres)
    call EchoVec_Real(vKinRateK ,  size(vKinRateK),   TRTAG='T', Unit=fres)
    Write(fres,*) "---------------------------------"
    call EchoVec_Real(vX,          size(vX),TRTAG='T', unit=fres)
    call EchoVec_Real(vEqSystem,   size(vEqSystem),TRTAG='T', unit=fres)
    Write(fres,*) " "
    end if   

  end subroutine Box_Info_Residual
  
end module M_Box_Residual_Dynam
