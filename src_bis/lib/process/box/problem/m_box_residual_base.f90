module M_Box_Residual_Base

  use M_Kinds
  use M_Trace

  use M_Box_Residual_Vars
 
  use M_Box_Vars  
  use M_Box_RefState_Vars, only : vG0rt_Aq, vLnK_As
  use M_Box_Kinetics_Vars, only : vKinRateK, Reaction_Limiter
  use M_Box_Solver_Vars,   only : ThetaRate
  use M_Box_System_Vars
  
  implicit none
  private

  !-- Public Functions
  public:: Box_Residual_Base

contains

  !-------

  function Box_Residual_Base(vX)
    !==================================================
    ! Purpose : compute the residual for the box 
    ! problem used in reaction transport coupled mode
    !==================================================
    use M_Box_Utils
    use M_Box_Prepare
    !--
    
    implicit none

    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(size(vX))     :: Box_Residual_Base
    
    !-- local variables
    integer :: iPr,iAs, iMk, iAq
    real(dp):: dni
    real(dp):: Lnai
    real(dp):: theta
    real(dp):: LnQ, LnK
    
    !---------------------------------------------------
    Box_Residual_Base = Zero

    !---- Prepare Terms --------------------------------
    
    call Box_Prepare_Residual(vX)
    
    !---- Mass balance equations for equilibrium components
    vEqBal(1:nCp) = Zero
    
    do iPr = 1,nCp
       iAq = iPr
       dni =  vXf(iAq) - vMolf(iAq)
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfPr(1:nCp,iPr)* &
            &  ( dni + dTime*Vout *vXf(iAq)  ) 
    end do
    
    do iAs = 1,nAs
       iAq = nPr + iAs
       dni =  vXf(iAq) - vMolf(iAq)
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfAs(1:nCp,iAs)* &
            &  ( dni + dTime*Vout *vXf(iAq) )
    end do
    
    do iMk = 1,nMk
       dni =  vXm(iMk) - vMolK(iMk)
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfMk(1:nCp,iMk)* dni
    end do
    
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
    
    theta = ThetaRate 
    do iMk=1,nMk
       ! iSpc = nPr + nAs + iMk
       vEqKin(iMk) = (vXm(iMk) - vMolK(iMk)) &
            - dTime* Reaction_Limiter *( theta*vVm(iMk) + (One - theta)*vKinRateK(iMk) )
    end do

    !---- Store the result in vEqSystem and form Residual 
    
    vEqSystem(1:nCp)         = vEqBal(1:nCp)
    vEqSystem(nCp+1:nCp+nAs) = vEqEquil(1:nAs)
    vEqSystem(nAq+1:nAq+nMk) = vEqKin(1:nMk)
    
    BOX_Residual_Base = vEqSystem
    
    call Box_Info_Residual (vX)   
    
  end function Box_Residual_Base

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
    
    if (LDebug_Residual) then
       !// Debug  
       call echovec_real(vX, size(vX), 'vX','T')
       call Box_Info_Variable(vX,'Variable vX')
       
       call echovec_real(vEqSystem, size(vX), 'vEqSystem','T')
       call Box_Info_Equation(vEqSystem,'Residual vEqSystem')
    end if
    
    if (LDebug_Residual) then
       call EchoVec_Real(vXf,         size(vXf),      TRTAG='T', Unit=33)
       call EchoVec_Real(vMolF - vXf, size(vXf),      TRTAG='T', Unit=35)
       call EchoVec_Real(vMolF ,      size(vXf),      TRTAG='T', Unit=36)
       call EchoVec_Real(vEqSystem,   size(vEqSystem),TRTAG='T', unit=34)
    end if   
    
  end subroutine Box_Info_Residual
  
end module M_Box_Residual_Base
