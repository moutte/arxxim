module M_Box_Jacobian_Dynam
  
  !==============================================
  ! Residual Computation Similar to Dynam Model
  !==============================================

  use M_Kinds
  use M_Trace
  
  use M_Box_Residual_Vars
  
  use M_Box_RefState_Vars,   only : vG0rt_Aq, vLnK_As
  use M_Box_Kinetics_Vars,   only : vKinRateK, Reaction_Limiter
  use M_Box_Solver_Vars,     only : ThetaRate
  use M_Box_System_Vars
  
  implicit none
  private
  
  !-- Public Functions
  public:: Box_Jacobian_Dynam
  
contains

  !---

  subroutine Box_Jacobian_Dynam(vX,tJac) 
    !==================================================
    ! Purpose   : compute the jacobian for the box 
    ! problem used in reaction transport coupled mode
    !==================================================
    use M_Matrix_Utils
    use M_Box_Utils
    !--
    use M_Box_Prepare
    use M_Box_Entries
    !--
    use M_Box_Vars     

    implicit none

    real(dp),dimension(:),                 intent(in) :: vX
    real(dp),dimension(size(vX),size(vX)), intent(out):: tJac

    !-- local variables
    
    integer :: iPr, iAs, iMk, iAq
    integer :: iEq, iEntry, nEq
    real(dp):: dni, dXf_dLnXf, dXm_dXm
    real(dp):: LnQ, LnK
    
    !---------------------------------------------------
    tJac=Zero

    !---- Prepare Terms --------------------------------
    call Box_Prepare_Jacobian(vX)
          
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
    
    do iPr = 1,nPr
       iAq = iPr
       vEqBal(1:nCp) = vEqBal(1:nCp) &
            &  + tAlfPr(1:nCp,iPr) + ( One + dTime*Vout  ) * vXf(iAq)
       !---
       iEntry = iPr
       dXf_dLnXf = vXf(iAq)
       tJac(1:nCp,iEntry ) =  tJac(1:nCp, iEntry) &
            & + tAlfPr(1:nCp,iPr) * ( One + dTime*Vout  ) * dXf_dLnXf
    end do

    do iAs = 1,nAs
       iAq = nPr + iAs
       vEqBal(1:nCp) = vEqBal(1:nCp) & 
           &  + tAlfAs(1:nCp,iAs)* ( One + dTime*Vout  ) * vXf(iAq)
       !---
       iEntry = nPr + iAs
       dXf_dLnXf = vXf(iAq)
       tJac(1:nCp,iEntry ) = tJac(1:nCp,iEntry ) &
            &    + tAlfAs(1:nCp,iAs) * ( One + dTime*Vout ) * dXf_dLnXf
    end do

    do iMk = 1,nMk
       
       dni = dTime*( ThetaRate*vVm(iMk) + (One-ThetaRate)* vKinRateK(iMk) ) 
       vEqBal(1:nCp) = vEqBal(1:nCp) + tAlfMk(1:nCp,iMk) * dni
       !---
              
       do iAq= 1, nAq
          tJac(1:nCp,iAq) = tJac(1:nCp,iAq) &
            &    + tAlfMk(1:nCp,iMk) * dTime * ThetaRate* dVm_dLnXf(iMk,iAq) 
       end do
       
       iEntry = nPr + nAs + iMk
       tJac(1:nCp,iEntry) = tJac(1:nCp,iEntry) &
            &    + tAlfMk(1:nCp,iMk) * dTime * ThetaRate* dVm_dXm(iMk,iMk)
    end do

     vEqBal(1:nCp) = vEqBal(1:nCp) - vTotf(1:nCp) + dTime* vQInj(1:nCp)

    !---- Mass action law equations for secondary species  
    
    vEqEquil = Zero

    do iAs = 1,nAs       

       iEq = nPr + iAs

       !-- activity contribution
       LnQ = Zero
       do iPr= 1,nPr
          iAq = iPr
          LnQ = LnQ + tNuAs(iAs,iPr) * vLnai(iAq)
          !--
          tJac(iEq,1:nAq) = tJac(iEq,1:nAq) + tNuAs(iAs, iPr) * dLnai_dLnXf(iAq,1:nAq)
       end do
       
       iAq = nPr + iAs
       LnQ = LnQ - One*vLnai(iAq)
       !--
       tJac(iEq,1:nAq) = tJac(iEq,1:nAq) - One*dLnai_dLnXf(iAq, 1:nAq)
       
       !-- reference state potential contribution
       LnK = vLnK_As(iAs)

       !-- add the contributions
       vEqEquil(iAs) = LnQ - LnK
       
    end do
    
    !---- Mass balance equations for kinetical components 

    vEqKin(1:nMk) = Zero

    do iMk=1,nMk       
       iEq = nPr + nAs + iMk
              
       vEqKin(iMk) = (vXm(iMk) - vMolK(iMk)) &
           & - dTime*( ThetaRate*vVm(iMk) + (One-ThetaRate)* vKinRateK(iMk) ) 
              
       dXm_dXm = One

       ! accumulation
       iEntry = nAq + iMk
       tJac(iEq,iEntry) =  tJac(iEq,iEntry) + dXm_dXm
       
       ! kinetic rate 
       tJac(iEq,1:nAq) =  tJac(iEq,1:nAq)  &
            &           - dTime *  ThetaRate* dVm_dLnXf(iMk,1:nAq) 
       

       iEntry = nAq + iMk
       tJac(iEq,iEntry) =   tJac(iEq,iEntry) &
            &           - dTime *  ThetaRate* dVm_dXm(iMk,iMk)
              
    end do

    !---- Store the result in vEqSystem and form Residual 
    
    vEqSystem(1:nCp)         = vEqBal(1:nCp)
    vEqSystem(nCp+1:nCp+nAs) = vEqEquil(1:nAs)
    vEqSystem(nAq+1:nAq+nMk) = vEqKin(1:nMk)
    
    !// Debug
    !call Matrix_Clear_Zeros(tJac)    
    call Box_Info_Jacobian_Pure(vX, tJac)

    !---- Apply Logarithm Scaling Options
    
    call Box_Jacobian_LogTransform(tJac)
    !--- Clear Zeros Coefficients
    !call Matrix_Clear_Zeros(tJac)  
 
    !====== end OF BOX JACOBIAN !! =====================

    !// Debug
    call Box_Info_Jacobian(vX, tJac) 
    call DEBUG_JACOBIAN(vX, tJac) 
    
    !//Debug
    !!call Box_Show_Jacobian_Error(vX,tJac)

  end subroutine Box_Jacobian_Dynam

  !---

  subroutine Box_Info_Jacobian_Pure(vX, tJac)
    use M_Box_Info
    use M_EchoMat
    use M_Box_Debug_Vars, only : LDebug_Jacobian
    implicit none
    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(:,:),intent(in) :: tJac
    
    if (LDebug_Jacobian) then
       write(*,*) '--------- Jacobian Before LogTransform -----------'
       call echomat_real(tJac, size(vX), size(vX), 'tJac')
    end if
    
  end subroutine Box_Info_Jacobian_Pure
  
  !---
  
  subroutine Box_Info_Jacobian(vX, tJac)
    !==================================================
    ! Purpose   : Info Debug Jacobian
    !==================================================
    use M_Box_Info
    use M_EchoMat
    use M_Box_Vars
    use M_Box_Debug_Vars, only : LDebug_Jacobian
    implicit none
    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(:,:),intent(in) :: tJac

    !// Debug  
    if (LDebug_Jacobian) then
       call Info_("")
       call Info_("========================================================")
       call Info_(" JACOBIAN DEBUG ")
       call Info_("========================================================")

       !// Physical vectors
       call Info_("----------------- Natural Vectors -------------------")
       call echovec_real(vMolF, nAq, 'vMolF','T')
       call echovec_real(vMolK, nMk, 'vMolK','T')
       call Box_Info_Amount(vMolF,vMolK,'Vectors vMolF, vMolK ')

       !// Debug  
       call Info_("----------------- Residual ------------------")
       call echovec_real(vX, size(vX), 'vX','T')
       call Box_Info_Variable(vX,'Variable vX')
       
       call echovec_real(vEqSystem, size(vX), 'vEqSystem','T')
       call Box_Info_Equation(vEqSystem,'Residual vEqSystem')
       
       write(*,*) '--------- Jacobian -----------'
       call echomat_real(tJac, size(vX), size(vX), 'tJac')
    end if
   
    if (LDebug_Jacobian) call Box_Show_Jacobian_Error(vX,tJac)
    
  end subroutine Box_Info_Jacobian

  !---

  subroutine Box_Show_Jacobian_Error(vX, tJac)
    !=====================================================
    ! Compute Errors Between tJac and tJacNum  
    !=====================================================
    use M_Numeric_Tools
    use M_Echomat
    use M_Box_Residual
    use M_Matrix_Utils
    use M_Info_Value
    !---
    implicit none
    !--    
    real(dp), intent(in) :: vX(:)    
    real(dp), intent(in):: tJac(:,:)
    !--
    real(dp) :: tJacNum(size(vX),size(vX))
    real(dp) :: tError(size(vX),size(vX))
    real(dp) :: tMaskError(size(vX),size(vX))
    real(dp) :: vResidualX(size(vX))
    real(dp) :: X(size(vX))
    real(dp) :: ErrorL1, ErrorL2, ErrorLinf
    integer :: N
    !--
    call Info_(" COMPUTE JACOBIAN ERRORS ")
    !// compute numerical jacobian
    X = vX
    N= size(vX)
    vResidualX = Box_Residual(X)
    call Jacobian_Numeric_Bis(Box_Residual,X,vResidualX,tJacNum)
    
    !// show differences
    write(*,*) "----------------- Entries ------------------"
    call echovec_real(vX, size(vX), 'vX','T')

    write(*,*) "----------------- Residual ------------------"
    call echovec_real(vResidualX, size(vX), 'vResidualX','T')

    write(*,*) '----------------- Jacobian Analytic ---------'
    call echomat_real(tJac, size(vX), size(vX), 'tJac')
    
    write(*,*) '----------------- Jacobian Numeric ----------'
    call echomat_real(tJacNum, size(vX), size(vX), 'tJacNum')
    
    tError = tJac-tJacNum
    ErrorL1   = Matrix_Norm_L1(tError)
    ErrorL2   = Matrix_Norm_L2(tError)
    ErrorLinf = Matrix_Norm_Linf(tError)

    write(*,*) '---------------- Error Jacobian -------------'
    call echomat_real(tError, size(vX), size(vX), 'tErrorJac tJac - tJacNum ')
    
    write(*,*) '---------------- Error Jacobian Filtered ----'
    call Matrix_Error(tError, ErrorL1/1000, tMaskError)
    call echomat_real(tMaskError, size(vX), size(vX), 'Mask Error')
    
    call Info_Value('ErrorL1  ', ErrorL1)
    call Info_Value('ErrorL2  ', ErrorL2)
    call Info_Value('ErrorLinf', ErrorLinf)
    
  end subroutine Box_Show_Jacobian_Error

  !---

  subroutine DEBUG_JACOBIAN(vX,tJac) 
    !=================================================================
    ! purpose : show the Jacobian
    !=================================================================
    use M_Trace
    use M_FILES
    use M_echomat
    use M_Box_Newton_Vars, only : DebJacob,  DebNewt
    implicit none
    real(dp),dimension(:), intent(in):: vX
    real(dp),dimension(:,:), intent(in):: tJac
    real(dp),dimension(size(tJac,1)):: dum
    !
    integer :: i,n, j
    n = size(tJac,1)
    !--
    DebJacob = .false.
    DebNewt  = .false.
    if(DebNewt .and. DebJacob) then !
       dum=MAXVAL(ABS(tJac),DIM=2) !Loop over rows to get the implicit scaling information
       write(fTrc,'(A)') "Box Jacobian"
       do i=1,n
          write(fTrc,'(E8.1,A1)',advance="no") dum(i), T_
       enddo
       write(fTrc,'(/,A,/)') "=MAX"
       do i=1,n
          do j=1,n
             if(tJac(i,j)/=Zero) then
                write(fTrc,'(E8.1,A1)',advance="no") tJac(i,j), T_
             else
                write(fTrc,'(A1,  A1)',advance="no") ".", T_
             end if
          enddo
          write(fTrc,'(I3)') i
       enddo
       !
    end if

  end subroutine DEBUG_JACOBIAN

end module M_Box_Jacobian_Dynam
