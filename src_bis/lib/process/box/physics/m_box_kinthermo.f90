module M_Box_KinThermo

  use M_Kinds
  use M_Trace

  use M_Box_Vars
  use M_Box_Thermo_Vars

  implicit none
  private

  character(len=10), parameter :: cQsKModel = "EXPLICIT"
  
  !// Public functions
  public :: BOX_Compute_Omega
  public :: BOX_Compute_QsK
  
contains

  subroutine BOX_Compute_Omega
     
    !==============================================================
    ! Purpose   : compute QsK OmegaK values for export
    !==============================================================
    use M_Box_System_Vars
    use M_Box_RefState_Vars
    use M_Safe_Functions

    implicit none
    !---    
    real(dp):: LnK, LnQ
    real(dp):: QsK, LnQsK   
    !---
    integer :: iMk, iAq, iPr
    real(dp) :: Lnai, Lnmi, LnActMk
    
    call Info_("Box_Compute_Affinity")
    
    do iMk = 1, nMk
       !-- activity contribution
       LnQ = Zero
       do iPr= 1,nPr
          iAq = iPr
          LnQ = LnQ + tNuMk(iMk,iPr) * vLnAct(iAq)
       end do
       LnActMk = Zero
       LnQ = LnQ - One*LnActMk
       
       !-- reference state potential contribution
       LnK =  vLnK_Mk(iMk)

       !-- add the two contributions
       LnQsK = LnQ - LnK
       QsK = FSafe_Exp(LnQsK)
       vOmegaK(iMk) = QsK
    end do

  end subroutine BOX_Compute_Omega

  !---
  
  subroutine BOX_Compute_QsK( &
       vLnXf,       & ! IN
       vLnai,       & ! IN
       dLnai_dLnXf, & ! IN
       vLnQsK,       & ! OUT
       dLnQsK_dLnXf, & ! OUT 
       vQsK,         & ! OUT
       dQsK_dLnXf )    ! OUT 

    !==============================================================
    ! Purpose   : compute QsK and its derivatives 
    !==============================================================
    implicit none
    !---
    real(dp), intent(in) :: vLnXf(:)
    !---
    real(dp), intent(in) :: vLnai(:)
    real(dp), intent(in) :: dLnai_dLnXf(:,:)
    !---
    real(dp), intent(out) :: vQsK(:)
    real(dp), intent(out) :: dQsK_dLnXf(:,:)
    !---
    real(dp), intent(out) :: vLnQsK(:)
    real(dp), intent(out) :: dLnQsK_dLnXf(:,:)

    call Info_("Box_Compute_QsK")
    
    select case ( cQsKModel )

    case ("EXPLICIT")
       call BOX_Compute_QsK_Explicit(  &
             vLnXf,        & ! IN
             vLnQsK,       & ! OUT
             dLnQsK_dLnXf, & ! OUT 
             vQsK,         & ! OUT
             dQsK_dLnXf )    ! OUT 

    case("implicit")
       call BOX_Compute_QsK_Implicit( &
             vLnai,        & ! IN
             dLnai_dLnXf,  & ! IN
             vLnQsK,       & ! OUT
             dLnQsK_dLnXf, & ! OUT 
             vQsK,         & ! OUT
             dQsK_dLnXf  )   ! OUT 
       
    end select

    call BOX_Info_QsK( vLnQsK, vQsK )

  end subroutine BOX_Compute_QsK

  !---

  subroutine BOX_Compute_QsK_Explicit( &
       vLnXf,        & ! IN
       vLnQsK,       & ! OUT
       dLnQsK_dLnXf, & ! OUT 
       vQsK,         & ! OUT
       dQsK_dLnXf )    ! OUT 

    !================================================================================
    ! Purpose   : compute QsK and its derivatives with specific activities
    !             Assume LnActw and vLnGam are explicit
    !================================================================================ 
    ! LnQ= < tNuMk(iMk,1:nCp)* vLnAct(1:nCp) > - vLnActMk(iMk) 
    ! LnK= < tNuMk(iMk,1:nCp), vSpc(1:nCp)%G0rt>  - vG0RtMK(iMk)
    ! LnQsK = LnQ - LnK  
    !================================================================================
    use M_Basis_Vars,  only : MWSv
    use M_Box_System_Vars
    use M_Box_RefState_Vars
    use M_Safe_Functions

    implicit none
    real(dp), intent(in) :: vLnXf(:)
    !---
    real(dp), intent(out) :: vLnQsK(:)
    real(dp), intent(out) :: dLnQsK_dLnXf(:,:)
    !---
    real(dp), intent(out) :: vQsK(:)
    real(dp), intent(out) :: dQsK_dLnXf(:,:)
    !---    
    real(dp):: LnK, LnQ
    real(dp):: QsK, LnQsK   
    !---
    integer :: iMk, iAq, iPr
    real(dp) :: Lnai, Lnmi, LnActMk, dLnActMk_dLnXf

    dLnQsK_dLnXf= Zero
    
    do iMk = 1, nMk

       !-- activity contribution
       LnQ = Zero
       do iPr= 1,nPr
          iAq = iPr
          if (iAq.eq.iAqH2O) then
             Lnai = vLnAct(iAqH2O)
             LnQ = LnQ + tNuMk(iMk,iPr) * Lnai
          else
             Lnmi = vLnXf(iAq) - vLnXf(iAqH2O) - log(MWSv) 
             Lnai = Lnmi + vLnGam(iAq) 
             LnQ = LnQ + tNuMk(iMk,iPr) * Lnai
             !---
             dLnQsK_dLnXf(iMk,iAq)    = dLnQsK_dLnXf(iMk,iAq)    + tNuMk(iMk,iPr)*One
             dLnQsK_dLnXf(iMk,iAqH2O) = dLnQsK_dLnXf(iMk,iAqH2O) - tNuMk(iMk,iPr)*One
          end if
       end do

       LnActMk = Zero
       dLnActMk_dLnXf = Zero
       LnQ = LnQ - One*LnActMk
       !--
       dLnQsK_dLnXf(iMk,iAq)    = dLnQsK_dLnXf(iMk,iAq) - One*dLnActMk_dLnXf

       !-- reference state potential contribution
       LnK =  vLnK_Mk(iMk)

       !-- add the two contributions
       LnQsK = LnQ - LnK
       vLnQsK(iMk) = LnQsK

       !-- recover QsK from LnQsK
       QsK = FSafe_Exp(LnQsK)
       vQsK(iMk) =  QsK
       dQsK_dLnXf(iMk,1:nAq) = dLnQsK_dLnXf(iMk,1:nAq) * QsK

    end do

  end subroutine BOX_Compute_QsK_Explicit

  !---

  subroutine BOX_Compute_QsK_Implicit( &
       vLnai,        & ! IN
       dLnai_dLnXf,  & ! IN
       vLnQsK,       & ! OUT
       dLnQsK_dLnXf, & ! OUT 
       vQsK,         & ! OUT
       dQsK_dLnXf  )   ! OUT 

    !=====================================================================
    ! Purpose   : compute QsK and its derivatives with implicit activity
    ! Idem Mass Action Law for Aqueous Secondary Species
    !=====================================================================
    ! LnQ= < tNuMk(iMk,1:nCp)* vLnAct(1:nCp) > - vLnActMk(iMk) 
    ! LnK= < tNuMk(iMk,1:nCp), vSpc(1:nCp)%G0rt>  - vG0RtMK(iMk)
    ! LnQsK = LnQ - LnK  
    !=====================================================================
    use M_Box_System_Vars
    use M_Box_RefState_Vars
    use M_Safe_Functions
    implicit none

    real(dp), intent(in) :: vLnai(:)
    real(dp), intent(in) :: dLnai_dLnXf(:,:)
    !---
    real(dp), intent(out) :: vLnQsK(:)
    real(dp), intent(out) :: dLnQsK_dLnXf(:,:)
    !---
    real(dp), intent(out) :: vQsK(:)
    real(dp), intent(out) :: dQsK_dLnXf(:,:)
    !---
    real(dp):: LnK, LnQ
    real(dp):: QsK, LnQsK   
    !---
    integer :: iMk, iAq, iPr
    real(dp) :: LnActMk, dLnActMk_dLnXf
    !---

    dLnQsK_dLnXf= Zero
     
    do iMk = 1, nMk

       !-- activity contribution
       LnQ = Zero
       do iPr= 1,nPr
          iAq = iPr
          LnQ = LnQ + tNuMk(iMk,iPr) * vLnai(iAq)
          !--
          dLnQsK_dLnXf(iMk, 1:nAq) = dLnQsK_dLnXf(iMk, 1:nAq) &
               &                     + tNuMk(iMk, iPr) * dLnai_dLnXf(iAq,1:nAq)
       end do
       LnActMk = Zero
       dLnActMk_dLnXf = Zero
       LnQ = LnQ - One*LnActMk
       !--
       dLnQsK_dLnXf(iMk,1:nAq) = dLnQsK_dLnXf(iMk, 1:nAq) - One*dLnActMk_dLnXf

       !-- reference state potential contribution
       LnK =  vLnK_Mk(iMk)

       !-- add the two contributions
       LnQsK = LnQ - LnK
       vLnQsK(iMk) = LnQsK

       !-- recover QsK from LnQsK
       QsK = FSafe_Exp(LnQsK)
       vQsK(iMk) =  QsK
       dQsK_dLnXf(iMk,1:nAq) = dLnQsK_dLnXf(iMk,1:nAq) * QsK

    end do

  end subroutine BOX_Compute_QsK_Implicit

  !---

  subroutine BOX_Info_QsK( vLnQsK, vQsK )

    !==================================================
    ! Purpose   : Info Debug Saturation Index QsK
    !==================================================
    use M_Global_Vars,    only : vSpc
    use M_Numeric_Const,  only : Ln10
    use M_Box_System_Vars,    only : vOrdMk
    use M_Box_Debug_Vars, only : LDebug_KinThermo
    implicit none

    real(dp), intent(in) :: vLnQsK(:)
    real(dp), intent(in) :: vQsK(:)

    integer :: iMk, iSpc    

    call Info_("BOX_Info_QsK")

    if ( LDebug_KinThermo ) then
       write(*,'(A6,1X,A20, 3A16)') &
            '  iMk ', 'iMkSpc%Name ', 'vLnQsK(iMk)', 'vLogQsK(iMk)', 'vQsK(iMk)'

       write(*,'(A6,1X,A20, 3A16)') &
            '------', '--------------------', '-----------', '-----------', '-----------'
       do iMk = 1, nMk
          iSpc = vOrdMk(iMk)
          write(*,'(I6,1X,A20, 3G16.6)') &
               iMk, vSpc(iSpc)%NamSp, vLnQsK(iMk),  vLnQsK(iMk)/Ln10, vQsK(iMk) 
       end do
    end if

  end subroutine BOX_Info_QsK

end module M_Box_KinThermo
