module M_Box_RefState

  use M_Kinds,only: dp
  use M_Box_RefState_Vars
  use M_Box_Vars
  
  implicit none
  private
  
  !// Public Functions
  public :: BOX_Compute_RefStatePotential
  public :: BOX_Compute_LogK

  !// Private Functions
  private :: BOX_Compute_RefStatePotential_Aq
  private :: BOX_Compute_RefStatePotential_Mk

  private :: BOX_Compute_LogK_As
  private :: BOX_Compute_LogK_Mk
  
contains

  subroutine BOX_Compute_RefStatePotential
    !==================================================
    ! Purpose   : compute G0Rt for Aq and Mk species
    !==================================================
    call Info_("BOX_Compute_RefStatePotential")
    call BOX_Compute_RefStatePotential_Aq
    call BOX_Compute_RefStatePotential_Mk

  end subroutine BOX_Compute_RefStatePotential

  !---

  subroutine BOX_Compute_RefStatePotential_Aq
    !==================================================
    ! Purpose   : compute G0Rt for Aqueous Species
    !==================================================
    use M_Global_Vars, only : vSpc
    use M_Box_System_Vars
    implicit none
    
    integer :: iAq, iSpc
    do iAq = 1, nAq
       iSpc = vOrdAq(iAq)
       vG0rt_Aq(iAq) = vSpc(iSpc)%G0rt
    end do

  end subroutine BOX_Compute_RefStatePotential_Aq

  !---
  
  subroutine BOX_Compute_RefStatePotential_Mk
    !==========================================================
    ! Purpose   : compute G0RtMk for Kinetic Mineral Species
    !==========================================================
    use M_Global_Vars, only : vSpc
    use M_Box_System_Vars
    implicit none
    integer :: iMk, iSpc
    
    do iMk = 1, nMk
       iSpc = vOrdMk(iMk)
       vG0rt_Mk(iMk)  = vSpc(iSpc)%G0rt
    end do
    
  end subroutine BOX_Compute_RefStatePotential_Mk

  !---
  
  subroutine BOX_Compute_LogK
    !===========================================================
    ! Purpose   : compute logK for As and Mk secondary species
    !             vlogKAs(iAs) and vlogKMk(iMk) 
    !===========================================================
    implicit none
    call Info_("BOX_Compute_LogK")
    call BOX_Compute_LogK_As
    call BOX_Compute_LogK_Mk

    call BOX_Info_RefState
    
  end subroutine BOX_Compute_LogK


  !---
  
  subroutine BOX_Compute_LogK_As
    !=======================================================
    ! Purpose   : compute LogK of secondary aqueous species
    !=======================================================
    use M_Box_System_Vars
    implicit none
    real(dp) :: DG0rt
    integer :: iAs, iAq, iPr
    !---
    
    do iAs = 1,nAs
       !-- reference state potential contribution
       DG0rt = Zero
       do iPr= 1,nPr
          iAq = iPr
          DG0rt = DG0rt + tNuAs(iAs,iPr) * vG0rt_Aq(iAq)
       end do
       
       iAq = nPr + iAs
       DG0rt = DG0rt - One*vG0rt_Aq(iAq)

       vLnK_As(iAs) = -DG0rt

    end do

  end subroutine BOX_Compute_LogK_As

  !---
  
  subroutine BOX_Compute_LogK_Mk
    !==================================================
    ! Purpose   : compute LogK of mineral species
    !==================================================
    use M_Box_System_Vars
    implicit none
    real(dp) :: DG0rt
    integer :: iAq, iMk, iPr
    !---
    do iMk = 1, nMk
       DG0rt = Zero
       do iPr= 1,nPr
          iAq = iPr
          DG0rt = DG0rt + tNuMk(iMk,iPr) * vG0rt_Aq(iAq)
       end do
       
       DG0rt = DG0rt - One * vG0rt_Mk(iMk)

       vLnK_Mk(iMk) = -DG0rt

    end do

  end subroutine BOX_Compute_LogK_Mk

  !---

  subroutine BOX_Info_RefState
    !==================================================
    ! Purpose   : Info Debug RefState
    !==================================================
    use M_Global_Vars,    only : vSpc
    use M_Numeric_Const,  only : Ln10
    use M_Box_Debug_Vars, only : LDebug_RefState
    use M_Box_System_Vars
    implicit none
    integer :: iAq, iSpc, iPr, iMk, iAs   
    
    call Info_("BOX_Info_RefState")
    !write(*,*) "LDebug_RefState=", LDebug_RefState
    if ( LDebug_RefState ) then
       
       write(*,'(2A6,1X,A20, 4A16)') &
       & ' Type ' ,' Index ', '  Name ',&
       & 'G0rt ', 'G0rt_Ln10 ',  'LnK' , 'Log10K'
       write(*,'(2A6,1X,A20, 4A16)') &
       & '------', '------', '--------------',&
       & '-------------', '-------------' , '----------', '------------'
       
       do iPr = 1,nPr
          iAq= iPr
          iSpc= vOrdAq(iAq)
          write(*,'(A6, I6,1X,A20, 4G16.6)') &
          & 'AQU Pr', iPr, vSpc(iSpc)%NamSp, &
          & vG0rt_Aq(iAq), vG0rt_Aq(iAq)/Ln10, Zero, Zero
       end do
       
       write(*,'(A)') "--"
       do iAs = 1,nAs
          iAq= nPr + iAs
          iSpc= vOrdAq(iAq)
          write(*,'(A6, I6,1X,A20, 4G16.6)') &
          & 'AQU As', iAs, vSpc(iSpc)%NamSp, &
          & vG0rt_Aq(iAq), vG0rt_Aq(iAq)/Ln10,  vLnK_As(iAs), vLnK_As(iAs)/Ln10
       end do
       
       write(*,'(A)') "--"
       do iMk = 1,nMk
          iSpc= vOrdMk(iMk)
          write(*,'(A6, I6,1X,A20, 4G16.6)') &
          & 'MIN Mk ', iMk, vSpc(iSpc)%NamSp, &
          & vG0rt_Mk(iMk), vG0rt_Mk(iMk)/Ln10, vLnK_Mk(iMk), vLnK_Mk(iMk)/Ln10
       end do
               
    end if

    
  end subroutine BOX_Info_RefState


end module M_Box_RefState
