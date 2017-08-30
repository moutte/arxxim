module M_Box_Volume

  use M_Kinds
  use M_Trace
  use M_Box_Vars
  use M_Box_Thermo_Vars
  implicit none
  private

  public :: BOX_Compute_Volume

contains

  !---

  subroutine BOX_Compute_Volume
    
    !----------------------------------------------------
    ! Compute the volume of fluid in the box 
    !----------------------------------------------------
    
    use M_Global_Vars,   only : vSpc, vKinFas
    use M_Box_System_Vars,   only : vOrdAq, vOrdMk
    
    implicit none
    real(dp):: MassF
    integer :: iAq,iSpc,iMk
    
    call Info_('BOX_Compute_Volume')
    
    MassF= Zero
    do iAq = 1,nAq
       iSpc = vOrdAq(iAq)
       MassF = MassF + vMolF(iAq)*vSpc(iSpc)%WeitKg
    end do
    
    VolF= MassF / RhoF    
    PhiF= VolF / Vbox   

    do iMk = 1,nMk
       iSpc= vOrdMk(iMk)
       vPhiK(iMk)= vMolK(iMk)*vSpc(iSpc)%WeitKg / vRhoK(iMk) / Vbox  
    end do

    !// Debug Info
    call Box_Info_Volume

  end subroutine BOX_Compute_Volume
  
  !---

 subroutine BOX_Info_Volume
    !==================================================
    ! Purpose   : Info Debug Volume
    !==================================================
    use M_Global_Vars,   only : vSpc, vKinFas
    use M_Box_System_Vars,   only : vOrdAq, vOrdMk
    use M_Box_Debug_Vars, only : LDebug_Volume

    implicit none
    integer :: iMk    
    
    call Info_("BOX_Info_Volume")

    if ( LDebug_Volume ) then
       write(*,'(5A16)') '  VBox      ', 'VolF       ', 'PhiF       ',    'PhiS       ', 'PhiF + PhiS'
       write(*,'(5A16)') '------------', '-----------', '-----------',    '-----------',  '-----------'        
       write(*,'(5G16.6)') VBox, VolF, PhiF,  sum(vPhiK(1:nMk)) ,    PhiF+ sum(vPhiK(1:nMk))
       write(*,*)
              
       !// Mineral Phase Volumes
       write(*,'(A6,1X,A12, 2A14)') '  iMk ', ' iMkFas%Name', ' vMolK     ',    ' vPhiK     '
       write(*,'(A6,1X,A12, 2A14)') '------', '------------', '-----------','-----------'
       
       do iMk = 1, nMk
          write(*,'(I6,1X,A12,2G14.5)') iMk, vKinFas(iMk)%NamKF,  vMolK(iMk), vPhiK(iMk)
       end do
       
       !// Fluid Phase Volumes
       write(*,*)
       write(*,'(A6,1X,A12, 2A14)') '  iFk ', ' iFkFas%Name', ' vMolF     ',    ' vPhiF     '
       write(*,'(A6,1X,A12, 2A14)') '------', '------------', '-----------','-----------'
       write(*,'(I6,1X,A12,2G14.5)') 1 , 'AQUEOUS     ',  sum(vMolF(1:nAq)), PhiF
       
    end if
    
  end subroutine BOX_Info_Volume
  
end module M_Box_Volume
