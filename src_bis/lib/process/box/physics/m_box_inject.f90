module M_Box_Inject

  use M_Kinds
  use M_Trace
  
  use M_Box_Vars
    
  implicit none
  private

  public :: BOX_Compute_Inject

contains

  subroutine BOX_Compute_Inject
    
    use M_Box_Coupler_Vars
    implicit none
    
    call Info_('BOX_Compute_Inject')
    
    if (LCouplerInject) then
       call BOX_Compute_Inject_Coupler
    else
       call BOX_Compute_Inject_System
    end if

  end subroutine BOX_Compute_Inject
    
  !---

  subroutine BOX_Compute_Inject_Coupler
    
    use M_Global_Vars,   only : vEle, vSpc
    use M_Box_System_Vars
    use M_Box_Debug_Vars, only : LDebug_Inject 
    implicit none
    integer :: iEle, iCp

    call Info_('BOX_Compute_Inject_Coupler')
    
    !// Qinj is given directly 
     if (LDebug_Inject) then
        write(*,'(A6, 1X, 3A16, 1A12)') ' iCp ' , ' vCpn(iCp)%Name ', ' Ele_Name ',  &
            '    Statut    ', 'vQInj    '
       write(*,'(A6, 1X, 3A16, 1A12)') ' ----' , ' ---------------', ' ---------',  &
            ' -------------', ' --------'
       do iCp = 1, nCp
          iEle = vCpnBox(iCp)%iEle
          write(*,'(I6, 1X, 2A16,1X, A15, G12.4)') iCp, vCpnBox(iCp)%NamCp,  vEle(iEle)%NamEl,&
               vCpnbox(iCp)%Statut, vQInj(iCp)
       end do
    end if
  end subroutine BOX_Compute_Inject_Coupler
  
  !---

  subroutine BOX_Compute_Inject_System
    
    use M_Global_Vars,   only : vEle, vSpc
    use M_Box_System_Vars
    use M_Box_Debug_Vars, only : LDebug_Inject 
    use M_Box_Thermo_Vars
    implicit none

    real(dp) :: MassF
    integer :: iAq, iEle, iCp

    call Info_('BOX_Compute_Inject_System')
    
    ! Compute Injection Rate in Term of Component 
    vCInj = Zero
    do iCp = 1, nCp
       do iAq = 1, nAq
          vCInj(iCp) = vCInj(iCp) + tAlfAq(iCp,iAq)* vMolFInj(iAq)
       end do
    end do
    
    !// Debug Check 
    MassF = Zero
    do iCp = 1, nCp
       iEle = vCpnBox(iCp)%iEle
       MassF = MassF + vCInj(iCp) * vEle(iEle)%WeitKg 
    end do
    
    !// Multiply the Concentration by Finj
    vQInj = ( vCInj/MassF * RhoF ) * Finj 

    if (LDebug_Inject) then
       write(*,*) "Mass Concentration Injection ( kg.m-3 ) =", MassF       
       write(*,'(A6, 1X, 3A16, 2A12)') ' iCp ' , ' vCpn(iCp)%Name ', ' Ele_Name ',  &
            '    Statut    ', '   vCInj  ' , 'vQInj    '
       write(*,'(A6, 1X, 3A16, 2A12)') ' ----' , ' ---------------', ' ---------',  &
            ' -------------', ' --------',  ' -------------'
       do iCp = 1, nCp
          iEle = vCpnBox(iCp)%iEle
          write(*,'(I6, 1X, 2A16,1X, A15, 2G12.4)') iCp, vCpnBox(iCp)%NamCp,  vEle(iEle)%NamEl,&
               vCpnbox(iCp)%Statut, vCInj(iCp), vQInj(iCp)
       end do
    end if
    
  end subroutine BOX_Compute_Inject_System
  

end module M_Box_Inject
