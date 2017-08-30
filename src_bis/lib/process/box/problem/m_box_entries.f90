module M_Box_Entries

  !====================================================================
  ! Purpose : apply operations on entries
  !--------------------------------------------------------------------
  !   - Recover variables form entries
  !   - Transform Jacobian depending from entry choice
  !====================================================================

  use M_Kinds
  use M_Box_Residual_Vars
  use M_Box_Debug_Vars
  use M_Box_Solver_Vars
  use M_Box_Utils
  implicit none
  private

  public :: BOX_Compute_Variables
  public :: BOX_Check_Variables
  public :: Box_Jacobian_LogTransform
  public :: BOX_Compute_Entries
  public :: BOX_Recover_Variables

contains

  subroutine BOX_Compute_Entries(vX)
    !======================================================
    ! purpose : Compute entries from natural variables
    !======================================================
    implicit none
    real(dp), intent(out) :: vX(:)
    integer :: iAq, iMk
    !--
    call Info_("BOX_Compute_Entries")

    !//Debug
    if (LDebug_Entries) then
       if (LogForFluid)     call Warning_("LogForFluid")
       if(.not.(LogForFluid)) call Warning_("Not LogForFluid")

       if(LogForMin) call Warning_("LogForMin")
       if(.not.(LogForMin))  call Warning_("Not LogForMin")
    end if

    !// Compute
    if(LogForFluid) then
       ! entry = log (aqueous species amount)
       vX(1:nAq) = vLnXf(1:nAq)
    else
       ! entry = aqueous species amount
       vX(1:nAq) = Limited_Exp(vLnXf(1:nAq))
    end if

    if(LogForMin) then
       ! entry = log (mineral species amount)
       vX(nAq+1:nAq+nMk) = Limited_Log(vXm(1:nMk))
    else
       ! entry = mineral species amount
       vX(nAq+1:nAq+nMk) = vXm(1:nMk)
    end if

    !// Debug Results
    if (LDebug_Entries) then
       write(*,'(A6,1X,A20, A16)') '  iAq ', 'vXf(iAq)', 'vX(iAq)'
       write(*,'(A6,1X,A20, A16)') '------', '------------', '-----------'
       do iAq = 1,nAq
          write(*,'(I6,1X,G20.10, G16.6)') iAq, vXf(iAq), vX(iAq)
       end do

       write(*,*)
       write(*,'(A6,1X,A20, A16)') '  iMk ', 'vXm(iMk)', 'vX(nAq+iMk)'
       write(*,'(A6,1X,A20, A16)') '------', '------------', '-----------'
       do iMk = 1,nMk
          write(*,'(I6,1X,G20.10, G16.6)') iMk, vXm(iMk), vX(nAq+iMk)
       end do
    end if

  end subroutine BOX_Compute_Entries

  !---

  subroutine BOX_Check_Variables(vX, OKSolution)
    !======================================================
    ! purpose : Check that the solution is valid
    !======================================================
    use M_Box_Vars,        only : EpsMineral

    implicit none
    logical :: OKSolution
    real(dp), intent(inout) :: vX(:)
    integer :: iAq, iMk

    call BOX_Compute_Variables(vX)

    OKSolution = .true.
    do iAq = 1,nAq
       if ( vXf(iAq) <= Zero ) then
          OKSolution = .false.
       end if
    end do

    do iMk = 1,nMk
       if ( vXm(iMk) < Zero ) then
          OKSolution = .false.
       end if
    end do

  end subroutine BOX_Check_Variables

  !---

  subroutine BOX_Recover_Variables(vX)

    use M_Box_Vars,        only : vMolK, vMolF
    implicit none
    real(dp), intent(out) :: vX(:)

    if(LogForFluid) then
       ! entry = log (aqueous species amount)
       vX(1:nAq) = Limited_Log(vMolF(1:nAq))
    else
       ! entry = aqueous species amount
       vX(1:nAq) = vMolF(1:nAq)
    end if

    if(LogForMin) then
       ! entry = log (mineral species amount)
       vX(nAq+1:nAq+nMk) = Limited_Log(vMolK(1:nMk))
    else
       ! entry = mineral species amount
       vX(nAq+1:nAq+nMk) = vMolK(1:nMk)
    end if

  end subroutine BOX_Recover_Variables

  subroutine BOX_Compute_Variables(vX)
    !======================================================
    ! purpose : Recover natural variables from entries
    !======================================================
    implicit none
    real(dp), intent(in) :: vX(:)
    integer :: iAq, iMk
    !--
    call Info_("BOX_Compute_Variables")

    !//Debug
    if (LDebug_Variables) then
       if (LogForFluid)     call Warning_("LogForFluid")
       if(.not.(LogForFluid)) call Warning_("Not LogForFluid")

       if(LogForMin) call Warning_("LogForMin")
       if(.not.(LogForMin))  call Warning_("Not LogForMin")
    end if

    !// Compute
    if(LogForFluid) then
       ! entry = log (aqueous species amount)
       vLnXf(1:nAq) = vX(1:nAq)
       vXf(1:nAq) = Limited_Exp(vX(1:nAq))
    else
       ! entry = aqueous species amount
       vXf(1:nAq) = vX(1:nAq)
       vLnXf(1:nAq) = Limited_Log(vX(1:nAq))
    end if

    if(LogForMin) then
       ! entry = log (mineral species amount)
       vXm(1:nMk) = Limited_Exp(vX(nAq+1:nAq+nMk))
    else
       ! entry = mineral species amount
       vXm(1:nMk) = vX(nAq+1:nAq+nMk)
    end if

    !//Debug Results
    if (LDebug_Variables) then
       write(*,'(A6,1X,A20, 2A16)') '------', '------------', '-----------', '-----------'
       write(*,'(A6,1X,A20, 2A16)') '  iAq ', 'vX(iAq)', 'vXf(iAq)', 'vLnXf(iAq)'
       write(*,'(A6,1X,A20, 2A16)') '------', '------------', '-----------', '-----------'
       do iAq = 1,nAq
          write(*,'(I6,1X,G20.10, 2G16.6)') iAq, vX(iAq), vXf(iAq), vLnXf(iAq)
       end do
       write(*,*)
       write(*,'(A6,1X,A20, A16)') '  iMk ', 'vX(nAq+iMk)', 'vXm(iMk)'
       write(*,'(A6,1X,A20, A16)') '------', '------------', '-----------'
       do iMk = 1,nMk
          write(*,'(I6,1X,G20.10, G16.6)') iMk, vX(nAq+iMk), vXm(iMk)
       end do
    end if

  end subroutine BOX_Compute_Variables

  !---

  subroutine Box_Jacobian_LogTransform(tJac)
    !=================================================================
    ! purpose : Apply options LogForFluid and LogForMin
    !=================================================================

    implicit none
    real(dp), intent(inout) :: tJac(:,:)
    !---
    real(dp) :: dXm_dLnXm_iMk
    real(dp) :: dLnXf_dXf_iAq
    integer  :: iEq, iEntry, iAq, iMk

    !--

    if (.not.(LogForFluid)) then
       do iEq = 1, size(tJac,1)
          do iAq = 1,nAq
             iEntry = iAq
             dLnXf_dXf_iAq = One/vXf(iAq)
             tJac(iEq,iEntry) = tJac(iEq,iEntry)* dLnXf_dXf_iAq
          end do
       end do
    end if

    !---

    if (LogForMin) then
       do iEq = 1, size(tJac,1)
          do iMk = 1,nMk
             iEntry = nAq+iMk
             dXm_dLnXm_iMk = vXm(iMk)
             tJac(iEq,iEntry) = tJac(iEq,iEntry)* dXm_dLnXm_iMk
          end do
       end do
    end if

  end subroutine Box_Jacobian_LogTransform

end module M_Box_Entries
