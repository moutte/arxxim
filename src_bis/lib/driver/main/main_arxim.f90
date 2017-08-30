subroutine main_arxim
!lt !DEC$ ATTRIBUTES DLLEXPORT :: main_arxim
!---------------
! Arxim Main
!---------------

!// Interactive Session
call Driver_Arxim

!// Static + Dynamic Compute Tmax = 0
! call preKINXIM

!// Static + Dynamic Compute Tmax = Tfinal
! call proKINXIM

end subroutine main_arxim

