module M_Box_TP_Update

  !=====================================================
  ! BOX MANAGER FOR TP CONDITIONS UPDATING 
  !=====================================================


  use M_Kinds
  use M_Trace
  implicit none
  private

  real(dp) :: TolTdgK = 0.1D0
  real(dp) :: TolPbar = 0.1D0

  public :: BOX_TP_Update
  !--
  private :: Box_Global_TP_Update
  private :: Box_KinFas_DG_UpDate
  private :: Box_KinModel_T_Update
  private :: TP_Changed

contains 

  !---

  subroutine BOX_TP_Update
    use M_Info_Value
    use M_Box_Vars, only : TdgK, tDgKBox, Pbar, PbarBox
    implicit none
    
    if ( TP_Changed( TdgK, TdgKBox, Pbar, PbarBox ) ) then
       call Box_TP_Info( HasChanged = .true. )        
       TdgK = TdgKBox
       Pbar = PbarBox
    else
       call Box_TP_Info( HasChanged = .false. )        
    end if

    call Box_Global_TP_Update
    !
    call Box_KinFas_DG_UpDate
    !
    call Box_KinModel_T_Update
    !
  end subroutine BOX_TP_Update

  !---

  logical function TP_Changed(T1, T2, P1, P2)

    implicit none
    real(dp) :: T1, P1
    real(dp) :: T2, P2
    !-----
    TP_Changed = ( ABS(T1 - T2) > TolTdgK) .or. ( ABS(P1 - P2) > TolPbar )

  end function TP_Changed
  
  !---

  subroutine BOX_Global_TP_Update
    !=================================================================
    ! Update Static TP Parameters for KinFas for Global System  
    !=================================================================
    use M_Box_Vars, only : TdgK, Pbar
    use M_Global_Tools,only: Global_TP_Update
    use M_Global_Vars, only: vSpcDtb, vSpc, vMixModel, vDiscretModel, vDiscretParam, vMixFas, vFas
    implicit none
    !
    integer:: nCp
    call Global_TP_Update( &
         & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
         & vSpc,vMixModel,vMixFas,vFas)                     !inout

  end subroutine BOX_Global_TP_Update

  !---

  subroutine Box_KinFas_DG_UpDate
    !=================================================================
    ! Update Static TP Parameters for KinFas
    !=================================================================

    ! nothing to do

  end subroutine Box_KinFas_DG_UpDate

  !---

  subroutine Box_KinModel_T_Update
    !=================================================================
    ! Update Static T Parameters for KinModel
    !=================================================================
    use M_Global_Vars, only : vKinModel
    use M_Box_Vars, only : TdgK, TimeFactor
    use M_T_KinModel
    !--
    
    integer :: iModel 
    do iModel = 1, size(vKinModel)
        call KinModel_CalcCoeffs(vKinModel(iModel),TimeFactor,TdgK)
     end do
     
  end subroutine Box_KinModel_T_Update

  !--- 

  subroutine Box_TP_Info ( HasChanged )
    !=================================================================
    ! Update Static T Parameters for KinModel
    !=================================================================

    use M_Box_Debug_Vars
    use M_Box_Vars, only : TdgK, tDgKBox, Pbar, PbarBox
    implicit none
    logical :: HasChanged

    if (LDebug_TPUpdate) then

       if (HasChanged) then 
          write(*,*) '---------------------------------------'
          write(*,*) "Warning !!! TP State Has Changed !!"
          write(*,*) "TdgK = ", TdgK, TdgKBox,  TdgKBox - TdgK
          write(*,*) "Pbar = ", Pbar, PbarBox,  PbarBox - Pbar
          write(*,*) '---------------------------------------'

       else
          write(*,*) '---------------------------------------'
          write(*,*) "Cool !! TP State Has Not Changed !!"
          write(*,*) "TdgK = ", TdgK, TdgKBox,  TdgKBox - TdgK
          write(*,*) "Pbar = ", Pbar, PbarBox,  PbarBox - Pbar
          write(*,*) '---------------------------------------'
       end if
    end if

  end subroutine Box_TP_Info

end module M_Box_TP_Update
