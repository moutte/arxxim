 module M_Box_Gamma

  use M_Kinds
  use M_Trace

  use M_Box_Vars
  use M_Box_Thermo_Vars

  implicit none
  private

  public :: BOX_Compute_Gamma
  public :: BOX_Compute_WaterProperties
  public :: BOX_Compute_Activity
  
contains

   !---

  subroutine BOX_Compute_WaterProperties
    !==================================================
    ! Purpose   : compute Lai0, LnGam and dLnGam_dvLnX
    !==================================================

    use M_Numeric_Const, only : Ln10
    use M_Global_Vars, only : vSpc
    use M_Basis_Vars,  only : isW, MWSv, vOrdAq
    use M_T_Species
    implicit none
    integer :: iAq, iSpc, iAqH
    real(dp) :: Z
    
    call Info_('BOX_Compute_Gamma')

    Ionic = Zero
    do iAq=1, nAq
       iSpc = vOrdAq(iAq)
       Z = vSpc(iSpc)%Z
       Ionic = Ionic + exp(vLnMolal(iAq))*Z*Z
    end do
    Ionic = Ionic * 0.5
    
    iAqH = Species_Index("H+", vSpc(vOrdAq(1:nAq)))
    if (iAqH>0) pH = -vLnAct(iAqH)/Ln10

    call Box_Info_Gamma

  end subroutine BOX_Compute_WaterProperties


  !---

  subroutine BOX_Compute_Gamma
    !==================================================
    ! Purpose   : compute Lai0, LnGam and dLnGam_dvLnX
    !==================================================

    use M_Global_Vars, only : vSpc, Solmodel
    use M_Basis_Vars,  only : isW, MWSv, vOrdAq
    use M_Solmodel_Calc_System
    
    implicit none
    integer :: iAq
    real(dp) :: OsmoSv
    real(dp) :: MsW
    call Info_('BOX_Compute_Gamma')

    call Solmodel_Calc_Physics( &
         & TdgK,Pbar,        & !IN  ! temperature and pressure
         & Solmodel,         & !IN  ! solvent phase model
         & vSpc(vOrdAq(1:nAq)), & !IN  ! vector of Aqueous Species
         & isW,              & !IN  ! index of solvent
         & vMolF,            & !IN  ! mole numbers
                                !--------------------------
         & vLnGam,           & !OUT ! gamma coefficient
         & vLnAct,           & !OUT ! activity
         & OsmoSv )            !OUT ! osmotic coefficient 

    do iAq=1, nAq
       vLnMolal(iAq) = log(vMolF(iAq)) -  log(vMolF(isW)) - log(MWSv) 
    end do

    call Box_Info_Gamma

  end subroutine BOX_Compute_Gamma

  !---

  subroutine Box_Compute_Activity( &
       &  vXf ,  &      ! IN
       &  vLnXf, &      ! IN
       &  vLnai, &      ! OUT
       &  dLnai_dLnXf ) ! OUT
    
    !==================================================
    ! Purpose   : compute Lnai and its derivatives
    !==================================================
    
    use M_Trace   
    use M_Basis_Vars,  only : isW, MWSv
    implicit none
    
    real(dp), intent(in)  :: vXf(:)
    real(dp), intent(in)  :: vLnXf(:)
    real(dp), intent(out) :: vLnai(:)
    real(dp), intent(out) :: dLnai_dLnXf(:, :)
    
    !---
    integer :: iAq
    real(dp) :: XfSum, Lnxw, Lnmi
    logical, parameter :: Limplicit_ACTIVITY_WATER= .false.
    !---
    call Info_('BOX_Compute_Activity')
    !---
    vLnai = Zero
    dLnai_dLnXf = Zero    
    !---
    
    if (Limplicit_ACTIVITY_WATER) then    
       ! Implicit Water Activity Model
    XfSum = SUM(vXf(1:nAq)) 
    Lnxw = vLnXf(isW) - log(XfSum)
    vLnai(isW) = Lnxw + vLnGam(isW) 
       
    dLnai_dLnXf(isW, 1:nAq) = - One/XfSum 
    dLnai_dLnXf(isW, isW ) = dLnai_dLnXf(isW, isW) + One 
    
    else       
       ! Explicit Water Activity Model
       vLnai(isW) = vLnAct(isW)
       dLnai_dLnXf(isW, 1:nAq) = Zero
    
    end if

    !---
    ! Solute Species iAq neq isW
    do iAq = 1, nAq
       if (.not. ( iAq .EQ. isW )) then
          Lnmi = vLnXf(iAq) - vLnXf(isW) - log(MWSv)
          vLnai( iAq ) = Lnmi + vLnGam(iAq) 

          dLnai_dLnXf(iAq,iAq) = One
          dLnai_dLnXf(iAq,isW) = - One           
       end if
    end do

    call BOX_Info_Activity(vXf, vLnXf, vLnai)

  end subroutine Box_Compute_Activity

  !---

  subroutine BOX_Info_Gamma
    !==================================================
    ! Purpose   : Info Debug Gamma
    !==================================================
    use M_Global_Vars,    only : vSpc
    use M_Basis_Vars,     only : vOrdAq
    use M_Box_Debug_Vars, only : LDebug_Gamma
    use M_Box_RefState_Vars
    implicit none
    integer :: iAq, iSpc    
    
    call Info_("BOX_Info_Gamma")

    if ( LDebug_Gamma ) then
       write(*,'(A6,1X,A20, 7A16)') &
            &  '  iAq ', 'iAqSpc%Name ', 'vMolF(iAq)', 'vMolal(iAq)', 'vG0rt(iAq)', &
            & 'vLnMolal(iAq)' ,'vLnGam(iAq)', 'vLnAct(iAq)', 'vAct(iAq)'

       write(*,'(A6,1X,A20, 7A16)') &
            & '------', '------------', '-----------', '-----------', '-----------', &
            & '-----------', '-----------', '-----------','-----------'

       do iAq = 1,nAq
          iSpc = vOrdAq(iAq)
          write(*,'(I6,1X,A20, 7G16.6)') &
            &  iAq, vSpc(iSpc)%NamSp, vMolF(iAq), exp(vLnMolal(iAq)), vG0rt_Aq(iAq), &
            &  vLnMolal(iAq), vLnGam(iAq), vLnAct(iAq), exp(vLnAct(iAq))
       end do
    end if

  end subroutine BOX_Info_Gamma

  !---

  subroutine BOX_Info_Activity(vXf, vLnXf, vLnai)
    !==================================================
    ! Purpose   : Info Debug Gamma
    !==================================================
    use M_Global_Vars,    only : vSpc
    use M_Box_System_Vars
    use M_Box_Debug_Vars, only : LDebug_Gamma
    
    implicit none
    real(dp), intent(in) :: vXf(:), vLnXf(:), vLnai(:)

    integer :: iAq, iSpc    

    if ( LDebug_Gamma ) then
       write(*,'(A6,1X,A20, 5A16)') &
            '  iAq ', 'iAqSpc%Name ', 'vXf(iAq)', 'vLnXf(iAq)', 'vLnGam(iAq)', 'vLnai(iAq)', 'vai(iAq)'

       write(*,'(A6,1X,A20, 5A16)') &
            '------', '------------', '-----------', '-----------', '-----------', '-----------', '-----------'

       do iAq = 1,nAq
          iSpc = vOrdAq(iAq)
          write(*,'(I6,1X,A20, 5G16.6)') &
               iAq, vSpc(iSpc)%NamSp, vXf(iAq), vLnXf(iAq), vLnGam(iAq), vLnai(iAq), exp(vLnai(iAq))
       end do
    end if

  end subroutine BOX_Info_Activity

end module M_Box_Gamma
