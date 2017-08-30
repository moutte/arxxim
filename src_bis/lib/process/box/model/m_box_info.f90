!! *******************************************************************
!! * File name:    M_Box_Info
!! * Purpose:      Show System State Information
!! * Author:       Anthony Michel (anthony.michel@ifp.fr)
!! * Created:      2008
!! ********************************************************************

module M_Box_Info

  use M_Kinds
  use M_Trace
  !
  implicit none
  private

  public :: Box_Info
  public :: Box_System_Info
  public :: Box_Info_TotF
  public :: Box_Info_Equation
  public :: Box_Info_Variable
  public :: Box_Info_Amount

contains

   subroutine Box_Info_TotF(vTot, legend, file)
    !==================================================
    ! Purpose   : Print Total Fluid Composition
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2009
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,      only : vSpc
    use M_Box_System_Vars,  only : vCpnBox
    !--
    implicit none
    real(dp), intent(in) :: vTot(:)
    integer, intent(in), optional :: file
    character(len=*), intent(in) :: legend 
    !--
    integer :: iCp
    character(len=20) :: nameVar
    integer :: f
    !--
    f= 6
    if (present(file)) f= file 
    !--
    write(f,'(3A)') '-------- ', trim(legend) ,' -------- '
    !--
    do iCp= 1, nCp 
       nameVar= vCpnBox(iCp)% NamCp
       write(f,'(A, I6, A10, A20, E16.6)') &
       & 'vTot ', iCp, ' Component ', nameVar, vTot(iCp)
    end do    

  end subroutine Box_Info_TotF

  !---

  subroutine Box_Info_Equation(vR, legend, file)
    !==================================================
    ! Purpose   : Print Equation Residual Value Information 
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2009
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,      only : vSpc
    use M_Box_System_Vars,  only : vOrdAq, vOrdMk, vCpnBox
    !--
    implicit none
    real(dp), intent(in) :: vR(:)
    integer, intent(in), optional :: file
    character(len=*), intent(in)              :: legend 
    !--
    integer :: iAs, iMk, iCp, iEq, iSpc
    character(len=20) :: nameVar
    integer :: f
    !--
    f= 6
    if (present(file)) f= file 
    !--
    write(f,'(3A)') '-------- ', trim(legend) ,' -------- '
    !--
    iEq= 0
    do iCp= 1, nCp 
       iEq= iEq + 1
       nameVar= vCpnBox(iCp)% NamCp
       write(f,'(A, I6, A10, A20, E16.8)') &
       & 'Eq ', iEq, 'BALANCE ', nameVar, vR(iEq)
    end do
    !--
    write(f,'(A)') '-- '
    do iAs= 1, nAs 
       iEq= iEq + 1
       iSpc= vOrdAq(nPr+iAs)
       nameVar= vSpc(iSpc)% NamSp
       write(f,'(A, I6, A10, A20, E16.8)') &
       & 'Eq ',  iEq,'EQUIL    ', nameVar, vR(iEq)
    end do
    !--
    write(f,'(A)') '-- '
    do iMk= 1, nMk 
      iEq= iEq + 1
      iSpc= vOrdMk(iMk)
      nameVar= vSpc(iSpc)% NamSp
      write(f,'(A, I6, A10, A20, E16.8)') &
      & 'Eq ', iEq, 'KINETIC  ', nameVar, vR(iEq)
    end do

  end subroutine Box_Info_Equation

  !---

  subroutine Box_Info_Variable(vX, legend, file)
    !==================================================
    ! Purpose   : Print Variable Value Information 
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2009
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,     only : vSpc
    use M_Box_System_Vars, only : vOrdAq, vOrdMk

    !--
    implicit none
    real(dp), intent(in) :: vX(:)
    integer, intent(in), optional :: file
    character(len=*), intent(in)              :: legend 
    !--
    integer :: iAq, iMk, iSpc, iVar
    character(len=20) :: nameVar
    integer :: f
    !--
    f= 6
    if (present(file)) f= file 
    !--
    write(f,'(3A)') '-------- ', trim(legend) ,' -------- '
    !--
    iVar= 0
    do iAq= 1, nAq 
       iVar= iVar + 1
       iSpc= vOrdAq(iAq)
       nameVar= vSpc(iSpc)%NamSp
       write(f,'(A, I6, 1X, A20, E16.8)') 'Var ', iVar, nameVar, vX(iVar)
    end do
    !--
    write(f,'(A)') '-- '
    do iMk= 1, nMk 
      iVar= iVar + 1
      iSpc= vOrdMk(iMk)
      nameVar= vSpc(iSpc)%NamSp
      write(f,'(A, I6, 1X, A20, E16.8)') 'Var ', iVar, nameVar, vX(iVar)
    end do

  end subroutine Box_Info_Variable

  !---

  subroutine Box_Info_Amount(vXf, vXMk , legend, file)
    !==================================================
    ! Purpose   : Print System Information in a File
    !
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2009
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,     only : vSpc
    use M_Box_System_Vars, only : vOrdAq, vOrdMk
    !--
    implicit none
    real(dp), intent(in) :: vXf(:), vXMk(:)
    integer, intent(in), optional :: file
    character(len=*), intent(in)  :: legend 
    !--
    integer :: iAq, iMk, iSpc
    !--
    character(len=20) :: nameVar
    integer :: f
    !--
    f= 6
    if (present(file)) f= file 
    !--
    write(f,'(3A)') '-------- ', trim(legend) ,' -------- '
    !--
    do iAq= 1, nAq 
       iSpc= vOrdAq(iAq)
       nameVar= vSpc(iSpc)%NamSp
       write(f,'(A, I6, A20, E16.8)') 'AQ ', iAq, ' '//nameVar, vXf(iAq)
    end do
    !--
    write(f,'(A)') '-- '
    do iMk= 1, nMk 
       iSpc= vOrdMk(iMk)
       nameVar= vSpc(iSpc)% NamSp
       write(f,'(A, I6, A20, F16.8)') 'MK ', iMk, ' '//nameVar, vXMk(iMk)
    end do

  end subroutine Box_Info_Amount

  !---

  subroutine Box_System_Info(file)
    !==================================================
    ! Purpose   : Print System Information in a File
    !
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2009
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,     only : vSpc, vEle
    use M_Box_System_Vars
    use M_Box_Debug_Vars, only : LDebug_System
    !--
    implicit none
    integer, intent(in), optional :: file
    !--
    integer :: iPr, iAs, iMk, iAq, iSpc
    integer :: iVar, nVar, nEle
    integer, allocatable :: iAqSpc(:)
    integer, allocatable :: iPrSpc(:)
    integer, allocatable :: iAsSpc(:)
    integer, allocatable :: iMkSpc(:)
    character(len=23) :: nameVar, typeVar
    character(len=23), allocatable :: vNameVar(:)
    integer :: f
    character(len=6) :: I1, I2, I3, I4
    !--
    nEle= size(vEle)

    f= 6
    if (present(file)) f= file 

    call Info_("Box_System_Info")  
    
    if (LDebug_System) then
       nSpc= size(vSpc)

       allocate(iAqSpc(nSpc))
       allocate(iPrSpc(nSpc))
       allocate(iAsSpc(nSpc))
       allocate(iMkSpc(nSpc))

       nVar= nPr + nAs + nMk
       allocate(vNameVar(nVar))

       !--

       iAqSpc= 0
       iPrSpc= 0
       iAsSpc= 0
       iMkSpc= 0
       do iVar= 1, nVar
          vNameVar(iVar)= '-----'
       end do

       !--

       do iAq= 1, nAq 
          iSpc= vOrdAq(iAq)
          iAqSpc ( iSpc )= iAq 
       end do

       do iPr= 1, nPr 
          iSpc= vOrdAq(iPr)
          iPrSpc ( iSpc )= iPr 
       end do

       do iAs= 1, nAs 
          iSpc= vOrdAq(nPr+iAs)
          iAsSpc ( iSpc )= iAs 
       end do

       do iMk= 1, nMk 
          iSpc= vOrdMk(iMk)
          iMkSpc ( iSpc )= iMk 
       end do

       write(f,'(A)') &
       & "----------------------------------------------------------------------------"
       write(f,'(A6, 1X,A20,A6,A6,A6,A6)') &
       & '  iSpc', &
       & '     Species          ', &
       & '   iAq', &
       & '   iPr', &
       & '   iAs', &
       & '   iMk'
       write(f,'(A6, 1X,A20,A6,A6,A6,A6)') &
       & ' -----', ' ---------------------',&
       & ' -----', ' -----', ' -----', ' -----'
       do iSpc= 1, nSpc 
          I1= Int2Str(iAqSpc(iSpc))
          I2= Int2Str(iPrSpc(iSpc))
          I3= Int2Str(iAsSpc(iSpc))
          I4= Int2Str(iMkSpc(iSpc))
          write(f,'(I6, 1X, A20, A6,A6,A6, A6)') &
          & iSpc, vSpc(iSpc)%NamSp, I1, I2, I3, I4          
       end do

       write(f,'(A)') ''

       write(f,'(A)') &
       & "----------------------------------------------------------------------------"
       write(f,'(A6, 1X, A20,A20,A6,A6)') &
       & '  iVar',&
       & ' Variable             ', &
       & '    TypeVar     ', &
       & ' iType', &
       & '  iSpc'
       write(f,'(A6, 1X, A20,A20,A6,A6)') &
       & ' -----',&
       & '----------------------',&
       & ' ---------------', &
       & ' -----', &
       & ' -----'
       iVar= 0
       typeVar= '  AQU PRIMARY'
       do iPr= 1, nPr 
          iVar= iVar + 1
          iSpc= vOrdAq(iPr)
          nameVar= vSpc(iSpc)% NamSp
          vNameVar(iVar)= nameVar
          write(f,'(I6, 1X, A20, A20, I6, I6)') &
          & iVar, nameVar, typeVar, iPr, iSpc
       end do

       write(f,'(A)') "-"
       typeVar= '  AQU SECONDARY'
       do iAs= 1, nAs 
          iVar= iVar + 1
          iSpc= vOrdAq(nPr+iAs)
          nameVar= vSpc(iSpc)%NamSp
          vNameVar(iVar)= nameVar
          write(f,'(I6, 1X, A20, A20, I6, I6)') &
          & iVar, nameVar, typeVar, iAs, iSpc
       end do

       write(f,'(A)') "-"
       typeVar= '  MIN KINETIC'
       do iMk= 1, nMk 
          iVar= iVar + 1
          iSpc= vOrdMk(iMk)
          nameVar= vSpc(iSpc)%NamSp
          vNameVar(iVar)= nameVar
          write(f,'(I6, 1X, A20, A20, I6, I6)') &
          & iVar, nameVar, typeVar, iMk, iSpc
       end do

       write(f,'(A)') " "
       write(f,'(A)') "------------------------------------------------------------------------------"
       call Print_Compo(transpose(tAlfPr), nPr, nEle, vNameVar(1:nPr),                 vEle(1:nEle)%NamEl, 'Formula Matrix tAlfPr' )
       call Print_Compo(transpose(tAlfAs), nAs, nEle, vNameVar(nPr+1:nPr+nAs),         vEle(1:nEle)%NamEl, 'Formula Matrix tAlfAs' )
       call Print_Compo(transpose(tAlfMk), nMk, nEle, vNameVar(nPr+nAs+1:nPr+nAs+nMk), vEle(1:nEle)%NamEl, 'Formula Matrix tAlfMk' )
       call Print_Compo(transpose(tAlfAq), nAq, nEle, vNameVar(1:nAq),                 vEle(1:nEle)%NamEl, 'Formula Matrix tAlfAq' )

       write(f,'(A)') " "
       write(f,'(A)') "------------------------------------------------------------------------------"
       call Print_Compo(tNuAs, nAs, nPr, vNameVar(nPr+1:nPr+nAs),         vNameVar(1:nPr), 'Stokio Matrix tNuAs' )
       call Print_Compo(tNuMk, nMk, nPr, vNameVar(nPr+nAs+1:nPr+nAs+nMk), vNameVar(1:nPr), 'Stokio Matrix tNuMk' )
       write(f,'(A)') " "
       write(f,'(A)') "------------------------------------------------------------------------------"
       call Info_("End Box_System_Info")  

       deallocate(iAqSpc)
       deallocate(iPrSpc)
       deallocate(iAsSpc)
       deallocate(iMkSpc)
       deallocate(vNameVar)

    end if

  end subroutine Box_System_Info

  !---

  subroutine Box_Info(file, legend)
    !==================================================
    ! Purpose   : Print Readeable Information about  
    !             the System State
    !
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_Box_Vars
    use M_Global_Vars,     only : vSpc, vKinFas
    use M_Box_System_Vars    
    use M_Box_Thermo_Vars
    !
    implicit none
    integer, intent(in), optional :: file
    character(len=*), intent(in), optional :: legend
    !--
    real(dp):: MolBox,    MolPhase
    real(dp):: RhoBox,    RhoPhase
    real(dp):: MassBox,   MassPhase
    real(dp):: VolBox,    VolPhase
    integer :: iAq, iMk
    integer :: f
    !---
    f= 6
    if (present(file)) f= file 

    call Info_("Box_Info")    

    write(f,"(A)")
    if (present(legend)) then
       write(f,"(A)") "----------------------------------------"
       write(f,"(A)") trim(Legend)
       write(f,"(A)") "----------------------------------------"
    end if
    write(f,"(A)")
    !// TP Conditions 
    write(f,"(2A16)") "--------------", "--------------"
    write(f,'(A16,E16.6)') "Date[s]",      Time
    write(f,'(A16,E16.6)') "Temperature[K]",      TdgK
    write(f,'(A16,E16.6)') "Pressure[bar]",       Pbar
    write(f,'(A16,E16.6)') "Volume[m3]",          VBox
    write(f,"(2A16)") "--------------", "--------------"

    !// Box
    MolBox    = SUM(vMolF(1:nAq)) + SUM(vMolK(1:nMk)) 
    MassBox   = SUM(vMolF(1:nAq)*vSpc(vOrdAq(1:nAq))%WeitKg) + SUM(vMolK(1:nMk)*vSpc(vOrdMk(1:nMk))%WeitKg)
    VolBox    = SUM(vMolF(1:nAq)*vSpc(vOrdAq(1:nAq))%WeitKg/RhoF) + SUM(vMolK(1:nMk)*vSpc(vOrdMk(1:nMk))%WeitKg/vRhoK(1:nMk))
    RhoBox    = MassBox/VolBox

    call Print_Macro_Header(f)

    !//-------------
    call Print_Macro_Line(f)
    call Print_Macro( f, "Model", "ARXIM", &
         "System " , 1, "BOX" , &
         MolBox, MolBox, VolBox, RhoBox  )

    !//-------------
    call Print_Macro_Line(f)

    !// Aqueous Phase
    MolPhase = SUM(vMolF(1:nAq))    
    MassPhase= SUM(vMolF(1:nAq)*vSpc(vOrdAq(1:nAq))%WeitKg)
    VolPhase = SUM(vMolF(1:nAq)*vSpc(vOrdAq(1:nAq))%WeitKg/RhoF)
    RhoPhase = MassPhase/VolPhase

    call Print_Macro( f, "System", "BOX", &
         "Phase " , 1, "AQUEOUS" , &
         MolPhase, MolBox, VolPhase, RhoPhase  )

    !// Mineral Phases 
    do iMk= 1, nMk
       MolPhase= vMolK(iMk)
       MassPhase= vMolK(iMk)*vSpc(vOrdMk(iMk))%WeitKg
       VolPhase = vMolK(iMk)*vSpc(vOrdMk(iMk))%WeitKg/vRhoK(iMk)
       RhoPhase = MassPhase/VolPhase       

       call Print_Macro( f, "System", "BOX", &
            "Phase" , 1+iMk, vKinFas(iMk)%NamKF , &
            MolPhase, MolBox, VolPhase, RhoPhase )
    end do

    !//-------------
    call Print_Macro_Line(f)

    !// Aqueous Phase
    MolPhase= SUM(vMolF(1:nAq))

    do iAq= 1, nAq
       call Print_Macro( f, "Phase", "AQUEOUS", &
            "Species " , iAq, vSpc(vOrdAq(iAq))%NamSp , &
            vMolF(iAq), MolPhase)
    end do

    !// Mineral Phases 
    do iMk= 1, nMk
       call Print_Macro_Line(f)
       MolPhase= vMolK(iMk)

       call Print_Macro( f, "Phase", vKinFas(iMk)%NamKF, &
            "Species",  1, vSpc(vOrdMk(iMk))%NamSp , &
            vMolK(iMk), MolPhase)
    end do

    call Print_Macro_Line(f)

  end subroutine Box_Info

  !// Utils

  subroutine Print_Macro_Header(f)
    implicit none
    integer :: f
    !---
    write(f,'(3A16,A6, 5A16)') &
    & "ParentType", "ParentName", "Type" , "Index", "Name", &
    & "n[Mol]", "x[%Mol]", "Volume[m3]", "Rho[kg.m-3]"

  end subroutine Print_Macro_Header

  !--

  subroutine Print_Macro_Line(f)
    implicit none
    integer :: f
    !---
    write(f,'(3A16,A6, 5A16)') &
    & "--------------", "--------------", "--------------" , "-----", &
    & "--------------", &
    & "--------------", "--------------", "--------------",  "--------------"
  end subroutine Print_Macro_Line

  !--

  subroutine Print_Macro(f, ParentTypeName, ParentName, TypeName, Index, Name, Value, Total, Vol, Rho)
    implicit none
    integer,          intent(in) :: f
    character(len=*), intent(in) :: TypeName, Name, ParentName, ParentTypeName
    real(dp),         intent(in) :: Value, Total
    real(dp),         intent(in), optional  :: Vol, Rho
    integer,          intent(in) :: Index
    !---
    if (present(Rho)) then
       write(f,'(3A16,I6, A16, 4E16.8)') &
       & trim(ParentTypeName), trim(ParentName), trim(TypeName), &
       & Index, trim(Name), Value, Value/Total, Vol, Rho
    else
       write(f,'(3A16,I6, A16, 2E16.8)') &
       & trim(ParentTypeName), trim(ParentName), trim(TypeName), &
       & Index, trim(Name), Value, Value/Total
    end if
  end subroutine Print_Macro

  !---

  subroutine Print_Compo(tCompo, nbSpc, nbEle, nameSpc, nameEle, legend, file)
    !============================================================
    ! purpose : sparse composition matrix informations
    !===========================================================
    implicit none

    real(dp),         intent(in)              :: tCompo(:,:)
    integer,          intent(in)              :: nbSpc, nbEle
    character(len=*), intent(in)              :: nameSpc(:), nameEle(:)

    character(len=*), intent(in)              :: legend 
    integer,          intent(in), optional    :: file
    !--- local variables
    integer :: f
    integer :: iEle, iSpc
    real(dp):: compo
    !
    f= 6
    if (present(file)) f= file

    write(f,'(3A)') '-------- ', trim(legend) ,' -------- '

    do iSpc= 1, nbSpc   

       do iEle= 1, nbEle
          compo= tCompo(iSpc,iEle)
          if (.not.(compo==0.D0)) &
          & write(f,'(A20,A20,G14.4)') nameSpc(iSpc), nameEle(iEle), compo
       end do
       write(f,'(A)') ' '
    end do

  end subroutine Print_Compo

  !--

  function Int2Str(I) result (S)
    implicit none
    character(len=6) :: S
    character(len=6) :: Str
    integer, intent(in) :: I
    !---
    if (I.eq.0) then 
       Str= '     -'
    else
       write(Str,'(I6)') I
    end if
    S= Str
  end function INT2STR

end module M_Box_Info

