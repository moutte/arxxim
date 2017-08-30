module M_Box_Test

  use M_Kinds
  use M_Trace
  implicit none
  public

contains

  !---

  subroutine Box_Start_Test
    !==================================================
    ! Purpose   : Start init
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_TimeMng
    use M_T_TimeLine
    use M_T_TimeStepMng
    implicit none

    type(T_TimeLine)    :: TimeLine
    type(T_TimeStepMng) :: TimeStepMng

    call Info_("Box_Start_Test")

    call TimeLine_Init(TimeLine, 0._dp, 10._dp)
    call TimeMng_Set_TimeLine(TimeLine)

    call TimeStepMng_Init(TimeStepMng, 1._dp, 0.01_dp, 5._dp)
    call TimeMng_Set_TimeStepMng(TimeStepMng)

    call TimeMng_Reset_TimeStep
    call TimeMng_Info

  end subroutine Box_Start_Test

  !----

  subroutine Box_Restart_Test
    !==================================================
    ! Purpose   : Restart init
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================

    use M_TimeMng
    use M_T_TimeLine
    use M_T_TimeStepMng
    implicit none

    type(T_TimeLine)    :: TimeLine
    type(T_TimeStepMng) :: TimeStepMng

    call Info_("Box_Restart")
    call TimeLine_Init(TimeLine, 10._dp, 20._dp)
    call TimeMng_Set_TimeLine(TimeLine)

    call TimeMng_Reset_TimeStep
    call TimeMng_Info

  end subroutine Box_Restart_Test

  !---

  subroutine Box_Calc_Init_Info
    !=============================================
    ! Info About The Problem Variables 
    !=============================================
    use M_Global_Vars, only : vSpc, vEle, vKinFas, vKinModel, vFas
    use M_Basis_Vars,  only : vOrdPr, vOrdAq, vOrdAs, vOrdMs, tAlfMs, tAlfAs,tNuAs,tAlfPr
    use M_System_Vars, only : vCpn
    use M_Box_Debug_Vars, only : LDebug_Calc
    
    implicit none
    integer :: nSpc, nAq, nPr, nAs, nMk, nCp
    integer :: N
    
    nSpc= size(vSpc)
    nAq=  size(vOrdAq)
    nPr=  size(vOrdPr)
    nAs=  size(vOrdAs)
    nMk = size(vKinFas)
    nCp = size(vCpn)
    
    !---
    if (LDebug_Calc) then
       call Info_("Box_Calc_Init_Check")

       N= size(vEle)
       call InfoList("Element",                   N,  vEle(:)%NamEl)

       call InfoList("Species",                   nSpc,  vSpc(:)%NamSp)
       call InfoList("Component",                 nCp,   vCpn(:)%NamCp)
       call InfoList("Aqueous-Species",           nAq,   vSpc(vOrdAq(:))%NamSp)
       call InfoList("Primary-Aqueous-Species",   nPr,   vSpc(vOrdPr(:))%NamSp)
       call InfoList("Secondary-Aqueous-Species", nAs,   vSpc(vOrdAs(:))%NamSp)

       N= size(vOrdMs)
       call InfoList("Mineral-Species",  N,  vSpc(vOrdMs(:))%NamSp)

       N= size(vKinModel)
       call InfoList("Mineral-Kinetics-Models",   N,   vKinModel(:)%Name )        

       call InfoList("Mineral-Kinetic-Phases",    nMk,  vKinFas(:)%NamKF )
       call InfoList("Mineral-Kinetic-Species",   nMk,  vSpc(vFas(vKinFas(1:nMk)%iFas)%iSpc)%NamSp)

       N= size(vFas)
       call InfoList("Phases",    N, vFas(:)%NamFs )

    end if

  end subroutine Box_Calc_Init_Info

  !---

  subroutine Box_Calc_System_Info
    use M_Box_System_Vars
    use M_Global_Vars, only : vSpc
    use M_Box_Debug_Vars, only : LDebug_Calc
    implicit none
    
    if (LDebug_Calc) then
       write(*,*) "================================"
       write(*,*) "nSpc =", size(vSpc)
       write(*,*) "================================"
       write(*,*) "nAq  =", nAq
       write(*,*) "nMk  =", nMk
       write(*,*) "================================"
       write(*,*) "nCp  =", nCp
       write(*,*) "nPr  =", nPr, shape(tAlfPr)    
       write(*,*) "nAs  =", nAs, shape(tAlfAs), shape(tNuAs)
       write(*,*) "nMk  =", nMk, shape(tAlfMk), shape(tNuMk)
       write(*,*) "================================"
       write(*,*) "shape(tAlfPr)  =", shape(tAlfPr)
       write(*,*) "shape(tAlfAs)  =", shape(tAlfAs)
       write(*,*) "shape(tAlfMk)  =", shape(tAlfMk)
       write(*,*) "================================"
       write(*,*) "shape(tNuAs)  =", shape(tNuAs)
       write(*,*) "shape(tNuMk)  =", shape(tNuMk)
       write(*,*) "================================"
    end if
    
  end subroutine Box_Calc_System_Info
       
  !---
  
  subroutine Box_Load_Check_Residual
    !=========================================================
    ! Box_Load_Check_Residual
    !=========================================================
    use M_Box_Info
    use M_Box_Utils
    use M_echomat
    use M_Box_Vars
    use M_Box_Debug_Vars, only : LDebug_Check

    use M_Box_Residual_Vars
    use M_Box_Residual
    use M_Box_Prepare
    use M_Box_Entries
    use M_Box_Jacobian
    use M_Box_Load_Vars
    
    implicit none
    !---
    call Info_("=====================================================")
    call Info_("Box_Load_Check_Residual")
    
    if (LDebug_Check) then
       call Box_Compute_Entries(vX)
       call Box_Prepare_TimeStep
       call Box_Prepare_GammaStep(vX)

       !// Physical vectors
       call Info_("----------------- Natural Vectors -------------------")
       call echovec_real(vMolF, nAq, 'vMolF','T')
       call echovec_real(vMolK, nMk, 'vMolK','T')
       call Box_Info_Amount(vMolF,vMolK,'Vectors vMolF, vMolK ')
       
       !// Compute residual   
       call Info_("----------------- Compute Residual ------------------")
       vR = Zero
       vR = Box_Residual(vX)
       
       !// Debug  
       call echovec_real(vX, size(vX), 'vX','T')
       call Box_Info_Variable(vX,'Variable vX')
       
       call echovec_real(vR, size(vX), 'vR','T')
       call Box_Info_Equation(vR,'Residual vR')
    
       !// Compute jacobian 
       call Info_("Compute Jacobian")
       tJac = Zero
       call Box_Jacobian(vX, tJac)

       !// Debug  
       call echomat_real(tJac, size(vX), size(vX), 'tJac')

    end if

    call Info_("Box_Load_Check_Residual")
    call Info_("=====================================================")
    

  end subroutine Box_Load_Check_Residual

  !---

  subroutine InfoList(Tag, n, vStrName)
    implicit none
    character(len=*), intent(in) :: Tag
    integer, intent(in) :: n
    character(len=*), intent(in) :: vStrName(:)
    !
    integer :: I
    !
    write(*,*) '-----------'
    call Info_("Info List : "// trim(Tag))
    write(*,'(A,I4)') "Size = ", n 
    do I = 1, n
       write(*,'(A,A,I4,A,A)') Tag, '[', I, '] = ', vStrName(I)
    end do

  end subroutine InfoList
end module M_Box_Test
