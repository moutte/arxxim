module M_Box_System
  
  use M_Kinds
  use M_Trace
  
  use M_Box_System_Vars
  
  public ::  Box_System_Init
  public ::  Box_System_Build
  public ::  Box_System_Info
  public ::  Box_System_Reset_Basis

contains

  !---
  
!!$  subroutine Box_System_Set_vOrdAq
!!$    !==================================================
!!$    ! set vOrdAq by using input vOrdAq
!!$    ! use M_Basis_Vars with renaming options 
!!$    !==================================================
!!$    use M_Basis_Vars,        only: vOrdAqBasis => vOrdAq ! current in basis
!!$    use M_Box_System_Vars,   only: vOrdAqBox   => vOrdAq ! saved in box
!!$    !---
!!$    vOrdAqBox = vOrdAqBasis 
!!$    
!!$  end subroutine Box_System_Set_vOrdAq
  
  subroutine Box_System_Reset_Basis
    !===============================================
    ! Set vCpnBox to be the Current Component system
    !===============================================
    use M_System_Vars,  only: vCpn  
    use M_Basis,        only: Basis_Build
     
    implicit none
    logical :: LBuffer = .false.
    call Info_("Box_System_Init")
  !---
    !// Debug Info
    call Box_System_Info_Component
    
    !// Set vCpnBox Current Component system
    vCpn= vCpnbox
    call Basis_Build(LBuffer, vCpn)   

	!// Debug Info
    call Box_System_Info_Component 
  
  end subroutine Box_System_Reset_Basis


  !---
  
  subroutine Box_System_Init
    !===============================================
    ! Init the Chemical System and store properties 
    ! For the box problem
    !===============================================
    use M_System_Vars,  only: vCpn
    use M_System_Tools, only: System_Build
    use M_Basis,        only: Basis_Change
    use M_Equil 
    
    implicit none
    call Info_("Box_System_Init")
    !---
    !// Build Global Chemical System 
    call System_Build 

    !// Compute Speciation of the bloc System
    call Equil_Calc("SPC")       
    
    !// Set Components Inert
    call Basis_Change("DYN",vCpn)
    
    !// Allocate Box Sytstem Structure
    call Box_System_Build
    
  end subroutine Box_System_Init

  !---
  
  subroutine Box_System_Build
    !===============================================
    ! Build Box System from Basis System
    !===============================================
    use M_Box_Read
    use M_Global_Vars,  only: vSpc, vEle, vFas, vKinFas
    use M_Basis_Vars,   only: vOrdAqBasis => vOrdAq, vOrdAs, vOrdPr, tAlfSp, tNuSp
    use M_System_Vars,  only: vCpn
    use M_Box_Param_Vars
    
    use M_T_Phase
    use M_T_Species
    implicit none

    Type (T_Phase) :: MkFas
    integer :: iMk, iSpc
    logical :: Ok
    
    call Info_("Box_System_Build")
    
    !//============ System Size Parameters 
    
    nSpc= size(vSpc)

    !// Aqueous Species 
    nAq=  size(vOrdAqBasis)
       
    !// Components
    nCp= size(vCpn)  

    !// Primary and Secondary Aqueous Species 
    nPr= size(vOrdPr)
    if (.not.(nPr==nCp)) call Fatal_("Wrong number of Primary Species")
    nAs= nAq - nPr
          
    call Box_Read_Rock(Ok)

    !//========= Kinetic Mineral Phases
    !// Kinetic Mineral Species 
    nMk=  size(vKinFas) 
    
    !// Check kinetic phases are pure minerals
    
    do iMk = 1,nMk
      MkFas = vFas(vKinFas(iMk)%iFas)
      !! select case(MkFas%Typ)
      !! case("PURE") 
      !!   call Info_("Add A New Kinetic Mineral Species : "//vSpc(MkFas%iSpc)%NamSp)
      !! case default
      !!   call Fatal_("Box KinFas type not supported = "//vKinFas(iMk)%NamKF//MkFas%Typ)
      !! end select
      if(MkFas%iSpc>0) then
        call Info_( &
        & "Add A New Kinetic Mineral Species : "//vSpc(MkFas%iSpc)%NamSp)
      else
        call Fatal_( &
        & "Box KinFas type not supported = "//vKinFas(iMk)%NamKF//" is not pure")
      end if
    end do
    
    !// Allocate System
    call Box_System_Vars_New
    
    !//========== Component
    vOrdAq = vOrdAqBasis 
    vCpnBox = vCpn
    
    !//========== Aqueous Phase System
    iAqH2O= Species_Index("H2O",vSpc(vOrdAq(1:nAq)))
    if (iAqH2O.eq.Zero) call Fatal_("No Index for Solvent H2O")
       
    tAlfPr(1:nCp, 1:nPr) = tAlfSp(1:nCp, vOrdPr(:))
    tAlfAs(1:nCp, 1:nAs) = tAlfSp(1:nCp, vOrdAs(:))
    
    tAlfAq(1:nCp, 1:nPr ) = tAlfPr
    tAlfAq(1:nCp, nPr+1:nAq) = tAlfAs
    
    tNuAs(1:nAs, 1:nPr) = tNuSp(vOrdAs(:),:)
    
    !//========== Mineral Phase System
    do iMk = 1,nMk
       MkFas = vFas(vKinFas(iMk)%iFas)
       vOrdMk(iMk) = MkFas%iSpc    
    end do
    
    do iMk = 1,nMk 
       iSpc = vOrdMk(iMk) 
       tNuMk(iMk,1:nPr)=  tNuSp(iSpc,1:nPr)
       tAlfMk(1:nCp,iMk)= tAlfSp(1:nCp,iSpc)
    end do
    
    call Box_System_Info
    
  end subroutine Box_System_Build
  
  !---
  
  subroutine Box_System_Info
    use M_EchoMat
    use M_Box_Debug_Vars
    use M_Global_Vars, only : vSpc, vEle
    use M_T_Component, only : Component_Print

    implicit none
    integer :: icp
    
    call Info_("Box_System_Info")

    if (LDebug_System) then
    
       write(*,*) "=================="
       write(*,*) "nAq = ", nAq
       write(*,*) "nMk = ", nMk
       write(*,*) "-------------"
       write(*,*) "nCp = ", nCp
       write(*,*) "nPr = ", nPr
       write(*,*) "nAs = ", nAs
       
       write(*,*) "=================="
       do iCp=1,nCp 
          call Component_Print(6,vEle,vSpc,vCpnBox(iCp)) 
       end do

       write(*,*) "==================="
       call echomat_real(tAlfPr, nCp, nPr, 'Formul Pr tAlfPr', 'T')
       call echomat_real(tAlfAs, nCp, nAs, 'Formul As tAlfAs', 'T')
       write(*,*) "==================="
       call echomat_real(tAlfMk, nCp, nMk, 'Formul Mk tAlfMk', 'T')
       call echomat_real(tAlfAq, nCp, nAq, 'Formul Aq tAlfAq', 'T')       
       write(*,*) "==================="
       call echomat_real(tNuAs,  nAs, nPr, 'Stoik  As tNuAs')
       call echomat_real(tNuMk,  nMk, nPr, 'Stoik  Mk tNuMk')
    end if
    
  end subroutine Box_System_Info

  !---

  subroutine Box_System_Info_Component
    !===============================================
    !  Show Components Differences
    !===============================================
    use M_System_Vars,  only: vCpn
    use M_Global_Vars,  only: vSpc, vEle
    use M_Box_Debug_Vars
    implicit none
    !
    integer :: iCp
    call Info_("Box_System_Info_Component")
    
    if (LDebug_System) then
       write(*,*) "======================================================"
       write(*,*) " Compare Components vCpn | vCpnBox                    "
       write(*,*) "======================================================"
       
       do iCp=1,nCp 
          write(*,*) "============= Component [", iCp,"]"
          call Component_Print(6,vEle,vSpc,vCpn(iCp)) 
          call Component_Print(6,vEle,vSpc,vCpnBox(iCp))
       end do
       write(*,*) "======================================================"
    end if
    
  end subroutine Box_System_Info_Component

end module M_Box_System
