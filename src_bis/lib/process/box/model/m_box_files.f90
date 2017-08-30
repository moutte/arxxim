module M_Box_Files

  use M_Kinds
  use M_Trace
  use M_Box_Files_Vars

  implicit none
  !
  private
  !
  public :: Box_Files_EnTeteFMnk
  public :: Box_Files_EnTeteFMol   
  public :: Box_Files_EnTeteFEle
  !
  public :: Box_Files_ShoMin
  public :: Box_Files_ShoAqu
  public :: Box_Files_ShoEle
  !---
  public :: Box_Files_Init
  public :: Box_Files_Close

contains

  !---

  subroutine Box_Files_EnTeteFMnk    
    use M_Box_Param_Vars, only : nMk
    use M_Global_Vars,    only : vKinFas
    use M_Box_TimeLoop_Vars,   only : TUnit
    implicit none
    !
    integer::i,N
    integer:: f
    !
    f= fBoxMnK
    N= nMk
    !
    !// Time Step
    !write(f,'(2(A15,A1))',advance="no") "iStep   ", T_,"Time    ",T_, "Time_"//TUnit, T_
    write(f,'(2(A15,A1))',advance="no") "iStep   ", T_, "Time_"//TUnit, T_
    
    !// General Properties
    write(f,'(2(A15,A1))',advance="no") "pH    ",T_,"PhiFluid    ",T_

    !// PhiM Minerals
    do i=1,N
       write(f,'(A15,A1)',advance="no") "PhiM_"//trim(vKinFas(i)%NamKF),T_
    enddo

    !// LogQsk
    do i=1,N
       write(f,'(A15,A1)',advance="no") "LogQsK_"//trim(vKinFas(i)%NamKF),T_
    enddo

    !// Mole Number
    do i=1,N
       write(f,'(A15,A1)',advance="no") "MolNr"//trim(vKinFas(i)%NamKF),T_

    enddo

    write(f,*)

  end subroutine Box_Files_EnTeteFMnk

  !---

  subroutine Box_Files_EnTeteFMol   
    use M_Box_Param_Vars,  only : nAq
    use M_Global_Vars,     only : vSpc
    use M_Box_System_Vars, only : vOrdAq, vOrdMk
    use M_Box_TimeLoop_Vars,   only : TUnit
    implicit none
    !
    integer::iAq, iSpc, N
    integer:: f
    !
    f= fBoxMol
    N= nAq
    !    
    !// Time Step
    !write(f,'(3(A15,A1))',advance="no") "iStep   ", T_,"Time    ",T_, "Time_"//TUnit, T_
    write(f,'(2(A15,A1))',advance="no") "iStep   ",T_, "Time_"//TUnit, T_
    
    !// General Properties
    write(f,'(2(A15,A1))',advance="no") "pH    ",T_,"PhiFluid    ",T_

    !// Mole Number  Aqueous Species
    do iAq=1,N
       iSpc = vOrdAq(iAq)
       write(f,'(A15,A1)',advance="no") 'n_'//trim(vSpc(iSpc)%NamSp),T_
    enddo

    write(f,*)

  end subroutine Box_Files_EnTeteFMol

  !---

  subroutine Box_Files_EnTeteFEle  
    use M_Box_Param_Vars,  only : nCp
    use M_Box_TimeLoop_Vars,   only : TUnit
    use M_Box_System_Vars, only : vCpnBox
    implicit none
    !
    integer::iCp, N
    integer:: f
    !
    f= fBoxEle
    N= nCp
    !    
    !// Time Step
    !write(f,'(3(A15,A1))',advance="no") "iStep   ", T_,"Time    ",T_, "Time_"//TUnit, T_
    write(f,'(2(A15,A1))',advance="no") "iStep   ",T_, "Time_"//TUnit, T_
    
    !// General Properties
    write(f,'(2(A15,A1))',advance="no") "pH    ",T_,"PhiFluid    ",T_

    !// Mole Number Aqueous Element
    do iCp=1,N
       write(f,'(A15,A1)',advance="no") 'nF_'//trim(vCpnBox(iCp)%NamCp),T_
    enddo

   !// Mole Number Mineral Element
    do iCp=1,N
       write(f,'(A15,A1)',advance="no") 'nS_'//trim(vCpnBox(iCp)%NamCp),T_
    enddo

    write(f,*)

  end subroutine Box_Files_EnTeteFEle

  !---

  subroutine Box_Files_ShoMin
    use M_Box_Vars,          only : Time, PhiF, vMolK, vPhiK, nMk
    use M_Box_Thermo_Vars,   only : pH, vOmegaK
    use M_Box_TimeLoop_Vars, only : iTimeStep, TUnit
    use M_SetOption
    implicit none
    !
    integer::I,N
    integer:: f
    real(dp) :: TFact
    !
    N= nMk
    f= fBoxMnK
    call ComputeTimeFactor( TFact, Tunit)

    !// Time Step
    !write(f,'(I15,A1,2(G15.8,A1))',advance="no") iTimeStep,T_,Time, T_, Time/TFact,T_
    write(f,'(I15,A1,1(G15.8,A1))',advance="no") iTimeStep, T_, Time/TFact,T_
    
    !// General Properties
    write(f,'(2(G15.8,A1))',advance="no") pH,T_, PhiF,T_
    
    !// PhiM Minerals
    do i=1,N
       write(f,'(G15.8,A1)',advance="no") vPhiK(i),T_
    enddo

     !// LogQsk
    do i=1,N
       write(f,'(G15.8,A1)',advance="no") log10(vOmegaK(i)),T_
    enddo

     !// Mole Number
    do i=1,N
       write(f,'(G15.8,A1)',advance="no") vMolK(i),T_
    enddo

    write(f,*)
    
  end subroutine Box_Files_ShoMin

 !---

  subroutine Box_Files_ShoAqu
    use M_Box_Vars,          only : Time, PhiF, vMolF, nAq
    use M_Box_Thermo_Vars,   only : pH
    use M_Box_TimeLoop_Vars, only : iTimeStep, TUnit
    use M_SetOption   
    implicit none
    !
    integer::iAq,N
    integer:: f
    real(dp) :: TFact
    !
    N= nAq
    f= fBoxMol
    call ComputeTimeFactor( TFact, Tunit)
   
    !// Time Step
    !write(f,'(I15,A1,2(G15.8,A1))',advance="no") iTimeStep,T_,Time, T_, Time/TFact,T_
    write(f,'(I15,A1,1(G15.8,A1))',advance="no") iTimeStep, T_, Time/TFact,T_
    
    !// General Properties
    write(f,'(2(G15.8,A1))',advance="no") pH,T_, PhiF,T_
    
    !// Mole Number  Aqueous Species
    do iAq=1,N
       write(f,'(G15.8,A1)',advance="no") vMolF(iAq),T_
    enddo

    write(f,*)
    
  end subroutine Box_Files_ShoAqu

  !---

  subroutine Box_Files_ShoEle
    use M_Box_Vars,          only : Time, PhiF, vMolF, vMolK, nAq, nMk
    use M_Box_Thermo_Vars,   only : pH
    use M_Box_TimeLoop_Vars, only : iTimeStep, TUnit
    use M_Box_System_Vars,   only : nCp, tAlfAq, tAlfMk
    use M_SetOption   
    implicit none
    !
    integer:: iAq, iMk, iCp
    integer:: f
    real(dp) :: TFact
    real(dp) :: vTotFCp, vTotMkCp
    !    
    f= fBoxEle
    call ComputeTimeFactor( TFact, Tunit)
   
    !// Time Step
    !write(f,'(I15,A1,2(G15.8,A1))',advance="no") iTimeStep,T_,Time, T_, Time/TFact,T_
    write(f,'(I15,A1,1(G15.8,A1))',advance="no") iTimeStep, T_, Time/TFact,T_
    
    !// General Properties
    write(f,'(2(G15.8,A1))',advance="no") pH,T_, PhiF,T_
    
    !// Mole Number  Aqueous Species
    do iCp=1,nCp 
       vTotFCp = 0.D0
       do iAq = 1,nAq
          vTotFCp = vTotFCp + tAlfAq(iCp,iAq)*vMolF(iAq)
       end do
       write(f,'(G15.8,A1)',advance="no") vTotFCp, T_
    end do

    do iCp=1,nCp 
       vTotMkCp = 0.D0
       do iMk = 1,nMk
          vTotMkCp = vTotMkCp + tAlfMk(iCp,iMk)*vMolK(iMk)
       end do
       write(f,'(G15.8,A1)',advance="no") vTotMkCp, T_
    end do
    
    write(f,*)
    
  end subroutine Box_Files_ShoEle

  !--- 

  subroutine Box_Files_Init
    use M_Files,        only: DirOut, Files_Index_Write
    use M_IOtools,      only: GetUnit
    implicit none
    integer :: IOS
    !---    
    if(fBoxMnK==0) then
       call GetUnit(fBoxMnK)
       open(UNIT=fBoxMnK,iostat=IOS,file=trim(DirOut)//"_minmol.restab")
       if (IOS>0) call Stop_("Cannot open output File ["//trim(DirOut)//"][_minmol.restab]")
       !
       call Files_Index_Write(fHtm, &
            & trim(DirOut)//"_minmol.restab", &
            & "DYNAMIC: at each time step, pH, vol'fractions fluid and minerals, logQsK, etc")
       call Box_Files_EnTeteFMnk
    end if
    
    if(fBoxMol==0) then
       call GetUnit(fBoxMol)
       open(UNIT=fBoxMol,iostat=IOS, file=trim(DirOut)//"_aqueous.restab")
       if (IOS>0) call Stop_("Cannot open output File ["//trim(DirOut)//"][_aqueous.restab]")
       !
       call Files_Index_Write(fHtm, &
            & trim(DirOut)//"_aqueous.restab", &
            & "DYNAMIC: at each time step, pH, aqueous mol numbers")
       call Box_Files_EnTeteFMol
    end if

    if(fBoxEle==0) then
       call GetUnit(fBoxEle)
       open(UNIT=fBoxEle,iostat=IOS, file=trim(DirOut)//"_element.restab")
       if (IOS>0) call Stop_("Cannot open output File ["//trim(DirOut)//"][_element.restab]")
       !
       call Files_Index_Write(fHtm, &
            & trim(DirOut)//"_elmement.restab", &
            & "DYNAMIC: at each time step, pH, element mol numbers")
       call Box_Files_EnTeteFEle
    end if

  end subroutine Box_Files_Init
  
  !---
  
  subroutine Box_Files_Close
    implicit none
    
    if(fBoxMnK>0)  then
       close(fBoxMnK)
       fBoxMnK=    0
    end if

    if(fBoxMol>0)  then
       close(fBoxMol)
       fBoxMol=    0
    end if

    if(fBoxEle>0)  then
       close(fBoxEle)
       fBoxEle=    0
    end if

  end subroutine Box_Files_Close

end module M_Box_Files

