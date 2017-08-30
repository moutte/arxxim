module M_Numeric_Show

  use M_Kinds
  use M_Trace
  implicit none
  private
  
  public :: TestJacob_Sho
  public :: Newton_Sho
  public :: Jacobian_Sho
  public :: NewtonError_Sho

contains 

  subroutine NewtonError_Sho(iErr)
    implicit none
    integer :: iErr
    character(len=128) :: Msg
    call NewtonError_Decode(iErr, Msg)
    call Message_(Msg)    
  end subroutine NewtonError_Sho

  !---
  
  subroutine NewtonError_Decode(iErr, Str)
    !-----------------------------------------------------------------------------------------
    ! Decode Error Code 
    !-----------------------------------------------------------------------------------------
    !iErr= 0 : convergence reached -> OK
    !iErr=-1 : MaxIts reached without convergence 
    !iErr=-2 : singular jacobian
    !iErr=-3 : roundoff problem in linesearch
    !iErr=-4 : no convergence on vFunc, problem with gradient ?
    !iErr=-5 : no convergence on vFunc, vX very close to previous vX -> stationary point ??
    !-----------------------------------------------------------------------------------------
    implicit none
    integer :: iErr
    character(len=*) :: Str
    select case (iErr)
    case(0)  ; Str = "Convergence reached -> OK"
    case(-1) ; Str = "iErr = -1 : MaxIts reached without convergence"
    case(-2) ; Str = "iErr = -2 : Singular jacobian"
    case(-3) ; Str = "iErr = -3 : Roundoff problem in linesearch"
    case(-4) ; Str = "iErr = -4 : No convergence on vFunc, problem with gradient ?"
    case(-5) ; Str = "iErr = -5 : no convergence on vFunc, vX very close to previous vX -> stationary point ??"
    case default 
       call Fatal_("Unknown Error code")
    end select
    
  end subroutine NewtonError_Decode

  !---

  subroutine TestJacob_Sho(tJacob,tJacobNum)
    use M_IoTools,only: GetUnit, OutStrVec
    implicit none
    !
    real(dp),dimension(:,:),intent(in):: tJacob,tJacobNum
    !
    real(dp),dimension(size(tJacob,1))  :: xx,xr
    integer:: ff,I
    !
    write(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN"
    write(fTrc,'(A)') "results in zzDynam_TestJacob.log"
    call GetUnit(ff)
    open(ff,file="debug_dynam_jacob.log")
    do I=1,size(tJacob,1)
       call OutStrVec(ff,tJacob(I,:),   Opt_I=I,Opt_S="JacobAna",Opt_C="G")
       call OutStrVec(ff,tJacobNum(I,:),Opt_I=I,Opt_S="JacobDif",Opt_C="G")
       !
       xx(:)= tJacob(I,:) -tJacobNum(I,:)
       call OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltJac",   Opt_C="G")
       !
       xr(:)= tJacob(I,:) +tJacobNum(I,:)
       where(xr(:)>1.D-16) ; xx(:)= xx(:) / xr(:)
       elsewhere           ; xx(:)= Zero
       endwhere
       call OutStrVec(ff,xx(:),         Opt_I=I,Opt_S="DeltRel",   Opt_C="G")
    end do
    write(ff,*)
    close(ff)
    write(fTrc,'(/,A,/)') "!!!!!! TEST JACOBIAN end"
  end subroutine TestJacob_Sho

  !---
  
  subroutine Newton_Sho(vX,vFunc,Its)
    use M_Numeric_Const,only: Ln10
    use M_Numeric_Tools,only: fNewtF,fNewtR,fNewt_I
    use M_IoTools,      only: OutStrVec
    implicit none
    !
    real(dp),intent(in):: vX(:),vFunc(:)
    integer, intent(in):: Its
    !
    fNewt_I= fNewt_I +1
    if(fNewtF>0) call OutStrVec(fNewtF,vX(:)/Ln10,   Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
    if(fNewtR>0) call OutStrVec(fNewtR,ABS(vFunc(:)),Opt_I=fNewt_I,Opt_J=Its,Opt_C="G")
    !
  end subroutine Newton_Sho

  !---

  subroutine Jacobian_Sho(t)
    implicit none
    real(dp),intent(in):: t(:,:)
    integer:: I,J
    do I=1,size(t,1)
       do J=1,size(t,2)
          if(t(I,J)/=0) then; write(fTrc,'(F7.1,A1)',advance="no") t(I,J), T_
          else              ; write(fTrc,'(A1,  A1)',advance="no") "0",       T_
          end if
       end do
       write(fTrc,*)
    end do
    return
  end subroutine Jacobian_Sho


end module M_Numeric_Show
