module M_Solmodel_Calc_Davies
!-----------------------------------------------------------------------
! Davies models of activity for Aqueous Species
!
! DAV_1 :  Davies 1
!   after Grenthe et al, Modelling In Aquatic Chemistry, OECD,Chap.IX,p.330
!   (same as in Bethke,'96, p110)
!   Davies, with Gamma(Neutral)=1.0
!
! DAV_2 :  Davies 2
!   after Grenthe et al, Modelling In Aquatic Chemistry, OECD,p.330
!   Davies, with Gamma(Neutral) a function of IonStrength
!
! SAMSON : Samson et al. 1999
!   Modeling chemical activity effects in strong ionic solutions.
!   Computational Materials Science 15: 285-294
!
!-----------------------------------------------------------------------

  use M_Numeric_Const
  use M_Kinds
  use M_Trace,only: fTrc,iDebug,T_

  implicit none
  private

  public:: Solmodel_Calc_Davies

contains

  subroutine Solmodel_Calc_Davies( &
  & vZSp, vMolal, &
  & dhA, iModel, MWSv, &
  & vLnGam, LnActSv, OsmoSv) !,ERR)

    integer, intent(in)::  vZSp(:)
    real(dp),intent(in)::  vMolal(:)
    real(dp),intent(in)::  dhA
    real(dp),intent(in)::  MWSv
    integer, intent(in)::  iModel

    real(dp),intent(out):: vLnGam(:)
    real(dp),intent(out):: LnActSv
    real(dp),intent(out):: OsmoSv

    real(dp):: LnGam_Charged, LnGam_Neutral
    real(dp):: Sum_Solute
    real(dp):: IonStr, SqrtIoSt
    integer :: vZ2(size(vZSp))
    integer :: nSolut,iAq
    !-------------------------------------------------------------------

    nSolut= size(vMolal)
    vZ2(:)= vZSp(:)*vZSp(:)

    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= dot_product(vZ2(:),vMolal(:)) /Two

    !-- "TRUNCATED DAVIES"
    !-- if(IonStr>0.3) IonStr- 0.3D0

    SqrtIoSt= sqrt(IonStr)

    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12

    !---------------------------------------------------------< Solute--
    select case(iModel)
    !
    case(6,7) !"DAV_1","DAV_2")
      LnGam_Charged = -dhA &
      &               *( SqrtIoSt/(One + SqrtIoSt) -0.3D0*IonStr ) &
      &               *Ln10
      ! there's also a version of Davies equation
      ! where coeff 0.2 is used instead of 0.3 !?!
    !
    case(9) !("SAMSON")
      LnGam_Charged = -dhA &
      &               *( SqrtIoSt/(One + SqrtIoSt) &
      &                 -(0.2D0 -4.17D-5*IonStr) *IonStr ) &
      &               *Ln10
    !
    end select
    
    select case(iModel)
    !
    case(6,9) !("DAV_1","SAMSON")
      LnGam_Neutral = Zero
      !-> official Davies equation
    !
    case(7) !"DAV_2")
      LnGam_Neutral = 0.064D0*IonStr*Ln10
      !-> non-official equation (where is it from ??)
    !
    end select
    
    do iAq=1,nSolut
      if (vZSp(iAq)/=0) then
        !------------------Charged Species--
        ! e.g. Bethke, eq.7.4, p110
        vLnGam(iAq)= vZ2(iAq) *LnGam_Charged
      else
        !------------------Neutral Species--
        vLnGam(iAq)= LnGam_Neutral
      end if
    end do
    !---------------------------------------------------------</Solute--
    
    !-------------------------------------------------------- Solvent --
    Sum_Solute = sum(vMolal(:))
    !___________LewisRandall,23-39,p347; EQ3NR-Doc, p.39
    LnActSv= Ln10 *Two/3.0D0 *dhA  *SqrtIost**3 *Sigma(SqrtIoSt) &
    &      - Ln10            *dhA  *IonStr*IonStr &
    &      - Sum_Solute
    LnActSv=  LnActSv * MWSv
    OsmoSv= - LnActSv /MWsv / Sum_Solute
    !--------------------------------------------------------/Solvent --

  end subroutine Solmodel_Calc_Davies

  real(dp) function Sigma(x)
  ! for calcul. Water Acivity,
  ! cf equ.23-40 in LewisRandall'61, used also in EQ3NR
    real(dp),intent(in):: x
    
    Sigma=  3.0D0*(One + x - One/(One+x) - Two*log(One+x))/x/x/x
  
  end function Sigma

end module M_Solmodel_Calc_Davies
