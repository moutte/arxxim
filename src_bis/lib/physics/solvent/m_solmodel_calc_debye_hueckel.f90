module M_Solmodel_Calc_Debye_Hueckel
!-----------------------------------------------------------------------
! Debye-Hueckel models of activity for Aqueous Species
!
! DH1    : DebyeHueckel
! DH2    : DebyeHueckel + Bdot
!
! DH1EQ3 : DebyeHueckel + EQ3 for neutral species
! DH2EQ3 : DebyeHueckel + Bdot + EQ3 for neutral species
!
!-----------------------------------------------------------------------

  use M_Kinds
  use M_Numeric_Const
  use M_Trace,only: fTrc,iDebug,T_
  implicit none
  private

  public:: Solmodel_Calc_Debye_Hueckel

contains

  subroutine Solmodel_Calc_Debye_Hueckel( &
  & TdgK, Pbar, vZSp, vMolal,          & !IN
  & dhA, dhB, vA0, bDot, iModel, MWSv, & !IN
  & vLnGam, LnActSv, OsmoSv)

    integer, intent(in) :: vZSp(:)
    real(dp),intent(in) :: vMolal(:)
    real(dp),intent(in) :: TdgK, Pbar
    real(dp),intent(in) :: vA0(:)
    real(dp),intent(in) :: dhA
    real(dp),intent(in) :: dhB
    real(dp),intent(in) :: bDot
    real(dp),intent(in) :: MWSv
    integer, intent(in) :: iModel

    real(dp),intent(out):: vLnGam(:)
    real(dp),intent(out):: LnActSv
    real(dp),intent(out):: OsmoSv

    real(dp):: Sum_Solute
    real(dp):: LnGam_Neutral
    real(dp):: A0,X
    real(dp):: bDot_Coeff
    real(dp):: IonStr, SqrtIoSt
    integer :: vZ2(size(vZSp))
    integer :: nSolut,iAq

    nSolut= size(vMolal)
    Sum_Solute= SUM(vMolal(:))
    vZ2(:)= vZSp(:)*vZSp(:)
    !
    !Ion Strength=  sum(m(i)*z2(i)) /2
    IonStr= dot_product(vZ2(:),vMolal(:)) /Two
    SqrtIoSt= SQRT(IonStr)
    !
    !-------------------------------------------------------------Solute
    !------------------------------------bDot Models Specific Parameters
    bDot_Coeff= Zero
    !
    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12
    !
    select case(iModel)

    case(4,5) !("DH2","DH2EQ3")
      bDot_Coeff= bDot

    case(10) !("HKF81") ! provisional
      bDot_Coeff= 0.03D0 !bDot
      ! value taken by Dolejs & Wagner,GCA, 2008

    end select
    !
    !-------------------------------------Eq3 Models Specific Parameters
    LnGam_Neutral= Zero
    select case(iModel)

    case(3,5) !("DH1EQ3", "DH2EQ3")
       ! retrieved from EQ36, v7.0, user's manual (Wolery, 1992)
       ! activ'coeff' of CO2(aq) in a (otherwise pure) NaCl solution
       ! as a function of ionic strength
       ! uses an expression proposed by Drummond, 1981
       LnGam_Neutral= (-1.0312D0 + 0.0012806D0*TdgK +255.9D0/TdgK)*IonStr &
       &            - (0.4445D0  - 0.001606D0 *TdgK              )*IonStr/(One+IonStr)
       !from Wolery'92, p41:
       !this should be applied only to non-polar species (O2(aq), H2(aq), N2(aq), ...),
       !for which a salting out effect is expected

    case(10) !("HKF81")
      LnGam_Neutral= 0.03D0 *IonStr *Ln10
      ! value taken by Dolejs & Wagner,GCA, 2008

    end select

    do iAq=1,nSolut
      if (vZSp(iAq)/=0) then
        !------------------------------------------------Charged Species
        select case(iModel)

        case(2,3,4,5) !("DH1","DH1EQ3","DH2","DH2EQ3")
          X= - dhA *vZ2(iAq) *SqrtIoSt/(One+vA0(iAq)*dhB*SqrtIoSt) &
          &  + bDot_Coeff*IonStr
          vLnGam(iAq)= X *Ln10

        case(10) !("HKF81")
          A0= 3.72D0
          ! Y=
          X= - dhA*vZ2(iAq)*SqrtIoSt/(One+A0*dhB*SqrtIoSt) &
          &  + log(One +MWsv *Sum_Solute) &
          &  + bDot_Coeff*IonStr
          vLnGam(iAq)= X *Ln10

        end select
      else
        !------------------------------------------------Neutral species
        vLnGam(iAq)= LnGam_Neutral
      end if
    end do
    !-----------------------------------------------------------/ Solute

    !----------------------------------------------------------- solvent
    !DebyeHueckel + BDot, cf EQ3NR-Doc, p.41
    LnActSv= Ln10 *Two/3.0D0 *dhA *SqrtIost**3 *Sigma(4.0D0*dhB*SqrtIoSt) &
    &      - Ln10            *bDot*IonStr*IonStr &
    &      - Sum_Solute
    !vLnAct(isW)= Zero
    LnActSv=  LnActSv *MWSv
    OsmoSv= - LnActSv /MWsv /Sum_Solute !osmotic coeff'
    !----------------------------------------------------------/ solvent

  end subroutine Solmodel_Calc_Debye_Hueckel

  function Sigma(x)
    ! for calcul. Water Acivity,
    ! cf equ.23-40 in LewisRandall'61, used also in EQ3NR
    real(dp):: Sigma
    real(dp),intent(in):: x

    Sigma= 3.0D0 *(One + x - One/(One+x) - Two*log(One+x))/x/x/x

  end function Sigma

end module M_Solmodel_Calc_Debye_Hueckel
