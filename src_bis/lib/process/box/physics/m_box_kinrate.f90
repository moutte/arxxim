module M_Box_KinRate 

  !====================================================
  ! Rate calculations on kinetic phases
  !====================================================

  use M_Kinds
  use M_Trace

  implicit none
  private
  !

  public:: Box_KinRate_CalcQsK
  public:: Box_KinRate_SatState
  
  public:: Box_KinRate_CalcQsKFactor   
  public:: Box_KinRate_CalcActivFactor

contains

  !---

  subroutine Box_KinRate_CalcQsK( &
       & nCi,       & !IN
       & nCx,       & !IN
       & DG0rt,     & !IN
       & vNu,       & !IN
       & vLnX,      & !IN
       & vLnGam,    & !IN
       & vLnActBuf, & !IN
       & LnActW,    & !IN
       & QsK,       & !OUT
       & dQsK_dLnXi)  !OUT

    !========================================================================================
    !  purpose : calc. Q/K and d(Q/K)/dLnX of a phase of properties DG,vNu 
    !            for a fluid composition vLnX composition input in affinity 
    !            (log(mole numbers))
    !========================================================================================

    !----------------------------------------------------------------------------------------
    !
    ! lnAct= lnMolal + lnGamma=  lnMol + lnGam - lnMolH2O - lnMWH2O
    ! QsK=exp(dot_product(tNuMk(iMk,1:nCp),vLnAct(1:nCp)) - vDG_Mk(iMk))
    ! or,directly from G0,
    ! QsK=exp(dot_product(tNuMk(iMk,1:nCp),vLnAct(1:nCp) +vSpc(1:nCp)%G0) - vKinFas(iMk)%G0)
    !
    ! caveat: Solvent should follow mole fraction scale, not molality ...
    ! here, we assume "constant" Solvent activity:
    !
    !----------------------------------------------------------------------------------------

    use M_Numeric_Const, only: MinExpDP, MaxExpDP
    use M_Basis_Vars,    only: isW, MWSv
    !
    implicit none  
    integer,              intent(in) :: nCi        !nr inert comp'nt
    integer,              intent(in) :: nCx        !nr mobile comp'nt
    real(dp),             intent(in) :: DG0rt      !deltaG/RT of formation reaction
    real(dp),dimension(:),intent(in) :: vNu        !stoikio of formation reaction
    real(dp),dimension(:),intent(in) :: vLnX       !log(Mole Nrs. Prim.Species) !1:nCi
    real(dp),dimension(:),intent(in) :: vLnGam     !log(Act.Coeff. Prim.Species) !1:nCi
    real(dp),dimension(:),intent(in) :: vLnActBuf  !log(ActivityBufferSpecies) !1:nCx
    real(dp),             intent(in) :: LnActW     !log(Activity Solvent)
    real(dp),             intent(out):: QsK
    real(dp),dimension(:),intent(out):: dQsK_dLnXi  !1:nCi
    !
    real(dp)::LnQsK

    !>>> Warning :: Strong asumption : isW = 1 !!!

    !// Inert Primary Species Contribution 
    LnQsK=     vNu(isW)  *LnActW &
         + dot_product(vNu(2:nCi),vLnX(2:nCi))           & !log(mole numbers)
         + dot_product(vNu(2:nCi),vLnGam(2:nCi))         & !log(gammas)
         -         SUM(vNu(2:nCi))*(vLnX(isW)+log(MWSv)) & !log(mass Solvent)
         - DG0rt                                           !equiv -logK

    !// Mobile Primary Species Contribution 
    if (nCx>0) &
         & LnQsK= LnQsK & !in dynamic mode, always nAx=0 -> nCx=nMx
         &  + dot_product(vNu(nCi+1:nCi+nCx),vLnActBuf(1:nCx)) !mobile species

    !// Tranform QsK 
    if (LnQsK>MinExpDP .and. LnQsK<MaxExpDP) QsK= exp(LnQsK)
    if (LnQsK<=MinExpDP) QsK= exp(MinExpDP)
    if (LnQsK>=MaxExpDP) QsK= exp(MaxExpDP)

    !// Compute Derivatives 
    dQsK_dLnXi(:)=  Zero
    dQsK_dLnXi(isW)= - SUM(vNu(2:nCi))*QsK
    dQsK_dLnXi(2:nCi)=  vNu(2:nCi) *QsK

  end subroutine Box_KinRate_CalcQsK

  !---

  subroutine Box_KinRate_SatState( &       
  & M,         & !IN
  & nMol,      & !IN
  & nMolMinim, & !IN
  & QsK,       & !IN
  & Iota,      & !IN 
  & cSatur)      !OUT

    !==========================================================================
    ! purpose : from QsK deduce cSatur, for branching to dissol/precip
    !--------------------------------------------------------------------------
    !  DISSOLU      | MINIMAL      || INERT                   | PRECIPI        |
    !  QsK < 1-Iota | QsK < 1-Iota || 1-Iota < QsK <Iota+Qss  | QsK > Iota+Qss |
    !  nMol> Min    | nMol < Min   ||       -                 |       -        |
    !==========================================================================

    use M_T_KinFas,only: T_KinFas
    implicit none

    type(T_KinFas),  intent(in) :: M     
    real(dp),        intent(in) :: nMol   
    real(dp),        intent(in) :: nMolMinim 
    real(dp),        intent(in) :: QsK   
    real(dp),        intent(in) :: Iota   

    character,intent(out):: cSatur
    !--
    real(dp) :: QsKseuil
    !--
    QsKseuil = M%QsKSeuil

    if ( QsK < One - Iota ) then
       if ( nMol > nMolMinim ) then
          cSatur="D" !"DISSOLU"
       else
          cSatur="M" !"MINIMAL" 
       end if
    elseif ( QsK > QsKseuil + Iota ) then
       cSatur="P" !"PRECIPI" 
    else
       cSatur="I" !"INERT"  
    end if

  end subroutine Box_KinRate_SatState

  !---

  subroutine Box_KinRate_CalcQsKFactor(&
  & cSatur,    & !IN
  & M,         & !IN: Kinetic model
  & QsK,       & !IN
  & dQsK_dLnXCi,  & !IN
  & VmQsK,     & !OUT
  & dVmQsK_dLnX)  !OUT

    !==================================================================
    ! purpose : Compute VmQsK = (1 - QsK^a)^b or (QsK^a-1)^b
    !==================================================================

    use M_T_KinModel,only: T_KinModel
    implicit none

    character,             intent(in) :: cSatur
    type(T_KinModel),      intent(in) :: M
    real(dp),              intent(in) :: QsK
    real(dp),dimension(:),intent(in) :: dQsK_dLnXCi    ! (1:nCi)
    real(dp),              intent(out):: VmQsK
    real(dp),dimension(:), intent(out):: dVmQsK_dLnX ! (1:nAq)
    !
    real(dp)::dVmQsk_dQsK
    integer ::nCi

    !// Partial derivative relative to QsK
    select case(trim(cSatur))

    case("D") !"DISSOLU"
       VmQsK= - (One - QsK**M%AlfaD)**M%BetaD
       dVmQsK_dQsK= M%BetaD &
            &         *(One - QsK**M%AlfaD)**(M%BetaD-One) &
            &         *M%AlfaD *QsK**(M%AlfaD-One)

    case("P") !"PRECIPI"
       VmQsK= (QsK**M%AlfaP - One)**M%BetaP
       dVmQsk_dQsK= M%BetaP &
            &       *(QsK**M%AlfaP - One)**(M%BetaP-One) &
            &       *M%AlfaP *QsK**(M%AlfaP-One)

    case default ! I(NERT), M(INIMAL), etc...
       VmQsK= Zero
       dVmQsK_dQsK= Zero

    end select

    !// Compose Derivative  
    nCi= size(dQsK_dLnXCi)
    dVmQsK_dLnX= Zero
    dVmQsK_dLnX(1:nCi) = dVmQsK_dQsK * dQsK_dLnXCi(1:nCi)

  end subroutine Box_KinRate_CalcQsKFactor

  !---

  subroutine Box_KinRate_CalcActivFactor(&
  & cSatur, &  !IN:  "DISSOLU"/"PRECIPI"
  & M,      &  !IN:  kinetic model of mineral
  & vLnAct, &  !IN:  fluid species activities
  & VmAct,  &  !OUT: factor depending on fluid chemistry only
  & dVmAct_dLnX ) !OUT: its derivative vs ln(X)

    !==================================================================
    ! purpose : Compute VmAct = k(T) * ai^ni *aj^nj
    !           this factor depends on fluid chemistry only
    !==================================================================

    use M_T_KinModel
    use M_T_Species,only: T_Species
    use M_Box_Vars, only : nAq
    use M_Global_Vars, only : vSpc 
    use M_T_KinFas,only: T_KinFas
    use M_Basis_Vars, only : vPrmBk
    implicit none

    character,            intent(in) :: cSatur
    type(T_KinModel),     intent(in) :: M
    real(dp),dimension(:),intent(in) :: vLnAct

    real(dp),             intent(out):: VmAct
    real(dp),dimension(:),intent(out):: dVmAct_dLnX 

    integer ::I
    integer :: iAq, jAq
    real(dp):: kT, X

    select case(trim(cSatur))

    case("D") !("DISSOLU")
       VmAct= Zero
       do I=1,M%NTermD         
          kT= M%Kd(I)
          X= One
          if (M%ISpcDiss(I)>0) then
             iAq = vPrmBk(M%ISpcDiss(I))
             if (iAq>nAq) call Fatal_("Wrong Aqueous Species Index in KinModel")
             X= X * exp( M%N_d(I) * vLnAct(iAq) )
          end if
          if (M%JSpcDiss(I)>0) then
             jAq = vPrmBk(M%JSpcDiss(I))
             if (jAq>nAq) call Fatal_("Wrong Aqueous Species Index in KinModel")
             X= X * exp( M%NJ_d(I)* vLnAct(jAq) )
          end if
          VmAct= VmAct + kT * X
       enddo

    case("P") ! ("PRECIPI")
       VmAct= Zero
       do I=1,M%NTermP
          kT= M%Kp(I)
          X= One
          if (M%ISpcDiss(I)>0) then
             iAq = vPrmBk(M%ISpcPrec(I))
             if (iAq>nAq) call Fatal_("Wrong Aqueous Species Index in KinModel")
             X= X * exp( M%N_p(I)  * vLnAct(iAq) )
          end if
          if (M%JSpcPrec(I)>0) then
             jAq = vPrmBk(M%JSpcPrec(I))
             if (jAq>nAq) call Fatal_("Wrong Aqueous Species Index in KinModel")
             X= X * exp( M%NJ_p(I) * vLnAct(jAq) )
          end if
          VmAct= VmAct + kT * X
       enddo

    case default ! M=MINIMAL, I=INERT, etc ...
       VmAct= Zero

    end select

    !-------------------------------------------------------------------
    ! provisionally, dVmAct_dLnX is not computed in this version  
    ! we assume that VLnAct is an explicit term         
    !-------------------------------------------------------------------
    dVmAct_dLnX= Zero  


  end subroutine Box_KinRate_CalcActivFactor

end module M_Box_KinRate
