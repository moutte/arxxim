module M_Solmodel_Calc_System

  !---------------------------------------------------------------
  ! Test for a New organization to compute Gamma
  ! Solvent Activity Calc Interface for all Models   
  ! Compute intermediate properties and call subroutines 
  !---------------------------------------------------------------

  use M_Kinds
  use M_Trace

  implicit none
  private 

  public :: Solmodel_Calc_System
  public :: Solmodel_Calc_Physics

contains
  
  !-- WARNING
  !-- when implementing a new activity model,
  !-- add its name to the list of activity model names
  
  subroutine Solmodel_Calc_System(&
  & TdgK,Pbar,        & !IN
  & SolModel,         & !IN  ! solvent phase model
  & isW,vSpcAq,vLAx,  & !IN
  & vMolF,            & !IN
  & vLnAct,vLnGam,    & !INOUT
  & vTooLow_,OsmoSv)    !OUT
  
    ! IN = vXPi,vXAs
    !    = MolNumbers of INERT aqu.species in system
    !      (i.e. for 1 kg water or 1 box)
    ! OUT= vLnAct & vLnGam
    
    use M_IOTools
    use M_Numeric_Const,only: TinyDP,Ln10,MinExpDP
    use M_Dtb_Const,only: T_CK
    use M_T_Tpcond, only: T_TPCond
    use M_T_Species,only: T_Species
    use M_T_SolModel
    use M_Safe_Functions
    !--
    implicit none
    !
    type(T_Species), intent(in):: vSpcAq(:)
    type(T_SolModel),intent(in):: SolModel
    real(dp),intent(in)   :: TdgK,Pbar
    integer, intent(in)   :: isW
    logical, intent(in)   :: vLAx(:)  !mobile aqu'species
    real(dp),intent(inout):: vMolF(:)  !mole numbers
    real(dp),intent(inout):: vLnAct(:) 
    real(dp),intent(inout):: vLnGam(:)
    logical, intent(out)  :: vTooLow_(:)
    real(dp),intent(out)  :: OsmoSv !solvent's osmotic coeff'
    !vMolF=  input for inert species, output for mobile species
    !vLnAct= activities, input for mobile species, output for others
    !activity coeff', input for mobile species, output for all (includ'g mobile)
    !
    !! type(T_Species),dimension(:),allocatable:: vSpcSolut
    !! real(dp),       dimension(:),allocatable:: vMolal,vLnGamSolut,vA0
    !! integer,        dimension(:),allocatable:: vZSp
    !! real(dp):: TdgK
    !! real(dp):: LnActSv
    real(dp):: MolWeitSv
    integer :: nAq,iAq,J
    !  
    !--
    ! compute molality for external species and test vTooLow_
    !--
    call Solmodel_Calc_Physics( &
    & TdgK,Pbar,        & !IN  ! temperature and pressure
    & SolModel,         & !IN  ! solvent phase model
    & vSpcAq,           & !IN  ! vector of Aqueous Species
    & isW,              & !IN  ! index of solvent
    & vMolF,            & !IN  ! mole numbers
    !--------------------------
    & vLnGam,           & !OUT ! gamma coefficient
    & vLnAct,           & !OUT ! activity
    & OsmoSv )            !OUT ! osmotic coefficient 
    !---
    nAq = size(vMolF)
    MolWeitSv = vSpcAq(isW)%WeitKg
    !
    !------------ vMolF of External Species
    !
    J= 0 
    do iAq=1,nAq
      if (.not.(iAq==isW)) then
        J= J+1
        if(vLAx(iAq)) then 
          vMolF(iAq)= FSafe_Exp(vLnAct(iAq)-vLnGam(iAq)) &
          &         * vMolF(isW)*MolWeitSv
        end if
      end if
    end do
    !
    !------------ vTooLow Logical
    !
    vTooLow_ = .false.
    J= 0 
    do iAq=1,nAq
      if (.not.(iAq==isW)) then
        if( ( vLnAct(iAq) - vLnGam(iAq) ) <= MinExpDP) then
          vTooLow_(iAq)=.true.
        else 
          vTooLow_(iAq)=.false.
        end if
      end if
    end do
    
  end subroutine Solmodel_Calc_System
  !---
  subroutine Solmodel_Calc_Physics( &
  & TdgK,Pbar,        & !IN  ! temperature and pressure
  & SolModel,         & !IN  ! sol'model activ'model
  & vSpcAq,           & !IN  ! vector of Aqueous Species
  & isW,              & !IN  ! index of solvent
  & vMolF,            & !IN  ! mole numbers
  !--------------------------
  & vLnGam,           & !OUT ! gamma coefficient
  & vLnAct,           & !OUT ! activity
  & OsmoSv )            !OUT ! osmotic coefficient 
  !
  use M_IOTools
  use M_Safe_Functions
  !--
  use M_Solmodel_Calc_Pitzer
  use M_Solmodel_Calc_Debye_Hueckel
  use M_Solmodel_Calc_Davies
  !~ use M_Solmodel_Calc_SIT
  use M_Solmodel_Calc_Dilute
  !--
  use M_T_Species
  use M_T_SolModel
  !----
  implicit none
  !
  real(dp),        intent(in):: TdgK,Pbar
  type(T_SolModel),intent(in):: SolModel
  type(T_Species), intent(in):: vSpcAq(:)
  integer,         intent(in):: isW       
  real(dp),        intent(in):: vMolF(:)  ! mole numbers
  real(dp),intent(out)  :: vLnGam(:) ! gamma coefficient
  real(dp),intent(out)  :: vLnAct(:) ! activity
  real(dp),intent(out)  :: OsmoSv    ! osmotic coefficient
  !---
  integer :: nAq, nSolute, iAq, iSolute
  real(dp), dimension(:),allocatable:: vMolal
  real(dp), dimension(:),allocatable:: vLnGamSolute
  real(dp), dimension(:),allocatable:: vA0
  integer,  dimension(:),allocatable:: IdxSolute
  integer,  dimension(:),allocatable:: vZSp
  !---
  real(dp):: MolWeitSv, MassSv
  real(dp):: LnActSv
  real(dp):: dhA, dhB, bDot, Rho, Eps
    !
    !---
    !// Retrieve Debye-Hueckel Coeffcients
    dhA=  SolModel%Dat%dhA
    dhB=  SolModel%Dat%dhB
    bDot= SolModel%Dat%bDot
    
    Rho=  SolModel%Dat%Rho
    Eps=  SolModel%Dat%Eps
    !
    !// Build solutes indices   
    nAq = size(vMolF)
    nSolute = nAq - 1 
    allocate(IdxSolute(nSolute))
    !
    iSolute = 0
    do iAq = 1, nAq
      if (.not.(iAq==isW)) then
        iSolute = iSolute+1
        IdxSolute(iSolute) = iAq
      end if
    end do
    
    !// Compute solvent properties
    MolWeitSv = vSpcAq(isW)%WeitKg
    MassSv = MolWeitSv * vMolF(isW)
    OsmoSv = One ! default value

    !// Compute solute properties
    allocate(vMolal(1:nSolute))
    vMolal(1:nSolute) = vMolF(IdxSolute(:))/MassSv
    !
    allocate(vZSp(1:nSolute))
    vZSp(1:nSolute)= vSpcAq(IdxSolute(:))%Z
    !
    allocate(vA0(1:nSolute))
    vA0(:)= vSpcAq(IdxSolute(:))%AquSize
    where(vA0(:)==Zero .and. vZSp(:)/=0) vA0(:)= 3.72D0 ! provisional !!!
    !
    allocate(vLnGamSolute(1:nSolute))
    vLnGamSolute = Zero

    !// Compute activities
    
    ! "IDEAL  ",              & ! 1
    ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
    ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
    ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
    ! "PITZER ", "SAMSON ",   & ! 8, 9
    ! "HKF81  ", "SIT    ",   & ! 10,11
    ! "name12 "               & ! 12
    
    ! select case(SolModel%ActModel)
    select case(SolModel%iActModel)

    case(1) ! ("IDEAL")
      call Solmodel_Calc_Dilute( &
      !! & vMolal,MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv) 

    case(8) !("PITZER")
      call Solmodel_Calc_Pitzer( &
      & vSpcAq(IdxSolute), &
      & MolWeitSv, &
      & Rho,Eps,TdgK,vMolal,&
      & vLnGamSolute, LnActSv, OsmoSv)

    case(2,3,4,5) !("DH1","DH2","DH1EQ3","DH2EQ3")
      call Solmodel_Calc_Debye_Hueckel( &
      & TdgK, Pbar, vZSp, vMolal, &
      & dhA, dhB, vA0, bDot, SolModel%iActModel, MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv)

    case(6,7,9) !("DAV_1","DAV_2","SAMSON")
      call Solmodel_Calc_Davies( &
      & vZSp, vMolal, &
      & dhA, SolModel%iActModel, MolWeitSv, &
      & vLnGamSolute, LnActSv, OsmoSv)

    !~ case(11) !("SIT")
      !~ call Solmodel_Calc_SIT( &
      !~ & vZSp, vMolal, &
      !~ & dhA, SolModel%ActModel, MolWeitSv, &
      !~ & vLnGamSolute, LnActSv, OsmoSv)
    
    end select

    !// Store activity and gamma for solute species 
    vLnGam(IdxSolute(:)) = vLnGamSolute(1:nSolute)
    vLnAct(IdxSolute(:)) = FSafe_vLog(vMolal(1:nSolute)) + vLnGamSolute(1:nSolute) 

    !// Store activity and gamma for solvent 
    vLnAct(isW)= LnActSv
    vLnGam(isW)= LnActSv + log( One + MolWeitSv *SUM(vMolal(:)) )

    deallocate(vMolal)
    deallocate(vLnGamSolute)
    deallocate(vZSp)
    deallocate(vA0)

  end subroutine Solmodel_Calc_Physics

end module M_Solmodel_Calc_System
