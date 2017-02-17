module M_SolModel_Calc
!--
!-- calculations of activity coefficients
!-- of solutes and solvent in aqueous solution
!-- to be applicable to any molality based solution model
!--
  use M_Kinds
  use M_Trace,only: Stop_,fTrc,T_,iDebug
  use M_T_DtbLogKTbl,only: DimLogK_Max
  !
  implicit none
  !
  private
  !
  public:: Solmodel_CalcGamma
  
contains

subroutine Solmodel_CalcGamma( &
& TdgK,Pbar,        & !IN
& SolModel,         & !IN
& isW,vSpcAq,vLAx,  & !IN
& vMolF,            & !IN
& vLnAct,vLnGam,    & !INOUT
& vTooLow_,OsmoSv) !OUT

  use M_Solmodel_Calc_System
  
  use M_T_Species,only: T_Species
  use M_T_SolModel,only: T_SolModel
  !--
  type(T_Species),intent(in):: vSpcAq(:)
  type(T_SolModel),intent(in):: SolModel
  real(dp),intent(in)   :: TdgK,Pbar
  integer, intent(in)   :: isW
  logical, intent(in)   :: vLAx(:)  
  real(dp),intent(inout):: vMolF(:) 
  real(dp),intent(inout):: vLnAct(:) 
  real(dp),intent(inout):: vLnGam(:)
  logical, intent(out)  :: vTooLow_(:)
  real(dp),intent(out)  :: OsmoSv
  !---
  logical :: OkNew = .false. ! .true. => switch to new version
  !-----------

  if (OkNew) then 
    call Solmodel_Calc_System(&
    & TdgK,Pbar,        & !IN
    & SolModel,         & !IN
    & isW,vSpcAq,vLAx,  & !IN
    & vMolF,            & !IN
    & vLnAct,vLnGam,    & !INOUT
    & vTooLow_,OsmoSv) !OUT
  else
    call Solmodel_CalcGamma_Old(&
    & TdgK,Pbar,        & !IN
    & SolModel,         & !IN
    & isW,vSpcAq,vLAx,  & !IN
    & vMolF,            & !IN
    & vLnAct,vLnGam,    & !INOUT
    & vTooLow_,OsmoSv) !OUT
  end if

end subroutine Solmodel_CalcGamma

subroutine Solmodel_CalcGamma_Old( &
!--
!-- WARNING
!-- when implementing a new activity model,
!-- add its name in the list of activity models
!--
& TdgK,Pbar,        & !IN
& SolModel,         & !IN
& isW,vSpcAq,vLAx,  & !IN
& vMolF,            & !IN
& vLnAct,vLnGam,    & !INOUT
& vTooLow_,OsmoSv) !OUT

!.IN = vXPi,vXAs=  MolNumbers of INERT aqu.species in system
!                  (i.e. for 1 kg water or 1 box)
!.OUT= vLnAct & vLnGam

  use M_IOTools
  use M_Numeric_Const,only: TinyDP,Ln10,MinExpDP
  use M_Dtb_Const, only: T_CK
  use M_T_Tpcond,  only: T_TPCond
  use M_T_Species, only: T_Species
  use M_T_SolModel,only: T_SolModel
  
  use M_Solmodel_Calc_Pitzer
  use M_Solmodel_Calc_Debye_Hueckel
  use M_Solmodel_Calc_Davies
  !~ use M_Solmodel_Calc_SIT
  use M_Solmodel_Calc_Dilute
  
  type(T_Species), intent(in):: vSpcAq(:)
  type(T_SolModel),intent(in):: SolModel
  
  real(dp),intent(in)   :: TdgK,Pbar
  integer, intent(in)   :: isW
  logical, intent(in)   :: vLAx(:)   !mobile aqu'species
  
  real(dp),intent(inout):: vMolF(:)  !mole numbers
  real(dp),intent(inout):: vLnAct(:) 
  real(dp),intent(inout):: vLnGam(:)
  
  logical, intent(out)  :: vTooLow_(:)
  real(dp),intent(out)  :: OsmoSv !SolModel's osmotic coeff'
  
  ! vMolF=  input for inert species, output for mobile species
  ! vLnAct= activities, input for mobile species, output for others
  ! activity coeff', input for mobile species, output for all (includ'g mobile)
  
  type(T_Species),dimension(:),allocatable:: vSpcSolut
  real(dp),       dimension(:),allocatable:: vMolal,vLnGamSolut,vA0
  integer,        dimension(:),allocatable:: vZSp
  
  real(dp):: LnActSv
  real(dp):: MolWeitSv,Rho,Eps,dhA,dhB,bDot !,Pbar
  integer :: nAq,iAq,N,J !,nAi
  integer :: iModel
  
  MolWeitSv= vSpcAq(isW)%WeitKg
  
  nAq= size(vMolF)
  !nAi= nAq-nAx !nr of non mobile au'species
  vTooLow_=.false.
  
  dhA= SolModel%Dat%dhA
  dhB= SolModel%Dat%dhB
  bDot=SolModel%Dat%bDot
  Rho= SolModel%Dat%Rho
  Eps= SolModel%Dat%Eps
  
  N= count(vSpcAq(:)%Typ=="AQU") - 1 !mod 19/06/2008 09:38
  !-> solute species only
  allocate(vSpcSolut(N))
  allocate(vMolal(N))
  allocate(vLnGamSolut(N))
  allocate(vA0(N),vZsp(N))
  
  ! build arrays for solute species:
  ! vMolal (= molalities), vSpcSolut, vLnGamSolut
  ! Molality(iAq)=MoleNumber(iAq)/MassSolvent; vMolal(isW)=1/MolWeitSv=55.51
  ! vSpcSolut is used mainly as input for Pitzer calculations
  ! (separation solute / SolModel)
  J= 0
  do iAq=1,size(vSpcAq)
    if(vSpcAq(iAq)%Typ=="AQU" .and. vSpcAq(iAq)%NamSp/="H2O") then
      J= J+1
      vSpcSolut(J)= vSpcAq(iAq)
      vLnGamSolut(J)= vLnGam(iAq)
      !! vMolalSolut(J)= vMolal(iAq)
      if(vLAx(iAq)) then 
        vMolal(J)= & !for mobile aqu'species, activity=IN
        & exp(vLnAct(iAq) - vLnGam(iAq))
        vMolF(iAq)= vMolal(J) *vMolF(isW) *MolWeitSv
      else !for inert aqu'species; MolWeitSv=MolWeitH2O
        vMolal(J)= vMolF(iAq) &
        & /vMolF(isW) &
        & /MolWeitSv
      end if
    end if
  end do
  
  vZSp(:)= vSpcSolut(:)%Z
  vA0(:)=  vSpcSolut(:)%AquSize
  where(vA0(:)==Zero .and. vZsp(:)/=0) vA0(:)= 3.72D0 ! provisional !!!!
  where(vZsp(:)==0) vA0(:)= Zero
  
  vLnGam=Zero
  OsmoSv=One
  
  ! "IDEAL  ",              & ! 1
  ! "DH1    ", "DH1EQ3 ",   & ! 2, 3
  ! "DH2    ", "DH2EQ3 ",   & ! 4, 5
  ! "DAV_1  ", "DAV_2  ",   & ! 6, 7
  ! "PITZER ", "SAMSON ",   & ! 8, 9
  ! "HKF81  ", "SIT    ",   & ! 10,11
  ! "name12 "               & ! 12
  
  !~ sModel= trim(SolModel%ActModel)
  iModel=SolModel%iActModel 
  
  select case(iModel)

  case(1) !("IDEAL", actually should be "DILUTE" !!)
    call Solmodel_Calc_Dilute( &
    & vLnGamSolut, LnActSv, OsmoSv)
     
  case(2,3,4,5,10) !("DH1","DH2","DH1EQ3","DH2EQ3","HKF81")
    call Solmodel_Calc_Debye_Hueckel( &
    & TdgK, Pbar, vZSp, vMolal, &
    & dhA, dhB, vA0, bDot, iModel, MolWeitSv, &
    & vLnGamSolut, LnActSv, OsmoSv)
    
  case(8) !("PITZER")
    call Solmodel_Calc_Pitzer( &
    & vSpcSolut, &
    & MolWeitSv, &
    & Rho,Eps,TdgK,vMolal,&
    & vLnGamSolut,LnActSv,OsmoSv)
    
  case(6,7,9) !("DAV_1","DAV_2","SAMSON")
    call Solmodel_Calc_Davies( &
    & vZSp, vMolal, &
    & dhA, iModel, MolWeitSv, &
    & vLnGamSolut, LnActSv, OsmoSv)
  
  !~ case("SIT")
    !~ call Solmodel_Calc_SIT( &
    !~ & vZSp, vMolal, &
    !~ & dhA, sModel, MolWeitSv, &
    !~ & vLnGamSolut, LnActSv, OsmoSv)
    
  end select
  
  !-------- back substitute solute species from vLnGamSolut in vLnGam --
  J= 0
  do iAq=1,size(vSpcAq)
    if(vSpcAq(iAq)%Typ=="AQU" .and. vSpcAq(iAq)%NamSp/="H2O") then
      J= J+1
      vLnGam(iAq)= vLnGamSolut(J)
      !
      if(vLAx(iAq)) then 
        !-- MOBILE SPECIES --
        !-- -> apply act'coeff' to calculate #mole from activity
        if(vLnAct(iAq) -vLnGam(iAq) <=MinExpDP) then
          vTooLow_(iAq)=.true.
          vMolal(J)= exp(MinExpDP)
        else 
          vTooLow_(iAq)=.false.
          vMolal(J)= exp(vLnAct(iAq) -vLnGamSolut(J))
        end if
        vMolF(iAq)= vMolal(J) *vMolF(isW) *MolWeitSv
      else
        !-- INERT SPECIES --
        !-- -> compute activity
        if(vMolal(J) <=TinyDP) then
          vTooLow_(iAq)=.true.
          vLnAct(iAq)= log(TinyDP)
        else
          vTooLow_(iAq)=.false.
          vLnAct(iAq)= log(vMolal(J)) + vLnGamSolut(J)
        end if
      end if
    end if
  end do
  !---------------------------------------------------/solute species--
  
  !--------------------------------------------------------- SolModel --
  vLnAct(isW)= LnActSv
  !if(sModel /= "IDEAL") then
  if(iModel /= 1) then
    vLnGam(isW)= LnActSv + log( One + MolWeitSv *SUM(vMolal(:)) )
  else
    vLnGam(isW)= Zero
  end if
  ! for SolModel, Gamma is the rational activity coeff', generally called Lambda
  !--------------------------------------------------------/ SolModel --
  
  deallocate(vSpcSolut)
  deallocate(vMolal)
  deallocate(vLnGamSolut)
  deallocate(vA0,vZsp)

end subroutine Solmodel_CalcGamma_Old

end module M_SolModel_Calc

