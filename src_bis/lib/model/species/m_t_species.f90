module M_T_Species
  use M_Kinds
  !
  implicit none
  private
  !
  public:: T_Species, T_SpcData, T_SpeciesDtb
  public:: Species_Zero
  public:: Species_Index
  public:: Species_Rename
  public:: Species_Stoikio
  public:: Species_Stoikio_Calc
  public:: Species_EntropyZero
  !
  public:: nElMax
  !
  type:: T_SpcData !container for storing a species' data
    real(dp):: &   !store current values
    & Mole, &      !mole number in reference system
    !              (either in 1 kg H2O, either in the box volume, etc.)
    & LAct, &      !ln(activity)
    & LGam         !ln(activ'coeff)
  end type T_SpcData
  !
  type:: T_SpeciesThermo
  ! container for storing current thermo data of a SPECIES
  ! for future use ??
    real(dp):: G0rt,V0,H0,S0,Cp0,LnPhi,logK
  end type T_SpeciesThermo
  !
  integer,parameter:: nElMax= 64
  !
  type T_SpeciesDtb
    ! pointer to a species in a database-specific array
    character(len=7) :: DtbModel 
    !selector: "AQU_HKF","MIN_HKF","LOGKTBL","LOGKANL",etc
    integer          :: Indx     
    !index of species in the corresponding vDtbEoS
    character(len=15):: DtbTrace 
    !retrace origin of data, herited from vDtb*%Num
  end type T_SpeciesDtb
  !
  type:: T_Species !container for storing a SPECIES' structure
    character(len=3 ):: Typ   !AQU/MIN/GAS/XCH
    character(len=23):: NamSp !
    !
    integer:: iDtb
    integer:: iDiscret
    ! iDtb/=0
    ! -> species pointing to a pure species in the thermodynamic database
    ! iDtb gives index of species in vSpcDtb
    ! iDiscret/= 0 -> species generated by discretization of a mixture
    ! iDiscret gives index of species in vDiscretParam
    !
    character(len=71):: Formula
    integer:: vStoikio(0:nElMax)
    integer:: Div 
    !"formula factor", i.e. common denominator of all element numbers
    !
    integer:: Z   !charge
    ! Note: actual charge is Z/real(Div) ???
    !
    real(dp):: &  !
    & WeitKg,  &  !molecule weight (formula weight), in Kg
    & AquSize, &  !for aqu.species, size parameter for Debye-Huckel relations
    !! & Vs,Ss,Cps,Hs,Gs, &    !solvation terms, for aqu'species, at current (T,P) 
    & G0rt,V0,H0,S0,Cp0,LnFug
    !-> values at current (T,P) !! (NOT those at ref'conditions)
    !G0rt=G/RT of formatION 
    !-> will produce deltaG's of ASSOCIATION/PRECIPITATION
    !caveat:
    !logK database contain, in lieu of G0, the log10K 
    !-> G0/RT= - log10K *ln10
    !
    type(T_SpcData):: Dat !store values of vMolF, vLnAct, vLnGam, ..
    !
  end type T_Species
  !
  ! the field "phase" could be added in order to deal with 
  ! vMulti-solution systems, 
  ! e.g. equilibrium of aqueous phase with a gas mixture or a solid solution
  !
  ! currently, the phase is associated through the iPhase index
  !
  ! in simple cases, originally considered, of equil. between aqu. solution
  ! and pure minerals or gases, the field "typ" is sufficient
  !
  ! NB1: what is referred to as "Gibbs energies" is actually G/RT, 
  ! i.e. -LnK
  ! NB2: if the database is in "logK", it is generally log10, 
  ! then LnK=Ln10*logK
  ! deltaG's for sec'species and minerals are for formatION reactions:
  !   Prim'Species -> Second'Species
  ! -> DG0(iSc)=  vSpc(iSc)%G0 - dot_product(tStNu(1:nCp,iSc),vSpc(1:nCp)%G0)
  ! each time basis is changed, DG0 must be re-calculated
  ! -> thus, DG0 not stored in vSpc
  !
  logical,parameter:: Use_Synonym= .true.
  
contains

integer function Species_Index(Str,V) !,nSpc_)
!--
!-- position of species with %NamSp--Str in vSpc_(1:nSpc_) --
!--
  character(len=*),intent(in):: Str
  type(T_Species), intent(in):: V(:)
  !
  integer::I
  !
  Species_Index=0
  if(size(V)==0) return
  I=0
  do
    I=I+1
    !if(trim(Str)==trim(V(I)%NamSp)) then
    if(Synonym(trim(Str),trim(V(I)%NamSp))) then
      Species_Index=I
      exit
    end if
    if(I==size(V)) exit
  end do
  !->> if Str not found, I=0
  !
  return
end function Species_Index

logical function Synonym(S1,S2)
  character(len=*),intent(in):: S1,S2
  !
  integer:: J1,J2
  character(len=4):: S
  character(len=23):: SS1,SS2
  !
  Synonym= (trim(S1)==trim(S2))
  if(Synonym) return
  if(.not. Use_Synonym) return
  !
  SS1= trim(S1)
  SS2= trim(S2)
  !
  !--- to consider ",AQ" and "(AQ)" as synonyms
  S= ",AQ"
  J1= INDEX(SS1,trim(S))  ;  if(J1>1) SS1= SS1(1:J1-1)//"_AQ"
  J2= INDEX(SS2,trim(S))  ;  if(J2>1) SS2= SS2(1:J2-1)//"_AQ"
  S= "(AQ)"
  J1= INDEX(SS1,trim(S))  ;  if(J1>1) SS1= SS1(1:J1-1)//"_AQ"
  J2= INDEX(SS2,trim(S))  ;  if(J2>1) SS2= SS2(1:J2-1)//"_AQ"
  if(trim(SS1)==trim(SS2)) then
    Synonym= .true.
    return
  end if
  !
  !--- to consider ",G" and "(G)" as synonyms
  S= ",G"
  J1= INDEX(SS1,trim(S))  ;  if(J1>1) SS1= SS1(1:J1-1)//"_G"
  J2= INDEX(SS2,trim(S))  ;  if(J2>1) SS2= SS2(1:J2-1)//"_G"
  S= "(G)"
  J1= INDEX(SS1,trim(S))  ;  if(J1>1) SS1= SS1(1:J1-1)//"_G"
  J2= INDEX(SS2,trim(S))  ;  if(J2>1) SS2= SS2(1:J2-1)//"_G"
  if(trim(SS1)==trim(SS2)) then
    Synonym= .true.
    return
  end if
  !
  return
end function Synonym

subroutine Species_Rename(sName)
  character(len=*),intent(inout):: sName
  !
  integer:: J
  character(len=4):: S
  
  if(Use_Synonym) return
  !
  !--- to consider ",AQ" and "(AQ)" as synonyms
  S= ",AQ"
  J= INDEX(sName,trim(S))
  if(J>1) sName= sName(1:J-1)//"_AQ"
  S= "(AQ)"
  J= INDEX(sName,trim(S))
  if(J>1) sName= sName(1:J-1)//"_AQ"
  
  !--- to consider ",G" and "(G)" as synonyms
  S= ",G"
  J= INDEX(sName,trim(S))
  if(J>1) sName= sName(1:J-1)//"_G"
  S= "(G)"
  J= INDEX(sName,trim(S))
  if(J>1) sName= sName(1:J-1)//"_G"
  
end subroutine Species_Rename

subroutine SpeciesDtb_Zero(S)
  type(T_SpeciesDtb),intent(out):: S
  !
  S%DtbModel=  "none"
  S%Indx=      0
  S%DtbTrace=  "none"
  !
  return
end subroutine SpeciesDtb_Zero

subroutine Species_Zero(S)
  type(T_Species),intent(out):: S
  !
  S%Typ=     "Z"
  S%NamSp=   "Z"
  S%Formula= "Z"
  !! S%DtbTrace="Z"
  !! S%DtbModel="Z"
  S%iDtb=     0
  S%iDiscret= 0
  S%Z=        0
  S%Div=      1
  S%WeitKg=   Zero
  S%AquSize=  3.72D0 !! default aqu'size
  !
  S%vStoikio(:)= 0 != S%vStoikio(1:nEl), stoichiometry in integer values
  S%vStoikio(1)= 1 != the common divisor of the chemical formula
  !
  S%G0rt=  Zero
  S%V0=    One
  S%H0=    Zero
  S%S0=    Zero
  S%Cp0=   Zero
  S%LnFug= Zero
  !
  !! S%Rho=      One
  !!S%S0Ele= Zero
  !S%Vs,S%Ss,S%Cps,S%Hs,S%Gs
  !S%G0rt,S%V0,S%H0,S%S0,S%Cp0,S%logK 
end subroutine Species_Zero

subroutine Species_Stoikio_Calc(vEle,vSpc,Ok)
  use M_T_Element,  only: T_Element,Element_Index
  !
  type(T_Element), intent(in)   :: vEle(:)
  type(T_Species), intent(inout):: vSpc(:) !modifs in nDiv,Z,Oxy
  logical,         intent(out)  :: Ok
  !
  integer:: I,ieOx
  logical:: fOk
  !
  ieOx= Element_Index("OX_",vEle)
  Ok= .true.
  !
  do I=1,size(vSpc)
    !if(vSpc(i)%iDiscret==0) then
    if(vSpc(i)%iDtb>0) then
      ! for discretized species,
      ! use DiscretSpecies_Stoikio_Calc
      call Species_Stoikio(vEle,ieOx,vSpc(I),fOk)
      if(.not. fOk) then
        Ok= .false.
        !! Msg= "Stoikiometry problem in species "//trim(vSpc(I)%NamSp)
        return
      end if
    end if
  end do
  !
  return
end subroutine Species_Stoikio_Calc

subroutine Species_Stoikio(vEle,ieOx,S,fOk)
!--
!-- build stoichiometry vector of species S according to element base vEle
!-- NEW: stoichiometry saved as S%vStoikio
!--
  use M_T_Element,only: T_Element,Formula_Read
  !
  type(T_Element),intent(in)::    vEle(:)
  integer,        intent(in)::    ieOx
  type(T_Species),intent(inout):: S
  logical,        intent(out)::   fOk
  !
  integer:: vStoik(1:size(vEle))
  integer:: ZSp,nDiv,nEl,ZExc
  !
  nEl= size(vEle)
  S%vStoikio(:)= 0
  !
  !write(71,'(2A)') S%NamSp,S%Formula
  call Formula_Read(S%Formula,vEle,ZSp,nDiv,fOk,vStoik)
  !
  if(.not. fOk) return
  !
  S%vStoikio(1:nEl)= vStoik(1:nEl)
  S%vStoikio(0)=     nDiv
  S%vStoikio(nEl+1)= ZSp
  !
  S%Div= nDiv !-> "formula factor"
  S%Z=   ZSp  !-> species charge
  !
  ZExc=dot_product(vStoik(1:nEl),vEle(1:nEl)%Z) !-> oxidation state
  !
  if(ieOx/=0) then
    S%vStoikio(ieOx)= ZSp - ZExc
    fOk=.true.
  !<added 2010May
  else
    fOk= (Zsp==ZExc)
  !</
  end if
  !
end subroutine Species_Stoikio

real(dp) function Species_EntropyZero(vEle,S)
  use M_Dtb_Const,only: S0_Hydrogen
  use M_T_Element,only: T_Element
  !
  type(T_Element),intent(in):: vEle(:)
  type(T_Species),intent(in):: S
  !
  integer:: N
  !
  N= size(vEle)
  !
  Species_EntropyZero= &
  & ( dot_product(S%vStoikio(1:N),vEle(1:N)%S0)  &
  &   - S%vStoikio(N+1) *S0_Hydrogen           ) &
  & /real(S%vStoikio(0))
  !
end function Species_EntropyZero

real(dp) function Species_WeitKg(vEle,S)
  use M_T_Element,only: T_Element
  type(T_Element),intent(in):: vEle(:)
  type(T_Species),intent(in):: S
  !
  integer:: nEl
  !
  nEl= size(vEle)
  Species_WeitKg= &
  & sum(vEle(1:nEl)%WeitKg * S%vStoikio(1:nEl)) /real(S%vStoikio(0))
  !
end function Species_WeitKg

end module M_T_Species
