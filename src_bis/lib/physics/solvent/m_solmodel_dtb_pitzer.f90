module M_Solmodel_Pitzer_Dtb
!--
!-- Module routines associées à Pitzer
!-- Livraison 10/2006: Nicolas Ferrando
!-- Remoduling : Anthony Michel
!--
!-- modified : J.Moutte, 11/2006
!-- -> M_Pitzer_Dtb-  database reading, and TP computations
!-- -> M_Pitzer_Calc- calculations
!--
  use M_Kinds
  use M_Trace,    only: fTrc,iDebug,T_,Stop_
  use M_T_Species,only: T_Species
  implicit none
  !
  private
  !
  !-- public --
  public:: Solmodel_Pitzer_Dtb_Init
  public:: Solmodel_Pitzer_Dtb_TPUpdate
  public:: Solmodel_Pitzer_Dtb_Clean
  public:: Solmodel_Pitzer_Dtb_TPtest
  !
  real(dp),dimension(:,:),  allocatable,public:: tBeta0,tBeta1,tBeta2
  real(dp),dimension(:,:),  allocatable,public:: tCPhi,tTheta,tLamda
  real(dp),dimension(:,:),  allocatable,public:: tAlfa1,tAlfa2
  real(dp),dimension(:,:,:),allocatable,public:: tZeta,tPsi
  !
  !-- private --
  type(T_Species),allocatable:: vSolut(:)
  !
  integer,dimension(:,:),  allocatable:: tI_Beta0,tI_Beta1,tI_Beta2
  integer,dimension(:,:),  allocatable:: tI_CPhi,tI_Theta,tI_Lamda
  integer,dimension(:,:,:),allocatable:: tI_Zeta,tI_Psi
  !
  real(dp),dimension(:),allocatable:: vBeta0,vBeta1,vBeta2
  real(dp),dimension(:),allocatable:: vCPhi,vTheta,vLamda
  real(dp),dimension(:),allocatable:: vZeta,vPsi
  real(dp),dimension(:),allocatable:: vAlfa1,vAlfa2
  !
  ! usage:
  !  k= tI_Beta0(i,j)
  !  if(k/=0) then
  !    tBeta0(i,j)= vBeta0(k)

  ! function CoeffAtTP()
  !    Res= T%X0
  !    &  + T%X3 *(TdgK -Tref) &
  !    &  + T%X1 *(One/TdgK -One/Tref) &
  !    &  + T%x2 *log(TdgK/Tref)       &
  !    &  + T%x4 *(TdgK*TdgK -Tref*Tref)

  !~ !! former version: static fitting parameter
  !~ type:: T_Fitt
    !~ integer:: iModel
    !~ real(dp):: vX(1:8)
  !~ end type T_Fitt

  type:: T_Fitt
    character(len=30):: namModel
    integer:: iModel ! model selector
    integer:: N      ! dimension of vX (model- dependent)
    real(dp),allocatable:: vX(:)
  end type T_Fitt
  !
  type(T_Fitt),dimension(:),allocatable:: vF_Beta0,vF_Beta1,vF_Beta2
  type(T_Fitt),dimension(:),allocatable:: vF_CPhi,vF_Theta,vF_Lamda
  type(T_Fitt),dimension(:),allocatable:: vF_Zeta,vF_Psi

contains

subroutine Solmodel_Pitzer_Dtb_Init(vSpc)
  use M_T_Species, only: T_Species
  type(T_Species),intent(in):: vSpc(:)
  !
  logical:: Error
  character(len=80):: ErrorMsg
  !
  call Solmodel_Pitzer_Dtb_Clean
  !
  call Solmodel_Pitzer_Dtb_Read(vSpc,Error,ErrorMsg)
  !
  if(Error) &
  & call Stop_("Error in Solmodel_PitzerDtb_Read:"//trim(ErrorMsg))
  !
  !call Solmodel_PitzerDtb_TPUpdate(TdgK)
  !
end subroutine Solmodel_Pitzer_Dtb_Init

subroutine Solmodel_Pitzer_Dtb_Read( &
& vSpc,           &
& Error,ErrorMsg  &
& )
!--
!-- read interaction coeffs from the Pitzer database
!-- must be called AFTER vSpc has been built
!-- -> vPitz will be a subset of vSpc
!--
  use M_IoTools
  use M_Files,        only: NamFPtz
  use M_Numeric_Const,only: ln10
  use M_T_Species,    only: T_Species,Species_Index
  !
  type(T_Species),intent(in):: vSpc(:)
  !
  logical,         intent(out):: Error
  character(len=*),intent(out):: ErrorMsg
  !
  real(dp),parameter:: Tref= 298.15D0
  !
  logical:: L_ScanSpecies
  !
  character(len=512):: L,W,W2
  logical           :: EoL
  !
  ! character(len=7):: FittFormat
  character(len=30):: namModel
  !
  logical :: Ok
  integer :: vISpc(3)
  integer :: iSp1,iSp2,iSp3,nSl,I,J,N
  integer :: F,ios,ff
  integer :: iModel
  integer :: iBeta0,iBeta1,iBeta2
  integer :: iCPhi,iTheta,iLamda
  integer :: iPsi,iZeta
  real(dp):: vXread(dimV) !,X0,vX(dimV)
  type(T_Fitt):: Coeff
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Solmodel_PitzerDtb_Read"
  !
  call Solmodel_Pitzer_Dtb_Clean
  ! Pitzer_Dtb_Init_Done=.false.
  !
  L_ScanSpecies= .false.
  !
  Error= .false.
  ErrorMsg= "OK"
  !
  !!  if(present(TdgK)) then; T=TdgK
  !!  else;                   T=298.15D0
  !!  end if
  !!  if(present(Pbar)) then; P=Pbar
  !!  else;                   P=1.01325D0 !1 atm in bars
  !!                          !caveat: some functions may use Pascal ?...
  !!  end if
  !!
  !! nSl= size(vSpc)
  !
  ! nSl is number of solute species
  nSl= count(vSpc(:)%Typ=="AQU") - 1 !mod 19/06/2008 09:38
  !
  allocate(vSolut(1:nSl))
  !
  J= 0
  do I=1,size(vSpc)
    if(vSpc(I)%Typ=="AQU" .and. vSpc(I)%NamSp/="H2O") then
      J= J+1
      vSolut(J)= vSpc(I)
    end if
  enddo
  !
  iModel= 0
  namModel= "none"
  !
  allocate(tBeta0(1:nSl,1:nSl))      ;  tBeta0=Zero
  allocate(tBeta1(1:nSl,1:nSl))      ;  tBeta1=Zero
  allocate(tBeta2(1:nSl,1:nSl))      ;  tBeta2=Zero
  allocate(tCPhi (1:nSl,1:nSl))      ;  tCPhi= Zero
  allocate(tTheta(1:nSl,1:nSl))      ;  tTheta=Zero
  allocate(tLamda(1:nSl,1:nSl))      ;  tLamda=Zero
  allocate(tAlfa1(1:nSl,1:nSl))      ;  tAlfa1=Zero
  allocate(tAlfa2(1:nSl,1:nSl))      ;  tAlfa2=Zero
  allocate(tZeta (1:nSl,1:nSl,1:nSl));  tZeta= Zero
  allocate(tPsi  (1:nSl,1:nSl,1:nSl));  tPsi=  Zero
  !
  allocate(tI_Beta0(1:nSl,1:nSl))    ;  tI_Beta0= 0
  allocate(tI_Beta1(1:nSl,1:nSl))    ;  tI_Beta1= 0
  allocate(tI_Beta2(1:nSl,1:nSl))    ;  tI_Beta2= 0
  allocate(tI_CPhi (1:nSl,1:nSl))    ;  tI_CPhi=  0
  allocate(tI_Theta(1:nSl,1:nSl))    ;  tI_Theta= 0
  allocate(tI_Lamda(1:nSl,1:nSl))    ;  tI_Lamda= 0

  !allocate(tI_Alfa1(1:nSl,1:nSl));        tI_Alfa1=.false.
  !allocate(tI_Alfa2(1:nSl,1:nSl));        tI_Alfa2=.false.

  allocate(tI_Zeta (1:nSl,1:nSl,1:nSl)) ; tI_Zeta= 0
  allocate(tI_Psi  (1:nSl,1:nSl,1:nSl)) ; tI_Psi=  0
  !
  call GetUnit(F)
  !
  !----- reading database, FIRST pass -> read tI_Beta0, tI_Beta1, ... --
  open(F,file=trim(NamFPtz),STATUS='OLD')

  DoFile1: do

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile1
    !
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle doFile1 !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit DoFile1
    !
    case("PITZER")
      !
      doLine1: do
        !
        read(F,'(A)',iostat=ios) L
        if(ios/=0) exit doFile1
        !
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle doLine1 !skip comment lines
        call AppendToEnd(L,W,EoL)

        select case(W)

        case("ENDINPUT");        exit doFile1
        case("ENDPITZER","END"); exit doLine1

        case("format")
          call LinToWrd(L,W2,EoL)
          select case(trim(W2))
          case("NEW")  ;  L_ScanSpecies= .true.
          case("OLD")  ;  L_ScanSpecies= .false.
          case default
            Error= .true.
            ErrorMsg= trim(W2)//"= not valid as format"
            return
          end select

        case("FITTING")
          !call LinToWrd(L,W2,EoL)
          cycle doLine1
        !
        !--------------------------------------- 2-species parameters --
        case("BETA0","BETA1", "BETA2","CPHI", "LAMBDA","THETA")
          if(L_ScanSpecies) then
            call LinToWrd(L,W2,EoL)
            call Species_Scan(W2,2,vISpc,Error)
            if(Error) then
              ErrorMsg= trim(W2)//" -> Error in species list ?"
              return
            else
              iSp1= vISpc(1)
              iSp2= vISpc(2)
            end if
          else
            call LinToWrd(L,W2,EoL)  ;  iSp1= Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL)  ;  iSp2= Species_Index(W2,vSolut)
          end if
          !
          !--- skip the lines with species not in current species set --
          if(iSp1*iSp2==0) cycle doLine1
          !---/
          !
          if(trim(W)=="LAMBDA") then
            if(vSolut(iSp2)%Z/=0) then
              ErrorMsg= "LAMBDA coeff: species 2 must be neutral"
              if(iDebug>2) write(fTrc,'(A,/,A)') &
              & trim(ErrorMsg), &
              & trim(vSolut(iSp1)%NamSp)//"="// &
              & trim(vSolut(iSp2)%NamSp)
              Error= .true.
              return
            end if
          else
            call Sort2sp(iSp1,iSp2,Error)
          end if
          !
          if(Error) &
          & call Stop_("Pitzer database: error in line "//trim(W)//trim(W2))
          !
          !if(iSp1*iSp2>1) then !-> the 2 species are in database
          select case(trim(W))
            case("BETA0")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta0)
            case("BETA1")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta1)
            case("BETA2")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Beta2)
            case("CPHI")      ; Ok= Ok_2Sp(iSp1,iSp2,tI_CPhi )
            case("THETA")     ; Ok= Ok_2Sp(iSp1,iSp2,tI_Theta)
            case("LAMBDA")    ; Ok= Ok_2Sp(iSp1,iSp2,tI_Lamda)
            !! case("Alfa") !! not used
            !!! Alfa is calculated according to Harvie-etal-84
            !!! see also Pitzer, ReviewInMineralogy-17, p104
          end select
            !
          !end if
        !endcase Beta...
        !--------------------------------------/ 2-species parameters --
        !
        !--------------------------------------- 3-species parameters --
        case("PSI","ZETA")
          !
          if(L_ScanSpecies) then
            call LinToWrd(L,W2,EoL)
            call Species_Scan(W2,3,vISpc,Error)
            if(Error) then
              ErrorMsg= trim(W2)//" -> Error in species list ?"
              return
            else
              iSp1= vISpc(1)
              iSp2= vISpc(2)
              iSp3= vISpc(3)
            end if
          else
            call LinToWrd(L,W2,EoL); iSp1=Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL); iSp2=Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL); iSp3=Species_Index(W2,vSolut)
          end if
          !
          if(iSp1*iSp2*iSp3==0) cycle doLine1
          !
          if(trim(W)=="PSI") then
            if(vSolut(iSp1)%Z*vSolut(iSp2)%Z<0) then
              ErrorMsg= "PSI coeff: Species 1 & 2 must be ions same charge"
              if(iDebug>2) write(fTrc,'(A,/,A)') &
              & trim(ErrorMsg), &
              & trim(vSolut(iSp1)%NamSp)//"="// &
              & trim(vSolut(iSp2)%NamSp)//"="// &
              & trim(vSolut(iSp3)%NamSp)
              Error= .true.
              return
            end if
          end if
          !
          if(trim(W)=="ZETA") then
            if(vSolut(iSp1)%Z*vSolut(iSp2)%Z>0) then
              ErrorMsg= "PSI coeff: Species 1 & 2 must be ions opposite charge"
              if(iDebug>2) write(fTrc,'(A,/,A)') &
              & trim(ErrorMsg), &
              & trim(vSolut(iSp1)%NamSp)//"="// &
              & trim(vSolut(iSp2)%NamSp)//"="// &
              & trim(vSolut(iSp3)%NamSp)
              Error= .true.
              return
            end if
          end if
          !
          ! call Sort3sp(iSp1,iSp2,iSp3,Error)
          !
          call Sort2sp(iSp1,iSp2,Error)
          !
          if(Error) &
          & call Stop_("Pitzer database: error in line "//trim(W)//trim(W2))
          !
          select case(W)
            !case("PSI");   tI_Psi (iSp1,iSp2,iSp3)= 1
            !case("ZETA");  tI_Zeta(iSp1,iSp2,iSp3)= 1
            case("PSI")  ;  Ok= Ok_3Sp(iSp1,iSp2,iSp3,tI_Psi )
            case("ZETA") ;  Ok= Ok_3Sp(iSp1,iSp2,iSp3,tI_Zeta)
          end select
          !
        !endcase PSI...
        !--------------------------------------/ 3-species parameters --
        case default
          Error= .true.
          ErrorMsg= trim(W)//" is unknown keyword"
          return
        !
        end select

      enddo doLine1
    !end case("PITZER")
    
    end select

  enddo DoFile1
  !
  close(F)
  !----/ reading database, FIRST pass -> read tI_Beta0, tI_Beta1, ... --

  !------------------------------------------------------ allocations --
  iBeta0= count(tI_Beta0(:,:)==1)
  iBeta1= count(tI_Beta1(:,:)==1)
  iBeta2= count(tI_Beta2(:,:)==1)
  iCPhi=  count(tI_CPhi (:,:)==1)
  iTheta= count(tI_Theta(:,:)==1)
  iLamda= count(tI_Lamda(:,:)==1)
  !
  allocate(vBeta0(iBeta0))  ;  allocate(vF_Beta0(iBeta0))
  allocate(vBeta1(iBeta1))  ;  allocate(vF_Beta1(iBeta1))
  allocate(vBeta2(iBeta2))  ;  allocate(vF_Beta2(iBeta2))
  !
  allocate(vAlfa1(iBeta1))
  allocate(vAlfa2(iBeta2))
  !
  allocate(vCPhi (iCPhi) )  ;  allocate(vF_CPhi (iCPhi) )
  allocate(vTheta(iTheta))  ;  allocate(vF_Theta(iTheta))
  allocate(vLamda(iLamda))  ;  allocate(vF_Lamda(iLamda))
  !
  iBeta0= 0
  iBeta1= 0
  iBeta2= 0
  iCPhi=  0
  iTheta= 0
  iLamda= 0
  !
  iPsi=  count(tI_Psi(:,:,:)==1)
  iZeta= count(tI_Zeta(:,:,:)==1)
  allocate(vPsi(iPsi))     ;  allocate(vF_Psi(iPsi))
  allocate(vZeta(iZeta))   ;  allocate(vF_Zeta(iZeta))
  !
  iPsi=  0
  iZeta= 0
  !------------------------------------------------------/allocations --

  !-------- reading database, SECOND pass -> read vBeta0, vBeta1, ... --
  open(F,file=trim(NamFPtz),STATUS='OLD')
  !
  DoFile2: do

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit doFile2
    !
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle doFile2 !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit doFile2
    !
    case("PITZER")
      !
      doLine2: do
        !
        read(F,'(A)',iostat=ios) L
        if(ios/=0) exit doFile2
        !
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle doLine2 !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        select case(W)
        !
        case("ENDINPUT");        exit doFile2
        case("ENDPITZER","END"); exit doLine2
        !
        case("FITTING")
          call LinToWrd(L,W2,EoL)
          namModel= trim(W2)
          select case(trim(W2))
            case("none")               ;  iModel= 0
            case("PHREEQC")            ;  iModel= 1
            case("KUEHN")              ;  iModel= 2
            case("CHRISTOV-2004")      ;  iModel= 3
            case("PITZER-1984")        ;  iModel= 4
            case("PABALAN-1987")       ;  iModel= 5
            case("POLYA-2007")         ;  iModel= 6
            case("LI-DUAN-2007")       ;  iModel= 7
            case default
              Error= .true.
              ErrorMsg= trim(W2)//"= unknown keyword for FITTING"
              return
          end select
        !
        !------------------------------------- 2-species parameters --
        case("BETA0","BETA1", "BETA2","CPHI", "LAMBDA","THETA")
          !------------------------------- retrieve species indexes --
          if(L_ScanSpecies) then
            call LinToWrd(L,W2,EoL)
            call Species_Scan(W2,2,vISpc,Error)
            if(Error) then
              ErrorMsg= trim(W2)//" -> Error in species list ?"
            else
              iSp1= vISpc(1)
              iSp2= vISpc(2)
            end if
          else
            call LinToWrd(L,W2,EoL)  ;  iSp1= Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL)  ;  iSp2= Species_Index(W2,vSolut)
          end if
          !
          if(iSp1*iSp2==0) cycle doLine2
          !
          call Sort2sp(iSp1,iSp2,Error)
          !--------------------------------/ retrieve species indexes --
          !
        !--------------------------------------- 3-species parameters --
        case("PSI","ZETA")
          !--------------------------------- retrieve species indexes --
          if(L_ScanSpecies) then
            call LinToWrd(L,W2,EoL)
            call Species_Scan(W2,3,vISpc,Error)
            if(Error) then
              ErrorMsg= trim(W2)//" -> Error in species list ?"
              return
            else
              iSp1= vISpc(1)
              iSp2= vISpc(2)
              iSp3= vISpc(3)
            end if
          else
            call LinToWrd(L,W2,EoL); iSp1=Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL); iSp2=Species_Index(W2,vSolut)
            call LinToWrd(L,W2,EoL); iSp3=Species_Index(W2,vSolut)
          end if
          !
          if(iSp1*iSp2*iSp3==0) cycle doLine2
          !
          !! call Sort3sp(iSp1,iSp2,iSp3,Error)
          call Sort2sp(iSp1,iSp2,Error)
          !
          !--------------------------------/ retrieve species indexes --
        case default
          Error= .true.
          ErrorMsg= trim(W)//" is unknown keyword"
        !
        end select
        !
        !-------------------------------------------- read parameters --
        !
        call ReadRValsV(L,N,vXread)
        !
        select case(iModel)
          case(0)  ;  N= 1  ! "none"
          case(1)  ;  N= 5  ! "PHREEQ"
          case(2)  ;  N= 5  ! "KUEHN"
          case(3)  ;  N= 9  ! "CHRISTOV-2004"
          case(4)  ;  N= 21 ! "PITZER-1984", Na+=Cl-
          case(5)  ;  N= 12 ! "PABALAN-1987"
          case(6)  ;  N= 8  ! "POLYA-2007"
          case(7)  ;  N= 11 ! "LI-DUAN-2007"
          case default
            call Stop_("Error on Pitz_Dtb / iModel")
        end select
        !
        allocate(Coeff%vX(N))
        !
        !!--- for static fitting model
        !Coeff%vX(1:8)= vXread(1:8)
        !Coeff%iModel= iModel
        !!---/
        Coeff%N= N
        Coeff%namModel= trim(namModel)
        Coeff%iModel= iModel
        Coeff%vX(1:N)= vXread(1:N)
        !
        select case(W)
      
        !--------------------------------------- 2-species parameters --
        case("BETA0")
          iBeta0= iBeta0 +1
          tI_Beta0(iSp1,iSp2)= iBeta0
          allocate(vF_Beta0(iBeta0)%vX(Coeff%N))
          vF_Beta0(iBeta0)= Coeff
        case("BETA1")
          iBeta1= iBeta1 +1
          tI_Beta1(iSp1,iSp2)= iBeta1
          allocate(vF_Beta1(iBeta1)%vX(Coeff%N))
          vF_Beta1(iBeta1)= Coeff
        case("BETA2")
          iBeta2= iBeta2 +1
          tI_Beta2(iSp1,iSp2)= iBeta2
          allocate(vF_Beta2(iBeta2)%vX(Coeff%N))
          vF_Beta2(iBeta2)= Coeff
        case("CPHI")
          iCPhi= iCPhi +1
          tI_CPhi(iSp1,iSp2)= iCPhi
          allocate(vF_CPhi(iCPhi)%vX(Coeff%N))
          vF_CPhi(iCPhi)= Coeff
        case("THETA")
          iTheta= iTheta +1
          tI_Theta(iSp1,iSp2)= iTheta
          allocate(vF_Theta(iTheta)%vX(Coeff%N))
          vF_Theta(iTheta)= Coeff
        case("LAMBDA")
          iLamda= iLamda +1
          tI_Lamda(iSp1,iSp2)= iLamda
          allocate(vF_Lamda(iLamda)%vX(Coeff%N))
          vF_Lamda(iLamda)= Coeff

        !! case("Alfa") !! not used
        !!! Alfa is calculated according to Harvie-etal-84
        !!! see also Pitzer, ReviewInMineralogy-17, p104

        !--------------------------------------- 3-species parameters --
        case("PSI")
          iPsi= iPsi+1
          !vPsi(iPsi)= X
          tI_Psi(iSp1,iSp2,iSp3)= iPsi
          allocate(vF_Psi(iPsi)%vX(Coeff%N))
          vF_Psi(iPsi)= Coeff
        case("ZETA")
          iZeta= iZeta +1
          tI_Zeta(iSp1,iSp2,iSp3)= iZeta
          allocate(vF_Zeta(iZeta)%vX(Coeff%N))
          vF_Zeta(iZeta)= Coeff

        end select
        !
        deallocate(Coeff%vX)
        !-------------------------------------------/ read parameters --
        !
      enddo doLine2
    !endcase("block")
    
    end select

  enddo doFile2
  !
  close(F)
  !--------/reading database, SECOND pass -> read vBeta0, vBeta1, ... --
  !
  !----------------------------------------------------- Alfa1,tAlfa2 --
  !--- following Pitzer-RevMin17-p104, Harvie-etal-84
  !--- cf CHRISTOV-2004-Moeller-GCA68-3718
  do iSp1=1,nSl
    do iSp2=1,nSl
      if(tI_Beta1(iSp1,iSp2)>0) then
        if( ABS(vSolut(iSp1)%Z)<2 .or. ABS(vSolut(iSp2)%Z)<2 ) then
          vAlfa1(tI_Beta1(iSp1,iSp2))= 2.0D0
        else
          vAlfa1(tI_Beta1(iSp1,iSp2))= 1.4D0
        end if
      end if

      if(tI_Beta2(iSp1,iSp2)>0) then
        if(ABS(vSolut(iSp1)%Z)<2 .or. ABS(vSolut(iSp2)%Z)<2) then
          !-> ignore the tBeta2 term in B_MX expression
          vAlfa2(tI_Beta2(iSp1,iSp2))= Zero
        else
          vAlfa2(tI_Beta2(iSp1,iSp2))= 12.0D0
        end if
      end if
      !
    enddo
  enddo
  !----------------------------------------------------/ Alfa1,tAlfa2 --
  !
  call Solmodel_Pitzer_Dtb_TPUpdate(298.15D0,1.0D0)
  !
  if(iDebug>0) then !=========================================< trace ==
    call GetUnit(ff)
    open(ff,file="debug_pitzer_data.log")
    !
    call Solmodel_Pitzer_Dtb_Check(ff,vSolut)
    call Solmodel_Pitzer_Dtb_Check_Old(ff,vSolut)
    !
    close(ff)
  end if !=====================================================< trace ==
  ! Solmodel_Pitzer
  if(iDebug>0) write(fTrc,'(A,/)') "</ Solmodel_Pitzer_ Solmodel_PitzerDtb_Read"
  !
contains

logical function Ok_2Sp(I1,I2,V)
!--
!-- check that a species pair is not already read,
!-- update the table V
!--
  integer,intent(in)   :: I1,I2
  integer,intent(inout):: V(:,:)
  !
  Ok_2Sp= (V(I1,I2)==0)
  if(Ok_2Sp) then
    V(I1,I2)= 1
  else
    print '(2(A,1X))',trim(vSolut(I1)%NamSp),trim(vSolut(I2)%NamSp)
    call Stop_("Pitz_Dtb.Error: already allocated")
  end if
  !
end function Ok_2Sp

logical function Ok_3Sp(I1,I2,I3,V)
!--
!-- check that a species triplet is not already read,
!-- update the table V
!--
  integer,intent(in)   :: I1,I2,I3
  integer,intent(inout):: V(:,:,:)
  !
  Ok_3Sp= (V(I1,I2,I3)==0)
  if(Ok_3Sp) then
    V(I1,I2,I3)= 1
  else
    call Stop_("Error in Pitz_Dtb: triplet already allocated")
  end if
  !
end function Ok_3Sp

end subroutine Solmodel_Pitzer_Dtb_Read

subroutine Sort2sp(I1,I2,Error)
!--
!-- order the I1,I2 pair,
!-- to put the data in the upper triangle
!--
  integer,intent(inout):: I1,I2
  logical,intent(out)  :: Error
  !
  integer:: J1,J2
  !
  Error= (I1==I2)
  !
  if(Error) return
  !
  J1= MIN(I1,I2)  ;  J2= MAX(I1,I2)
  I1= J1          ;  I2= J2
  !
end subroutine Sort2sp

subroutine Sort3sp(I1,I2,I3,Error)
!-- order I1,I2,I3 to I1-I2-I3
  integer,intent(inout):: I1,I2,I3
  logical,intent(out)  :: Error
  !
  integer:: J1,J2,J3
  !
  Error= (I1==I2 .or. I1==I3 .or. I2==I3)
  if(Error) return
  !
  J1= MIN(I1,I2,I3)
  J3= MAX(I1,I2,I3)
  !
  if(I1/=J1 .and. I1/=J3) J2= I1
  if(I2/=J1 .and. I2/=J3) J2= I2
  if(I3/=J1 .and. I3/=J3) J2= I3
  !
  I1= J1  ;  I2= J2  ;  I3= J3
  !
end subroutine Sort3sp

subroutine Species_Scan(Str,N,iSpc,Error)
!--
!-- scan a species string, where species are separated by "-",
!-- and compute corresponding indexes in vSolut
!--
  use M_T_Species,only: Species_Index
  character(len=*),intent(in):: Str
  integer,intent(in) :: N
  integer,intent(out):: iSpc(:)
  logical,intent(out):: Error
  !
  character(len=23):: Str1,Str2
  integer:: I,J,K
  !
  Error= .false.
  !
  Str1= trim(Str)//"=" !append a separator at end
  Str2= ""
  K= 0
  do
    !write(91,'(2A)') "str1=",trim(Str1)
    I= INDEX(Str1,'=') ! first occurence of separator
    J= len_trim(Str1)
    if(I==0) exit
    write(Str2,'(A)') trim(str1(1:I-1)) ! Str2= left of first separator
    write(Str1,'(A)') trim(str1(I+1:J)) ! Str1= right of first separator
    K= K+1
    iSpc(K)= Species_Index(Str2,vSolut)
    !write(91,'(A,I3,1X,A)') "    I,str2=",iSpc(K),trim(Str2)
    if(K==N) exit
  enddo
  !pause
  !write(91,*)
  !
  Error= (K<N)
  !
  return
end subroutine Species_Scan

subroutine Solmodel_Pitzer_Dtb_TPUpdate(TdgK,Pbar)
  real(dp),intent(in):: TdgK,Pbar
  !
  call Fitt_TPUpdate(TdgK,Pbar, vF_Beta0, vBeta0)
  call Fitt_TPUpdate(TdgK,Pbar, vF_Beta1, vBeta1)
  call Fitt_TPUpdate(TdgK,Pbar, vF_Beta2, vBeta2)
  !
  call Fitt_TPUpdate(TdgK,Pbar, vF_CPhi,  vCPhi)
  call Fitt_TPUpdate(TdgK,Pbar, vF_Theta, vTheta)
  call Fitt_TPUpdate(TdgK,Pbar, vF_Lamda, vLamda)
  !
  call Fitt_TPUpdate(TdgK,Pbar, vF_Psi,  vPsi)
  call Fitt_TPUpdate(TdgK,Pbar, vF_Zeta, vZeta)
  !
  call Pitzer_TableFill
  !
  return
end subroutine  Solmodel_Pitzer_Dtb_TPUpdate

subroutine Pitzer_TableFill
  !
  call TableFill(tI_Beta0, vBeta0, tBeta0)
  call TableFill(tI_Beta1, vBeta1, tBeta1)
  call TableFill(tI_Beta2, vBeta2, tBeta2)
  !
  call TableFill(tI_Beta1, vAlfa1, tAlfa1)
  call TableFill(tI_Beta2, vAlfa2, tAlfa2)
  !
  call TableFill(tI_CPhi,  vCPhi,  tCPhi )
  call TableFill(tI_Theta, vTheta, tTheta)
  call TableFill(tI_Lamda, vLamda, tLamda)
  !
  call TableFill3(tI_Psi,  vPsi,  tPsi)
  call TableFill3(tI_Zeta, vZeta, tZeta)
  !
  return
end subroutine Pitzer_TableFill

subroutine TableFill(tI_In,vCoef,tR_Out)
  integer, intent(in) :: tI_In(:,:)
  real(dp),intent(in) :: vCoef(:)
  real(dp),intent(out):: tR_Out(:,:)
  !
  integer:: I,J,K
  !
  do I=1,size(tI_In,1)
    do J=1,size(tI_In,2)
      K= tI_In(I,J)
      if(K>0) tR_Out(I,J)= vCoef(K)
    enddo
  enddo
  !
end subroutine TableFill

subroutine TableFill3(tI_In,vCoef,tR_Out)
  integer, intent(in) :: tI_In(:,:,:)
  real(dp),intent(in) :: vCoef(:)
  real(dp),intent(out):: tR_Out(:,:,:)
  !
  integer:: I,J,K,L
  !
  do I=1,size(tI_In,1)
    do J=1,size(tI_In,2)
      do K=1,size(tI_In,3)
        L= tI_In(I,J,K)
        if(L>0) tR_Out(I,J,K)= vCoef(L)
      enddo
    enddo
  enddo
  !
end subroutine TableFill3

subroutine Fitt_TPUpdate(T,P,vFitt,vCoef)
  real(dp),intent(in):: T,P !TdgK,Pbar
  type(T_Fitt),intent(in):: vFitt(:)
  real(dp),intent(out):: vCoef(:)
  !
  integer:: I
  real(dp),parameter:: Tref= 298.15D0
  type(T_Fitt):: Fit
  !
  do I=1,size(vCoef)
    !
    Fit= vFitt(I)
    !
    select case(Fit%iModel)

    case(0) !"none"
      vCoef(I)= Fit%vX(1) !"none"

    case(1) !"PHREEQ"
      vCoef(I)= Fit%vX(1) &
      &       + Fit%vX(2) *(One/T -One/Tref) &
      &       + Fit%vX(3) *log(T/Tref)       &
      &       + Fit%vX(4) *(T     -Tref)     &
      &       + Fit%vX(5) *(T*T   -Tref*Tref)
      !
    case(2) !"KUEHN"
      vCoef(I)= Fit%vX(1) &
      &       + Fit%vX(2) *(One/T -One/Tref) &
      &       + Fit%vX(3) *log(T/Tref)       &
      &       + Fit%vX(4) *(T     -Tref)     &
      &       + Fit%vX(5) *(T*T   -Tref*Tref)
      !
    case(3) !"CHRISTOV-2004"
      vCoef(I)= Fit%vX(1) &              ! a1
      &       + Fit%vX(2) *T           & ! a2*T
      &       + Fit%vX(3) *T*T         & ! a3*T**2
      &       + Fit%vX(4) *T**3        & ! a4*T**3
      &       + Fit%vX(5) *(One/T)     & ! a5/T
      &       + Fit%vX(6) *log(T)      & ! a6*ln(T)
      &       + Fit%vX(7) /(T-263.0D0) & ! a7/(T-263)
      &       + Fit%vX(8) /(680.0D0-T) & ! a8/(680-T)
      &       + Fit%vX(9) /(T-227.0D0)   ! a9/(T-227)
      !
      ! AZAROUAL, Chemical Geology,
      ! f(T)= a1 + a2 T  + a3 /T         + a4 lnT       + a5 /(T- 263)
      !          + a6 T² + a7 /(680 - T) + a8 /(T- 227)
      ! nearly the same equation as Christov, 2004,
      ! with a different order and an additional term in /(T- 227)
      ! a1 = 1
      ! a2 = 2
      ! a3 = 5
      ! a4 = 6
      ! a5 = 7
      ! a6 = 3
      ! a7 = 8
      ! a8 = 9
      !
    case(4) !"PITZER-1984", Na+=Cl-
      ! K.S. Pitzer, J.C. Peiper and R.H. Busey, 1984
      ! Thermodynamic properties of aqueous sodium chloride solutions,
      ! Journal of Physical and Chemical Reference Data 13 (1) (1984), pp. 1-102.
      vCoef(I)= & !
      &  Fit%vX(1)/T &
      +  Fit%vX(2)   + Fit%vX(3) *P + Fit%vX(4) *P**2 + Fit%vX(5)*P**3                   &
      +  Fit%vX(6)*log(T)                                                                &
      + (Fit%vX(7)   + Fit%vX(8) *P + Fit%vX(9) *P**2 + Fit%vX(10)*P**3)*T               &
      + (Fit%vX(11)  + Fit%vX(12)*P + Fit%vX(13)*P**2                  )*T**2            &
      + (Fit%vX(14)  + Fit%vX(15)*P + Fit%vX(16)*P**2 + Fit%vX(17)*P**3)/(T-227.0D0)     &
      + (Fit%vX(18)  + Fit%vX(19)*P + Fit%vX(20)*P**2 + Fit%vX(21)*P**3)/(680.0D0 - T)
      !
    case(5) !"PABALAN-1987"
      vCoef(I)= & !
      &   Fit%vX(1)  + Fit%vX(2)*P            &
      + ( Fit%vX(3)  + Fit%vX(4)*P )  /T      &
      +   Fit%vX(5)                   *log(T) &
      + ( Fit%vX(6)  + Fit%vX(7)*P )  *T      &
      + ( Fit%vX(8)  + Fit%vX(9)*P )  *T**2   &
      +   Fit%vX(10)                  /(T-227.0D0) &
      + ( Fit%vX(11) + Fit%vX(12)*P ) /(647.0D0-T)
      !
    case(6) !"POLYA-2007"
      vCoef(I)= & !
      & Fit%vX(1)    &
      + Fit%vX(2) *T &
      + Fit%vX(3) /T &
      + Fit%vX(4) /(T-210.0D0) &
      + Fit%vX(5) /(647.0D0-T) &
      + Fit%vX(6) *(T-443.0D0)**3/3.0D0 &
      + Fit%vX(7) *(P-One) &
      + Fit%vX(8) *(P-One)**2/2.0D0

    case(7) ! "LI-DUAN-2007"
      vCoef(I)= & !
      & Fit%vX(1)             &
      + Fit%vX(2) *T          &
      + Fit%vX(3) /T          &
      + Fit%vX(4)*T**2        &
      + Fit%vX(5)/(630.0D0 - T)   &
      + Fit%vX(6) *P          &
      + Fit%vX(7)*P*log(T)    &
      + Fit%vX(8)*P/T         &
      + Fit%vX(9) *P/(630.0D0 - T) &
      + Fit%vX(10)*P**2/(630.0D0 - T)**2  &
      + Fit%vX(11)*T*log(P)

    case default
      call Stop_("in Fitt_TPUpdate, iModel out of range")

    end select
    !
  end do
  !
  return
end subroutine Fitt_TPUpdate

subroutine FGENERIQUE (INDIC,NOM,T,P,A,AV,TREF1,TREF2,TREF3,FF)
!!! UNuseD !!!
! Calcul des parametres de Pitzer en fonction de T et P
!
! Arguments d'entree :
! --------------------
! INDIC : si 0 : fonction par defaut (Monnin, Chemical Geology, 153(1999)187-209
!   sinon : l'user specifie sa fonction
! NOM : nom de la fonction si l'user la specifie
! T : Temperature (K)
! P : Pression (Pa)
! A : parametres de la fonction par defaut
! AV: parametres de la fonction par defaut
! TREF1 : temperature de reference 1 pour la fonction par defaut (K)
! TREF2 : temperature de reference 2 pour la fonction par defaut (K)
! TREF3 : temperature de reference 3 pour la fonction par defaut (K)
!
! Arguments de sortie :
! ---------------------
! FF : valeur calculee par la fonction
  !
  integer,      intent(in):: INDIC
  character*200,intent(in):: NOM
  real(dp),     intent(in):: A(10),AV(10)
  real(dp),     intent(in):: T,P,TREF1,TREF2,TREF3
  !
  real(dp),     intent(out):: FF
  !
  real(dp):: F0,FV,PBAR
  !
  PBAR=P*1.0d-5
  !
  if (INDIC==0) then
    !----------------------------- Fonction par defaut (Monnin, 1999) --
    F0= A(1)                 &
    & + A(2) *T             &
    & + A(3) /T             &
    & + A(4) *log(T)        &
    & + A(5) /(T - TREF1)   &
    & + A(6) *T**2          &
    & + A(7) /(TREF2 - T)   &
    & + A(8) /(T - TREF3)   &
    & + A(9) *T**3          &
    & + A(10)*T**4
    !
    !-- Volume function (pressure correction)
    FV= AV(1)              &
    & + AV(2) *T            &
    & + AV(3) /T            &
    & + AV(4) *log(T)       &
    & + AV(5) /(T - TREF1)  &
    & + AV(6) *T**2         &
    & + AV(7) /(TREF2 - T)  &
    & + AV(8) /(T - TREF3)  &
    & + AV(9) *T**3         &
    & + AV(10)*T**4
    !
    FF=  F0 + PBAR * FV
    !----------------------------/ Fonction par defaut (Monnin, 1999) --
  else
    select case (NOM)
      !  Ecrire ici case ('Nom de la nouvelle fonction')
      !  Ecrire ici FF=
      !  case (
    end select
  end if
  !
  return
end subroutine FGENERIQUE

!~ subroutine Fitt_TPUpdate_(TdgK,vFitt,vCoef)
  !~ real(dp),intent(in):: TdgK
  !~ type(T_Fitt),intent(in):: vFitt(:)
  !~ real(dp),intent(out):: vCoef(:)
  !~ !
  !~ integer:: I
  !~ real(dp),parameter:: Tref= 298.15D0
  !~ !
  !~ do I=1,size(vCoef)
    !~ vCoef(I)= vFitt(I)%X0 &
    !~ &       + vFitt(I)%X1 *(One/TdgK -One/Tref) &
    !~ &       + vFitt(I)%X2 *log(TdgK/Tref)       &
    !~ &       + vFitt(I)%X3 *(TdgK     -Tref)     &
    !~ &       + vFitt(I)%X4 *(TdgK*TdgK -Tref*Tref)
  !~ end do
  !~ !
  !~ return
!~ end subroutine Fitt_TPUpdate_

subroutine Solmodel_Pitzer_Dtb_TPtest
  use M_IoTools
  use M_Files,only: DirOut
  !
  integer :: f1,f2
  ! integer :: I,J,K
  real(dp):: vT(1:9),Pbar
  ! type(T_Fitt):: Fitt
  ! real(dp),allocatable:: Res(:,:)
  ! character(len=30):: cFormat
  !
  vT(1:9)=(/ 0.0D0,  25.0D0, 50.0D0, 75.0D0, &
  &        100.0D0, 125.0D0,150.0D0,175.0D0, &
  &        200.0D0 /)
  Pbar= 1.0D0
  !
  !write(cFormat,'(A,I3,A)') '(3(A,1X),',n,'(G12.3,1X))'
  !write(cFormat,'(A)') '(3(A,1X),5(G12.3,1X))'
  !
  !write(ff,cFormat,advance="no") "X=",(x(k),k=1,n)
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< Pitzer_Dtb_TPtest"
  !
  call GetUnit(f1)
  open(f1,file=trim(DirOut)//"_pitzer_tptest.tab")
  call GetUnit(f2)
  open(f2,file=trim(DirOut)//"_pitzer_tpcoeffs.tab")
  !
  call TPtest(f1,f2, "BETA0", vT,Pbar, tI_Beta0, vF_Beta0,vBeta0)
  call TPtest(f1,f2, "BETA1", vT,Pbar, tI_Beta1, vF_Beta1,vBeta1)
  call TPtest(f1,f2, "BETA2", vT,Pbar, tI_Beta2, vF_Beta2,vBeta2)
  call TPtest(f1,f2, "CPHI",  vT,Pbar, tI_CPhi,  vF_CPhi, vCPhi )
  call TPtest(f1,f2, "THETA", vT,Pbar, tI_Theta, vF_Theta,vTheta)
  call TPtest(f1,f2, "LAMBDA",vT,Pbar, tI_Lamda, vF_Lamda,vLamda)
  !~ allocate(Res(9,size(vBeta0)))
  !~ do I=1,9
    !~ call Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_Beta0, vBeta0)
    !~ Res(I,:)= vBeta0(:)
  !~ enddo
  !~ do I=1,size(tI_Beta0,1)
    !~ do J=1,size(tI_Beta0,2)

      !~ if(tI_Beta0(I,J)>0) then
        !~ print *, &
        !~ & trim(vSolut(I)%NamSp), " ", &
        !~ & trim(vSolut(J)%NamSp)
        !~ Fitt= vF_Beta0(tI_Beta0(I,J))
        !~ write(f1,'(3(A12,1X),9(G15.6,1X))') &
        !~ & "BETA0", &
        !~ & trim(vSolut(I)%NamSp), &
        !~ & trim(vSolut(J)%NamSp), &
        !~ & (Res(K,tI_Beta0(I,J)), K=1,9)
        !~ write(f2,'(4(A12,1X),24(G15.6,1X))') &
        !~ & trim(Fitt%namModel), &
        !~ & "BETA0", &
        !~ & trim(vSolut(I)%NamSp), &
        !~ & trim(vSolut(J)%NamSp), &
        !~ & (Fitt%vX(K),           K=1,Fitt%N)
      !~ end if

    !~ enddo
  !~ enddo
  !~ deallocate(Res)
  !
  !~ allocate(Res(9,size(vBeta1)))
  !~ do I=1,9
    !~ call Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_beta1, vbeta1)
    !~ Res(I,:)= vbeta1(:)
  !~ enddo
  !~ do I=1,size(tI_beta1,1)

    !~ do J=1,size(tI_beta1,2)

      !~ if(tI_beta1(I,J)>0) then
        !~ Fitt= vF_Beta1(tI_Beta1(I,J))
        !~ write(f1,'(3(A12,1X),9(G15.6,1X))') &
        !~ & "BETA1", &
        !~ & trim(vSolut(I)%NamSp), &
        !~ & trim(vSolut(J)%NamSp), &
        !~ & (Res(K,tI_beta1(I,J)), K=1,9)
        !~ write(f2,'(4(A12,1X),24(G15.6,1X))') &
        !~ & trim(Fitt%namModel), &
        !~ & "BETA1", &
        !~ & trim(vSolut(I)%NamSp), &
        !~ & trim(vSolut(J)%NamSp), &
        !~ & (Fitt%vX(K),           K=1,Fitt%N)
      !~ end if

    !~ enddo

  !~ enddo
  !~ deallocate(Res)
  !
  close(f1)
  close(f2)
  !
  if(iDebug>0) print '(/,A,/)',"=!= results in pitzer_tptest.tab =!="
  !
  if(iDebug>0) write(fTrc,'(A,/)') "<  Solmodel_Pitzer_Pitzer_Dtb_TPtest"
  !
end subroutine Solmodel_Pitzer_Dtb_TPtest

subroutine TPtest(f1,f2,Str,vT,Pbar,tI_Beta,vF_Beta,vBeta)
  !
  integer,     intent(in):: f1,f2
  character(*),intent(in):: Str
  real(dp),    intent(in):: vT(:)
  real(dp),    intent(in):: Pbar
  integer,     intent(in):: tI_Beta(:,:)
  type(T_Fitt),intent(in):: vF_Beta(:)
  real(dp),    intent(in):: vBeta(:)
  !
  type(T_Fitt):: Fitt
  integer:: I,J,K,nT
  real(dp),allocatable:: Res(:,:)
  real(dp),allocatable:: vX(:)
  
  nT= size(vT)
  !
  allocate(vX(size(vBeta)))
  allocate(Res(9,size(vBeta)))
  !
  do I=1,nT
    call Fitt_TPUpdate(vT(I)+273.15D0,Pbar, vF_Beta, vX)
    Res(I,:)= vX(:)
  enddo
  !
  do I=1,size(tI_Beta,1)
    do J=1,size(tI_Beta,2)

      if(tI_Beta(I,J)>0) then
        print *, &
        & trim(vSolut(I)%NamSp), " ", &
        & trim(vSolut(J)%NamSp)
        Fitt= vF_Beta(tI_Beta(I,J))
        write(f1,'(2(A,1X),9(G15.6,1X))') &
        & trim(Str), &
        & trim(vSolut(I)%NamSp)//"="//trim(vSolut(J)%NamSp), &
        & (Res(K,tI_Beta(I,J)), K=1,nT)
        write(f2,'(3(A,1X),24(G15.6,1X))') &
        & trim(Fitt%namModel), &
        & trim(Str), &
        & trim(vSolut(I)%NamSp)//"="//trim(vSolut(J)%NamSp), &
        & (Fitt%vX(K), K=1,Fitt%N)
      end if

    enddo
  enddo
  !
  deallocate(Res)
  deallocate(vX)

end subroutine TPtest

subroutine ReadRVals4(Line,x1,x2,x3,x4)
  use M_IOTools,only:LinToWrd,WrdToReal
  character(len=*) Line
  character(255)   Word
  real(dp)::x1,x2,x3,x4
  integer  ::i
  logical  ::EoL
  real(dp),dimension(1:4)::vX
  vX=0.0
  EoL=.true.
  do i=1,4
    call LinToWrd(Line,Word,EoL)
    call WrdToReal(Word,vX(i))
    if(EoL) exit
  enddo
  x1=vX(1) ; x2=vX(2) ; x3=vX(3) ; x4=vX(4)
end subroutine ReadRVals4

subroutine Solmodel_Pitzer_Dtb_Check(ff,vSpc)
  use M_T_Species,only: T_Species
  
  type(T_Species),intent(in):: vSpc(:)
  integer,intent(in):: ff
  
  call Dtb_Check(ff,"BETA0",    vSpc, tI_Beta0 , vBeta0)
  call Dtb_Check(ff,"BETA1",    vSpc, tI_Beta1 , vBeta1)
  call Dtb_Check(ff,"BETA2",    vSpc, tI_Beta2 , vBeta2)
  call Dtb_Check(ff,"ALFA1",    vSpc, tI_Beta1 , vAlfa1)
  call Dtb_Check(ff,"ALFA2",    vSpc, tI_Beta2 , vAlfa2)
  call Dtb_Check(ff,"CPHI",     vSpc, tI_CPhi  , vCPhi)
  call Dtb_Check(ff,"THETA",    vSpc, tI_Theta , vTheta)
  call Dtb_Check(ff,"LAMBDA",   vSpc, tI_Lamda , vLamda)
  call Dtb_Check_Dim3(ff,"ZETA",vSpc, tI_Zeta  , vZeta)
  call Dtb_Check_Dim3(ff,"PSI", vSpc, tI_Psi   , vPsi)
  
  return
end subroutine Solmodel_Pitzer_Dtb_Check

subroutine Solmodel_Pitzer_Dtb_Check_Old(ff,vSpc)
  use M_T_Species,only: T_Species
  
  type(T_Species),intent(in):: vSpc(:)
  integer,intent(in):: ff
  
  call Dtb_Check_Old(ff,"BETA0",vSpc,tBeta0)
  call Dtb_Check_Old(ff,"BETA1",vSpc,tBeta1)
  call Dtb_Check_Old(ff,"BETA2",vSpc,tBeta2)
  call Dtb_Check_Old(ff,"ALFA1",vSpc,tAlfa1)
  call Dtb_Check_Old(ff,"ALFA2",vSpc,tAlfa2)
  call Dtb_Check_Old(ff,"CPHI", vSpc,tCPhi)
  call Dtb_Check_Old(ff,"THETA",vSpc,tTheta)
  call Dtb_Check_Old(ff,"LAMBDA",vSpc,tLamda)
  call Dtb_Check_Dim3_Old(ff,"ZETA",vSpc,tZeta)
  call Dtb_Check_Dim3_Old(ff,"PSI", vSpc,tPsi)
  
  return
end subroutine Solmodel_Pitzer_Dtb_Check_Old

subroutine Dtb_Check(f,Str,vSpc,tI_Coef,vCoef)

  use M_T_Species,only: T_Species
  !
  integer,        intent(in):: f
  character(*),   intent(in):: Str
  type(T_Species),intent(in):: vSpc(:)
  integer,        intent(in):: tI_Coef(:,:)
  real(dp),       intent(in):: vCoef(:)
  !
  integer :: I,J
  real(dp):: X
  !
  write(f,'(/,A,/)') trim(Str)
  do I=1,size(tI_Coef,1)
    do J=1,size(tI_Coef,2)
      if(tI_Coef(I,J)>0) then
        X= vCoef(tI_Coef(I,J))
        write(f,'(G15.6,A1,2(A,A1))') &
        & X,t_,trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
      end if
    enddo
  enddo
  !
end subroutine Dtb_Check

subroutine Dtb_Check_Dim3(f,Str,vSpc,tI_Coef,vCoef)
!--
!-- same as Dtb_Check, for 3D tables
!--
  use M_T_Species,only: T_Species
  !
  integer,        intent(in):: f
  character(*),   intent(in):: Str
  type(T_Species),intent(in):: vSpc(:)
  integer,        intent(in):: tI_Coef(:,:,:)
  real(dp),       intent(in):: vCoef(:)
  !---------------------------------------------------------------------
  integer :: I,J,K
  real(dp):: X
  !
  write(f,'(/,A,/)') trim(Str)
  do I=1,size(tI_Coef,1)
    do J=1,size(tI_Coef,2)
      do K=1,size(tI_Coef,3)
        if(tI_Coef(I,J,K)>0) then
          X= vCoef(tI_Coef(I,J,K))
          write(f,'(G15.6,A1,3(A,A1))') &
          & X,t_, &
          & trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_,trim(vSpc(K)%NamSp),t_
        end if
      enddo
    enddo
  enddo
  
  return
end subroutine Dtb_Check_Dim3

subroutine Dtb_Check_Old(f,Str,vSpc,vCoef)
!-- 
  use M_T_Species,only: T_Species
  !---------------------------------------------------------------------
  integer,        intent(in):: f
  character(*),   intent(in):: Str
  type(T_Species),intent(in):: vSpc(:)
  real(dp),       intent(in):: vCoef(:,:)
  !
  integer :: I,J
  real(dp):: X
  !---------------------------------------------------------------------
  
  write(f,'(/,A,/)') trim(Str)
  do I=1,size(vCoef,1)
    do J=1,size(vCoef,2)
      X= vCoef(I,J)
      if(X/=Zero) &
      & write(f,'(G15.6,A1,2(A,A1))') &
      & X,t_,trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_
    enddo
  enddo
  
  return
end subroutine Dtb_Check_Old

subroutine Dtb_Check_Dim3_Old(f,Str,vSpc,vCoef)
! same as Dtb_Check, for 3D tables
  use M_T_Species,only: T_Species
  !
  integer,        intent(in):: f
  character(*),   intent(in):: Str
  type(T_Species),intent(in):: vSpc(:)
  real(dp),       intent(in):: vCoef(:,:,:)
  !
  integer :: I,J,K
  real(dp):: X
  !
  write(f,'(/,A,/)') trim(Str)
  do I=1,size(vCoef,1)
    do J=1,size(vCoef,2)
      do K=1,size(vCoef,3)
        X= vCoef(I,J,K)
        if(X/=Zero) &
        & write(f,'(G15.6,A1,3(A,A1))') &
        & X,t_,trim(vSpc(I)%NamSp),t_,trim(vSpc(J)%NamSp),t_,trim(vSpc(K)%NamSp),t_
      enddo
    enddo
  enddo
  !
end subroutine Dtb_Check_Dim3_Old

subroutine Solmodel_Pitzer_Dtb_Clean
  !
  if (allocated(vSolut)) deallocate(vSolut)
  !
  if (allocated(tBeta0))   deallocate(tBeta0)
  if (allocated(tBeta1))   deallocate(tBeta1)
  if (allocated(tBeta2))   deallocate(tBeta2)
  if (allocated(tCPhi))    deallocate(tCPhi)
  if (allocated(tTheta))   deallocate(tTheta)
  if (allocated(tLamda))   deallocate(tLamda)
  if (allocated(tAlfa1))   deallocate(tAlfa1)
  if (allocated(tAlfa2))   deallocate(tAlfa2)
  if (allocated(tZeta))    deallocate(tZeta)
  if (allocated(tPsi))     deallocate(tPsi)
  !
  if (allocated(tI_Beta0))   deallocate(tI_Beta0)
  if (allocated(tI_Beta1))   deallocate(tI_Beta1)
  if (allocated(tI_Beta2))   deallocate(tI_Beta2)
  if (allocated(tI_CPhi))    deallocate(tI_CPhi)
  if (allocated(tI_Theta))   deallocate(tI_Theta)
  if (allocated(tI_Lamda))   deallocate(tI_Lamda)
  if (allocated(tI_Zeta))    deallocate(tI_Zeta)
  if (allocated(tI_Psi))     deallocate(tI_Psi)
  !
  if (allocated(vF_Beta0))   deallocate(vF_Beta0)
  if (allocated(vF_Beta1))   deallocate(vF_Beta1)
  if (allocated(vF_Beta2))   deallocate(vF_Beta2)
  if (allocated(vF_CPhi))    deallocate(vF_CPhi)
  if (allocated(vF_Theta))   deallocate(vF_Theta)
  if (allocated(vF_Lamda))   deallocate(vF_Lamda)
  if (allocated(vF_Zeta))    deallocate(vF_Zeta)
  if (allocated(vF_Psi))     deallocate(vF_Psi)
  !
  if (allocated(vBeta0))   deallocate(vBeta0)
  if (allocated(vBeta1))   deallocate(vBeta1)
  if (allocated(vBeta2))   deallocate(vBeta2)
  if (allocated(vAlfa1))   deallocate(vAlfa1)
  if (allocated(vAlfa2))   deallocate(vAlfa2)
  if (allocated(vCPhi))    deallocate(vCPhi)
  if (allocated(vTheta))   deallocate(vTheta)
  if (allocated(vLamda))   deallocate(vLamda)
  if (allocated(vZeta))    deallocate(vZeta)
  if (allocated(vPsi))     deallocate(vPsi)
  !
end subroutine Solmodel_Pitzer_Dtb_Clean

end module M_Solmodel_Pitzer_Dtb
