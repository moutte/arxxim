module M_Solmodel_Read
!--
!-- reading solvent parameters from database
!--
  use M_Trace,only: iDebug,fTrc,Stop_,T_
  use M_Kinds

  implicit none

  private

  public:: Solmodel_Solvent_Read
  public:: SolModel_Read

contains

!--subroutine Solmodel_Solvent_Read-------------------------------------
!< reads Solvent parameters from database
!! density, dielectric, Debye-Hueckel param's
!-----------------------------------------------------------------------
subroutine Solmodel_Solvent_Read( &
& nTP,              & ! IN
& vTdgC,            & ! IN
& SolModel,         & ! INOUT
& Rho_Ok, Eps_Ok, DHA_Ok, DHB_Ok, BDot_Ok, & !OUT
& Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)  !OUT

  use M_IOTools !, only:dimV,LinToWrd
  use M_Dtb_Const, only: T_CK
  use M_Files,     only: NamFLogK,NamFPtz,NamFInn
  use M_T_SolModel,only: T_SolModel,T_SolModelDat
  use M_T_SolModel,only: vSolModelAct
  use M_Solmodel_Vars,only: T_Spline

  integer,           intent(in)   :: nTP       !< dim'n of vTdgC
  real(dp),          intent(in)   :: vTdgC(:)  !< temperature array for spline
  type(T_SolModel),  intent(inout):: SolModel  !< solution model
  logical,           intent(out)  :: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  type(T_Spline),    intent(inout):: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl

  character(len=512):: L,W
  real(dp)::vX(dimV)
  logical :: Ok,EoL
  !logical:: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  integer:: ios,mDum,f,I

  type(T_SolModelDat):: vSolvDat(nTP)
  ! real(dp),allocatable:: B(:), C(:), D(:)

  !--

  if(iDebug>0) write(fTrc,'(/,A)') "< Solmodel_Solvent_Read"

  Ok=.false.
  !
  Rho_Ok= .false.
  Eps_Ok= .false.
  DHA_Ok= .false.
  DHB_Ok= .false.
  BDot_Ok=.false.

  call GetUnit(f)
  open(f,file=trim(NamFLogK),STATUS='OLD')
  !
  DoFile: do

    read(f,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile

    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile

    call AppendToEnd(L,W,EoL)

    if(W=="ENDINPUT") exit DoFile

    if(W=="SOLVENT") then

      Ok=.true.

      DoSolvent: do

        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)

        if(W(1:1)=='!')               cycle DoSolvent
        call AppendToEnd(L,W,EoL)
        if(W=="ENDINPUT")             exit  DoFile
        if(W=="END") exit DoSolvent
        if(W=="ENDSOLVENT")           exit DoSolvent

        if(INDEX("RHO_EPS_BdoT_DHA_DHB",trim(W))>0) then
          call ReadRValsV(L,mDum,vX)
          if(mDum<nTP) call Stop_("!!! in SolModel block, data nr < dim' TP.table !!!")
        end if

        select case(trim(W))

        case("NAME")
          call LinToWrd(L,W,EoL)
          SolModel%Name=trim(W)

        case("MODEL")
          call LinToWrd(L,W,EoL)
          !if(len_trim(W)>15) call Stop_(trim(W)//" = SolModel model name too long !!!")
          !
          SolModel%iActModel= 0
          do i=1,size(vSolModelAct)
            if(trim(W)==trim(vSolModelAct(i))) then
              SolModel%iActModel= i
              exit
            end if
          enddo
          !~ if(INDEX(Solmodel_ModelList,trim(W)//"_")==0) &
          !~ & call Stop_(trim(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          if(SolModel%iActModel==0) &
          & call Stop_(trim(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          !~ SolModel%ActModel=trim(W)
          !~ if(SolModel%ActModel=="PITZER") NamFPtz= trim(NamFInn)
          if(SolModel%iActModel==8) NamFPtz= trim(NamFInn)

        case("RHO")
          if(nTP>0) then
            vSolvDat(1:nTP)%Rho=  vX(1:nTP)  ;   Rho_Ok= .true.
          end if

        case("EPS")
          if(nTP>0) then
            vSolvDat(1:nTP)%Eps=  vX(1:nTP)  ;   Eps_Ok= .true.
          end if

        case("DHA")
          if(nTP>0) then
            vSolvDat(1:nTP)%DHA=  vX(1:nTP)  ;   DHA_Ok= .true.
          end if

        case("DHB")
          if(nTP>0) then
            vSolvDat(1:nTP)%DHB=  vX(1:nTP)  ;   DHB_Ok= .true.
          end if

        case("BdoT")
          if(nTP>0) then
            vSolvDat(1:nTP)%BDot= vX(1:nTP)  ;   BDot_Ok=.true.
          end if

        end select

      enddo DoSolvent

    end if !W=="SolModel"

  enddo DoFile

  close(f)

  !~ if(nTP>0) then

    if(Rho_Ok) &
    & call Spline_Init(nTP,vTdgC(:),vSolvDat(:)%Rho,Rho_spl)
    if(Eps_Ok) &
    & call Spline_Init(nTP,vTdgC(:),vSolvDat(:)%Eps,Eps_spl)
    if(DHA_Ok) &
    & call Spline_Init(nTP,vTdgC(:),vSolvDat(:)%DHA,DHA_spl)
    if(DHB_Ok) &
    & call Spline_Init(nTP,vTdgC(:),vSolvDat(:)%DHB,DHB_spl)
    if(BDot_Ok) &
    & call Spline_Init(nTP,vTdgC(:),vSolvDat(:)%BDot,BDot_spl)

  !~ end if

  ! A.M Error Detected Here
  !!if(iDebug>0) write(fTrc,'(A,E15.8)') "DHA",vSolvDat(1)%DHA
  !!if(iDebug>0) write(fTrc,'(A,E15.8)') "DHB",vSolvDat(1)%DHB

  if(iDebug>0) write(fTrc,'(A,/)') "</ Solmodel_Solvent_Read"

end subroutine Solmodel_Solvent_Read

subroutine SolModel_Read(vSpc,Ok,MsgError)
!--
!-- initialize vSolModel --
!--

  use M_IOTools

  use M_Files,     only: NamFInn
  use M_T_Species, only: T_Species,Species_Index
  use M_T_SolModel,only: T_SolModel,T_SolModelDat
  use M_T_SolModel,only: vSolModelAct
  !
  use M_Global_Vars,only: vSolModel
  !---------------------------------------------------------------------
  type(T_Species), intent(in) :: vSpc(:)
  logical,         intent(out):: Ok
  character(len=*),intent(out):: MsgError
  !---------------------------------------------------------------------
  type(T_SolModel):: S
  type(T_SolModel):: vSolTmp(10)
  integer:: vIsolTmp(300)
  !
  character(len=512):: L,W
  logical :: sEoL
  integer:: ios
  integer:: F
  !logical:: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  integer:: I,K,N,nS,M
  !---------------------------------------------------------------------
  if(iDebug>0) write(fTrc,'(A,/)') "< SolModel_Read"

  call GetUnit(F)
  open(F,file=trim(NamFInn))

  Ok= .true.
  MsgError= "OK"

  N= 0

  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,sEol)
    if(W(1:1)=='!')   cycle DoFile
    call AppendToEnd(L,W,sEoL)
    !
    if(W=="ENDINPUT") exit DoFile
    !
    !-----------------------------------------build a new solution model
    if(W=="ELECTROLYTE.MODEL") then

      !~ ELECTROLYTE.MODEL MYAQUEOUS-MODEL
        !~ MODEL IDEAL
        !~ SOLVENT H2O
        !~ POLE
          !~ H2O
          !~ OH-
          !~ H+
          !~ CA+2
          !~ NA+
          !~ NACL(AQ)
        !~ end
      !~ end

      !~ call MixModel_Zero(S)
      !
      call LinToWrd(L,W,sEol) 
      S%Name=trim(W)
      !
      if(iDebug>1) print '(/,A,/)',"Reading electrolyte model "//trim(S%Name)
      !
      DoSol: do
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        !
        call LinToWrd(L,W,sEol) !; if(iDebug>2) write(fTrc,'(A2,A)') "W=",trim(W)
        !
        if(W(1:1)=='!') cycle DoSol !skip comment lines
        call AppendToEnd(L,W,sEoL)
        !
        select case(W)
        !
        case("ENDINPUT"); exit DoFile
        !
        case("END","ENDELECTROLYTE.MODEL")
        !---------------------------------- end of one ELECTROLYTE block
        !-------------------------- save the T_SolModel value in vSolTmp
          if(S%nSpecies>0 .and. S%iSolvent>0) then
            N=N+1
            vSolTmp(N)%Name=      S%Name
            vSolTmp(N)%iActModel= S%iActModel
            vSolTmp(N)%iSolvent=  S%iSolvent
            vSolTmp(N)%nSpecies=  S%nSpecies
            allocate(vSolTmp(N)%vISpecies(S%nSpecies))
            do I=1,S%nSpecies
              vSolTmp(N)%vISpecies(I)= vIsolTmp(I)
            enddo
          end if
          !
          exit DoSol
        !----------------------/save the T_MixModel value in linked list
        !
        case("MODEL") !MODEL IDEAL / POLE / SITE / SPECIAL
        !---------------------------------------------------- read MODEL
          call LinToWrd(L,W,sEol)

          M= 0
          do i=1,size(vSolModelAct)
            if(trim(W)==trim(vSolModelAct(i))) then
              !SolModel%iActModel= i
              M= I
              exit
            end if
          enddo
          !~ if(INDEX(Solmodel_ModelList,trim(W)//"_")==0) &
          !~ & call Stop_(trim(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          if(M==0) then
            Ok= .false.
            MsgError=         trim(W)//" = UNKNOWN ACTIVITY MODEL !!!"
            return !----------------------------------------------return
          end if
          !
          !S%ActModel=trim(W)
          S%iActModel=M
          !
          if(iDebug>0) write(fTrc,'(2A)') "MODEL=",trim(W)
          !
          cycle DoSol
          !not really necessary,
          !... unless the current value of W has changed to POLE or MARGULES ...
        !----------------------------------------------------/read MODEL 
        !
        case("SOLVENT")
        !----------------------------------------------read SOLVENT name
          call LinToWrd(L,W,sEol)
          K= Species_Index(W,vSpc)
          if(K==0) then
            Ok= .false.
            MsgError=                            "error in SOLVENT name"
            return !----------------------------------------------------
          end if
          S%iSolvent= K
          !
          if(iDebug>0) write(fTrc,'(2A)') "SOLVENT=",vSpc(S%iSolvent)%NamSp
          !
        !-------------------------------------------/read SOLVENT name--
        !
        !------------------------------------------------read POLE block
        case("POLE")
          !
          nS= 0
          !
          if(iDebug>0) write(fTrc,'(A)') "<< Read_POLE_block"
          !
          DoPole: do
            !
            read(F,'(A)') L; call LinToWrd(L,W,sEol)
            if(W(1:1)=='!') cycle DoPole
            call AppendToEnd(L,W,sEoL)
            !
            if(trim(W)=="END") exit DoPole
            if(trim(W)=="ENDPOLE") exit DoPole
            !
            K= Species_Index(W,vSpc)
            !
            if(K<1) then
              Ok= .false.
              MsgError=                trim(W)//"species not in SPECIES"
              return !--------------------------------------------return
            end if
            if(vSpc(K)%Typ/="AQU") then
              Ok= .false.
              MsgError=                       trim(W)//" is not aqueous"
              return !--------------------------------------------return
            end if
            !
            nS= nS+1
            vIsolTmp(nS)= K
            !~ print *,"vIsolTmp(nS)=",vIsolTmp(nS)
            !
          enddo DoPole

          S%nSpecies= nS

        end select
        !
      end do DoSol
      !
    end if !if(W=="MIXTURE.MODEL" .or. W=="SOLUTION.MODEL")
  end do DoFile
  !
  close(F)

  if(N>0) then

    if(allocated(vSolModel)) deallocate(vSolModel)
    allocate(vSolModel(N))

    do I=1,N
      vSolModel(I)%Name=      vSolTmp(I)%Name
      vSolModel(I)%iActModel= vSolTmp(I)%iActModel
      vSolModel(I)%iSolvent=  vSolTmp(I)%iSolvent
      vSolModel(I)%nSpecies=  vSolTmp(I)%nSpecies
      nS= vSolModel(I)%nSpecies
      allocate(vSolModel(I)%vISpecies(nS))
      vSolModel(I)%vISpecies(1:nS)= vSolTmp(I)%vISpecies(1:nS)
    enddo

    if(iDebug>2) then
      print *,"ELECTROLYTE MODELS"
      do I=1,N
        S= vSolModel(I)
        print *,S%Name,vSolModelAct(S%iActModel),S%nSpecies
        print *,S%vISpecies(1:nS)
      enddo
    end if

  end if

  if(iDebug>0) write(fTrc,'(A,/)') "</ SolModel_Read"

end subroutine SolModel_Read

subroutine Spline_Init( &
& N, vX, vY, &
& S)
  use M_Solmodel_Vars,only: T_Spline
  use M_CMM_Spline

  integer,       intent(in) :: N
  real(dp),      intent(in) :: vX(N)
  real(dp),      intent(in) :: vY(N)
  !
  type(T_Spline),intent(out):: S

  real(dp):: B(N),C(N),D(N)

  call CMM_Spline_Compute(N,vX(1:N),vY(1:N),B,C,D)

  S%Dimm= N
  S%vX(1:N)= vX(1:N)
  S%vY(1:N)= vY(1:N)
  S%vSplineB(1:N)= B(1:N)
  S%vSplineC(1:N)= C(1:N)
  S%vSplineD(1:N)= D(1:N)
  
  return
end subroutine Spline_Init

end module M_Solmodel_Read

