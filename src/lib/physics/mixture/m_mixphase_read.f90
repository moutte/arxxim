module M_MixPhase_Read
!--
!-- tools for reading mixture phases (phase of variable composition)
!-- from text files,
!-- and build an array of T_MixPhase (container for a mixture phase)
!--
  use M_Trace,only: iDebug,fTrc,T_,Warning_
  use M_Kinds
  use M_T_MixPhase,only: T_MixPhase
  !
  implicit none
  !
  private
  !
  public:: MixPhase_BuildLnk
  public:: MixPhase_LnkToVec
  !
  public:: T_LnkFas
  !
  type:: T_LnkFas
    type(T_MixPhase):: Value
    type(T_LnkFas),pointer ::Next
  end type T_LnkFas

contains

subroutine LnkFas_Build(First,Input_,L,P)
  logical               :: First
  type(T_MixPhase)      :: Input_
  type(T_LnkFas),pointer:: L
  type(T_LnkFas),pointer:: P
  !
  if(First) nullify(L)
  !
  if(First) then
    allocate(L)     ; nullify(L%next)     ; L%Value=Input_     ; P=>L
  else
    allocate(P%next); nullify(P%next%next); P%next%Value=Input_; P=>P%next
  end if
  !
end subroutine LnkFas_Build

subroutine MixPhase_BuildLnk(vSpc,vMixModel,N,LnkFas,Ok,MsgError)
!--
!-- scan the MIXTURE block(s)
!--
  use M_IoTools
  use M_Files,     only: NamFInn
  use M_T_Species, only: T_Species,Species_Index
  use M_T_MixModel !,only: T_MixModel,MixModel_Index,MixModel_XPoleToXSite
  use M_T_MixPhase,only: T_MixPhase !,MixPhase_ActPole
  !
  type(T_Species), intent(in) :: vSpc(:)
  type(T_MixModel),intent(in) :: vMixModel(:)
  integer,         intent(out):: N
  type(T_LnkFas),  pointer    :: LnkFas
  logical,         intent(out):: Ok
  character(*),    intent(out):: MsgError
  !
  type(T_LnkFas),pointer:: pFas
  !
  type(T_MixPhase)  :: MixFas
  type(T_MixModel)  :: MixModel
  character(len=255):: L !,sList0 !=elements in Database
  character(len=80) :: W !,V1,V2
  logical :: EoL !,Ok
  integer :: iS,iEl,I,J,K !,M,iPhas,iCp,iSp,I,J,K,M,nCp,fInn !,iEl !,J
  real(dp):: X1 !,Act
  integer :: ios,fInn
  real(dp),allocatable:: vXAtom(:)
  !
  if(idebug>1) write(fTrc,'(/,A)') "< MixPhase_BuildLnk"
  !
  Ok= .true.
  MsgError= "Ok"
  !
  call GetUnit(fInn)
  open(fInn,file=trim(NamFInn))
  !Ok=.false.
  !
  N=0
  !
  DoFile: do
    !
    read(fInn,'(A)',iostat=ios) L ; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)
    !
    case("ENDINPUT"); exit DoFile
    !
    case("SOLUTION","MIXTURE") !mixture name and composition
      !
      if(trim(W)=="SOLUTION") &
      & call Warning_(trim(W)//" soon Obsolete, better use MIXTURE !!!")
      !
      !SOLUTION SOL1
      !  MODEL FELDSPAR_SS
      !  COMPOSITION
      !    POLE2 0.5
      !    POLE3 0.5
      !  endCOMPOSITION
      !endSOLUTION
      !../..
      call LinToWrd(L,W,EoL) !-> phase name, SOL1

      if(Species_Index(W,vSpc)>0) then
        Ok= .false.
        MsgError="Mixture name should not match with any species name !!!"
        return !--------------------------------------------------return
      end if

      MixFas%Name=trim(W)
      !
      DoMixFas: do
        !
        read(fInn,'(A)',iostat=ios) L  ;  if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoMixFas
        call AppendToEnd(L,W,EoL)
        !
        select case(W)
        !
        case("ENDINPUT")                       ; exit DoFile
        !
        case("END","ENDSOLUTION","ENDMIXTURE") ; exit DoMixFas
        !
        !------------------------------------------------read model name
        case("MODEL")
        ! read model name and check it is among the model base vMixModel
        ! -> return its index, iS, in array vMixModel_
          call LinToWrd(L,W,EoL)
          !
          !find the mixture name in vMixModel%Name
          iS= MixModel_Index(vMixModel,W)

          if(iS==0) then
            Ok= .false.
            MsgError=                 trim(W)//"= MIXTURE.MODEL UNKNOWN"
            return !----------------------------------------------return
          end if

          if(idebug>1) write(fTrc,'(2A)') "MODEL=", trim(W)
          MixFas%iModel=   iS
          !
          !=MixModel
        !endcase("MODEL")
        !-----------------------------------------------/read model name
        !
        !-----------------------------------------/read phase compositon
        case("COMPOSITION")
          !
          if(idebug>1) write(fTrc,'(A)') "!!! !!!Read_COMPOSITION"
          !
          if(MixFas%iModel==0) then
            Ok= .false.
            MsgError=             "must define the model beforehand ..."
            return !----------------------------------------------return
          end if
          !
          MixModel=vMixModel(MixFas%iModel)
          !
          !! select case(trim(MixModel%Model))

          !! case("IDEAL","POLE","MOLECULAR","FELSPAR")
          !-------------------input on several lines, one per end-member
            do
              !
              read(fInn,'(A)') L; call LinToWrd(L,W,EoL)
              !-> end-member name, or "!", or "END", or "ENDCOMPOSITION"
              if(W=="!") cycle !-> comment line
              call AppendToEnd(L,W,EoL)
              if(W=="END") exit
              if(W=="ENDCOMPOSITION") exit
              !
              I= Species_Index(trim(W),vSpc)
              !
              if(I==0) then
                Ok= .false.
                MsgError=   trim(W)//" <-end-member NOT FOUND in species list"
                return !------------------------------------------return
              end if
              !
              !--- search for the end-member with index I
              !--- in the end-member list, vIPole(:), of mixing model MixModel
              !--- and read mole fraction
              K=0
              do
                K=K+1
                if(K>MixModel%NPole) exit
                !
                if(I==MixModel%vIPole(K)) then
                  !
                  call LinToWrd(L,W,EoL) !-> rest of line -> composition
                  call WrdToReal(W,X1)
                  MixFas%vLPole(K)= (X1>Zero)
                  MixFas%vXPole(K)=  X1 !save to MixFas%composition
                  !
                  if(idebug>1) write(fTrc,'(2(A,I3))') &
                  & "K=",K,", (MixFas%iModel)%vIPole(K)=",MixModel%vIPole(K)
                  !
                  exit
                end if
              end do
              !---/ search for the end-member named W
              !
              if(K>MixModel%NPole) then
                Ok= .false.
                MsgError=      trim(W)//" is NOT AN end-MEMBER of this model"
                return !------------------------------------------return
              end if
              !
            end do
          !endcase("IDEAL","POLE","MOLECULAR","FELSPAR")
          !
          !-------------------------for SITE models: read site fractions
          !!! case("SITE")
          !!! ! input should be in number of atoms (not fractions ...) (???)
          !!! ! one line per site, reproducing the structure
          !!! ! of the SITE block of the MIXTURE.MODEL
          !!!   M=0
          !!!   K=0
          !!!   DoSite: do !M=1,MixModel%NSite
          !!!     read(fInn,'(A)') L
          !!!     call LinToWrd(L,W,EoL)
          !!!     call AppendToEnd(L,W,EoL)
          !!!     !
          !!!     if(W=="!")              cycle DoSite !-> comment line
          !!!     if(W=="END")            exit DoSite
          !!!     if(W=="ENDCOMPOSITION") exit DoSite
          !!!     !
          !!!     M= M+1 !site index
          !!!     if(M > MixModel%NSite) exit DoSite
          !!!     !
          !!!     do
          !!!     !-- scan the rest of the line
          !!!     !-- -> retrieve composition data
          !!!       call LinToWrd(L,W,EoL)
          !!!       if(W=="!") cycle DoSite !end of line -> read a new line
          !!!       call WrdToReal(W,X1) !save to MixFas%composition
          !!!       !MixFas%vLPole(K)=(X1>Zero)
          !!!       !
          !!!       K=K+1
          !!!       if(K>MixModel%NAtom) exit DoSite
          !!!       ! save to vMixFas(J)%composition
          !!!       MixFas%vXAtom(K)= X1 / MixModel%vMulti(M)
          !!!       !
          !!!       if(idebug>1) &
          !!!       & write(fTrc,'(I3,1X,2A,G15.6)') K,MixModel%vNamAtom(K),"=",X1
          !!!       !
          !!!     end do
          !!!   end do DoSite
          !!! !endcase("SITE")
          !------------------------/for SITE models: read site fractions
          !
          !! end select !case (trim(MixModel%Model))
          !
          if(MixModel%Model==Mix_Site) then
          ! to initialize vLPole, calculate the ideal activities of the endmembers
            !
            allocate(vXAtom(MixModel%NAtom))
            call MixModel_XPoleToXSite( &
            & MixModel, &
            & MixFas%vXPole(1:MixModel%NPole), &
            & vXAtom(1:MixModel%NAtom))
            !
            !------------------------------------------------------trace
            if(idebug>1) then
              !
              write(fTrc,'(/,A,/)') "POLE FRACTIONS"
              do i=1,MixModel%NPole
                write(fTrc,'(G12.4,1X)',advance="NO") MixFas%vXPole(I)
              end do
              write(fTrc,*)
              !
              write(fTrc,'(/,A,/)') "SITE FRACTIONS"
              I= MixModel%vIAtomSite(1)
              do iEl=1,MixModel%NAtom
                J= MixModel%vIAtomSite(iEl)
                if(J/=I) then
                  I= J
                  write(fTrc,*)
                end if
                write(fTrc,'(A,1X,G12.4,1X)',advance="NO") &
                & MixModel%vNamAtom(iEl), vXAtom(iEl)
              end do
              write(fTrc,*)
              !
            end if
            !-----------------------------------------------------/trace

            if(idebug>1) &
            & write(fTrc,'(/,A,/)') "for site models, initialize vLPole"
            !
            do K=1,MixModel%NPole
              !if(MixModel%vHasPole(K)) then
                X1= MixModel_Site_ActivIdeal(MixModel,K,vXAtom)
                if(idebug>1) write(fTrc,'(A,1X,G15.6)') MixModel%vNamPole(K),X1
                MixFas%vLPole(K)= (X1>Zero)
              !else
              !  MixFas%vLPole(K)= .false.
              !end if
            end do
            !
            if(idebug>1) write(fTrc,'(/,A,/)') "END vLPole"
            !
          end if
          !
          !-------------------------------------------------------------
          !-- for the model of the phase, MixModel%Model, to be accepted,
          !-- must check whether all 'active' endmembers of MixModel%Model
          !-- are found in vSpc
          !-------------------------------------------------------------
          do K=1,MixModel%NPole
            if(MixFas%vLPole(K)) then !only for "active" endmembers
              !
              J=Species_Index(MixModel%vNamPole(K),vSpc)
              !
              if(J==0) then
                Ok= .false.
                MsgError= trim(MixFas%Name) &
                & //" SPECIES NOT IN BASE:"//trim(MixModel%vNamPole(K))
                return !------------------------------------------return
              end if
              !
              MixModel%vIPole(K)=J
              !
            end if
          end do
          N=N+1
          !
          if(idebug>1) &
          & write(fTrc,'(I3,1X,A24,A24)') N, MixFas%Name, MixModel%Name
          call LnkFas_Build(N==1,MixFas,LnkFas,pFas)
          !
        !endcase("COMPOSITION")
        !----------------------------------------/read phase composition
        !
        end select !case(W)!
        !
        if(allocated(vXAtom)) deallocate(vXAtom)
        !
      end do DoMixFas
      !endcase("MIXTURE")
    end select !case(W)
  end do DoFile
  !
  close(fInn)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ MixPhase_BuildLnk"
  !
  return
end subroutine MixPhase_BuildLnk
  !
subroutine MixPhase_LnkToVec(Lnk,V)
  use M_T_MixPhase,only: T_MixPhase
  !
  type(T_LnkFas),               pointer    :: Lnk
  type(T_MixPhase),dimension(:),intent(out):: V
  !
  type(T_LnkFas),pointer:: pCur,pPrev
  integer::I
  !
  pCur=>Lnk
  I=0
  do while (associateD(pCur))
    I=I+1
    V(I)=pCur%Value
    pPrev=>pCur; pCur=> pCur%next; deallocate(pPrev)
  end do
  !
  return
end subroutine MixPhase_LnkToVec

end module M_MixPhase_Read
