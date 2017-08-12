module M_MixModel_Read
  use M_Kinds
  use M_Trace,   only: iDebug,fTrc,T_
  use M_T_MixModel,only: T_MixModel
  implicit none
  private
  !
  public:: MixModel_BuildLnk
  public:: MixModel_LnkToVec
  public:: MixModel_Check
  !
  public:: T_LnkMixModel
  !
  type:: T_LnkMixModel
    type(T_MixModel)::Value
    type(T_LnkMixModel),pointer::Next
  end type T_LnkMixModel
  !
contains

subroutine MixModel_BuildLnk( &
& vSpc, &
& Lnk,N,Ok,MsgError)
!--
!-- reading MIXTURE.MODEL blocks from file
!-- -> build a linked list of T_MixModel
!--
  use M_IOTools
  use M_T_MixModel
  use M_Files,    only: NamFSol
  use M_T_Species,only: T_Species,Species_Index
  !---------------------------------------------------------------------
  type(T_Species),    intent(in) :: vSpc(:)
  type(T_LnkMixModel),pointer    :: Lnk
  integer,            intent(out):: N
  logical,            intent(out):: Ok
  character(*),       intent(out):: MsgError
  !---------------------------------------------------------------------
  type(T_LnkMixModel),pointer:: pCur
  !
  character(len=512):: L,W
  character(len=3)  :: Str
  character(len=80) :: SiteList
  !character(len=255):: sFormul
  logical            :: sEol
  logical            :: MargulTheriak= .false.
  integer            :: ios,F,K,J,iEl,iS,iP,offset_,iMargIndex !,M1,M2
  type(T_MixModel)   :: S
  type(T_Margul)     :: M
  logical            :: vHasPole(MaxPole)
  real(dp),dimension(dimV)::vX
  real(dp),allocatable:: vXpol(:),vXatom(:)
  real(dp):: Y
  !
  !----------------------------------------------for reading SITE models
  ! real(dp)::NormCoeff
  !
  integer,dimension(0:MaxSite):: vNEle
  ! vNEle = nr of possible elements on each site, vNEle<=MaxMulti
  !       -> gives structure of input composition vector
  ! vNEle(0)=0,
  !       -> to enable i=i+vNEle(iSite-1)
  !          for localization of i in composition vector ...
  !
  character(3*MaxMulti):: vSiteString(1:MaxSite)
  !--------------------------------------------------------------------/
  !
  if(idebug>1) write(fTrc,'(/,A)') "<-----------------MixModel_BuildLnk"
  !
  Ok= .true.
  MsgError= "Ok"
  !
  N=0
  nullify(Lnk)
  !
  call GetUnit(F)
  open(F,file=trim(NamFSol))

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
    if(W=="MIXTURE.MODEL" .or. W=="SOLUTION.MODEL") then
    !
    ! keyword SOLUTION.MODEL should become obsolete for symmetric models
    ! MIXTURE.MODEL is the preferred keyword for symmetric models,
    ! SOLUTION.MODEL should be for assymetric solvent / solute mixtures
    !
      call MixModel_Zero(S)
      vHasPole(:)= .false.
      !
      call LinToWrd(L,W,sEol) !MIXTURE.MODEL OLIVINE ->S%Name="OLIVINE"
      S%Name=trim(W)
      !
      if(idebug>1) write(fTrc,'(/,2A)') "MIXTURE=",trim(S%Name)
      if(iDebug>2) print '(/,A,/)',"Reading mixture model "//trim(S%Name)
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
        case("END","ENDMIXTURE.MODEL","ENDSOLUTION.MODEL")
        !------------------------------------- end of one MIXTURE block,
        !------------------- -> save the T_MixModel value in linked list
          if((S%nPole>1) .and. &
          &  (count(vHasPole(:) .EQV. .true.)>1)) then
            if(iDebug>2) print *,"------------------------ACCEPTED"
            !
            if(.not. all(vHasPole(1:S%nPole))) &
            & call MixModel_Clean(vHasPole, S)
            !
            N=N+1
            if(N==1) then !--fill the linked list
              allocate(Lnk); nullify(Lnk%next)
              Lnk%Value=S; pCur=>Lnk
            else
              allocate(pCur%next); nullify(pCur%next%next)
              pCur%next%Value=S; pCur=>pCur%next
            end if
          else
            if(iDebug>2) print *,"------------------------REJECTED"
          end if
          !
          exit DoSol
        !---------------------/ save the T_MixModel value in linked list
        !
        case("MODEL") !MODEL IDEAL / POLE / SITE / SPECIAL
        !-----------------------------------------------read MODEL index
          ! Mix_Molecular= 1, & !
          ! Mix_Felspar=   2, & !
          ! Mix_Site=      3, & !
          ! Mix_Special=   4, & !
          !
          call LinToWrd(L,W,sEol)
          !
          select case(trim(W))
          case("MOLECULAR","IDEAL","POLE")  ;  S%Model= Mix_Molecular ! "MOLECULAR"
          case("FELSPAR")                   ;  S%Model= Mix_Felspar   ! "FELSPAR"
          case("SITE")                      ;  S%Model= Mix_Site      ! "SITE"
          case("SPECIAL")                   ;  S%Model= Mix_Special   ! "SPECIAL"
          case("EXCHANGE")                  ;  S%Model= Mix_Exchange
          case("SURFACE")                   ;  S%Model= Mix_Surface
          case default
            Ok= .false.
            MsgError=trim(W)//"= Unknown MODEL in MIXTURE.MODEL"
          end select
          !
          if(idebug>1) write(fTrc,'(2A)') "MODEL=",trim(W)
          !
          S%NSite=1
          if(.not. sEol) then !=former version, read MULTIplicity
            call LinToWrd(L,W,sEol)
            call WrdToInt(W,S%vMulti(1))
          end if
          !
          cycle DoSol
          !not really necessary,
          !... unless the current value of W has changed to POLE or MARGULES ...
        !----------------------------------------------------/read MODEL
        !
        case("MULTI") ! needed only for MOLECULAR models with Multi/=1
        !----------------------------------------------read MULTIPLICITY
          call LinToWrd(L,W,sEol)
          call WrdToInt(W,K)
          if(K<0 .or. K>MaxMulti) then
            Ok= .false.
            MsgError=                            "error in MULTIPLICITY"
            return !----------------------------------------------return
          end if
          !
          S%vMulti(1)= K
          !
          if(idebug>1) write(fTrc,'(A,I3)') "MULTI=",S%vMulti(1)
          !
        !---------------------------------------------/read MULTIPLICITY
        !
        case("SITE")
        !------------------------------------------------read SITE block
          !
          ! example
          !   SITE
          !     !SiteName / Multiplicty / Elements
          !     T   2  SI_AL_
          !     M   1  MG_FE_AL_VAC
          !     VA4 4  MG_FE_AL_
          !   end SITE
          !
          S%nSite= 0
          vNEle=   0
          S%NAtom= 0
          iEl=     0
          SiteList= ""
          !
          ! if(trim(S%Model)/="SITE") then
          if(S%Model /= Mix_Site) then
            Ok= .false.
            MsgError=             "only SITE mixtures have a SITE block"
            return !----------------------------------------------return
          end if
          !
          DoSite: do
            !
            read(F,'(A)') L; call LinToWrd(L,W,sEol)
            if(W(1:1)=='!') cycle DoSite
            call AppendToEnd(L,W,sEoL)
            if(W=="END") exit DoSite
            if(W=="ENDSITE") exit DoSite
            !
            S%nSite=S%nSite+1
            !
            if(S%nSite>MaxSite) then
              Ok= .false.
              MsgError=                  "in SITE block: TOO MANY SITES"
              return !--------------------------------------------return
            end if
            !
            !----------first word of line= name of site, up to 3 lettres
            call Str_Append(W,3)
            S%vNamSite(S%nSite)=trim(W)
            SiteList= trim(SiteList)//trim(W)
            !
            if(iDebug>2) write(fTrc,*) &
            & "NamSite,SiteList= ",S%vNamSite(S%nSite),trim(SiteList)
            !------------------------------------------------/first word
            !
            !------------second word of line= site multiplicity, integer
            call LinToWrd(L,W,sEol)
            call WrdToInt(W,S%vMulti(S%nSite))
            !
            if(S%vMulti(S%nSite)>MaxMulti) then
              Ok= .false.
              MsgError=trim(S%vNamSite(S%nSite))//"= Multiplicity too high"
              return !--------------------------------------------return
            end if
            !-----------------------------------------------/second word
            !
            !---------------------- list of elements potentially present
            if(S%Model == Mix_Site) then
              !
              if(sEol) then
                Ok= .false.
                MsgError= trim(S%vNamSite(S%nSite))//" NO ATOM LIST ???"
                return !------------------------------------------return
              end if
              !
              call LinToWrd(L,W,sEol)
              if(iDebug>2) write(fTrc,*) trim(W)
              !
              vNEle(S%nSite)= len_trim(W)/3 ! how many elements on the site
              vSiteString(S%nSite)= trim(W)
              !
              if(S%NAtom + vNEle(S%nSite) > MaxAtom) then
                Ok= .false.
                MsgError=         trim(W)//": Too Many Elements in SITE"
                return !------------------------------------------return
              end if
              !
              do K=1,vNEle(S%nSite)
                !
                S%vIAtomSite(K+S%NAtom)= S%nSite
                S%vAtomMulti(K+S%NAtom)= S%vMulti(S%nSite)
                !
                !----name of the site fraction
                S%vNamAtom(K +S%NAtom)= trim(W(3*K-2:3*K))//S%vNamSite(S%nSite)
                !
              end do
              !
              S%NAtom= S%NAtom + vNEle(S%nSite) ! total nr site fractions
              !
            end if
            !----------------------/list of elements potentially present
            !
          end do DoSite
          !check that sum(S%vMulti(1:S%nSite)) <= MaxElPole
        !-----------------------------------------------/read SITE block
        !
        !----------------------------------------------- read POLE block
        case("POLE")
          !
          S%nPole=     0
          S%tPoleAtom= 0
          !
          if(idebug>1) write(fTrc,'(A)') "<<------------Read_POLE_block"
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
            !old! if(K<1) then
            !old!   Ok= .false.
            !old!   MsgError=        trim(W)//" not in base-> model not valid"
            !old!   return !----------------------------------------return
            !old! end if
            !
            !if(idebug>1) &
            !& write(fTrc,'(A,I3,2A)') "Species_Index=",K," <-POLE= ",trim(W)
            !
            if(S%nPole<MaxPole) then
              S%nPole=S%nPole+1
            else
              Ok= .false.
              MsgError=                           "TOO MANY end-members"
              return !--------------------------------------------return
            end if
            !
            vHasPole(S%nPole)= K>0
            !
            if(vHasPole(S%nPole)) then
              !
              !----------------------------------verify the species type
              if(S%nPole==1) then
                S%Typ= trim(vSpc(K)%Typ)
                !
                select case(S%Model)
                case(Mix_Surface)
                  if(S%Typ/="SUR") then
                    Ok= .false.
                    MsgError= "should have a SUR(face) Species"
                    return !--------------------------------------return
                  end if
                case(Mix_Exchange)
                  if(S%Typ/="EXC") then
                    Ok= .false.
                    MsgError= "should have a EXC(hange) Species"
                    return !--------------------------------------return
                  end if
                end select
              else
                if(trim(vSpc(K)%Typ)/=trim(S%Typ)) then
                  Ok= .false.
                  MsgError= "END-members should be same type, either GAS or MIN"
                  return !----------------------------------------return
                end if
              end if
              !---------------------------------/verify the species type
              !
              S%vIPole(S%nPole)=K ! index of end-member name in species list
              !
            else
              !
              S%vIPole(S%nPole)=0
              !
            end if
            !
            S%vNamPole(S%nPole)=trim(W)
            !
            if(idebug>1) write(fTrc,'(A5,I3,A1,I3,A1,A)') &
            & "POLE=",S%nPole,T_,K,T_,trim(W)
            !
            !--------------for SITE models, reading pole / site relation
            if(S%Model == Mix_Site) then
              !
              ! scan the rest of the input line
              ! to retrieve the site repartition for the end-member
              !
              do iS=1,S%NSite
                if(sEol) then
                  Ok= .false.
                  MsgError=trim(S%vNamPole(S%nPole))//" NO ATOM LIST ???"
                  return !----------------------------------------return
                end if
                !
                call LinToWrd(L,W,sEol)
                !!!
                !!!M2 2 MG_FE_     -> vElSite(1),vElSite(2),            ( NEle(1)=2 )
                !!!M1 1 MG_FE_AL_  -> vElSite(3),vElSite(4),vElSite(5)  ( NEle(2)=3 )
                !!!T  4 AL_SI_     -> etc ...
                !!!A  1 K__
                !!!OH 2 OH_
                !!!
                !!!SIDEROPHYLLITE FE_FE_ AL_ AL_AL_SI_SI_ K__ OH_OH_ !->
                !!!
                if(iDebug>2) write(fTrc,'(2A)') "Site String= ",trim(W)
                !
                if(len_trim(W)/=3*S%vMulti(iS)) then
                  Ok= .false.
                  MsgError=trim(W)//" has length not consistent with vMulti"
                  return !----------------------------------------return
                end if
                !
                if(iS==1) then
                  offset_= 0
                else
                  offset_= offset_ +vNEle(iS-1)
                end if
                !
                do K=1,S%vMulti(iS)
                  !
                  str=trim(W(3*K-2:3*K)) ! extract the element name
                  !
                  iEl= index(vSiteString(iS),str)
                  if(iEl<=0) then
                    Ok= .false.
                    MsgError=       trim(str)//" <-Atom not Found   !!!"
                    return !--------------------------------------return
                  end if
                  iEl= (iEl+2)/3
                  S%tPoleAtom(S%nPole,iEl+offset_)= S%tPoleAtom(S%nPole,iEl+offset_) +1
                  !
                end do
                !
              end do !!do iS=1,S%NSite
              !
            end if
            !-------------/for SITE models, reading pole / site relation
            !
          end do DoPole
          !
          !----------------------------------------normalization factors
          if(S%Model == Mix_Site) then
            allocate(vXpol(S%NPole))
            allocate(vXatom(S%NAtom))
            do iP=1,S%NPole
              vXpol(:)=  Zero
              vXpol(iP)= One
              call MixModel_XPoleToXSite(S,vXpol(:),vXatom(:))
              Y= One
              do iEl= 1,S%NAtom
                if(S%tPoleAtom(iP,iEl)/=0) Y= Y *vXAtom(iEl)**S%vAtomMulti(iEl)
                !S%tPoleAtom(iP,iEl)
              end do
              S%vPoleCoeff(iP)= Y
            end do
            deallocate(vXpol)
            deallocate(vXatom)
            !
            !! print *,"normalization factors:"
            do iP=1,S%NPole
              if(S%vPoleCoeff(iP)==Zero) then
                Ok= .false.
                MsgError=        "Problem with normalization factor !!!"
                return !------------------------------------------return
              end if
              S%vPoleCoeff(iP)= One /S%vPoleCoeff(iP)
            end do
            !! pause
          end if
          !---------------------------------------/normalization factors
          !
          if(S%Model /= Mix_Site) S%NAtom= S%nPole
          !
          !--------------------------------------------------------trace
          if(iDebug>2 .and. S%Model==Mix_Site) then

            write(fTrc,'(/,A)') "tPoleAtom"
            do iS=1,S%NPole
              do iEl=1,S%NAtom
                write(fTrc,'(I7,1X)',advance="NO") S%tPoleAtom(iS,iEl)
              end do
              write(fTrc,*)
            end do

            write(fTrc,'(A)') "vAtomMulti"
            do iEl=1,S%NAtom
              write(fTrc,'(I7,1X)',advance="NO") S%vAtomMulti(iEl)
            end do
            write(fTrc,*)

            write(fTrc,'(A)') "tPoleAtom"
            do iS=1,S%NPole
              do iEl=1,S%NAtom
                write(fTrc,'(F7.2,1X)',advance="NO") &
                & S%tPoleAtom(iS,iEl) /real(S%vAtomMulti(iEl))
              end do
              write(fTrc,*)
            end do

            write(fTrc,'(/,A)') "Normalization:"
            do iP=1,S%NPole
              write(fTrc,'(A,1X,F7.2)') S%vNamPole(iP), S%vPoleCoeff(iP)
            end do

          end if
          !-------------------------------------------------------/trace
          !
          if(idebug>1) write(fTrc,'(A)') "<<-----------/Read_POLE_block"
        !-----------------------------------------------/read POLE block
        !
        case("MARGULES")
        !--------------------------------------------read MARGULES block
          if(S%Model == Mix_Special) then
            Ok= .false.
            MsgError=           "SPECIAL MODEL -> NO MARGULES block !!!"
            return !----------------------------------------------return
          end if
          if(idebug>1) write(fTrc,'(A)') "<<--------Read_MARGULES_block"
          !
          MargulTheriak= .false.
          if (.not. sEol) then
            call LinToWrd(L,W,sEol)
            if (trim(W)=="THERIAK") MargulTheriak= .true.
          end if
          !
          !--examples---------------------------------------------------
          !
          ! <case of 'molecular' mixture>
          !
          ! <old format, a la Capitani>
          !  MARGULES
          !    GROSSULAR PYROPE
          !    & 112     21560.    18.79    .10
          !    & 122     69200.    18.79    .10
          !    GROSSULAR ALMANDINE
          !    & 112     20320.     5.08    .17
          !    & 122      2620.     5.08    .09
          !  endMARGULES
          !
          ! <new format, a la Berman>
          ! <indices used in ijkl codes refer to order in the end-member list>
          !  MODEL MOLECULAR
          !  POLE
          !    GROSSULAR
          !    PYROPE
          !    ALMANDINE
          !  end
          !  MARGULES                          !                 vDegree(:)=
          !    112     21560.    18.79    .10  !GROSSULAR PYROPE       2,1,0
          !    122     69200.    18.79    .10  !GROSSULAR PYROPE       1,2,0
          !    113     20320.     5.08    .17  !GROSSULAR ALMANDINE    2,0,1
          !    133      2620.     5.08    .09  !GROSSULAR ALMANDINE    1,0,2
          !    ../..
          !  end MARGULES
          !
          ! <case of 'SITE-mixing' mixture, with Margules param's for each site>
          !  MODEL SITE
          !  SITE
          !    C 2 CA_MG_
          !  end SITE
          !  MARGULES SITE
          !    C 223     1307.40      0.00      0.01 !G_xs+= W223 *X_MG^2 *X_FE
          !    C 233     2092.40      0.00      0.06 !G_xs+= W233 *X_MG   *X_FE^2
          !    C 112    85529.00     18.79      0.21 !G_xs+= W112 *X_CA^2 *X_MG
          !    C 122    50874.90     18.79      0.02 !G_xs+= W122 *X_CA   *X_MG^2
          !    ../..
          !  end
          !
          if(MargulTheriak) then
          
          doMargulTheriak: do
          
            !type tMarg !Margules parameter for up to Ternary mixing
            !  integer::Dim !2 or 3
            !  integer,dimension(1:3)::iPole,Power
            !  real::WG,WH,WS,WV,WCp,WK
            !end type tMarg
            !
            read(F,'(A)') L ; call LinToWrd(L,W,sEol)
            !
            if(W(1:1)=="!") cycle doMargulTheriak
            call AppendToEnd(L,W,sEoL)
            if(W=="END" .or. W=="ENDMARGULES") exit doMargulTheriak
            !
            !M%nMargul=M%nMargul+1
            !
            ! <old format, a la Capitani>
            !  MARGULES THERIAK
            !    GROSSULAR PYROPE
            !    & 112     21560.    18.79    .10
            !    & 122     69200.    18.79    .10
            !    GROSSULAR ALMANDINE
            !    & 112     20320.     5.08    .17
            !    & 122      2620.     5.08    .09
            !  endMARGULES
            !
            ! wrk ! if(W(1:1)/='&') then
            ! wrk ! ! scan the GROSSULAR PYROPE line
            ! wrk ! ! -> dimension of M, indices of Poles
            ! wrk !   !
            ! wrk !   M%NPole=0
            ! wrk !   do
            ! wrk !     K=0
            ! wrk !     do 
            ! wrk !       K=K+1 
            ! wrk !       if(W==trim(S%vNamPole(K))) then
            ! wrk !         !find occurrence of "GROSSULAR" in end member list
            ! wrk !         M%NPole=M%NPole+1
            ! wrk !         M%iPole(M%NPole)=K
            ! wrk !       end if
            ! wrk !       if(K>S%nPole) exit
            ! wrk !     end do
            ! wrk !     if(sEol) exit
            ! wrk !     call LinToWrd(L,W,sEol)
            ! wrk !   end do !end scan the Names
            ! wrk !   !
            ! wrk !   if(M%Npole<2) call Stop_(trim(W)//" <-endMEMBERS NOT FOUND")
            ! wrk !   cycle doMargulTheriak
            ! wrk !   !
            ! wrk ! end if
            ! wrk ! !
            ! wrk ! if(W(1:1)=='&') then
            ! wrk ! !read Margules Parameters
            ! wrk !   !
            ! wrk !   S%nMarg=S%nMarg+1
            ! wrk !   call LinToWrd(L,W,sEol)
            ! wrk !   M%Power(1:M%NPole)=0 !initial values
            ! wrk !   do K=1,len_trim(W) !processing the k1k2... string
            ! wrk !     J=CarToInt(W(K:K))
            ! wrk !     if(J>0.and.J<=M%NPole) M%Power(J)=M%Power(J)+1
            ! wrk !   end do
            ! wrk !   !1222 1122 1112 1222 1223 1233
            ! wrk !   call ReadRValsV(L,K,vX)
            ! wrk !   M%WH= vX(1)  ;  M%WS=vX(2)  ;  M%WV=vX(3)  ;  M%WCp=vX(4)
            ! wrk !   M%Kohler=FLOOR(vX(5))
            ! wrk !   !M%WG= M%WH &
            ! wrk !   !&   + M%WCP*(T-T0) &
            ! wrk !   !&   -(M%WS +M%WCP*log(T/T0))*T &
            ! wrk !   !&   + M%WV *P
            ! wrk !   S%vMarg(S%nMarg)=M
            ! wrk !   !
            ! wrk ! end if

          end do doMargulTheriak
          
          else
          
          doMargul: do
            ! cf M_T_MixModel :
            !
            !<old tMargul type>
            ! type tMargul !Margules parameter for up to quaternary mixing
            !   integer::Dim !2 to 4
            !   integer,dimension(1:4)::iPole,Power
            !   real::WG,WH,WS,WV,WCp,WK
            ! end type tMargul
            !
            !<new tMargul type>
            !  better for using systematic evaluation & derivation
            !  of monomials & polynomials
            !
            ! type tMargul !Margules parameter
            !   integer :: vDegree(1:MaxAtom)
            !   real(dp):: WG,WH,WS,WV,WCp,WK
            ! end type tMargul
            !-> max' possible dimension of vDegree is
            !   max number of site fractions (for a site model)
            !   or max number of end-members (for a molecular model)
            !
            !-> vDegree is given a fixed dimension MaxAtom
            !
            !-> for a site model, vDegree informs
            !   the degree of the site fraction
            !   and also which site is the Margules for
            !
            read(F,'(A)') L; call LinToWrd(L,W,sEol)
            !
            if(W(1:1)=="!") cycle doMargul
            call AppendToEnd(L,W,sEoL)
            if(W=="END" .or. W=="ENDMARGULES") exit doMargul
            !
            iMargIndex= 1
            !
            !--------------------------------------------for SITE models
            !-----------------------find index of string W in S%vNamSite
            if(S%Model == Mix_Site) then
              !
              if(index(W,"-")==2) then
                iMargIndex= CarToInt(W(1:1))
                W= W(3:len_trim(W))
                !! print *,"I, W=", iMargIndex,trim(W)  ; call pause_
              else
                call Str_Append(W,3)
                J= index(SiteList,trim(W))
                !J= 0
                !do iS= 1,S%nSite
                !  if(trim(W) == trim(S%vNamSite(iS))) then
                !    J= iS
                !    exit
                !  end if
                !end do
                if(J==0) then
                  Ok= .false.
                  MsgError=      trim(W)//" is not a site of this model"
                  return !----------------------------------------return
                end if
                J= (J+2)/3
                !
                iMargIndex= J
                !
                call LinToWrd(L,W,sEol)
              end if
            end if
            !-------------------------------------------/for SITE models
            !
            !------------- retrieve the power coding string, e.g. 113 --
            !
            ! the "power coding string", examples:
            !
            !   for molecular model:
            !   indices refer to e-m order in POLE
            !          112 -> vDegree= 2,1,0,0,..,0
            !          133 ->          1,0,2,0,..,O
            !
            !   for site model:
            !   indices refer to site-fraction order in SITE
            !     T   2  SI_AL_
            !     M   1  MG_FE_AL_VAC
            !     VA4 4  MG_FE_AL_
            !  -> site fraction order: Si_T,Al_T,Mg_M,Fe_M,Al_M, ...
            !   MARGULES
            !   !input (12, 13, ..) refer to indices in site
            !   ! T-1=Si_T, T-2=Al_T, etc.
            !     T   12   -> vDegree= 1,1, 0,0,0,0, 0,0,0, 0,0,0,....
            !     M   13   -> vDegree= 0,0, 1,0,3,0, 0,0,0, 0,0,0,....
            !     VA4 13   -> vDegree= 0,0, 0,0,0,0, 1,0,3, 0,0,0,....
            !
            if(iMargIndex==1) then
              offset_= 0
            else
              offset_= sum(vNEle(1:iMargIndex-1))
            end if
            !
            M%vDegree(1:MaxAtom)= 0
            !
            do K=1,len_trim(W) !processing the k1k2... string
              !
              J= CarToInt(W(K:K))
              if(J>0 .and. J<=S%NPole) then
                M%vDegree(J +offset_)= M%vDegree(J +offset_) +1
              else
                Ok= .false.
                MsgError=   "problem with coding string Wijk= "//trim(W)
                return !------------------------------------------return
              end if
              !
            end do
            !
            !--------------------------------retrieve values of WH, etc.
            call ReadRValsV(L,K,vX)
            M%WH=vX(1)  ;  M%WS=vX(2)  ;  M%WV=vX(3)  ;  M%WCp=vX(4)
            M%Kohler=  FLOOR(vX(5))
            !----------------------------------------------------------/
            !M%WG= M%WH &
            !&   + M%WCP*(T-T0) &
            !&   -(M%WS +M%WCP*log(T/T0))*T &
            !&   + M%WV *P
            !
            if(S%nMarg==MaxMarg) then
              Ok= .false.
              MsgError=                              "TOO MANY Margules"
              return !--------------------------------------------return
            end if
            S%nMarg= S%nMarg+1
            !
            if(S%Model == Mix_Site) then
              S%vIMargSite(S%nMarg)= iMargIndex
            else
              S%vIMargSite(S%nMarg)= 1
            end if
            S%vMarg(S%nMarg)= M
            !
            !----------------------------------------------------- trace 
            if(idebug>1) then
              write(fTrc,'(A,I2,A)',advance="NO") "nMarg=",S%nMarg," vDegree="
              write(fTrc,'(*(I1,1X))') (S%vMarg(S%nMarg)%vDegree(J),J=1,S%NAtom)
            end if
            !-----------------------------------------------------/trace 
            !
          end do doMargul
          
          end if
          !
          if(idebug>1) write(fTrc,'(A)') "<<-------/Read_MARGULES_block"
        !endcase MARGULES
        !-------------------------------------------/read MARGULES block
        end select
        !
      end do DoSol
      !
    end if !if(W=="MIXTURE.MODEL" .or. W=="SOLUTION.MODEL")
  end do DoFile
  !
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "<----------------/MixModel_BuildLnk"
  !
end subroutine MixModel_BuildLnk

subroutine MixModel_Clean(vHasPole,MM)
!-----------------------------------------------------------------------
!-- for the mixing model MM, remove the end-members that are absent
!-- in the current species list, 
!-- and update the list of margules parameters
!-- CAVEAT: this routine is currently NOT implemented for SITE models !!
!-----------------------------------------------------------------------
  use M_T_MixModel,only: T_MixModel,MaxPole,MaxAtom
  use M_T_MixModel,only: Mix_Molecular,Mix_Felspar
  !---------------------------------------------------------------------
  logical,         intent(in)   :: vHasPole(:)
  type(T_MixModel),intent(inout):: MM
  !---------------------------------------------------------------------
  type(T_MixModel):: MNew
  integer:: i,j,N
  integer:: vPermut(MaxPole),vDegree(MaxAtom)
  logical:: Keep
  !---------------------------------------------------------------------
  if(idebug>1) write(fTrc,'(A)') "<<---------------------MixModel_Clean"
  !
  MNew= MM
  select case(MM%Model)
  case(Mix_Molecular,Mix_Felspar) ! ("MOLECULAR","FELSPAR")
    N=0
    do i=1,MM%nPole
      if(vHasPole(i)) then
        N=N+1
        MNew%vIPole(N)= MM%vIPole(i)
        MNew%vNamPole(N)= trim(MM%vNamPole(i))
        vPermut(N)= i
        ! MNew%vHasPole(N)= .true.
      end if
    end do
    MNew%nPole= N
    MNew%nSite= 0
    N=0
    do i=1,MM%nMarg
      vDegree(:)= 0
      ! vDegree(:) contains all degrees,
      !   listed in the same order as
      !   the e-m list in case of molecular model,
      !   or as the site fractions in case of site model
      ! if a e-m or a site fraction is not involved, its degree is zero
      Keep= .true.
      do j=1,MM%nPole
        if(MM%vMarg(i)%vDegree(j)>0) then
          if(vHasPole(j)) then
            vDegree(j)= MM%vMarg(i)%vDegree(vPermut(j))
          else
            Keep= .false.
          end if
        end if
      end do
      if(Keep) then
        N=N+1
        MNew%vMarg(N)= MM%vMarg(i)
        MNew%vMarg(N)%vDegree(:)= vDegree(:)
        !----------------------------------------------------------trace
        if(idebug>1) then
          write(fTrc,'(A,I1,A)',advance="NO") 'nMarg=',N," vDegree="
          write(fTrc,'(*(I1,1X))') (MNew%vMarg(N)%vDegree(J),J=1,MNew%nPole)
        end if
        !---------------------------------------------------------/trace
      end if
      MNew%NMarg= N
    end do
  end select
  MM= MNew
  !
  if(idebug>1) write(fTrc,'(A)') "<<--------------------/MixModel_Clean"
  !
  return
end subroutine MixModel_Clean

subroutine MixModel_LnkToVec(Lnk,vMixModel)
  use M_T_MixModel,only: T_MixModel
  !
  type(T_LnkMixModel),          pointer    :: Lnk
  type(T_MixModel),dimension(:),intent(out):: vMixModel
  !
  type(T_LnkMixModel),pointer:: pCur, pPrev
  integer:: K
  !
  !
  if(idebug>1) write(fTrc,'(/,A)') "< MixModel_LnkToVec"
  !
  pCur=>Lnk
  K=0
  do while (associateD(pCur))
    K=K+1
    vMixModel(K)=pCur%Value
    !
    if(idebug>1) write(fTrc,'(A,I3,3A,I3)') &
    & "k=",k, &
    & "  MixModel(k)%Name=",trim(vMixModel(K)%Name), &
    & "  MixModel(k)%nPole=",vMixModel(K)%nPole
    !
    pPrev=>pCur
    pCur=> pCur%next
    deallocate(pPrev)
  end do
  !
  if(idebug>1) write(fTrc,'(A,/)') "< MixModel_LnkToVec"
  !
end subroutine MixModel_LnkToVec

subroutine MixModel_Check(vSpc,vMixModel)
  use M_Dtb_Const,only: Tref,Pref
  use M_T_Species, only: T_Species
  use M_T_MixModel,only: T_MixModel, T_Margul, Mix_Site
  use M_T_MixModel,only: MixModel_Param_Update, MixModel_Name
  !
  type(T_Species), dimension(:),intent(in):: vSpc
  type(T_MixModel),dimension(:),intent(in):: vMixModel
  !
  type(T_MixModel)   :: SM
  type(T_Margul)     :: Marg
  integer:: iEl,I,iPol,K
  logical:: Err
  character(15):: Str
  !
  if(idebug>1) write(fTrc,'(/,A)') "< MixModel_Check"
  !
  do I=1,size(vMixModel)

    SM=vMixModel(I)

    write(fTrc,'(/,A,I3,/)') "Mixture ",I
    write(fTrc,'(2A)')   "Name=    ",SM%Name
    
    call MixModel_Name(SM,Str,Err)
    if(.not. Err) then
      write(fTrc,'(2A)')   "Model=   ",trim(Str)
    else
      write(fTrc,'(2A)')   "Model=   ","UNKNOWN !!"
    end if

    if(SM%Model==Mix_Site) then

      !iEl=0
      !do iSit=1,SM%NSite
      !  write(fTrc,*) "SITE=",trim(SM%vNamSite(iSit))," ,vMulti=",SM%vMulti(iSit)
      !  write(fTrc,*) "Atoms="
      !  iEl=iEl+1
      !  write(fTrc,*) "iEl ", iEl, " ", SM%vNamAtom(iEl)
      !end do

      write(fTrc,'(/,A,/)') "***vNamAtom"
      do iEl=1,SM%NAtom
        write(fTrc,'(A)',advance="NO") SM%vNamAtom(iEl)
      end do

      write(fTrc,'(/,A,/)') "***endmembers"
      do iPol=1,SM%NPole

        write(fTrc,'(2A)') "POLE=",trim(SM%vNamPole(iPol))

      end do

    end if !! SITE models
    !
    !compute values of SM%vMarg(:)%WG at given P,T
    call MixModel_Param_Update(Tref,Pref,SM)
    !
    write(fTrc,'(/,A,/)') "list of poles and their index in vSpc"

    do iPol=1,SM%nPole
      if(SM%vIPole(iPol)>0) &
      & write(fTrc,'(2X,A,I3,I3,1X,A23)') &
      & "Pole",          &
      & iPol,            &
      & SM%vIPole(iPol), &
      & vSpc(SM%vIPole(iPol))%NamSp
    end do

    if(SM%nMarg>0) then
      write(fTrc,'(/,A,/)') "Margules"
      do K=1,SM%nMarg
        Marg=SM%vMarg(K)
        do iPol=1,SM%nPole
          write(fTrc,'(I3,A1)',advance="no") Marg%vDegree(iPol),T_
        end do
        write(fTrc,'(G15.3,A1,G15.3)') Marg%WH,T_,Marg%WG
      end do
      write(fTrc,*)
    end if

    write(fTrc,'(/,A,I3,/)') "End Mixture ",I
  end do
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ MixModel_Check"
  return
end subroutine MixModel_Check

end module M_MixModel_Read

