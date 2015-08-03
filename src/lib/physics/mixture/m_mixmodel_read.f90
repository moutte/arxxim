MODULE M_MixModel_Read
  USE M_Kinds
  USE M_Trace,   ONLY: iDebug,fTrc,T_
  USE M_T_MixModel,ONLY: T_MixModel
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC:: MixModel_BuildLnk
  PUBLIC:: MixModel_LnkToVec
  PUBLIC:: MixModel_Check
  !
  PUBLIC:: T_LnkMixModel
  !
  TYPE:: T_LnkMixModel
    TYPE(T_MixModel)::Value
    TYPE(T_LnkMixModel),POINTER::Next
  ENDTYPE T_LnkMixModel
  !
CONTAINS

SUBROUTINE MixModel_BuildLnk( &
& vSpc, &
& Lnk,N,Ok,MsgError)
!--
!-- reading MIXTURE.MODEL blocks from file
!-- -> build a linked list of T_MixModel
!--
  USE M_IOTools
  USE M_T_MixModel
  USE M_Files,    ONLY: NamFSol
  USE M_T_Species,ONLY: T_Species,Species_Index
  !---------------------------------------------------------------------
  TYPE(T_Species),    INTENT(IN) :: vSpc(:)
  TYPE(T_LnkMixModel),POINTER    :: Lnk
  INTEGER,            INTENT(OUT):: N
  LOGICAL,            INTENT(OUT):: Ok
  CHARACTER(*),       INTENT(OUT):: MsgError
  !---------------------------------------------------------------------
  TYPE(T_LnkMixModel),POINTER:: pCur
  !
  CHARACTER(LEN=512):: L,W
  CHARACTER(LEN=3)  :: Str
  CHARACTER(LEN=80) :: SiteList
  !CHARACTER(LEN=255):: sFormul
  LOGICAL            :: sEol
  LOGICAL            :: MargulTheriak= .FALSE.
  INTEGER            :: ios,F,K,J,iEl,iS,iP,offset_,iMargIndex !,M1,M2
  TYPE(T_MixModel)   :: S
  TYPE(T_Margul)     :: M
  LOGICAL            :: vHasPole(MaxPole)
  REAL(dp),DIMENSION(dimV)::vX
  REAL(dp),ALLOCATABLE:: vXpol(:),vXatom(:)
  REAL(dp):: Y
  !
  !----------------------------------------------for reading SITE models
  ! REAL(dp)::NormCoeff
  !
  INTEGER,DIMENSION(0:MaxSite):: vNEle
  ! vNEle = nr of possible elements on each site, vNEle<=MaxMulti
  !       -> gives structure of input composition vector
  ! vNEle(0)=0,
  !       -> to enable i=i+vNEle(iSite-1)
  !          for localization of i in composition vector ...
  !
  CHARACTER(3*MaxMulti):: vSiteString(1:MaxSite)
  !--------------------------------------------------------------------/
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "<-----------------MixModel_BuildLnk"
  !
  Ok= .TRUE.
  MsgError= "Ok"
  !
  N=0
  NULLIFY(Lnk)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFSol))

  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,sEol)
    IF(W(1:1)=='!')   CYCLE DoFile
    CALL AppendToEnd(L,W,sEoL)
    !
    IF(W=="ENDINPUT") EXIT DoFile
    !
    !-----------------------------------------build a new solution model
    IF(W=="MIXTURE.MODEL" .OR. W=="SOLUTION.MODEL") THEN
    !
    ! keyword SOLUTION.MODEL should become obsolete for symmetric models
    ! MIXTURE.MODEL is the preferred keyword for symmetric models,
    ! SOLUTION.MODEL should be for assymetric solvent / solute mixtures
    !
      CALL MixModel_Zero(S)
      !
      CALL LinToWrd(L,W,sEol) !MIXTURE.MODEL OLIVINE ->S%Name="OLIVINE"
      S%Name=TRIM(W)
      !
      IF(iDebug>0) WRITE(fTrc,'(/,2A)') "MIXTURE=",TRIM(S%Name)
      IF(iDebug>2) PRINT '(/,A,/)',"Reading mixture model "//TRIM(S%Name)
      !
      DoSol: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        !
        CALL LinToWrd(L,W,sEol) !; IF(iDebug>2) WRITE(fTrc,'(A2,A)') "W=",TRIM(W)
        !
        IF(W(1:1)=='!') CYCLE DoSol !skip comment lines
        CALL AppendToEnd(L,W,sEoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT"); EXIT DoFile
        !
        CASE("END","ENDMIXTURE.MODEL","ENDSOLUTION.MODEL")
        !------------------------------------- end of one MIXTURE block,
        !------------------- -> save the T_MixModel value in linked list
          IF((S%nPole>1) .AND. &
          &  (COUNT(vHasPole(:) .EQV. .TRUE.)>1)) THEN
            IF(iDebug>2) PRINT *,"------------------------ACCEPTED"
            !
            IF(ANY(vHasPole(1:S%nPole) .EQV. .FALSE.)) &
            & CALL MixModel_Clean(vHasPole, S)
            !
            N=N+1
            IF(N==1) THEN !--fill the linked list
              ALLOCATE(Lnk); NULLIFY(Lnk%next)
              Lnk%Value=S; pCur=>Lnk
            ELSE
              ALLOCATE(pCur%next); NULLIFY(pCur%next%next)
              pCur%next%Value=S; pCur=>pCur%next
            END IF
          ELSE
            IF(iDebug>2) PRINT *,"------------------------REJECTED"
          END IF
          !
          EXIT DoSol
        !---------------------/ save the T_MixModel value in linked list
        !
        CASE("MODEL") !MODEL IDEAL / POLE / SITE / SPECIAL
        !-----------------------------------------------------read MODEL
          CALL LinToWrd(L,W,sEol)
          IF(TRIM(W)=="IDEAL") W="MOLECULAR"
          IF(TRIM(W)=="POLE") W="MOLECULAR"
          S%Model=TRIM(W)
          !
          IF(iDebug>0) WRITE(fTrc,'(2A)') "MODEL=",TRIM(W)
          !
          S%NSite=1
          IF(.NOT. sEol) THEN !=former version
            CALL LinToWrd(L,W,sEol)
            CALL WrdToInt(W,S%vMulti(1))
          ENDIF
          !
          CYCLE DoSol
          !not really necessary,
          !... unless the current value of W has changed to POLE or MARGULES ...
        !----------------------------------------------------/read MODEL
        !
        CASE("MULTI") ! needed only for MOLECULAR models with Multi/=1
        !----------------------------------------------read MULTIPLICITY
          CALL LinToWrd(L,W,sEol)
          CALL WrdToInt(W,K)
          IF(K<0 .or. K>MaxMulti) THEN
            Ok= .FALSE.
            MsgError=                            "error in MULTIPLICITY"
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          S%vMulti(1)= K
          !
          IF(iDebug>0) WRITE(fTrc,'(A,I3)') "MULTI=",S%vMulti(1)
          !
        !---------------------------------------------/read MULTIPLICITY
        !
        CASE("SITE")
        !------------------------------------------------read SITE block
          !
          ! example
          !   SITE
          !     !SiteName / Multiplicty / Elements
          !     T   2  SI_AL_
          !     M   1  MG_FE_AL_VAC
          !     VA4 4  MG_FE_AL_
          !   END SITE
          !
          S%nSite= 0
          vNEle=   0
          S%NAtom= 0
          iEl=     0
          SiteList= ""
          !
          IF(TRIM(S%Model)/="SITE") THEN
            Ok= .FALSE.
            MsgError=             "only SITE mixtures have a SITE block"
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          DoSite: DO
            !
            READ(F,'(A)') L; CALL LinToWrd(L,W,sEol)
            IF(W(1:1)=='!') CYCLE DoSite
            CALL AppendToEnd(L,W,sEoL)
            IF(W=="END") EXIT DoSite
            IF(W=="ENDSITE") EXIT DoSite
            !
            S%nSite=S%nSite+1
            !
            IF(S%nSite>MaxSite) THEN
              Ok= .FALSE.
              MsgError=                  "in SITE block: TOO MANY SITES"
              RETURN !--------------------------------------------RETURN
            ENDIF
            !
            !----------first word of line= name of site, up to 3 lettres
            CALL Str_Append(W,3)
            S%vNamSite(S%nSite)=TRIM(W)
            SiteList= TRIM(SiteList)//TRIM(W)
            !
            IF(iDebug>2) WRITE(fTrc,*) &
            & "NamSite,SiteList= ",S%vNamSite(S%nSite),TRIM(SiteList)
            !------------------------------------------------/first word
            !
            !------------second word of line= site multiplicity, integer
            CALL LinToWrd(L,W,sEol)
            CALL WrdToInt(W,S%vMulti(S%nSite))
            !
            IF(S%vMulti(S%nSite)>MaxMulti) THEN
              Ok= .FALSE.
              MsgError=TRIM(S%vNamSite(S%nSite))//"= Multiplicity too high"
              RETURN !--------------------------------------------RETURN
            ENDIF
            !-----------------------------------------------/second word
            !
            !---------------------- list of elements potentially present
            IF(TRIM(S%Model)=="SITE") THEN
              !
              IF(sEol) THEN
                Ok= .FALSE.
                MsgError= TRIM(S%vNamSite(S%nSite))//" NO ATOM LIST ???"
                RETURN !------------------------------------------RETURN
              ENDIF
              !
              CALL LinToWrd(L,W,sEol)
              IF(iDebug>2) WRITE(fTrc,*) TRIM(W)
              !
              vNEle(S%nSite)= LEN_TRIM(W)/3 ! how many elements on the site
              vSiteString(S%nSite)= TRIM(W)
              !
              IF(S%NAtom + vNEle(S%nSite) > MaxAtom) THEN
                Ok= .FALSE.
                MsgError=         TRIM(W)//": Too Many Elements in SITE"
                RETURN !------------------------------------------RETURN
              ENDIF
              !
              DO K=1,vNEle(S%nSite)
                !
                S%vIAtomSite(K+S%NAtom)= S%nSite
                S%vAtomMulti(K+S%NAtom)= S%vMulti(S%nSite)
                !
                !----name of the site fraction
                S%vNamAtom(K +S%NAtom)= TRIM(W(3*K-2:3*K))//S%vNamSite(S%nSite)
                !
              ENDDO
              !
              S%NAtom= S%NAtom + vNEle(S%nSite) ! total nr site fractions
              !
            ENDIF
            !----------------------/list of elements potentially present
            !
          ENDDO DoSite
          !check that SUM(S%vMulti(1:S%nSite)) <= MaxElPole
        !-----------------------------------------------/read SITE block
        !
        !----------------------------------------------- read POLE block
        CASE("POLE")
          !
          S%nPole=     0
          S%tPoleAtom= 0
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "<<------------Read_POLE_block"
          !
          DoPole: DO
            !
            READ(F,'(A)') L; CALL LinToWrd(L,W,sEol)
            IF(W(1:1)=='!') CYCLE DoPole
            CALL AppendToEnd(L,W,sEoL)
            !
            IF(TRIM(W)=="END") EXIT DoPole
            IF(TRIM(W)=="ENDPOLE") EXIT DoPole
            !
            K= Species_Index(W,vSpc)
            !
            !old! IF(K<1) THEN
            !old!   Ok= .FALSE.
            !old!   MsgError=        TRIM(W)//" not in base-> model not valid"
            !old!   RETURN !----------------------------------------RETURN
            !old! ENDIF
            !
            !IF(iDebug>0) &
            !& WRITE(fTrc,'(A,I3,2A)') "Species_Index=",K," <-POLE= ",TRIM(W)
            !
            IF(S%nPole<MaxPole) THEN
              S%nPole=S%nPole+1
            ELSE
              Ok= .FALSE.
              MsgError=                           "TOO MANY end-members"
              RETURN !--------------------------------------------RETURN
            ENDIF
            !
            vHasPole(S%nPole)= K>0
            !
            IF(vHasPole(S%nPole)) THEN
              !
              IF(S%nPole==1) THEN
                S%Typ= TRIM(vSpc(K)%Typ)
              ELSE
                IF(TRIM(vSpc(K)%Typ)/=TRIM(S%Typ)) THEN
                  Ok= .FALSE.
                  MsgError= "end-members should be same type, either GAS or MIN"
                  RETURN !----------------------------------------RETURN
                ENDIF
              ENDIF
              !
              S%vIPole(S%nPole)=K ! index of end-member name in species list
              !
            ELSE
              !
              S%vIPole(S%nPole)=0
              !
            ENDIF
            !
            S%vNamPole(S%nPole)=TRIM(W)
            !
            IF(iDebug>0) WRITE(fTrc,'(A5,I3,A1,I3,A1,A)') &
            & "POLE=",S%nPole,T_,K,T_,TRIM(W)
            !
            !--------------for SITE models, reading pole / site relation
            IF(TRIM(S%Model)=="SITE") THEN
              !
              ! scan the rest of the input line
              ! to retrieve the site repartition for the end-member
              !
              DO iS=1,S%NSite
                IF(sEol) THEN
                  Ok= .FALSE.
                  MsgError=TRIM(S%vNamPole(S%nPole))//" NO ATOM LIST ???"
                  RETURN !----------------------------------------RETURN
                ENDIF
                !
                CALL LinToWrd(L,W,sEol)
                !!!
                !!!M2 2 MG_FE_     -> vElSite(1),vElSite(2),            ( NEle(1)=2 )
                !!!M1 1 MG_FE_AL_  -> vElSite(3),vElSite(4),vElSite(5)  ( NEle(2)=3 )
                !!!T  4 AL_SI_     -> etc ...
                !!!A  1 K__
                !!!OH 2 OH_
                !!!
                !!!SIDEROPHYLLITE FE_FE_ AL_ AL_AL_SI_SI_ K__ OH_OH_ !->
                !!!
                IF(iDebug>2) WRITE(fTrc,'(2A)') "Site String= ",TRIM(W)
                !
                IF(LEN_TRIM(W)/=3*S%vMulti(iS)) THEN
                  Ok= .FALSE.
                  MsgError=TRIM(W)//" has length not consistent with vMulti"
                  RETURN !----------------------------------------RETURN
                ENDIF
                !
                IF(iS==1) THEN
                  offset_= 0
                ELSE
                  offset_= offset_ +vNEle(iS-1)
                ENDIF
                !
                DO K=1,S%vMulti(iS)
                  !
                  str=TRIM(W(3*K-2:3*K)) ! extract the element name
                  !
                  iEl= INDEX(vSiteString(iS),str)
                  IF(iEl<=0) THEN
                    Ok= .FALSE.
                    MsgError=       TRIM(str)//" <-Atom not Found   !!!"
                    RETURN !--------------------------------------RETURN
                  ENDIF
                  iEl= (iEl+2)/3
                  S%tPoleAtom(S%nPole,iEl+offset_)= S%tPoleAtom(S%nPole,iEl+offset_) +1
                  !
                ENDDO
                !
              ENDDO !!DO iS=1,S%NSite
              !
            ENDIF
            !-------------/for SITE models, reading pole / site relation
            !
          ENDDO DoPole
          !
          !----------------------------------------normalization factors
          IF(TRIM(S%Model)=="SITE") THEN
            ALLOCATE(vXpol(S%NPole))
            ALLOCATE(vXatom(S%NAtom))
            DO iP=1,S%NPole
              vXpol(:)=  Zero
              vXpol(iP)= One
              CALL MixModel_XPoleToXSite(S,vXpol(:),vXatom(:))
              Y= One
              DO iEl= 1,S%NAtom
                IF(S%tPoleAtom(iP,iEl)/=0) Y= Y *vXAtom(iEl)**S%vAtomMulti(iEl)
                !S%tPoleAtom(iP,iEl)
              ENDDO
              S%vPoleCoeff(iP)= Y
            ENDDO
            DEALLOCATE(vXpol)
            DEALLOCATE(vXatom)
            !
            !! PRINT *,"normalization factors:"
            DO iP=1,S%NPole
              IF(S%vPoleCoeff(iP)==Zero) THEN
                Ok= .FALSE.
                MsgError=        "Problem with normalization factor !!!"
                RETURN !------------------------------------------RETURN
              ENDIF
              S%vPoleCoeff(iP)= One /S%vPoleCoeff(iP)
            ENDDO
            !! PAUSE
          ENDIF
          !---------------------------------------/normalization factors
          !
          IF(TRIM(S%Model)/="SITE") S%NAtom= S%nPole
          !
          !--------------------------------------------------------trace
          IF(iDebug>2 .AND. TRIM(S%Model)=="SITE") THEN

            WRITE(fTrc,'(/,A)') "tPoleAtom"
            DO iS=1,S%NPole
              DO iEl=1,S%NAtom
                WRITE(fTrc,'(I7,1X)',ADVANCE="NO") S%tPoleAtom(iS,iEl)
              ENDDO
              WRITE(fTrc,*)
            ENDDO

            WRITE(fTrc,'(A)') "vAtomMulti"
            DO iEl=1,S%NAtom
              WRITE(fTrc,'(I7,1X)',ADVANCE="NO") S%vAtomMulti(iEl)
            ENDDO
            WRITE(fTrc,*)

            WRITE(fTrc,'(A)') "tPoleAtom"
            DO iS=1,S%NPole
              DO iEl=1,S%NAtom
                WRITE(fTrc,'(F7.2,1X)',ADVANCE="NO") &
                & S%tPoleAtom(iS,iEl) /REAL(S%vAtomMulti(iEl))
              ENDDO
              WRITE(fTrc,*)
            ENDDO

            WRITE(fTrc,'(/,A)') "Normalization:"
            DO iP=1,S%NPole
              WRITE(fTrc,'(A,1X,F7.2)') S%vNamPole(iP), S%vPoleCoeff(iP)
            ENDDO

          ENDIF
          !-------------------------------------------------------/trace
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "<<-----------/Read_POLE_block"
        !-----------------------------------------------/read POLE block
        !
        CASE("MARGULES")
        !--------------------------------------------read MARGULES block
          IF(TRIM(S%Model)=="SPECIAL") THEN
            Ok= .FALSE.
            MsgError=           "SPECIAL MODEL -> NO MARGULES BLOCK !!!"
            RETURN !----------------------------------------------RETURN
          END IF
          IF(iDebug>0) WRITE(fTrc,'(A)') "<<--------Read_MARGULES_block"
          !
          MargulTheriak= .FALSE.
          IF (.NOT. sEol) THEN
            CALL LinToWrd(L,W,sEol)
            IF (TRIM(W)=="THERIAK") MargulTheriak= .TRUE.
          END IF
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
          !  ENDMARGULES
          !
          ! <new format, a la Berman>
          ! <indices used in ijkl codes refer to order in the end-member list>
          !  MODEL MOLECULAR
          !  POLE
          !    GROSSULAR
          !    PYROPE
          !    ALMANDINE
          !  END
          !  MARGULES                          !                 vDegree(:)=
          !    112     21560.    18.79    .10  !GROSSULAR PYROPE       2,1,0
          !    122     69200.    18.79    .10  !GROSSULAR PYROPE       1,2,0
          !    113     20320.     5.08    .17  !GROSSULAR ALMANDINE    2,0,1
          !    133      2620.     5.08    .09  !GROSSULAR ALMANDINE    1,0,2
          !    ../..
          !  END MARGULES
          !
          ! <case of 'SITE-mixing' mixture, with Margules param's for each site>
          !  MODEL SITE
          !  SITE
          !    C 2 CA_MG_
          !  END SITE
          !  MARGULES SITE
          !    C 223     1307.40      0.00      0.01 !G_xs+= W223 *X_MG^2 *X_FE
          !    C 233     2092.40      0.00      0.06 !G_xs+= W233 *X_MG   *X_FE^2
          !    C 112    85529.00     18.79      0.21 !G_xs+= W112 *X_CA^2 *X_MG
          !    C 122    50874.90     18.79      0.02 !G_xs+= W122 *X_CA   *X_MG^2
          !    ../..
          !  END
          !
          IF(MargulTheriak) THEN
          
          doMargulTheriak: DO
          
            !TYPE tMarg !Margules parameter for up to Ternary mixing
            !  INTEGER::Dim !2 or 3
            !  INTEGER,DIMENSION(1:3)::iPole,Power
            !  REAL::WG,WH,WS,WV,WCp,WK
            !ENDTYPE tMarg
            !
            READ(F,'(A)') L ; CALL LinToWrd(L,W,sEol)
            !
            IF(W(1:1)=="!") CYCLE doMargulTheriak
            CALL AppendToEnd(L,W,sEoL)
            IF(W=="END" .OR. W=="ENDMARGULES") EXIT doMargulTheriak
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
            !  ENDMARGULES
            !
            ! wrk ! IF(W(1:1)/='&') THEN
            ! wrk ! ! scan the GROSSULAR PYROPE line
            ! wrk ! ! -> dimension of M, indices of Poles
            ! wrk !   !
            ! wrk !   M%NPole=0
            ! wrk !   DO
            ! wrk !     K=0
            ! wrk !     DO 
            ! wrk !       K=K+1 
            ! wrk !       IF(W==TRIM(S%vNamPole(K))) THEN
            ! wrk !         !find occurrence of "GROSSULAR" in end member list
            ! wrk !         M%NPole=M%NPole+1
            ! wrk !         M%iPole(M%NPole)=K
            ! wrk !       ENDIF
            ! wrk !       IF(K>S%nPole) EXIT
            ! wrk !     ENDDO
            ! wrk !     IF(sEol) EXIT
            ! wrk !     CALL LinToWrd(L,W,sEol)
            ! wrk !   ENDDO !end scan the Names
            ! wrk !   !
            ! wrk !   IF(M%Npole<2) CALL Stop_(TRIM(W)//" <-ENDMEMBERS NOT FOUND")
            ! wrk !   CYCLE doMargulTheriak
            ! wrk !   !
            ! wrk ! ENDIF
            ! wrk ! !
            ! wrk ! IF(W(1:1)=='&') THEN
            ! wrk ! !read Margules Parameters
            ! wrk !   !
            ! wrk !   S%nMarg=S%nMarg+1
            ! wrk !   CALL LinToWrd(L,W,sEol)
            ! wrk !   M%Power(1:M%NPole)=0 !initial values
            ! wrk !   DO K=1,LEN_TRIM(W) !processing the k1k2... string
            ! wrk !     J=CarToInt(W(K:K))
            ! wrk !     IF(J>0.AND.J<=M%NPole) M%Power(J)=M%Power(J)+1
            ! wrk !   ENDDO
            ! wrk !   !1222 1122 1112 1222 1223 1233
            ! wrk !   CALL ReadRValsV(L,K,vX)
            ! wrk !   M%WH= vX(1)  ;  M%WS=vX(2)  ;  M%WV=vX(3)  ;  M%WCp=vX(4)
            ! wrk !   M%Kohler=FLOOR(vX(5))
            ! wrk !   !M%WG= M%WH &
            ! wrk !   !&   + M%WCP*(T-T0) &
            ! wrk !   !&   -(M%WS +M%WCP*LOG(T/T0))*T &
            ! wrk !   !&   + M%WV *P
            ! wrk !   S%vMarg(S%nMarg)=M
            ! wrk !   !
            ! wrk ! ENDIF

          ENDDO doMargulTheriak
          
          ELSE
          
          doMargul: DO
            ! cf M_T_MixModel :
            !
            !<old tMargul type>
            ! TYPE tMargul !Margules parameter for up to quaternary mixing
            !   INTEGER::Dim !2 to 4
            !   INTEGER,DIMENSION(1:4)::iPole,Power
            !   REAL::WG,WH,WS,WV,WCp,WK
            ! ENDTYPE tMargul
            !
            !<new tMargul type>
            !  better for using systematic evaluation & derivation
            !  of monomials & polynomials
            !
            ! TYPE tMargul !Margules parameter
            !   INTEGER :: vDegree(1:MaxAtom)
            !   REAL(dp):: WG,WH,WS,WV,WCp,WK
            ! ENDTYPE tMargul
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
            READ(F,'(A)') L; CALL LinToWrd(L,W,sEol)
            !
            IF(W(1:1)=="!") CYCLE doMargul
            CALL AppendToEnd(L,W,sEoL)
            IF(W=="END" .OR. W=="ENDMARGULES") EXIT doMargul
            !
            iMargIndex= 1
            !
            !--------------------------------------------for SITE models
            !-----------------------find index of string W in S%vNamSite
            IF(TRIM(S%Model)=="SITE") THEN
              !
              IF(INDEX(W,"-")==2) THEN
                iMargIndex= CarToInt(W(1:1))
                W= W(3:LEN_TRIM(W))
                !! print *,"I, W=", iMargIndex,TRIM(W)  ; call pause_
              ELSE
                CALL Str_Append(W,3)
                J= INDEX(SiteList,TRIM(W))
                !J= 0
                !DO iS= 1,S%nSite
                !  IF(TRIM(W) == TRIM(S%vNamSite(iS))) THEN
                !    J= iS
                !    EXIT
                !  ENDIF
                !ENDDO
                IF(J==0) THEN
                  Ok= .FALSE.
                  MsgError=      TRIM(W)//" is not a site of this model"
                  RETURN !----------------------------------------RETURN
                ENDIF
                J= (J+2)/3
                !
                iMargIndex= J
                !
                CALL LinToWrd(L,W,sEol)
              ENDIF
            ENDIF
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
            IF(iMargIndex==1) THEN
              offset_= 0
            ELSE
              offset_= SUM(vNEle(1:iMargIndex-1))
            ENDIF
            !
            M%vDegree(1:MaxAtom)= 0
            !
            DO K=1,LEN_TRIM(W) !processing the k1k2... string
              !
              J= CarToInt(W(K:K))
              IF(J>0 .AND. J<=S%NPole) THEN
                M%vDegree(J +offset_)= M%vDegree(J +offset_) +1
              ELSE
                Ok= .FALSE.
                MsgError=   "problem with coding string Wijk= "//TRIM(W)
                RETURN !------------------------------------------RETURN
              ENDIF
              !
            ENDDO
            !
            !--------------------------------retrieve values of WH, etc.
            CALL ReadRValsV(L,K,vX)
            M%WH=vX(1)  ;  M%WS=vX(2)  ;  M%WV=vX(3)  ;  M%WCp=vX(4)
            M%Kohler=  FLOOR(vX(5))
            !----------------------------------------------------------/
            !M%WG= M%WH &
            !&   + M%WCP*(T-T0) &
            !&   -(M%WS +M%WCP*LOG(T/T0))*T &
            !&   + M%WV *P
            !
            IF(S%nMarg==MaxMarg) THEN
              Ok= .FALSE.
              MsgError=                              "TOO MANY Margules"
              RETURN !--------------------------------------------RETURN
            ENDIF
            S%nMarg= S%nMarg+1
            !
            IF(TRIM(S%Model)=="SITE") THEN
              S%vIMargSite(S%nMarg)= iMargIndex
            ELSE
              S%vIMargSite(S%nMarg)= 1
            ENDIF
            S%vMarg(S%nMarg)= M
            !
            !----------------------------------------------------- trace 
            IF(iDebug>0) THEN
              WRITE(fTrc,'(A,I2,A)',ADVANCE="NO") "nMarg=",S%nMarg," vDegree="
              WRITE(fTrc,'(*(I1,1X))') (S%vMarg(S%nMarg)%vDegree(J),J=1,S%NAtom)
            ENDIF
            !-----------------------------------------------------/trace 
            !
          END DO doMargul
          
          END IF
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "<<-------/Read_MARGULES_block"
        !ENDCASE MARGULES
        !-------------------------------------------/read MARGULES block
        END SELECT
        !
      END DO DoSol
      !
    ENDIF !IF(W=="MIXTURE.MODEL" .OR. W=="SOLUTION.MODEL")
  END DO DoFile
  !
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "<----------------/MixModel_BuildLnk"
  !
ENDSUBROUTINE MixModel_BuildLnk

SUBROUTINE MixModel_Clean(vHasPole,MM)
!-----------------------------------------------------------------------
!-- for the mixing model MM, remove the end-members that are absent
!-- in the current species list, 
!-- and update the list of margules parameters
!-- CAVEAT: this routine is currently NOT implemented for SITE models !!
!-----------------------------------------------------------------------
  USE M_T_MixModel,ONLY: T_MixModel,MaxPole,MaxAtom
  !---------------------------------------------------------------------
  LOGICAL,         INTENT(IN)   :: vHasPole(:)
  TYPE(T_MixModel),INTENT(INOUT):: MM
  !---------------------------------------------------------------------
  TYPE(T_MixModel):: MNew
  INTEGER:: i,j,N
  INTEGER:: vPermut(MaxPole),vDegree(MaxAtom)
  LOGICAL:: Keep
  !---------------------------------------------------------------------
  IF(iDebug>0) WRITE(fTrc,'(A)') "<<---------------------MixModel_Clean"
  !
  MNew= MM
  SELECT CASE(TRIM(MM%Model))
  CASE("MOLECULAR","FELSPAR")
    N=0
    DO i=1,MM%nPole
      IF(vHasPole(i)) THEN
        N=N+1
        MNew%vIPole(N)= MM%vIPole(i)
        MNew%vNamPole(N)= TRIM(MM%vNamPole(i))
        vPermut(N)= i
        ! MNew%vHasPole(N)= .TRUE.
      END IF
    END DO
    MNew%nPole= N
    MNew%nSite= 0
    N=0
    DO i=1,MM%nMarg
      vDegree(:)= 0
      ! vDegree(:) contains all degrees,
      !   listed in the same order as
      !   the e-m list in case of molecular model,
      !   or as the site fractions in case of site model
      ! if a e-m or a site fraction is not involved, its degree is zero
      Keep= .True.
      DO j=1,MM%nPole
        IF(MM%vMarg(i)%vDegree(j)>0) THEN
          IF(vHasPole(j)) THEN
            vDegree(j)= MM%vMarg(i)%vDegree(vPermut(j))
          ELSE
            Keep= .FALSE.
          END IF
        END IF
      END DO
      IF(Keep) THEN
        N=N+1
        MNew%vMarg(N)= MM%vMarg(i)
        MNew%vMarg(N)%vDegree(:)= vDegree(:)
        !----------------------------------------------------------trace
        IF(iDebug>0) THEN
          WRITE(fTrc,'(A,I1,A)',ADVANCE="NO") 'nMarg=',N," vDegree="
          WRITE(fTrc,'(*(I1,1X))') (MNew%vMarg(N)%vDegree(J),J=1,MNew%nPole)
        ENDIF
        !---------------------------------------------------------/trace
      END IF
      MNew%NMarg= N
    END DO
  END SELECT
  MM= MNew
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') "<<--------------------/MixModel_Clean"
  !
  RETURN
END SUBROUTINE MixModel_Clean

SUBROUTINE MixModel_LnkToVec(Lnk,vMixModel)
  USE M_T_MixModel,ONLY: T_MixModel
  !
  TYPE(T_LnkMixModel),          POINTER    :: Lnk
  TYPE(T_MixModel),DIMENSION(:),INTENT(OUT):: vMixModel
  !
  TYPE(T_LnkMixModel),POINTER:: pCur, pPrev
  INTEGER:: K
  !
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixModel_LnkToVec"
  !
  pCur=>Lnk
  K=0
  DO WHILE (ASSOCIATED(pCur))
    K=K+1
    vMixModel(K)=pCur%Value
    !
    IF(iDebug>0) WRITE(fTrc,'(A,I3,3A,I3)') &
    & "k=",k, &
    & "  MixModel(k)%Name=",TRIM(vMixModel(K)%Name), &
    & "  MixModel(k)%nPole=",vMixModel(K)%nPole
    !
    pPrev=>pCur
    pCur=> pCur%next
    DEALLOCATE(pPrev)
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "< MixModel_LnkToVec"
  !
ENDSUBROUTINE MixModel_LnkToVec

SUBROUTINE MixModel_Check(vSpc,vMixModel)
  USE M_Dtb_Const,ONLY: Tref,Pref
  USE M_T_Species, ONLY: T_Species
  USE M_T_MixModel,ONLY: T_MixModel, T_Margul, MixModel_Margul_Wg_Calc
  !
  TYPE(T_Species), DIMENSION(:),INTENT(IN):: vSpc
  TYPE(T_MixModel),DIMENSION(:),INTENT(IN):: vMixModel
  !
  TYPE(T_MixModel)   :: SM
  TYPE(T_Margul)     :: Marg
  INTEGER::iEl,I,iPol,K
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixModel_Check"
  !
  DO I=1,SIZE(vMixModel)

    SM=vMixModel(I)

    WRITE(fTrc,'(/,A,I3,/)') "Mixture ",I
    WRITE(fTrc,'(2A)')   "Name=    ",SM%Name
    WRITE(fTrc,'(2A)')   "Model=   ",SM%Model

    IF(TRIM(SM%Model)=="SITE") THEN

      !iEl=0
      !DO iSit=1,SM%NSite
      !  WRITE(fTrc,*) "SITE=",TRIM(SM%vNamSite(iSit))," ,vMulti=",SM%vMulti(iSit)
      !  WRITE(fTrc,*) "Atoms="
      !  iEl=iEl+1
      !  WRITE(fTrc,*) "iEl ", iEl, " ", SM%vNamAtom(iEl)
      !ENDDO

      WRITE(fTrc,'(/,A,/)') "***vNamAtom"
      DO iEl=1,SM%NAtom
        WRITE(fTrc,'(A)',ADVANCE="NO") SM%vNamAtom(iEl)
      ENDDO

      WRITE(fTrc,'(/,A,/)') "***endmembers"
      DO iPol=1,SM%NPole

        WRITE(fTrc,'(2A)') "POLE=",TRIM(SM%vNamPole(iPol))

      ENDDO

    ENDIF !! SITE models
    !
    !compute values of SM%vMarg(:)%WG at given P,T
    CALL MixModel_Margul_Wg_Calc(Tref,Pref,SM)
    !
    WRITE(fTrc,'(/,A,/)') "list of poles and their index in vSpc"

    DO iPol=1,SM%nPole
      IF(SM%vIPole(iPol)>0) &
      & WRITE(fTrc,'(2X,A,I3,I3,1X,A23)') &
      & "Pole",          &
      & iPol,            &
      & SM%vIPole(iPol), &
      & vSpc(SM%vIPole(iPol))%NamSp
    ENDDO

    IF(SM%nMarg>0) THEN
      WRITE(fTrc,'(/,A,/)') "Margules"
      DO K=1,SM%nMarg
        Marg=SM%vMarg(K)
        DO iPol=1,SM%nPole
          WRITE(fTrc,'(I3,A1)',ADVANCE='NO') Marg%vDegree(iPol),T_
        ENDDO
        WRITE(fTrc,'(G15.3,A1,G15.3)') Marg%WH,T_,Marg%WG
      ENDDO
      WRITE(fTrc,*)
    ENDIF

    WRITE(fTrc,'(/,A,I3,/)') "End Mixture ",I
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ MixModel_Check"
  RETURN
ENDSUBROUTINE MixModel_Check

ENDMODULE M_MixModel_Read

