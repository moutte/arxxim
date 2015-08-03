MODULE M_SpeciesDtb_Read
!--
!-- routines for reading data for species from the databases
!--
  USE M_Kinds
  USE M_Trace,    ONLY: iDebug,fTrc,fHtm,Stop_,T_,Warning_,Pause_
  USE M_T_Species,ONLY: T_SpeciesDtb
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: SpeciesDtb_BuildLnk
  PUBLIC:: SpeciesDtb_LnkToVec
  PUBLIC:: LnkSpc_Build
  PUBLIC:: Species_Read_AquSize
  PUBLIC:: Species_Write_AquSize
  !
  PUBLIC:: T_LnkSpc
  !
  TYPE:: T_LnkSpc
    TYPE(T_SpeciesDtb)      :: Value
    TYPE(T_LnkSpc), POINTER :: Next
  ENDTYPE T_LnkSpc
  !
CONTAINS

SUBROUTINE Species_RemoveDoublons
  USE M_T_Species,  ONLY: T_Species,Species_Index
  USE M_Global_Vars,ONLY: vSpc
  !
  ! TYPE(T_Species),ALLOCATABLE:: vSpcNew(:)
  ! LOGICAL,ALLOCATABLE:: vSkip(:)
  INTEGER:: N,NNew,I

  N= SIZE(vSpc)
  NNew= N
  DO WHILE(N>0)
    I= Species_Index(vSpc(N)%NamSp,vSpc)
    IF(I<N) THEN
      PRINT *,"replace species ",trim(vSpc(N)%NamSp), trim(vSpc(I)%NamSp)
      CALL Pause_
      vSpc(I)= vSpc(N)
      NNew= NNew -1
    ENDIF
    N= N-1
  ENDDO

  IF(NNew<N) THEN
    PRINT *,"Old/NewDim",N,NNew
    CALL Pause_
  ENDIF

  RETURN
END SUBROUTINE Species_RemoveDoublons

SUBROUTINE LnkSpc_Build(Init,Inputt,L,P)
  LOGICAL,           INTENT(IN):: init
  TYPE(t_SpeciesDtb),INTENT(IN):: inputt
  TYPE(t_Lnkspc),    POINTER   :: l,p

  IF(init) THEN
    NULLIFY(l)
    ALLOCATE(l)
    NULLIFY(l%next)
    l%Value=inputt
    p=>l
  ELSE
    ALLOCATE(p%next)
    NULLIFY(p%next%next)
    p%next%Value=inputt
    p=>p%next
  ENDIF

  RETURN
ENDSUBROUTINE LnkSpc_Build

SUBROUTINE SpeciesDtb_BuildLnk(LnkSpc,N)
!--
!-- build vSpc directly from HSV files
!-- implements the S%Model and S%Indx
!-- which link the species S with the database models
!--
  USE M_Dtb_Vars, ONLY: vDtbMinHkf,vDtbMinThr,vDtbAquHkf,vDtbLogKtbl,vDtbLogKAnl
  !
  TYPE(T_LnkSpc),POINTER    :: LnkSpc
  INTEGER,       INTENT(OUT):: N
  !
  TYPE(T_LnkSpc),POINTER:: pCur !, pPrev
  TYPE(T_SpeciesDtb):: S
  INTEGER:: I
  !
  N= 0
  IF(SIZE(vDtbLogKtbl)>0) THEN
    DO I=1,SIZE(vDtbLogKtbl)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbLogKtbl
      S%DtbModel= "LOGKTBL"
      S%DtbTrace= TRIM(vDtbLogKtbl(I)%Num)
      CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    ENDDO
    ! stop here,
    ! because logK discrete databases can't be mixed
    ! with bases of other type
    RETURN !------------------------------------------------------RETURN
  ENDIF
  !
  N= 0
  IF(SIZE(vDtbLogKAnl)>0) THEN
    DO I=1,SIZE(vDtbLogKAnl)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbLogKAnl
      S%DtbModel= "LOGKANL"
      S%DtbTrace= TRIM(vDtbLogKAnl(I)%Num)
      CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    ENDDO
    ! stop here,
    ! because logK analytic databases can't be mixed
    ! with bases of other type
    RETURN !------------------------------------------------------RETURN
  ENDIF
  !
  N= 0
  !------------------------------------------------------aqueous species 
  !-------------------------------------place H2O as species Nr1 in list
  IF(SIZE(vDtbAquHkf)>0) THEN
    ! if there are aqueous species following the HKF equation of state,
    ! then the eos for H2O is necessarily the Haar-Gallagher-Kell model
    N=N+1
    S%Indx= 1
    S%DtbModel= "H2O_HGK" != Haar-Gallagher-Kell
    S%DtbTrace= "HAAR_etal"
    CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    !
    DO I=1,SIZE(vDtbAquHkf)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbAquHkf
      S%DtbModel= "AQU_HKF"
      S%DtbTrace= TRIM(vDtbAquHkf(I)%Num)
      CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    ENDDO
  ENDIF
  !
  !---------------------------------------then mineral (& gas ?) species
  IF(SIZE(vDtbMinHkf)>0) THEN
    DO I=1,SIZE(vDtbMinHkf)
      N=N+1
      S%Indx=     I !-> index of species S in vDtbMinHkf
      S%DtbModel= TRIM(vDtbMinHkf(I)%Typ)//"_HKF"
      S%DtbTrace= TRIM(vDtbMinHkf(I)%Num)
      CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    ENDDO
  ENDIF
  !
  !---------------------------------------idem, for Theriak-like formats
  IF(SIZE(vDtbMinThr)>0) THEN
    DO I=1,SIZE(vDtbMinThr)
      N=N+1
      S%Indx=     I
      S%DtbModel= TRIM(vDtbMinThr(I)%Typ)//"_THR"
      S%DtbTrace= TRIM(vDtbMinThr(I)%Num)
      CALL LnkSpc_Build(N==1,S,LnkSpc,pCur)
    ENDDO
  ENDIF

  RETURN
ENDSUBROUTINE SpeciesDtb_BuildLnk

SUBROUTINE SpeciesDtb_LnkToVec(Lnk,vSpc)
!--
!-- transfer the content of a linked list of species
!-- to an array of species,
!-- then deallocate the list
!--
  TYPE(T_LnkSpc),    POINTER    :: Lnk
  TYPE(T_SpeciesDtb),INTENT(OUT):: vSpc(:)
  !
  TYPE(T_LnkSpc),POINTER:: pCur, pPrev
  TYPE(T_SpeciesDtb):: S
  INTEGER:: N

  N= 0
  pCur=> Lnk
  DO WHILE (ASSOCIATED(pCur))
    S= pCur%Value
    N= N+1
    vSpc(N)= S
    pPrev=>  pCur
    pCur=>   pCur%next
    DEALLOCATE(pPrev)
  ENDDO

  RETURN
ENDSUBROUTINE SpeciesDtb_LnkToVec

SUBROUTINE Species_Read_AquSize(vSpc)
!--
!-- read size parameters of (some) aqueous species
!-- -> overwrite current value in vSpc
!--
  USE M_IOTools
  USE M_Files,    ONLY: NamFInn
  USE M_T_Species,ONLY: T_Species,Species_Index
  !
  TYPE(T_Species),INTENT(INOUT):: vSpc(:)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios,i
  !
  LOGICAL,ALLOCATABLE:: vIsPresent(:)

  IF(COUNT(vSpc(:)%Typ=="AQU")==0) RETURN

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Species_Read_AquSize"

  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))

  ALLOCATE(vIsPresent(SIZE(vSpc)))
  vIsPresent(:)= vSpc(:)%Z==0

  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    !
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)

    SELECT CASE(W)
    !
    CASE("SPECIES.SIZE")

      DoBlock: DO

        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        CALL AppendToEnd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines

        SELECT CASE(W)
          CASE("END","ENDSPECIES.SIZE")
            EXIT DoFile
        ENDSELECT

        i= Species_Index(W,vSpc)

        IF(i>0) THEN
          IF(vSpc(i)%Typ/="AQU") THEN
            CALL Warning_("Species "//TRIM(W)//" is not aqueous")
          ELSE
            CALL LinToWrd(L,W,EoL)
            CALL WrdToReal(W,vSpc(i)%AquSize)
            vIsPresent(i)= .TRUE.
            !print '(A,I3,1X,F7.2,1X,A)', &
            !& "AquSize ... ",i,vSpc(i)%AquSize,vSpc(I)%NamSp
            !pause_
          ENDIF
        ELSE
          CALL Warning_("Species "//TRIM(W)//" not in base")
        ENDIF

      ENDDO DoBlock

    ENDSELECT

  ENDDO DoFile
  CLOSE(F)

  IF(iDebug>0) THEN
    WRITE(fTrc,'(A,/)') "aqu'size parameters"
    DO i=1,SIZE(vSpc)
      IF(vSpc(i)%Typ=="AQU" .and. vSpc(i)%Z/=0) THEN
        IF(vSpc(I)%AquSize >1.E-6) THEN
          WRITE(fTrc,"(A,A1,F7.2)") &
          & vSpc(I)%NamSp, T_, vSpc(I)%AquSize
        ELSE
          WRITE(fTrc,"(A,A1,A)") &
          & vSpc(I)%NamSp, T_, " <- SIZE PARAMETER not defined"
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE(vIsPresent)

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Species_Read_AquSize"

  RETURN
ENDSUBROUTINE Species_Read_AquSize

SUBROUTINE Species_Write_AquSize(vSpc)
!--
!-- write size parameters of charged aqueous species
!--
  USE M_IOTools
  USE M_T_Species,ONLY: T_Species
  !
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  !
  INTEGER:: F,i
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Species_Write_AquSize"
  !
  CALL GetUnit(F)
  OPEN(F,FILE="aqusize.tab")
  !
  WRITE(F,'(A)') "SPECIES.SIZE"
  !
  DO i=1,SIZE(vSpc)
    IF(vSpc(i)%Typ=="AQU" .AND. vSpc(i)%Z/=0) &
    & WRITE(F,'(A,A1,F12.2)') vSpc(I)%NamSp, T_, vSpc(i)%AquSize
  ENDDO
  !
  WRITE(F,'(A)') "END"
  !
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Species_Write_AquSize"

  RETURN
ENDSUBROUTINE Species_Write_AquSize

SUBROUTINE Canevas
!--
!-- template for file reading
!--
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Canevas"
  !
  CALL GetUnit(F)
  OPEN(F,FILE="XXX")
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT")
      EXIT DoFile
    !
    CASE("BLOCK") !!!!!!canevas for READing one "block"
      !... I=0
      DoBlock: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        CALL AppendToEnd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
        !
        SELECT CASE(W)
        CASE("ENDINPUT")
          EXIT DoFile
        CASE("END","ENDBLOCK")
          EXIT DoBlock
          !-> in case reading data from several BLOCK..END blocks
          !CASE("END"); EXIT DoFile
          !!-> in CASE reading only the first available BLOCK
        CASE("XXX")
          !....
        CASE("YYY")
          !....
        ENDSELECT
        !...
        !
      ENDDO DoBlock
    !
    ENDSELECT
    !
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Canevas"

  RETURN
ENDSUBROUTINE Canevas

ENDMODULE M_SpeciesDtb_Read

!! INTEGER FUNCTION SpeciesDtb_Index(Str,V)
!! !.position of species with %NamSp==Str in vSpc_(1:nSpc_)
!!   CHARACTER(LEN=*),INTENT(IN):: Str
!!   TYPE(T_SpeciesDtb), INTENT(IN):: V(:)
!!   !
!!   INTEGER     ::I
!!   SpeciesDtb_Index=0
!!   IF(SIZE(V)==0) RETURN
!!   I=0
!!   DO
!!     I=I+1 !; IF(iDebug>0) WRITE(fTrc,'(A)') vCpn(I)%SpName
!!     IF(TRIM(Str)==TRIM(V(I)%NamSp)) THEN
!!       SpeciesDtb_Index=I
!!       EXIT
!!     ENDIF
!!     IF(I==SIZE(V)) EXIT
!!   ENDDO !IF Str not found -> I=0
!!   RETURN
!! ENDFUNCTION SpeciesDtb_Index
!! 
!! LOGICAL FUNCTION LnkSpc_Found(L,W,Spc)
!! !.find index of string W in linked list LnkSp
!! !.RETURN corresponding species Spc
!!   !
!!   !INTEGER::I_LnkSpc
!!   TYPE(T_LnkSpc),    POINTER    :: L
!!   CHARACTER(LEN=*),  INTENT(IN) :: W
!!   TYPE(T_SpeciesDtb),INTENT(OUT):: Spc
!!   TYPE(T_LnkSpc),    POINTER    :: P,pPrev
!!   INTEGER::I
!!   !
!!   P=>L; I=0
!!   LnkSpc_Found=.FALSE.
!!   DO WHILE (ASSOCIATED(P))
!!     I=I+1
!!     IF(TRIM(w)==TRIM(P%Value%NamSp)) THEN
!!       LnkSpc_Found=.TRUE.
!!       Spc=P%Value
!!       EXIT;
!!     ENDIF
!!     pPrev=>P
!!     P=> P%next
!!   ENDDO
!! ENDFUNCTION LnkSpc_Found

