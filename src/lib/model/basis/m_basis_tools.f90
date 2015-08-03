MODULE M_Basis_Tools
!--
!-- tools used by M_Basis for operations on basis (stoichiometry)
!--
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,fHtm,T_,Stop_,iDebug,Pause_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Basis_Dimensions
  PUBLIC:: Basis_Alloc
  PUBLIC:: Basis_Calc_iO_iH_iOx
  PUBLIC:: Basis_MixPhase_CpnUpdate
  PUBLIC:: Basis_SpcTypes_Build
  PUBLIC:: Basis_Cpn_Order
  PUBLIC:: Basis_SpcOrdPrm_Build
  PUBLIC:: Basis_AlfaNu_Build
  PUBLIC:: Basis_FindIBal
  PUBLIC:: Basis_Cpn_ISpc_Select
  PUBLIC:: Basis_Species_Select_Test
  PUBLIC:: Basis_FreeSet_Select
  PUBLIC:: Basis_StoikTable_Sho

  !~ PUBLIC:: Basis_NuTable_Change
  !~ PUBLIC:: Basis_Calc_iW_iH_iOH_iO2
  !~ PUBLIC:: Basis_RedoxAssign

CONTAINS

SUBROUTINE Basis_Alloc(nSp,nCp,nFs,nAq,nAs,nMs)
  USE M_Basis_Vars,ONLY: tStoikio
  USE M_Basis_Vars,ONLY: Basis_SpcTypes_Alloc,Basis_SpcOrdPrm_Alloc,Basis_AlfaNu_Alloc
  !
  INTEGER,INTENT(IN):: nSp,nCp,nFs,nAq,nAs,nMs
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Alloc"
  !
  !---allocate logical tables vLxx>
  CALL Basis_SpcTypes_Alloc(nSp)
  !
  !---allocate vOrdXx,vPrmXx>
  CALL Basis_SpcOrdPrm_Alloc(nCp,nSp,nAq,nAs,nMs)
  !
  !---allocate stoichio tables tAlfXx,tNuXx>
  CALL Basis_AlfaNu_Alloc(nCp,nSp,nFs,nAs,nMs)
  !
  IF(ALLOCATED(tStoikio)) DEALLOCATE(tStoikio)
  ALLOCATE(tStoikio(1:nCp,1:nSp))
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Alloc"
  !
END SUBROUTINE Basis_Alloc

SUBROUTINE Basis_Calc_iO_iH_iOx( & !
& vEle,vCpn, &    !in
& icO_,icH_,icOx) !out
!--
!-- find indexes of elements O, H, OX in vCpn (may change when basis change)
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  !
  TYPE(T_Element),  DIMENSION(:),INTENT(IN) :: vEle
  TYPE(T_Component),DIMENSION(:),INTENT(IN) :: vCpn
  INTEGER,                       INTENT(OUT):: icO_,icH_,icOx
  ! ranks of elements O__,H__,OX_ in Cpn of current component list
  INTEGER:: I
  !
  icO_=0; icH_=0; icOx=0
  DO I=1,SIZE(vCpn)
    SELECT CASE(vEle(vCpn(I)%iEle)%NamEl)
    CASE("O__"); icO_= I
    CASE("H__"); icH_= I
    CASE("OX_"); icOx= I
    END SELECT
  ENDDO
  IF(iDebug==5) THEN
    PRINT '(A,3I3)',"icO_icH_icOx= ",icO_,icH_,icOx
    CALL Pause_
  ENDIF
ENDSUBROUTINE Basis_Calc_iO_iH_iOx

SUBROUTINE Basis_Dimensions( & !
& BufferIsExtern,        & !IN
& vCpn,vSpc,nAq,nMn,nGs, & !IN
& nCi,nMx,nAx,nMs,nAs)     !OUT
!--
!-- count number of Inert / Mobile components
!-- assign species' statut
!--
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Phase,    ONLY: T_Phase
  !
  LOGICAL,          INTENT(IN) :: BufferIsExtern
  TYPE(T_Component),INTENT(IN) :: vCpn(:)
  TYPE(T_Species),  INTENT(IN) :: vSpc(:)
  INTEGER,          INTENT(IN) :: nAq,nMn,nGs
  INTEGER,          INTENT(OUT):: nCi,nMx,nAx,nMs,nAs
  !
  INTEGER:: nCp,nSp,iSp,iCp,N
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Dimensions"
  !
  N=nGs !currently not used ...
  nCp= SIZE(vCpn)
  nSp= SIZE(vSpc)
  nMx= 0
  nAx= 0
  !
  DO iCp=1,nCp
    !
    IF(vCpn(iCp)%Statut=="INERT") CYCLE 
    !--species not always assigned at this stage--
    !
    iSp=vCpn(iCp)%iSpc !-> index of Prim'species in vSpc(1:nSp)
    !
    SELECT CASE(vSpc(iSp)%Typ)

      CASE("AQU")
        !!!<BUFFER MODIFIED
        IF(BufferIsExtern) THEN
          IF(vCpn(iCp)%Statut=="BUFFER") &
          & CALL Stop_("Buffer species should be NON AQUEOUS !!!...")
        ELSE
          IF(vCpn(iCp)%Statut=="BUFFER") nAx=nAx+1
        ENDIF
        !!!</BUFFER MODIFIED
        IF(vCpn(iCp)%Statut=="MOBILE") nAx=nAx+1

      CASE("MIN"); nMx=nMx+1 !nr of Mobile Min.Species

      CASE("GAS"); nMx=nMx+1 !nr of Mobile Gas.Species (currently not different from Min ...)

    END SELECT
    !
  ENDDO
  !
  nCi= nCp -nAx -nMx
  nMs= nMn -nMx
  nAs= nAq -nCi -nAx
  !
  IF(iDebug>2) WRITE(fTrc,'(7(A,I4))') &
  & "nCp=",  nCp,"/ nCi=",nCi,"/ nAx=",nAx,"/ nAs=",nAs, &
  & "/ nMx=",nMx,"/ nMs=",nMs,"/ nMn=",nMn
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Dimensions"

  RETURN
ENDSUBROUTINE Basis_Dimensions

SUBROUTINE Basis_MixPhase_CpnUpdate( &
& vMixFas,   & !in
& vSpc,vCpn)   !inout
!--
!-- update activities of end-members in non-aqueous solutions
!-- update activities of mobile components
!--
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component
  USE M_T_MixPhase, ONLY: T_MixPhase
  !
  TYPE(T_MixPhase), INTENT(IN)   :: vMixFas(:)
  TYPE(T_Species),  INTENT(INOUT):: vSpc(:)
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  INTEGER :: iCp,iSp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "<== Basis_MixPhase_CpnUpdate"
  !
  !------------------update activities of externally buffered components
  DO iCp=1,SIZE(vCpn)
    iSp= vCpn(iCp)%iSpc
    !
    IF(vCpn(iCp)%iMix/=0) THEN
    !--- IF prim'species is in non'aqu'sol', THEN component is mobile --
      vSpc(iSp)%Dat%LAct= vMixFas(vCpn(iCp)%iMix)%vLnAct(vCpn(iCp)%iPol)
      vCpn(iCp)%LnAct= vSpc(iSp)%Dat%LAct
    ENDIF
    !
    IF(iDebug>0) WRITE(fTrc,'(A,I3,A1,A,A1,G15.6)') &
    & "CPN=",iCp,T_,vSpc(iSp)%NamSp,T_,vSpc(iSp)%Dat%LAct
  ENDDO
  !-----------------/
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "<==/ Basis_MixPhase_CpnUpdate"
  !
  RETURN
ENDSUBROUTINE Basis_MixPhase_CpnUpdate

SUBROUTINE Basis_SpcTypes_Build( & !
& vSpc,vCpn, & ! IN
& vLCi,vLAx,vLMx,vLAs,vLMs) ! OUT
!--
!-- build the logical arrays vLCi,vLAx,vLMx,vLAs,vLMs,...
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component
  !
  TYPE(T_Species),     INTENT(IN) :: vSpc(:)
  TYPE(T_Component),   INTENT(IN) :: vCpn(:)
  LOGICAL,DIMENSION(:),INTENT(OUT):: vLCi,vLAx,vLMx,vLAs,vLMs
  !
  INTEGER:: iCp,iSp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_SpcTypes_Build"
  !
  !-------------------------------------------update logical tables vLxx
  vLCi=.FALSE.
  vLAx=.FALSE.
  vLMx=.FALSE.
  vLAs(:)= vSpc(:)%Typ(1:3)=="AQU"
  vLMs(:)= vSpc(:)%Typ(1:3)=="MIN" .OR. vSpc(:)%Typ(1:3)=="GAS"
  !vLGs(:)= vSpc(:)%Typ=="GAS"
  !
  DO iCp=1,SIZE(vCpn)
    !
    iSp= vCpn(iCp)%iSpc !-> iSp is primary species
    !
    vLAs(iSp)= .FALSE.
    vLMs(iSp)= .FALSE.
    !
    SELECT CASE(TRIM(vCpn(iCp)%Statut))
    !
    CASE("INERT","BALANCE") !Intern'Species,
      vLCi(iSp)=.TRUE.
      !includes "balance" species for charge adjustment (in case pH fixed)
    !
    CASE("MOBILE","BUFFER")
      SELECT CASE(vSpc(iSp)%Typ(1:3))
        CASE("AQU")       ; vLAx(iSp)=.TRUE. !eXtern'Aqu.Species
        CASE("MIN","GAS") ; vLMx(iSp)=.TRUE. !eXtern'Min.Species
      END SELECT
    !
    END SELECT
  ENDDO
  !-----------------------------------------------/update logical tables
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_SpcTypes_Build"
  RETURN
ENDSUBROUTINE Basis_SpcTypes_Build

SUBROUTINE Basis_Cpn_Order( & !
& vSpc,nCi,nAx, & !in
& vCpn)           !out
!--
!-- sort the components, inert first, mobile next --
!--
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component
  !~ USE M_Basis_Vars, ONLY: iEleCpn,iCpnEle
  !~ USE M_Basis_Vars, ONLY: Basis_CpnPrm_Alloc
  !
  TYPE(T_Species),  INTENT(INOUT):: vSpc(:) !vSpc(iSp)%Dat%LAct=vCpn(iCp)%LnAct
  INTEGER,          INTENT(IN)   :: nCi,nAx
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  TYPE(T_Component),DIMENSION(:),ALLOCATABLE::vC !-> re-ordered vCpn
  INTEGER:: iCp,iCi,iAx,iMx,iSp !,nSp,nCp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Cpn_Order"
  !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(A)') "Components->Species"
    DO iCp=1,SIZE(vCpn)
      WRITE(fTrc,'(A,I3,A1,A)') &
      & "Cpn",iCp,T_,TRIM(vSpc(vCpn(iCp)%iSpc)%NamSp)
    ENDDO
  ENDIF
  !
  iCi=0 !index for Internal Prim'Aqu'Species
  iAx=0 !index for mobile Prim'Aqu'species, begins after nCi in vCpn
  iMx=0 !index of mobile Min.species, begins after Ci+nAx in vEle
  !
  !--------------------------reorder components, Inert / MobAqu / MobMin
  ALLOCATE(vC(1:SIZE(vCpn)))
  !
  DO iCp=1,SIZE(vCpn)
    !
    iSp=vCpn(iCp)%iSpc !-> iSp is primary species
    !
    SELECT CASE(TRIM(vCpn(iCp)%Statut))

    CASE("INERT","BALANCE") !Intern'Species
      !--------------------------------------- inert components first --
      iCi=iCi+1
      vC(iCi)=vCpn(iCp)
      !
    CASE("MOBILE","BUFFER")
      vSpc(iSp)%Dat%LAct= vCpn(iCp)%LnAct
      SELECT CASE(vSpc(iSp)%Typ(1:3))

      CASE("AQU") !eXtern'Aqu.Species
        !---------------------------------------- then mobile aqueous --
        iAx=iAx+1
        vC(nCi+iAx)=vCpn(iCp)

      CASE("MIN","GAS") !eXtern'Min.Species
        !----------------------------------------- mobile non'aqueous --
        iMx=iMx+1
        vC(nCi+nAx+iMx)=vCpn(iCp)

      END SELECT

    END SELECT
    !
  ENDDO
  !
  vCpn(:)= vC(:)
  !
  DEALLOCATE(vC)
  !-------------------------/reorder components, Inert / MobAqu / MobMin
  !
  !-------------------------------------------- build iEleCpn,iCpnEle --
  !~ CALL Basis_CpnPrm_Alloc(SIZE(vCpn))
  !~ CALL Basis_CpnPrm_Calc(vCpn,iEleCpn,iCpnEle)
  !-------------------------------------------/ build iEleCpn,iCpnEle --
  !
  !----------------------------------------------------------------trace
  IF(iDebug>0) THEN
    WRITE(fTrc,'(/,A,/)') "Components, Re-Ordered"
    WRITE(fTrc,'(2(A,/))') &
    & "iCp,Cpn%Statut,vSpc(Cpn%iSpc)%NamSp,vCpn(iCp)%Mole,Cpn%LnAct/Ln10"
    !
    DO iCp=1,SIZE(vCpn)
      WRITE(fTrc,'(I3,A1,2(A,A1),2(G15.6,A1))') &
      & iCp,                        T_,&
      & vCpn(iCp)%Statut,           T_,&
      & vSpc(vCpn(iCp)%iSpc)%NamSp, T_,&
      & vCpn(iCp)%Mole,             T_,&
      & vCpn(iCp)%LnAct/Ln10,       T_
    ENDDO
    !CALL Pause_
  ENDIF
  !---------------------------------------------------------------/trace
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Cpn_Order"
  !
ENDSUBROUTINE Basis_Cpn_Order

!~ SUBROUTINE Basis_CpnPrm_Calc( & !
!~ & vCpn,             & ! IN
!~ & iEleCpn,iCpnEle)    ! OUT
  !~ USE M_T_Element,  ONLY: T_Element
  !~ USE M_T_Component,ONLY: T_Component
  !~ !
  !~ TYPE(T_Component),INTENT(IN) :: vCpn(:)
  !~ INTEGER,          INTENT(OUT):: iEleCpn(:),iCpnEle(:)
  !~ !
  !~ INTEGER:: I,J
  !~ !
  !~ iCpnEle(:)=vCpn(:)%iEle
  !~ !-> iEl=iCpnEle(iCp)= which Ele is associated with vCpn(iCp)
  !~ !
  !~ DO I=1,SIZE(vCpn)
    !~ DO J=1,SIZE(vCpn)
      !~ IF(J==iCpnEle(I)) iEleCpn(I)=J
    !~ ENDDO
    !~ !-> vCpn(iEleCpn(iEl)) is associated with element iEl
  !~ ENDDO
  !~ !
!~ ENDSUBROUTINE Basis_CpnPrm_Calc

SUBROUTINE Basis_SpcOrdPrm_Build( &
& vSpc,vCpn, &
& vLCi,vLAx,vLMx,vLAs,vLMs, &
& vOrdPr,vOrdAq,vOrdAs,vOrdMs, &
& vPrmFw,vPrmBk)
!--
!-- build vOrdPr, vOrdAq, vOrdAs, vOrdMs, vPrmFw, vPrmBk,
!-- variables used- vLAs,vLMs,.., nCi,nAx,nMx,nAs,nMs, ...
!-- vPrmFw(J) is index, in vSpc, of the J-th species in the current array
!-- vPrmBk(J) is index in the current array of the J-th species of vSpc
!--
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component
  !
  TYPE(T_Species),  DIMENSION(:),INTENT(IN) :: vSpc
  TYPE(T_Component),DIMENSION(:),INTENT(IN) :: vCpn
  LOGICAL,          DIMENSION(:),INTENT(IN) :: vLCi,vLAx,vLMx
  LOGICAL,          DIMENSION(:),INTENT(IN) :: vLAs,vLMs
  INTEGER,          DIMENSION(:),INTENT(OUT):: vOrdPr,vOrdAq,vOrdAs,vOrdMs
  INTEGER,          DIMENSION(:),INTENT(OUT):: vPrmFw,vPrmBk
  !
  INTEGER:: I,J,nCp,nSp
  INTEGER:: nAq,nCi,nAs,nAx,nMx,nMs
  ! INTEGER:: iCi,iAx,iMx
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_SpcOrdPrm_Build"
  !
  nSp= SIZE(vSpc)
  nCp= SIZE(vCpn)
  nAq= count(vSpc%Typ(1:3)=="AQU")
  !
  nCi= count(vLCi)
  nAs= count(vLAs)
  nAx= count(vLAx)
  nMx= count(vLMx)
  nMs= count(vLMs)
  !
  !--- Aqueous Inert Prim'Species --
  vOrdPr(1:nCi)= vCpn(1:nCi)%iSpc
  !--- Aqueous Mobile Species --
  IF(nAx>0) vOrdPr(nCi+1:nCi+nAx)= vCpn(nCi+1:nCi+nAx)%iSpc
  !--- non Aqu'Mobile Species --
  IF(nMx>0) vOrdPr(nCi+nAx+1:nCi+nAx+nMx)= vCpn(nCi+nAx+1:nCi+nAx+nMx)%iSpc
  !
  !--- Aqueous Inert Prim'Species --
  vOrdAq(1:nCi)= vCpn(1:nCi)%iSpc
  !--- Aqueous Mobile Species --
  IF(nAx>0) &
  & vOrdAq(nCi+nAs+1:nCi+nAs+nAx)= vCpn(nCi+1:nCi+nAx)%iSpc
  !
  !------------------ scan vSpc for aqu. second. species -> build vOrdAs 
  IF(nAs>0) THEN !
    J=0
    DO I=1,nSp
      IF(vLAs(I)) THEN !vSpc(I)%Typ=="AQU"
        J=J+1
        vOrdAs(J)=I
        !-> vOrdAs(J, J=1:nAs) points to index of Jth sec'aqu'species in vSpc
        vOrdAq(nCi+J)=I
      ENDIF
    ENDDO
  ENDIF
  !------------------/
  !
  !-------------------scan vSpc for min. second. species -> build vOrdMs
  IF(nMs>0) THEN
    J=0
    DO I=1,nSp
      IF(vLMs(I)) THEN !vSpc(I)%Typ=="MIN" .OR. vSpc(:)%Typ=="GAS"
        J=J+1
        vOrdMs(J)=I
        !-> vOrdMs(J, J=1:nMs) points to index of Jth sec'min'species in vSpc
      ENDIF
    ENDDO
  ENDIF
  !-------------------/
  !
  !nCi          //nAs           //nAx           //nMx       //nMs
  !Aqu.Inert.Spc//Aqu.Second.Spc//Aqu.Mobile.Spc//Min.Mobile//Min.notMobile
  !1            //nCi+1         //nCi+nAs+1     //nAq+1     //nAq+nMx+1
  vPrmFw(1:nCi)=                            vOrdPr(1:nCi)
  IF(nAx>0) vPrmFw(nCi+nAs+1:nCi+nAs+nAx)=  vOrdPr(nCi+1:nCi+nAx)
  IF(nMx>0) vPrmFw(nAq+1:nAq+nMx)=          vOrdPr(nCi+nAx+1:nCi+nAx+nMx)
  IF(nAs>0) vPrmFw(nCi+1:nCi+nAs)=          vOrdAs(1:nAs)
  IF(nMs>0) vPrmFw(nAq+nMx+1:nAq+nMx+nMs)=  vOrdMs(1:nMs)
  !
  DO I=1,nSp
    vPrmBk(vPrmFw(I))=I
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_SpcOrdPrm_Build"
  !
ENDSUBROUTINE Basis_SpcOrdPrm_Build

SUBROUTINE AlfaTable_Build( &
& nCp,nSp,tStoikio, &
& vIndexBuffer, &
& vTot, &
& tAlfSp)
!--
!-- build tAlfSp from tStoikio(1:nCp,1:nSp)
!--
  USE M_Basis_Ortho
  !
  INTEGER,  INTENT(IN) :: nCp,nSp
  REAL(dp), INTENT(IN) :: tStoikio(:,:)
  INTEGER,  INTENT(IN) :: vIndexBuffer(:)
  !LOGICAL,  INTENT(IN) :: HasBuffer
  REAL(dp), INTENT(INOUT):: vTot(:)
  REAL(dp), INTENT(OUT)  :: tAlfSp(:,:)
  !
  REAL(dp),ALLOCATABLE:: tXOrthoBasis(:,:)
  != indexes of buffer components in tStoikio
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< AlfaTable_Build"
  !
  !tStoikio:column iSp = stoikio of species iSp in terms of elements (1:nCp)
  !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !
  tAlfSp(1:nCp,1:nSp)= tStoikio(1:nCp,1:nSp)
  !
  IF(ANY(vIndexBuffer(:)>0)) THEN
    !
    ALLOCATE(tXOrthoBasis(nCp,nCp))
    !
    IF(iDebug>2) THEN
      CALL Basis_Ortho_Compute_Debug(vIndexBuffer,tAlfSp,tXOrthoBasis)
    ELSE
      CALL Basis_Ortho_Compute(vIndexBuffer,tAlfSp,tXOrthoBasis)
    ENDIF
    !
    tAlfSp= MATMUL(tXOrthoBasis,tAlfSp)
    vTot=   MATMUL(tXOrthoBasis,vTot)
    !
    DEALLOCATE(tXOrthoBasis)
    !
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ AlfaTable_Build"
  !
END SUBROUTINE AlfaTable_Build

SUBROUTINE AlfaTable_Extract( &
& vCpn,tAlfSp,vLCi,vLAx,vLMx,vLAs,vLMs, &
& tAlfPr,tAlfAs,tAlfMs)
!-----------------------------------------------------------------------
!.extract blocks from tAlfSp(1:nCp,1:nSp) to form tAlfPr,tAlfAs,tAlfMs
!.according to vLCi,vLAx,vLMx,vLAs,vLMs
!-----------------------------------------------------------------------
  USE M_T_Component,ONLY: T_Component
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  REAL(dp), INTENT(IN):: tAlfSp(:,:)
  LOGICAL,  INTENT(IN):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !
  REAL(dp), INTENT(OUT):: tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !
  INTEGER :: I,J
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< AlfaTable_Extract"
  !
  !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !
  !--extract blocks from tStoikio(1:nCp,1:nSp) to tAlfPr,tAlfAs,tAlfMs--
  J= 0
  !-------------------------------------subblock of inert prim'species--
  DO I=1,COUNT(vLCi)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  ENDDO
  !-------------------------------------subblock of mobile aqu'species--
  DO I=1,COUNT(vLAx)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  ENDDO
  !-------------------------------------subblock of mobile min'species--
  DO I=1,COUNT(vLMx)
    J= J+1
    tAlfPr(:,J)=tAlfSp(:,vCpn(J)%iSpc)
  ENDDO
  !
  tAlfAs(:,:)= Zero
  !---------------------------------extract columns of sec'aqu.species--
  IF(COUNT(vLAs)>0) CALL Extract(2,tAlfSp,vLAs,tAlfAs)
  !
  tAlfMs(:,:)= Zero
  !-----------------------------extract columns of sec'min.gas.species--
  IF(COUNT(vLMs)>0) CALL Extract(2,tAlfSp,vLMs,tAlfMs)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ AlfaTable_Extract"
  !
END SUBROUTINE AlfaTable_Extract

SUBROUTINE AlfaTable_CalcBuffered( & !
& vCpn,vSpc,tStoikio,vLCi,vLAs,vLMs, & !
& vTot, &
& tAlfPr,tAlfAs,tAlfMs) !!! NOT USED ???
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  USE M_Numeric_Mat,ONLY: LU_BakSub
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  TYPE(T_Species),  INTENT(IN):: vSpc(:)
  REAL(dp),         INTENT(IN):: tStoikio(:,:)
  LOGICAL,          INTENT(IN):: vLCi(:),vLAs(:),vLMs(:)
  REAL(dp),         INTENT(INOUT):: vTot(:)
  REAL(dp),         INTENT(INOUT):: tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !
  REAL(dp),ALLOCATABLE:: tTransform(:,:)
  INTEGER, ALLOCATABLE:: vIndx(:)
  REAL(dp),ALLOCATABLE:: vY(:)
  LOGICAL :: bSingul
  INTEGER :: I,nCp,nSp,nCi
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< AlfaTable_CalcBuffered"
  !
  nCp= SIZE(vCpn)
  nSp= SIZE(vSpc)
  nCi= count(vLCi)
  !
  !for mole balance operations, replace the elements by "standard" primary species
  != "INERT" component: single ion with default valency, e.g. Na->Na+,... -> no modIF'
  != "mobile" component: use species as component
  !  -> elements replaced by corresponding species, e.g. C by CO2
  !
  ALLOCATE(tTransform(nCp,nCp))
  ALLOCATE(vIndx(nCp))
  ALLOCATE(vY(nCp))
  !
  tTransform(1:nCp,1:nCp)=Zero
  DO I=1,nCp; tTransform(I,I)=One; ENDDO !identity matrix
  !
  IF(nCp>nCi) THEN
  ! for buffered species, replace element by species
  ! using corresponding column of formula matrix
    DO I=nCi+1,nCp
      tTransform(:,I)= tStoikio(:,vCpn(I)%iSpc)
    ENDDO
  ENDIF
  !
  CALL Compute_Transform(tTransform,vIndx,bSingul)
  IF(bSingul) CALL Stop_("SINGUL IN CalcAlfaTableBuffered")
  !-> tTransform gives Elements as functions of Int'Elements/Ext'Species
  !
  !-----------------------compute new total amounts of component species
  CALL LU_BakSub(tTransform,vIndx,vTot(1:nCp))
  !----------------------/
  !DO I=1,nSp
  !  vY=tAlfSp(1:nCp,I)
  !  CALL LU_BakSub(tTransform,vIndx,vY)
  !  tAlfSp(1:nCp,I)=vY(1:nCp)
  !ENDDO
  !
  DO I=1,nCp
    vY=tAlfPr(1:nCp,I)
    CALL LU_BakSub(tTransform,vIndx,vY)
    tAlfPr(1:nCp,I)=vY(1:nCp)
  ENDDO
  !
  DO I=1,COUNT(vLAs)
    vY=tAlfAs(1:nCp,I)
    CALL LU_BakSub(tTransform,vIndx,vY)
    tAlfAs(1:nCp,I)=vY(1:nCp)
  ENDDO
  !
  DO I=1,COUNT(vLMs)
    vY=tAlfMs(1:nCp,I)
    CALL LU_BakSub(tTransform,vIndx,vY)
    tAlfMs(1:nCp,I)=vY(1:nCp)
  ENDDO
  !
  DEALLOCATE(tTransform)
  DEALLOCATE(vIndx)
  DEALLOCATE(vY)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ AlfaTable_CalcBuffered"
ENDSUBROUTINE AlfaTable_CalcBuffered

SUBROUTINE NuTable_Build( & !
& vIndxPrimSpc,tStoikio,  & !IN
& tNuSp)                    !OUT
  USE M_Numeric_Mat,ONLY: LU_BakSub
  !
  INTEGER, INTENT(IN) :: vIndxPrimSpc(:)
  REAL(dp),INTENT(IN) :: tStoikio(:,:)
  REAL(dp),INTENT(OUT):: tNuSp(:,:)
  !
  REAL(dp),ALLOCATABLE:: tTransform(:,:)
  INTEGER, ALLOCATABLE:: vIndx(:)
  REAL(dp),ALLOCATABLE:: vY(:)
  LOGICAL :: bSingul
  INTEGER :: I,J,nCp
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< NuTable_Build"
  !
  nCp= SIZE(tStoikio,1)
  !
  ALLOCATE(tTransform(nCp,nCp))
  !
  !-------------------assemble tTransform from the columns of tStoikio--
  !-----------------------------------------that point to Prim'Species--
  J=0
  DO I=1,nCp
    J=J+1
    tTransform(:,J)= tStoikio(:,vIndxPrimSpc(I))
  ENDDO
  !
  !--------------------------> tTransform is tStoikio(nCi | nAx | nMx)--
  !! iCi= 0
  !! iAx= COUNT(vLCi)
  !! iMx= COUNT(vLCi) +COUNT(vLAx)
  !! DO I=1,SIZE(tNuSp,1)
  !!   IF(vLCi(I)) THEN
  !!     iCi= iCi+1
  !!     tTransform(:,iCi)= tStoikio(:,vCpn(iCi)%iSpc)
  !!   ELSEIF(vLAx(I)) THEN
  !!     iAx= iAx+1
  !!     tTransform(:,iAx)= tStoikio(:,vCpn(iAx)%iSpc)
  !!   ELSEIF(vLMx(I)) THEN
  !!     iMx= iMx+1
  !!     tTransform(:,iMx)= tStoikio(:,vCpn(iMx)%iSpc)
  !!   ENDIF
  !! ENDDO
  !
  ALLOCATE(vIndx(nCp))
  CALL Compute_Transform(tTransform,vIndx,bSingul)
  !
  IF(bSingul) CALL Stop_("SINGUL IN NuTable_Build")
  !
  ALLOCATE(vY(nCp))
  !
  !-- in tNuSp, species are in same order as in vSpc --
  DO I=1,SIZE(tNuSp,1)
    vY= tStoikio(1:nCp,I)
    CALL LU_BakSub(tTransform,vIndx,vY)
    !Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    tNuSp(I,1:nCp)= vY(1:nCp)
    !-> this includes transposition between tAlfa and tNu !!
  ENDDO
  !
  DEALLOCATE(tTransform)
  DEALLOCATE(vIndx)
  DEALLOCATE(vY)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ NuTable_Build"
  !
ENDSUBROUTINE NuTable_Build
  !
SUBROUTINE NuTable_Extract(  & !
& tNuSp,vLAs,vLMs, & !
& tNuAs,tNuMs)
  !
  REAL(dp),INTENT(IN) :: tNuSp(:,:)
  LOGICAL, INTENT(IN) :: vLAs(:),vLMs(:)
  REAL(dp),INTENT(OUT):: tNuAs(:,:),tNuMs(:,:)
  !
  !----extract blocks from tNuSp, for the different types of sec'species
  !----------------------------------------------- second'aqu'species --
  IF(COUNT(vLAs)>0) CALL Extract( & !
  & 1,tNuSp,vLAs, & !IN
  & tNuAs)          !OUT
  !----------------------------------------------- second'min'species --
  IF(COUNT(vLMs)>0) CALL Extract( & !
  & 1,tNuSp,vLMs, & !IN
  & tNuMs)          !OUT
  !------------------------------------------------------/extract blocks
  !
END SUBROUTINE NuTable_Extract

!~ SUBROUTINE Basis_NuTable_Change(  & !
!~ & vIndxPrimSpc, & !IN
!~ & tStoikio,     & !IN
!~ & tNuSp)          !OUT
!~ !--
!~ !-- NuTable_Build, simplified version for general use
!~ !-- - no extraction of tNuAs, tNuMs
!~ !--
  !~ USE M_Numeric_Mat,ONLY: LU_BakSub
  !~ !
  !~ INTEGER, INTENT(IN) :: vIndxPrimSpc(:)
  !~ REAL(dp),INTENT(IN) :: tStoikio(:,:)
  !~ REAL(dp),INTENT(OUT):: tNuSp(:,:)
  !~ !
  !~ REAL(dp),ALLOCATABLE:: tTransform(:,:)
  !~ INTEGER, ALLOCATABLE:: vIndx(:)
  !~ REAL(dp),ALLOCATABLE:: vY(:)
  !~ LOGICAL :: bSingul
  !~ INTEGER :: I,J,nCp
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_NuTable_Change"
  !~ !
  !~ nCp= SIZE(tStoikio,1)
  !~ !
  !~ ALLOCATE(tTransform(nCp,nCp))
  !~ !--------------------- assemble tTransform from columns of tStoikio --
  !~ J=0
  !~ DO I=1,nCp
    !~ J=J+1
    !~ tTransform(:,J)= tStoikio(:,vIndxPrimSpc(I))
  !~ ENDDO
  !~ !-------------------------> tTransform is tStoikio(nCi | nAx | nMx) --
  !~ !
  !~ ALLOCATE(vIndx(nCp))
  !~ CALL Compute_Transform(tTransform,vIndx,bSingul)
  !~ !
  !~ IF(bSingul) CALL Stop_("SINGUL IN Basis_NuTable_Change")
  !~ !
  !~ ALLOCATE(vY(nCp))
  !~ !
  !~ !-- in tNuSp, species will be in same order as in tStoikio, same as in vSpc --
  !~ DO I=1,SIZE(tNuSp,1)
    !~ vY= tStoikio(1:nCp,I)
    !~ CALL LU_BakSub(tTransform,vIndx,vY)
    !~ !Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    !~ tNuSp(I,1:nCp)= vY(1:nCp)
    !~ !-> this includes transposition between tAlfa and tNu !!
  !~ ENDDO
  !~ !
  !~ DEALLOCATE(tTransform)
  !~ DEALLOCATE(vIndx)
  !~ DEALLOCATE(vY)
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_NuTable_Change"
  !~ !
!~ ENDSUBROUTINE Basis_NuTable_Change

SUBROUTINE Compute_Transform(tTransform,vIndx,Error)
  USE M_Numeric_Mat,ONLY: LU_Decomp
  REAL(dp),INTENT(INOUT):: tTransform(:,:)
  INTEGER, INTENT(OUT)  :: vIndx(:)
  LOGICAL, INTENT(OUT)  :: Error
  !
  REAL(dp):: D
  !
  !--the transformation matrix, tTransform, must be invertible !!
  CALL LU_Decomp(tTransform,vIndx,D,Error)
  !
END SUBROUTINE Compute_Transform

SUBROUTINE Fas_NuTable_Build( & !
& vFas,vMixFas,vMixModel,tNuSp, & !in
& tNuFas)   !out
  USE M_T_Phase,    ONLY: T_Phase
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_T_MixModel, ONLY: T_MixModel
  !
  TYPE(T_Phase),     INTENT(IN) :: vFas(:)
  TYPE(T_MixPhase),  INTENT(IN) :: vMixFas(:)
  TYPE(T_MixModel),  INTENT(IN) :: vMixModel(:)
  REAL(dp),          INTENT(IN) :: tNuSp(:,:)
  !
  REAL(dp),          INTENT(OUT):: tNuFas(:,:)
  !
  TYPE(T_MixModel):: SM
  TYPE(T_MixPhase):: SF
  INTEGER:: I,J,K,iCp,NPole !,nFs
  INTEGER:: nCp
  !
  nCp=SIZE(tNuSp,2)
  !
  DO I=1,SIZE(vFas)
    !
    IF(vFas(I)%iSpc /= 0) THEN  ! "PURE","DISCRET"
      !
      K=vFas(I)%iSpc
      tNuFas(I,:)= tNuSp(K,:)
      !
    ELSE IF(vFas(I)%iMix /= 0) THEN ! "MIXT"
      !
      SF=vMixFas(vFas(I)%iMix) !-> the solution phase
      SM=vMixModel(SF%iModel)  !-> the corresponding model
      NPole=SM%NPole
      !
      ! WHERE(.NOT.SM%vHasPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
      ! WHERE(.NOT.SF%vLPole(1:NPole))   SF%vXPole(1:SM%NPole)=Zero
      !
      DO iCp=1,nCp
        tNuFas(I,iCp)= Zero
        !
        DO J=1,NPole
          IF(SF%vLPole(J)) &
          & tNuFas(I,iCp)= tNuFas(I,iCp) &
          &              + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
        ENDDO
        !
      ENDDO
      !
    ELSE IF (vFas(I)%iSol /= 0) THEN
      !
      K=vFas(I)%iSpc
      tNuFas(I,:)= tNuSp(K,:)
    END IF
  
    !! SELECT CASE(vFas(I)%Typ)
    !! !
    !! CASE("PURE") !,"DISCRET")
    !!   K=vFas(I)%iSpc
    !!   tNuFas(I,:)= tNuSp(K,:)
    !! !--
    !! CASE("MIXT")
    !!   !
    !!   SF=vMixFas(vFas(I)%iMix) !-> the solution phase
    !!   SM=vMixModel(SF%iModel)  !-> the corresponding model
    !!   NPole=SM%NPole
    !!   !
    !!   WHERE(.NOT.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
    !!   !
    !!   DO iCp=1,nCp
    !!     tNuFas(I,iCp)= Zero
    !!     !
    !!     DO J=1,NPole
    !!       IF(SF%vLPole(J)) &
    !!       & tNuFas(I,iCp)= tNuFas(I,iCp) &
    !!       &              + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
    !!     ENDDO
    !!     !
    !!   ENDDO
    !! !--
    !! CASE("SOLU")
    !!   K=vFas(I)%iSpc
    !!   tNuFas(I,:)= tNuSp(K,:)
    !!   
    !! !--
    !! END SELECT
    !
  ENDDO
  !
ENDSUBROUTINE Fas_NuTable_Build

SUBROUTINE Fas_AlfaTable_Build( & !
& vFas,vMixFas,vMixModel, tAlfSp, & !in
& tAlfFs)   !out
  USE M_T_Phase,    ONLY: T_Phase
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_T_MixModel, ONLY: T_MixModel
  !
  TYPE(T_Phase),     INTENT(IN) :: vFas(:)
  TYPE(T_MixPhase),  INTENT(IN) :: vMixFas(:)
  TYPE(T_MixModel),  INTENT(IN) :: vMixModel(:)
  REAL(dp),          INTENT(IN) :: tAlfSp(:,:)
  !
  REAL(dp),          INTENT(OUT):: tAlfFs(:,:)
  !
  TYPE(T_MixModel):: SM
  TYPE(T_MixPhase):: SF
  INTEGER:: I,J
  INTEGER:: iPur,iMix,iSol,iCp
  INTEGER:: nCp,NPole
  !
  nCp=SIZE(tAlfSp,1)
  !
  DO I=1,SIZE(vFas)
    !
    iPur= vFas(I)%iSpc
    iMix= vFas(I)%iMix
    !
    IF(iPur>0) THEN
      tAlfFs(:,I)= tAlfSp(:,iPur)
    ENDIF
    !
    IF(iMix>0) THEN
      !
      SF=vMixFas(iMix) !-> the solution phase
      SM=vMixModel(SF%iModel)  !-> the corresponding model
      NPole=SM%NPole
      !
      WHERE(.NOT.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
      !
      DO iCp=1,nCp
        tAlfFs(iCp,I)= Zero
        !
        DO J=1,NPole
          IF(SF%vLPole(J)) &
          & tAlfFs(iCp,I)= tAlfFs(iCp,I) &
          &              + tAlfSp(iCp,SM%vIPole(J))*SF%vXPole(J)
        ENDDO
     ENDDO
      !
    ENDIF
  ENDDO
  !
ENDSUBROUTINE Fas_AlfaTable_Build

SUBROUTINE Extract(cod,tIn,vL,tOut)
!--
!-- from matrix tIn,
!-- build the matrix tOut
!-- by extracting the rows (cod--1) or the columns (cod--2)
!-- that have the LOGICAL vL true
!--
  INTEGER,                INTENT(IN) :: cod
  REAL(dp),DIMENSION(:,:),INTENT(IN) :: tIn
  LOGICAL, DIMENSION(:),  INTENT(IN) :: vL
  REAL(dp),DIMENSION(:,:),INTENT(OUT):: tOut
  INTEGER::I,J
  J=0
  SELECT CASE(cod)
    CASE(1)
      DO I=1,SIZE(tIn,1)
        IF(vL(I)) THEN  ; J=J+1; tOut(J,:)=tIn(I,:)  ; ENDIF
      ENDDO
    CASE(2)
      DO I=1,SIZE(tIn,2)
        IF(vL(I)) THEN  ; J=J+1; tOut(:,J)=tIn(:,I)  ; ENDIF
      ENDDO
  END SELECT
ENDSUBROUTINE Extract

SUBROUTINE Basis_Cpn_ISpc_Select( &
& vSpc,isW,vCpn, & !in
& vIndex) !out
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component
  USE M_Numeric_Tools !-> sort, iMinLoc_R, iFirstLoc
  !
  TYPE(t_Species),  INTENT(IN) :: vSpc(:)
  INTEGER,          INTENT(IN) :: isW
  TYPE(t_Component),INTENT(IN) :: vCpn(:)
  INTEGER,          INTENT(OUT):: vIndex(:)
  !
  INTEGER :: nCp,nAq
  INTEGER :: iAq,iPr,I,J,K,JJ
  REAL(dp):: Y,Pivot
  REAL(dp),ALLOCATABLE:: tStoik(:,:)
  INTEGER, ALLOCATABLE:: vIPivot(:),vSpPrm(:)
  LOGICAL, ALLOCATABLE:: vSpDone(:)
  LOGICAL, ALLOCATABLE:: vCpDone(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Cpn_ISpc_Select"
  !
  nAq= count(vSpc%Typ=="AQU")
  nCp= SIZE(vCpn)
  !
  vIndex(:)= vCpn(:)%iSpc
  !
  ALLOCATE(vSpPrm(1:nAq))
  ALLOCATE(vSpDone(1:nAq)); vSpDone(1:nAq)=.FALSE.
  ALLOCATE(vCpDone(1:nCp)); vCpDone(1:nCp)=.FALSE.
  !
  CALL CalcPermut(vSpc(1:nAq)%Dat%Mole,vSpPrm)
  ! -> sort aqu.species in order of increasing abudance
  ! -> permutation array vPrm
  ! to apply permutation vSpPrm to array Arr: Arr=Arr(vSpPrm)
  !
  ALLOCATE(tStoik(nCp,nCp))   ;  tStoik(:,:)= Zero
  ALLOCATE(vIPivot(nCp))      ;  vIPivot(:)= 0
  !
  K=   0
  iPr= 0
  !
  !--------------------- species 1 is solvent, component 1 is solvent --
  K=   K+1
  iPr= iPr+1
  !vCpn(iPr)%iSpc= isW
  vIndex(iPr)= isW
  vSpDone(isW)= .true.
  vCpDone(iPr)= .true.
  tStoik(:,K)= vSpc(isW)%vStoikio(1:nCp) !tFormula(:,isW)
  vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
  !! vIPivot(K)= 1 !!???!!!
  Pivot= tStoik(vIPivot(K),K)
  tStoik(:,K)= tStoik(:,K) /Pivot
  !
  !--------------------------------scan the BUFFER and MOBILE components
  DO i=1,nCp
    !
    !-- read BUFFER and MOBILE stoikiometries in order to select,
    !-- in DO00 block,
    !-- primary aqueous species independent from these species
    !
    IF(vCpn(i)%Statut=="BUFFER" .OR. vCpn(i)%Statut=="MOBILE") THEN
      !
      K= K+1
      vIndex(i)= vCpn(i)%iSpc
      tStoik(:,K)= vSpc(vCpn(i)%iSpc)%vStoikio(1:nCp)
      DO J=1,K-1
        Y= tStoik(vIPivot(J),K)
        tStoik(:,K)= tStoik(:,K)- Y*tStoik(:,J)
        !-> substract column vIPivot(J) from column K
        !-> tStoik(vIPivot(J),K) will be 0
      ENDDO
      vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
      !-> index of first non zero element
      IF(iDebug>2) THEN
        DO JJ=1,nCp
          WRITE(fTrc,'(F7.2,1X)',ADVANCE="NO") tStoik(JJ,K)
        ENDDO
        WRITE(fTrc,'(2(A,I3),1X,A)') &
        & "K=",K, " /  iPr=",i,TRIM(vSpc(vIndex(i))%NamSp)
        !PAUSE_
      ENDIF
      Pivot= tStoik(vIPivot(K),K)
      tStoik(:,K)= tStoik(:,K) /Pivot
      !-> first non-zero element, tStoik(vIPivot(K),K), will be 1.
      !
      vCpDone(i)= .true.
      IF(vSpc(vCpn(i)%iSpc)%Typ=="AQU") vSpDone(vCpn(i)%iSpc)= .true.
      !
    ENDIF
    !
  ENDDO
  !-------------------------------/scan the BUFFER and MOBILE components
  !
  Do00: DO iPr=1,nCp
    !
    IF(vCpDone(iPr)) CYCLE
    !
    K= K+1
    IF(K>nCp) EXIT DO00
    !
    I= nAq+1
    !-----------------------------------------------scan aqueous species
    DoAqu: DO
      I= I-1
      IF(I==0) EXIT DoAqu
      !
      !--scan species beginning with most abundant----------------------
      iAq= vSpPrm(i)
      !print *,"Spc=",TRIM(vSpc(iAq)%NamSp),vSpc(iAq)%Dat%Mole  ;  pause_
      IF(vSpDone(iAq)) CYCLE DoAqu
      !
      tStoik(:,K)=  vSpc(iAq)%vStoikio(1:nCp) !tFormula(:,iAq)
      !
      IF(tStoik(vCpn(iPr)%iEle,K)==0) CYCLE DoAqu
      !
      DO J=1,K-1
        Y= tStoik(vIPivot(J),K)
        tStoik(:,K)= tStoik(:,K)- Y*tStoik(:,J)
        !-> tStoik(vIPivot(J),K) will be 0
      ENDDO
      ! vIPivot(K)= IFirstLoc(abs(tStoik(:,K)) >1.D-9)
      vIPivot(K)= iFirstLoc(tStoik(:,K) /= Zero)
      !-> index of first non zero element
      !
      IF(vIPivot(K)>nCp) THEN
      !--all coeff's are 0 -> this species is not independent-----------
      !---> CYCLE DOAqu --
        !
        IF(iDebug==4) THEN
          DO JJ=1,nCp; WRITE(fTrc,'(F7.2,1X)',ADVANCE="NO") tStoik(JJ,K); ENDDO
          WRITE(fTrc,'(A,I3,1X,2A)') &
          & "iFirst=",vIPivot(K),TRIM(vSpc(iAq)%NamSp),"-> removed"
        ENDIF
        !
        vSpDone(iAq)= .true.
        CYCLE DOAqu
      ENDIF
      !
      IF(iDebug==4) THEN
        DO JJ=1,nCp; WRITE(fTrc,'(F7.2,1X)',ADVANCE="NO") tStoik(JJ,K); ENDDO
        WRITE(fTrc,'(2(A,I3),1X,A)') &
        & "K=",K," /  iPr=",iPr,TRIM(vSpc(iAq)%NamSp)
      ENDIF
      !
      Pivot= tStoik(vIPivot(K),K)
      tStoik(:,K)= tStoik(:,K) /Pivot
      !
      !iPr= iPr + 1
      !vCpn(iPr)%iSpc= iAq
      vIndex(iPr)=  iAq
      vCpDone(iPr)= .true.
      vSpDone(iAq)= .true.
      !
      EXIT DoAqu
      !
    ENDDO DoAqu
    !----------------------------------------------/scan aqueous species
  ENDDO Do00
  !
  IF(iDebug==4) THEN
    WRITE(71,'(A,1X)',ADVANCE="NO") "BASIS="
    DO I=1,nCp
      WRITE(71,'(A15,1X)',ADVANCE="NO") vSpc(vIndex(I))%NamSp
    ENDDO
    WRITE(71,*)
  ENDIF
  !
  DEALLOCATE(vSpDone)
  DEALLOCATE(vCpDone)
  DEALLOCATE(vSpPrm)
  DEALLOCATE(tStoik)
  DEALLOCATE(vIPivot)
  !
  IF(iDebug==4) THEN
    PRINT *,"< Basis_Cpn_ISpc_Select"
    DO i=1,nCp
      WRITE(6,'(G15.6,1X,A)') &
      & vSpc(vIndex(I))%Dat%Mole, &
      & TRIM(vSpc(vIndex(I))%NamSp)
    ENDDO
    PRINT *,"</ Basis_Cpn_ISpc_Select"
    CALL Pause_
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Cpn_ISpc_Select"
  !
  RETURN
ENDSUBROUTINE Basis_Cpn_ISpc_Select

SUBROUTINE Basis_Species_Select_Test( & !
& vSpc,isW,vCpn,tFormul, & !in
& vIndex) !out
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  USE M_Numeric_Tools,ONLY: CalcPermut
  USE M_T_Species,    ONLY: T_Species
  USE M_T_Component,  ONLY: T_Component
  !
  TYPE(t_Species),  INTENT(IN) :: vSpc(:)
  INTEGER,          INTENT(IN) :: isW
  TYPE(t_Component),INTENT(IN) :: vCpn(:)
  REAL(dp),         INTENT(IN) :: tFormul(:,:)
  INTEGER,          INTENT(OUT):: vIndex(:)
  !
  INTEGER :: nCp,nAq,nCi
  INTEGER :: i
  INTEGER, ALLOCATABLE:: vOrder(:)
  !
  nAq= COUNT(vSpc%Typ=="AQU")
  nCp= SIZE(vCpn)
  nCi= COUNT(vCpn(:)%Statut=="INERT")
  !
  ALLOCATE(vOrder(nAq))
  !--------------------sort aqu.species in order of decreasing abudance,
  !------------------------------------- -> permutation array vSpvPrm --
  !-------- to apply permutation vSpPrm to array Arr: Arr-Arr(vSpPrm) --
  CALL CalcPermut(-vSpc(1:nAq)%Dat%Mole,vOrder)
  !--------------------/
  !
  vIndex(:)= 0
  vIndex(1)= isW
  vIndex(nCi+1:nCp)= vCpn(nCi+1:nCp)%iSpc
  !
  CALL Basis_Species_Select( &
  !& vSpc,isW,vCpn, & !in
  & tFormul,  & !in
  & vOrder,   & !in
  & vIndex)       !inout
  !! DO i=2,nCp
  !!   IF(vCpn(i)%Statut/="INERT") vIndex(i)= vCpn(i)%iSpc
  !! ENDDO
  !
  !! IF(iDebug==4) THEN
    PRINT *,"< Basis_Cpn_ISpc_Select"
    DO i=1,nAq
      WRITE(6,'(G15.6,1X,A)') &
      & vSpc(vOrder(i))%Dat%Mole,TRIM(vSpc(vOrder(i))%NamSp)
    ENDDO
    PRINT *,"== NEW BASIS =="
    DO i=1,nCp
      WRITE(6,'(G15.6,1X,A)') &
      & vSpc(vIndex(I))%Dat%Mole,TRIM(vSpc(vIndex(I))%NamSp)
    ENDDO
    PRINT *,"</Basis_Cpn_ISpc_Select"
    CALL Pause_
  !! ENDIF
  !
  DEALLOCATE(vOrder)
END SUBROUTINE Basis_Species_Select_Test

SUBROUTINE Basis_Species_Select( &
!& vSpc,isW,vCpn, & !in
& tFormulSpc, & !in
& vOrder,    &  !in
& vIndex)       !inout
!--
!-- build a set of independent species,
!-- comprising the solvent isW, the mobile or buffer species,
!-- and the most abundant aqueous species independent from these
!--
  !USE M_T_Species,    ONLY: T_Species
  !USE M_T_Component,  ONLY: T_Component
  !
  !TYPE(t_Species),  INTENT(IN) :: vSpc(:)
  !INTEGER,          INTENT(IN) :: isW
  !TYPE(t_Component),INTENT(IN) :: vCpn(:)
  REAL(dp),INTENT(IN)   :: tFormulSpc(:,:)
  INTEGER, INTENT(IN)   :: vOrder(:)
  INTEGER, INTENT(INOUT):: vIndex(:)
  !
  ! vIndex(1)= isW
  ! vIndex(nCi+1:nCp)= vCpn(nCi+1:nCp)%iSpc
  !
  INTEGER :: nCp,nAq,nCx
  INTEGER :: I,K
  LOGICAL :: Error
  REAL(dp),ALLOCATABLE:: tStoikSpc(:,:)
  REAL(dp),ALLOCATABLE:: tStoikCpn(:,:)
  !
  IF(iDebug==4) WRITE(fTrc,'(/,A)') "<-------------Basis_Species_Select"
  !
  nAq= SIZE(vOrder)
  nCp= SIZE(vIndex)
  nCx= COUNT(vIndex(:)>0)
  !
  ALLOCATE(tStoikCpn(nCp,nCx))
  ALLOCATE(tStoikSpc(nCp,nAq))
  !
  !------------------------build the stoikio table of imposed components
  K= 0
  DO i=1,nCp
    IF(vIndex(i)>0) THEN
      K= K+1
      tStoikCpn(1:nCp,K)= tFormulSpc(1:nCp,vIndex(i))
    ENDIF
  ENDDO
  !-----------------------/build the stoikio table of imposed components
  DO i=1,nAq
    tStoikSpc(1:nCp,i)= tFormulSpc(1:nCp,i)
  ENDDO
  !
  CALL Basis_FreeSet_Select( & !
  & tStoikSpc,    & !IN
  & vIndex,       & !OUT
  & Error,        & !OUT
  & vOrder=      vOrder,    & !IN
  & tStoikioCpn= tStoikCpn)   !IN
  !
  DEALLOCATE(tStoikSpc)
  DEALLOCATE(tStoikCpn)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</-------------Basis_Species_Select"
  !
  RETURN
ENDSUBROUTINE Basis_Species_Select

SUBROUTINE Basis_FreeSet_Select( & !
& tStoikio,    & !IN
& vIndex,      & !INOUT
& Error,       & !OUT
& vOrder,      & !IN
& tStoikioCpn)   !IN
!--
!-- using Gaussian elimination (without permutations),
!-- find a set of independent species (of stoikios tStoikio)
!-- independent also of those in tStoikioCpn
!--
  USE M_Numeric_Tools,ONLY: iMaxLoc_R
  !
  REAL(dp),INTENT(IN)   :: tStoikio(:,:)
  INTEGER, INTENT(INOUT):: vIndex(:)
  LOGICAL, INTENT(OUT)  :: Error
  INTEGER, INTENT(IN),OPTIONAL:: vOrder(:)
  REAL(dp),INTENT(IN),OPTIONAL:: tStoikioCpn(:,:)
  !
  INTEGER :: nCp,nSp,nX
  INTEGER :: iSp,I,J,K,N
  REAL(dp):: Y,Pivot
  REAL(dp),ALLOCATABLE:: tStoikTmp(:,:)
  INTEGER, ALLOCATABLE:: vIPivot(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "<--------------Basis_FreeSet_Select"
  !
  Error= .FALSE.
  !
  IF(PRESENT(tStoikioCpn)) THEN
    nX= SIZE(tStoikioCpn,2)
  ELSE
    nX= 0
  ENDIF
  nCp= SIZE(tStoikio,1)
  nSp= SIZE(tStoikio,2)
  !
  ALLOCATE(tStoikTmp(nCp,nCp))   ;  tStoikTmp(:,:)= Zero
  ALLOCATE(vIPivot(nCp))         ;  vIPivot(:)= 0
  !
  K= 0
  !----------------------------------------scan imposed basis components
  !--e.g. SOLVENT / BUFFER / MOBILE components--------------------------
  !--(they have vIndex(:)>0)--------------------------------------------
  DO i=1,nX
    !
    !-- read stoikiometries of species already imposed, in order to select,
    !-- in DO_0 block,
    !-- primary aqueous species independent from these species
    !
    K= K+1
    tStoikTmp(:,K)= tStoikioCpn(:,i)
    DO J=1,K-1
      Y= tStoikTmp(vIPivot(J),K)
      tStoikTmp(:,K)= tStoikTmp(:,K)- Y*tStoikTmp(:,J)
      !-> substract column vIPivot(J) from column K
      !-> tStoikTmp(vIPivot(J),K) will be 0
    ENDDO
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(:,K)))
    Pivot= tStoikTmp(vIPivot(K),K)
    IF(ABS(Pivot)<1.D-9) THEN
      != all coeff's are 0 -> this species is not independent -> EXIT ==
      Error= .TRUE.
      RETURN
    ENDIF
    !
    tStoikTmp(:,K)= tStoikTmp(:,K) /Pivot
    !-> first non-zero element, tStoikTmp(vIPivot(K),K), will be 1.
    !
  ENDDO
  !---------------------------------------/scan imposed basis components
  !
  N=1
  !-----------------------------------------------scan the other species
  DO_0: DO iSp=1,nSp
    !
    K= K+1
    IF(K>nCp) EXIT DO_0
    !
    IF(PRESENT(vOrder)) THEN ;  tStoikTmp(:,K)= tStoikio(:,vOrder(iSp))
    ELSE                     ;  tStoikTmp(:,K)= tStoikio(:,iSp)
    ENDIF
    !
    DO J=1,K-1
      Y= tStoikTmp(vIPivot(J),K)
      tStoikTmp(:,K)= tStoikTmp(:,K)- Y*tStoikTmp(:,J)
      !-> tStoikTmp(vIPivot(J),K) will be 0
    ENDDO
    vIPivot(K)= iMaxLoc_R(ABS(tStoikTmp(:,K)))
    Pivot= tStoikTmp(vIPivot(K),K)
    IF(ABS(Pivot)<1.D-9) THEN
      !-- all coeff's are 0 -> this species is not independent -> CYCLE DO_0 --
      ! IF(iDebug>2) THEN
      !   DO JJ=1,nCp; WRITE(6,'(F7.2,1X)',ADVANCE="NO") tStoikTmp(JJ,K); ENDDO
      !   WRITE(6,'(2(A,I3),1X,A)') "K=",K," /  iPr=",K,"Reject"
      ! ENDIF
      K= K-1
      CYCLE DO_0
    ENDIF
    !
    ! IF(iDebug>2) THEN
    !   DO JJ=1,nCp; WRITE(6,'(F7.2,1X)',ADVANCE="NO") tStoikTmp(JJ,K); ENDDO
    !   WRITE(6,'(2(A,I3),1X,A)') "K=",K," /  iPr=",K,"Accept"
    ! ENDIF
    !
    tStoikTmp(:,K)= tStoikTmp(:,K) /Pivot
    !
    !--TRICK: inert aqueous species are placed--------------------------
    !--AFTER solvent and BEFORE mobile species !!!----------------------
    ! vIndex(K-nX+1)= vOrder(iSp)
    DO WHILE(vIndex(N)>0)
      N= N+1
    ENDDO
    IF(PRESENT(vOrder)) THEN ;  vIndex(N)= vOrder(iSp)
    ELSE                     ;  vIndex(N)= iSp
    ENDIF
    !
  ENDDO DO_0
  !----------------------------------------------/scan the other species
  !
  DEALLOCATE(tStoikTmp)
  DEALLOCATE(vIPivot)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</-------------Basis_FreeSet_Select"
  !
  RETURN
ENDSUBROUTINE Basis_FreeSet_Select

SUBROUTINE Basis_AlfaNu_Build( & !
& vEle,vCpn,vSpc,tStoikio,     & !IN
& vFas,vMixFas,vMixModel,      & !
& vLCi,vLAx,vLMx,vLAs,vLMs,    & !
& LSho,LBuffer,                & !
& vOrdPr,vOrdAs,vOrdMs,        & !
& vTot,                        & !
& tAlfSp,                      & !OUT
& tAlfPr,tAlfAs,tAlfMs,        & !
& tNuSp,tNuAs,tNuMs,           & !
& tNuFas,tAlfFs)                 !

  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Phase,    ONLY: T_Phase
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_T_MixModel, ONLY: T_MixModel
  !
  TYPE(T_Element),  DIMENSION(:),  INTENT(IN):: vEle
  TYPE(T_Component),DIMENSION(:),  INTENT(IN):: vCpn
  TYPE(T_Species),  DIMENSION(:),  INTENT(IN):: vSpc
  REAL(dp),         DIMENSION(:,:),INTENT(IN):: tStoikio
  TYPE(T_Phase),    DIMENSION(:),  INTENT(IN):: vFas
  TYPE(T_MixPhase), DIMENSION(:),  INTENT(IN):: vMixFas
  TYPE(T_MixModel), DIMENSION(:),  INTENT(IN):: vMixModel
  LOGICAL,          DIMENSION(:),  INTENT(IN):: vLCi,vLAx,vLMx,vLAs,vLMs
  LOGICAL,                         INTENT(IN):: LSho,LBuffer
  INTEGER,          DIMENSION(:),  INTENT(IN):: vOrdPr,vOrdAs,vOrdMs
  !
  REAL(dp),         DIMENSION(:),  INTENT(INOUT):: vTot
  REAL(dp),         DIMENSION(:,:),INTENT(OUT)::   tAlfSp !,tAlfEle
  REAL(dp),         DIMENSION(:,:),INTENT(OUT)::   tAlfPr,tAlfAs,tAlfMs
  REAL(dp),         DIMENSION(:,:),INTENT(OUT)::   tNuSp,tNuAs,tNuMs
  REAL(dp),         DIMENSION(:,:),INTENT(OUT)::   tNuFas,tAlfFs
  !
  INTEGER:: nCi,nAs,nMs,nAx,nMx
  INTEGER:: I,J
  INTEGER:: vIndxPrimSpc(SIZE(vCpn))
  INTEGER:: vIndxBuffers(SIZE(vCpn))
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "<----------------Basis_AlfaNu_Build"
  !
  nCi= COUNT(vLCi)
  nAx= COUNT(vLAx)
  nMx= COUNT(vLMx)
  nAs= COUNT(vLAs)
  nMs= COUNT(vLMs)
  !
  vIndxBuffers(:)= 0  !
  IF(LBuffer) THEN
    DO I=1,SIZE(vCpn)
      IF(vCpn(I)%Statut=="BUFFER") vIndxBuffers(I)= vCpn(I)%iSpc
      !IF(vCpn(I)%Statut=="BUFFER" .OR. vCpn(I)%Statut=="MOBILE") &
      !& vIndxBuffers(I)= vCpn(I)%iSpc
    ENDDO
  ENDIF
  !
  CALL AlfaTable_Build( &
  & SIZE(vCpn),SIZE(vSpc),tStoikio, & !
  & vIndxBuffers, & !
  & vTot(:), & !
  & tAlfSp) !,tAlfEle)
  !
  CALL Fas_AlfaTable_Build( & !
  & vFas,vMixFas,vMixModel,tAlfSp, & !IN
  & tAlfFs)   !OUT
  !
  CALL AlfaTable_Extract( &
  & vCpn,tAlfSp,vLCi,vLAx,vLMx,vLAs,vLMs, &
  & tAlfPr,tAlfAs,tAlfMs)
  !
  !--------- component order in vCpn is Inert, Aqu'Mobile, Min'Mobile --
  J= 0
  !-------------------------------------- index of inert prim'species --
  DO I=1,nCi  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(I)%iSpc  ;  ENDDO
  !-------------------------------------- index of mobile aqu'species --
  DO I=1,nAx  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(J)%iSpc  ;  ENDDO
  !-------------------------------------- index of mobile min'species --
  DO I=1,nMx  ; J= J+1  ;  vIndxPrimSpc(J)= vCpn(J)%iSpc  ;  ENDDO
  !
  CALL NuTable_Build( & !
  & vIndxPrimSpc,tStoikio, & !IN
  & tNuSp) !OUT
  !
  CALL NuTable_Extract( & !
  & tNuSp,vLAs,vLMs, & !IN
  & tNuAs,tNuMs) !OUT
  !
  CALL Fas_NuTable_Build( & !
  & vFas,vMixFas,vMixModel,tNuSp, & !IN
  & tNuFas)   !OUT
  !
  IF(LSho) CALL AlfaTable_Sho( &
  & vEle,vCpn,vSpc,vFas,vOrdPr,vOrdAs,vOrdMs, &
  & tAlfPr,tAlfAs,tAlfMs,tAlfFs)
  !
  IF(LSho) CALL NuTable_Sho( &
  & vCpn,vSpc,vLCi,vLAx,vLMx,vLAs,vLMs, &
  & vOrdPr,vOrdAs,vOrdMs,tNuSp,tNuAs,tNuMs, &
  & .TRUE.)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</---------------Basis_AlfaNu_Build"
  !
ENDSUBROUTINE Basis_AlfaNu_Build

SUBROUTINE Basis_StoikTable_Sho( &
& vEle,vCpn,vSpc, &         !in
& tStoikio, &               !in
& vLCi,vLAx,vLMx,vLAs,vLMs) !in

  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut,Files_Index_Write
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  !
  TYPE(T_Element),  INTENT(IN):: vEle(:)
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  TYPE(T_Species),  INTENT(IN):: vSpc(:)
  REAL(dp),         INTENT(IN):: tStoikio(:,:)
  LOGICAL,          INTENT(IN):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !
  INTEGER:: F,I,iPr,nCp,nSp
  !
  nSp=SIZE(vSpc)
  nCp=SIZE(vCpn)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_stoik_for.log")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_stoik_for.log",&
  & "STOIKIO: stoichiometry matrix tStoikio, transposed")
  !
  !--- header
  DO iPr=1,nCp
    WRITE(F,'(A3,A1)',ADVANCE='NO') vEle(vCpn(iPr)%iEle)%NamEl, T_
  ENDDO
  WRITE(F,'(3(A,A1))') "Chg",T_,"_____species",T_,"STATUS",T_
  !---/
  !
  DO I=1,nSp
    DO iPr=1,nCp
      WRITE(F,'(F7.2,A1)',ADVANCE='NO') tStoikio(iPr,I), T_
    ENDDO
    WRITE(F,'(I3,A1,A12,A1)',ADVANCE='NO') vSpc(I)%Z, T_, TRIM(vSpc(I)%NamSp),T_
    IF(vLCi(I)) WRITE(F,'(A)') "INT_PRIMARY"
    IF(vLAx(I)) WRITE(F,'(A)') "EXT_AQUEOUS"
    IF(vLMx(I)) WRITE(F,'(A)') "EXT_MINERAL"
    IF(vLAs(I)) WRITE(F,'(A)') "INT_AQUEOUS"
    IF(vLMs(I)) WRITE(F,'(A)') "INT_MINERAL"
  ENDDO
  !
  CLOSE(F)
  !
ENDSUBROUTINE Basis_StoikTable_Sho

SUBROUTINE NuTable_Sho( &
& vCpn,vSpc, &
& vLCi,vLAx,vLMx,vLAs,vLMs, &
& vOrdPr,vOrdAs,vOrdMs, &
& tNuSp,tNuAs,tNuMs, &
& LogK_Sho)
  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut,Files_Index_Write
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  !
  TYPE(T_Component),DIMENSION(:),  INTENT(IN):: vCpn
  TYPE(T_Species),  DIMENSION(:),  INTENT(IN):: vSpc
  LOGICAL,          DIMENSION(:),  INTENT(IN):: vLCi,vLAx,vLMx,vLAs,vLMs
  INTEGER,          DIMENSION(:),  INTENT(IN):: vOrdPr,vOrdAs,vOrdMs
  REAL(dp),         DIMENSION(:,:),INTENT(IN):: tNuSp,tNuAs,tNuMs
  LOGICAL,                         INTENT(IN):: LogK_Sho
  !
  INTEGER :: F,iPr,I
  INTEGER :: nSp,nCp,nAs,nMs
  REAL(dp):: X,rScale
  !
  nSp=SIZE(vSpc)
  nCp=SIZE(vCpn)
  nAs=count(vLAs)
  nMs=count(vLMs)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_stoik_nu.log") !"OUT_STOIK2.TXT")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_stoik_nu.log",&
  & "STOIKIO: Nu Table: stoikio of second'species in terms of prim'species")
  !
  WRITE(F,'(2A)') &
  & "!.tNuSp, ALL SPECIES", &
  & "!.stoichiometry, logK_dtb, logK_reaction, logK_reaction/scale"
  !
  DO iPr=1,nCp
    WRITE(F,'(A7,A1)',ADVANCE='NO') TRIM(vSpc(vOrdPr(iPr))%NamSp),T_
  ENDDO
  !
  WRITE(F,'(A)') "Status and Name(iSpc)"
  !
  DO I=1,nSp
    !
    DO iPr=1,nCp
      IF(tNuSp(I,iPr)/=0) THEN; WRITE(F,'(F7.2,A1)',ADVANCE='NO') tNuSp(I,iPr), T_
      ELSE;                     WRITE(F,'(A7,  A1)',ADVANCE='NO')          ".", T_
      ENDIF
    ENDDO
    !
    IF(vLCi(I)) WRITE(F,'(A,A1)',ADVANCE='NO') "INT_PRIMARY",T_
    IF(vLAx(I)) WRITE(F,'(A,A1)',ADVANCE='NO') "EXT_AQUEOUS",T_
    IF(vLMx(I)) WRITE(F,'(A,A1)',ADVANCE='NO') "EXT_MINERAL",T_
    IF(vLAs(I)) WRITE(F,'(A,A1)',ADVANCE='NO') "INT_AQUEOUS",T_
    IF(vLMs(I)) WRITE(F,'(A,A1)',ADVANCE='NO') "INT_MINERAL",T_
    !
    IF(LogK_Sho) THEN
      X= vSpc(I)%G0rt - DOT_PRODUCT(tNuSp(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
      rScale= SUM(ABS(tNuSp(I,1:nCp))) !-> scaling factor !!!
      WRITE(F,'(G15.6,A1)',ADVANCE='NO') rScale,T_
      IF(rScale==Zero) rScale= One
      WRITE(F,'(3(F15.6,A1))',ADVANCE='NO') &
      & -vSpc(I)%G0rt /Ln10,T_,X/Ln10,T_,X/Ln10/rScale,T_
      != +log_10(K of Formation Reaction)
    ENDIF
    !
    WRITE(F,'(A)') TRIM(vSpc(I)%NamSp)
    !
  ENDDO
  !
  IF(nAs>0) THEN

    WRITE(F,'(/,A,/)') "!.tNuAs, SECOND.AQU, logK_dtb, logK_reaction"

    DO iPr=1,nCp
      WRITE(F,'(A7,A1)',ADVANCE='NO') TRIM(vSpc(vOrdPr(iPr))%NamSp),T_
    ENDDO

    WRITE(F,'(A)') "Name(iAs)"

    DO I=1,nAs
      DO iPr=1,nCp !---------------------------------------write stoikio
        WRITE(F,'(F7.2,A1)',ADVANCE='NO') tNuAs(I,iPr), T_
      ENDDO
      !
      IF(LogK_Sho) THEN
        X=  vSpc(vOrdAs(I))%G0rt &
        & - DOT_PRODUCT(tNuAs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
        !X=vG0As(I) - DOt_product(tNuAs(I,1:nCp),vG0Pr(1:nCp))
        WRITE(F,'(2(F15.6,A1))',ADVANCE='NO') &
        & -vSpc(vOrdAs(I))%G0rt /Ln10,T_,X/Ln10,T_ !=+log10K of Formation Reaction
      ENDIF
      !
      WRITE(F,'(A)') TRIM(vSpc(vOrdAs(I))%NamSp)
      !
    ENDDO

  ENDIF
  !
  IF(nMs>0) THEN

    WRITE(F,'(/,A,/)') "!.tNuMs, SECOND.MIN, logK_dtb, logK_reaction"

    DO iPr=1,nCp
      WRITE(F,'(A7,A1)',ADVANCE='NO') TRIM(vSpc(vOrdPr(iPr))%NamSp),T_
    ENDDO

    WRITE(F,'(A)') "Name(iMs)"

    DO I=1,nMs
      !
      DO iPr=1,nCp !===================================< write stoikio==
        WRITE(F,'(F7.2,A1)',ADVANCE='NO') tNuMs(I,iPr), T_
      ENDDO
      !
      IF(LogK_Sho) THEN
        X=  vSpc(vOrdMs(I))%G0rt &
        & - DOT_PRODUCT(tNuMs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
        !X=vG0As(iAs) - DOt_product(tNuAs(iAs,1:nCp),vG0Pr(1:nCp))
        WRITE(F,'(2(F15.6,A1))',ADVANCE='NO') &
        & -vSpc(vOrdMs(I))%G0rt /Ln10,T_,X/Ln10,T_ !=+log10K of Formation Reaction
      ENDIF
      !
      WRITE(F,'(A)') TRIM(vSpc(vOrdMs(I))%NamSp)
      !
    ENDDO

  ENDIF
  !
  CLOSE(F)
ENDSUBROUTINE NuTable_Sho

SUBROUTINE AlfaTable_Sho( &
& vEle,vCpn,vSpc,vFas,    &
& vOrdPr,vOrdAs,vOrdMs,   &
& tAlfPr,tAlfAs,tAlfMs,tAlfFs)
  !---------------------------------------------------------------------
  USE M_IoTools,    ONLY: GetUnit
  USE M_Files,      ONLY: DirOut,Files_Index_Write
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Phase,    ONLY: T_Phase
  !
  TYPE(T_Element),  DIMENSION(:),  INTENT(IN):: vEle
  TYPE(T_Component),DIMENSION(:),  INTENT(IN):: vCpn
  TYPE(T_Species),  DIMENSION(:),  INTENT(IN):: vSpc
  TYPE(T_Phase),    DIMENSION(:),  INTENT(IN):: vFas
  INTEGER,          DIMENSION(:),  INTENT(IN):: vOrdPr,vOrdAs,vOrdMs
  REAL(dp),         DIMENSION(:,:),INTENT(IN):: tAlfPr,tAlfAs,tAlfMs,tAlfFs
  !
  INTEGER:: nCp,nCi,nAs,nMs,nFs
  INTEGER:: F
  INTEGER:: I,J
  !---------------------------------------------------------------------
  nCp= SIZE(vCpn)
  nFs= SIZE(vFas)
  nCi= count(vCpn(:)%Statut=="INERT") + count(vCpn(:)%Statut=="BUFFER")
  nAs= SIZE(tAlfAs,2)
  nMs= SIZE(tAlfMs,2)
  !nAs=count(vLAs)
  !nMs=count(vLMs)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_stoik_alfa.log") !"out_stoik1.txt")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_stoik_alfa.log",&
  & "STOIKIO: Alfa Tables: stoikio of all species in terms of elements / prim'species")
  !
  WRITE(F,'(/,A,/)') "tAlfPr, Prim'Species"
  CALL NamCp_Sho
  DO I=1,nCp
    DO J=1,nCp; WRITE(F,'(F7.3,A1)',ADVANCE='NO') tAlfPr(J,I), T_; ENDDO
    WRITE(F,'(I3,A1,A12)') vSpc(vOrdPr(I))%Z, T_, TRIM(vSpc(vOrdPr(I))%NamSp)
  ENDDO
  !
  IF(nAs>0) THEN
    WRITE(F,'(/,A,/)') "tAlfAs, Aqu'Sec'Species"
    CALL NamCp_Sho
    DO I=1,nAs
      DO J=1,nCp; WRITE(F,'(F7.3,A1)',ADVANCE='NO') tAlfAs(J,I), T_; ENDDO
      WRITE(F,'(I3,A1,A12)') vSpc(vOrdAs(I))%Z, T_, TRIM(vSpc(vOrdAs(I))%NamSp)
    ENDDO
  ENDIF
  IF(nMs>0) THEN
    WRITE(F,'(/,A,/)') "tAlfMs, Min'Sec'Species"
    CALL NamCp_Sho
    DO I=1,nMs
      DO J=1,nCp; WRITE(F,'(F7.3,A1)',ADVANCE='NO') tAlfMs(J,I), T_; ENDDO
      WRITE(F,'(I3,A1,A12)') vSpc(vOrdMs(I))%Z, T_, TRIM(vSpc(vOrdMs(I))%NamSp)
    ENDDO
  ENDIF
  ! IF(nFs>0) THEN
  !   WRITE(F,'(/,A,/)') "tAlfFs, Phases"  
  !   DO J=1,nCi;     WRITE(F,'(A3,A1)',ADVANCE='NO') vEle(vCpn(J)%iEle)%NamEl, T_; ENDDO
  !   DO J=nCi+1,nCp; WRITE(F,'(A7,A1)',ADVANCE='NO') vEle(vCpn(J)%iEle)%NamEl, T_; ENDDO
  !   WRITE(F,*)
  !   DO I=1,nFs
  !     DO J=1,nCp; WRITE(F,'(F7.3,A1)',ADVANCE='NO') tAlfFs(J,I), T_; ENDDO
  !     WRITE(F,'(A)') TRIM(vFas(I)%NamFs)
  !   ENDDO
  ! ENDIF
  CLOSE(F)

CONTAINS

SUBROUTINE NamCp_Sho
  INTEGER:: K
  DO K=1,nCi;     WRITE(F,'(A3,A1)',ADVANCE='NO') vCpn(K)%NamCp, T_; ENDDO
  DO K=nCi+1,nCp; WRITE(F,'(A7,A1)',ADVANCE='NO') vCpn(K)%NamCp, T_; ENDDO
  WRITE(F,*)
END SUBROUTINE NamCp_Sho

ENDSUBROUTINE AlfaTable_Sho

SUBROUTINE Basis_FindIBal( & !
& vCpn,vSpc,               & ! IN
& iH_,iBal)                  ! OUT
!--
!-- find index of the BALANCE element -> iBal
!--
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  TYPE(T_Species),  INTENT(IN):: vSpc(:)
  INTEGER,          INTENT(IN) :: iH_
  INTEGER,          INTENT(OUT):: iBal
  !---------------------------------------------------------------------
  INTEGER::I
  !
  !-------------- finds for which element the mole balance constraint --
  !------------------------------ is replaced by a charge balance one --
  iBal= 0
  DO I=1,SIZE(vCpn)
    IF(TRIM(vCpn(I)%Statut)=="BALANCE") THEN
      iBal=I
      EXIT
    ENDIF
  ENDDO
  !
  !<new, deleted 200911>
  !~ IF(iBal>0) THEN
    !~ IF (vSpc(vCpn(iBal)%iSpc)%Z==0) &
    !~ & CALL Stop_("Balance Species MUST BE Charged Species")
  !~ ENDIF
  !</new>
  !
  IF(iBal==0 .AND. vCpn(iH_)%Statut=="INERT") iBal= iH_
  !
  !---------------------------- option ELECTRONEUTRALITY NOT ENFORCED --
  !=< used mainly for construction of pH diagram using a PATH CHANGE, ==
  !----------------------------------- not for speciation of a system --
  IF(iBal==0 &
  & .AND. vCpn(iH_)%Statut/="INERT" &
  & .AND. iDebug>0) &
  & PRINT '(A)', &
  & "CAVEAT: NO Element For Electron Balance -> ELECTRONEUTRALITY NOT ENFORCED"
  !! & CALL Stop_("NO Element For Electron Balance ????")
  !
ENDSUBROUTINE Basis_FindIBal

ENDMODULE M_Basis_Tools

!~ SUBROUTINE AlfaTable_Calc( &
!~ & vCpn,vSpc,tStoikio,vLCi,vLAx,vLMx,vLAs,vLMs,bOH2, &
!~ & tAlfSp,tAlfPr,tAlfAs,tAlfMs)
!~ !.called in Basis_AlfaNu_Build
  !~ USE M_T_Component,ONLY: T_Component
  !~ USE M_T_Species,  ONLY: T_Species
  !~ USE M_Numeric_Mat,      ONLY: LU_Decomp, LU_BakSub
  !~ !
  !~ TYPE(T_Component),INTENT(IN):: vCpn(:)
  !~ TYPE(T_Species),  INTENT(IN):: vSpc(:)
  !~ REAL(dp), INTENT(IN):: tStoikio(:,:)
  !~ LOGICAL,  INTENT(IN):: vLCi(:),vLAx(:),vLMx(:),vLAs(:),vLMs(:)
  !~ LOGICAL,  INTENT(IN):: bOH2
  !~ !
  !~ REAL(dp), INTENT(OUT):: tAlfSp(:,:),tAlfPr(:,:),tAlfAs(:,:),tAlfMs(:,:)
  !~ !
  !~ REAL(dp),ALLOCATABLE:: tTransform(:,:)
  !~ INTEGER, ALLOCATABLE:: Indx(:)
  !~ REAL(dp),ALLOCATABLE:: Y(:)
  !~ !
  !~ LOGICAL :: bSingul
  !~ !
  !~ REAL(dp):: D
  !~ INTEGER :: I,nCp,nSp,nCi,nAx,nMx,nAs,nMs
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< AlfaTable_Calc"
  !~ !
  !~ nCp= SIZE(vCpn)
  !~ nSp= SIZE(vSpc)
  !~ nCi= count(vLCi)
  !~ nAs= count(vLAs)
  !~ nAx= count(vLAx)
  !~ nMx= count(vLMx)
  !~ nMs= count(vLMs)
  !~ !
  !~ ALLOCATE(tTransform(nCp,nCp))
  !~ ALLOCATE(Indx(nCp))
  !~ ALLOCATE(Y(nCp))
  !~ !
  !~ !tStoikio:column iSp = stoikio of species iSp in terms of elements (1:nCp)
  !~ !tAlf_:column iSp = stoikio of species iSp in terms of nCi elements & nCm components
  !~ !
  !~ !extract blocks from tStoikio(1:nCp,1:nSp) to tAlfPr,tAlfAs,tAlfMs
  !~ !
  !~ tAlfSp(1:nCp,1:nSp)= tStoikio(1:nCp,1:nSp)
  !~ !
  !~ !DO I=1,nCp
  !~ !  tTransform(:,I)= vCpn(I)%vStoikio(1:nCp)
  !~ !ENDDO
  !~ !
  !~ IF(bOH2) THEN
    !~ tTransform(1:nCp,1:nCp)=Zero
    !~ DO I=1,nCp; tTransform(I,I)=One; ENDDO !-> identity matrix
    !~ tTransform(1:nCp,1)= tStoikio(1:nCp,1) !-> column 1, O replaced by OH2
    !~ CALL LU_Decomp(tTransform,Indx,D,bSingul)
    !~ IF(bSingul) CALL Stop_("SINGUL IN AlfaTableCalc")
    !~ !-> tTransform gives Elements as Functions of Int'Elements/Ext'Species
    !~ !
    !~ DO I=1,nSp
      !~ Y= tStoikio(1:nCp,I)
      !~ CALL LU_BakSub(tTransform,Indx,Y)
      !~ tAlfSp(1:nCp,I)=Y(1:nCp)
    !~ ENDDO
    !~ !
  !~ ENDIF
  !~ !
  !~ DEALLOCATE(tTransform)
  !~ DEALLOCATE(Indx)
  !~ DEALLOCATE(Y)
  !~ !
  !~ !__________________________________________subblock of inert prim'species
  !~ DO I=1,nCi
    !~ tAlfPr(:,I)=tAlfSp(:,vCpn(I)%iSpc)
  !~ ENDDO
  !~ !__________________________________________subblock of mobile aqu'species
  !~ DO I=1,nAx
    !~ tAlfPr(:,nCi+I)=tAlfSp(:,vCpn(nCi+I)%iSpc)
  !~ ENDDO
  !~ !__________________________________________subblock of mobile min'species
  !~ DO I=1,nMx
    !~ tAlfPr(:,nCi+nAx+I)=tAlfSp(:,vCpn(nCi+nAx+I)%iSpc)
  !~ ENDDO
  !~ !
  !~ IF(nAs>0) CALL Extract(2,tAlfSp,vLAs,tAlfAs) !extract columns of sec'aqu.species
  !~ IF(nMs>0) CALL Extract(2,tAlfSp,vLMs,tAlfMs) !extract columns of sec'min.gas.species
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ AlfaTable_Calc"
  !~ !
!~ ENDSUBROUTINE AlfaTable_Calc

!~ SUBROUTINE FasTable_Calc( & !
!~ & vCpn,vFas,vMixFas,vMixModel, & !in
!~ & tNuSp, tAlfSp, & !in
!~ & tNuFas,tAlfFs)   !out

  !~ USE M_T_Component,ONLY: T_Component
  !~ USE M_T_Phase,    ONLY: T_Phase
  !~ USE M_T_MixPhase, ONLY: T_MixPhase
  !~ USE M_T_MixModel, ONLY: T_MixModel
  !~ !
  !~ TYPE(T_Component), INTENT(IN) :: vCpn(:)
  !~ TYPE(T_Phase),     INTENT(IN) :: vFas(:)
  !~ TYPE(T_MixPhase),  INTENT(IN) :: vMixFas(:)
  !~ TYPE(T_MixModel),  INTENT(IN) :: vMixModel(:)
  !~ REAL(dp),          INTENT(IN) :: tNuSp(:,:), tAlfSp(:,:)
  !~ !
  !~ REAL(dp),          INTENT(OUT):: tNuFas(:,:),tAlfFs(:,:)
  !~ !
  !~ TYPE(T_MixModel):: SM
  !~ TYPE(T_MixPhase):: SF
  !~ INTEGER:: I,J,K,iCp,NPole !,nFs
  !~ INTEGER:: nCp
  !~ !
  !~ nCp=SIZE(vCpn)
  !~ !nFs=SIZE(vFas)
  !~ !
  !~ DO I=1,SIZE(vFas)
    !~ !
    !~ SELECT CASE(vFas(I)%Typ)
      !~ !
      !~ CASE("PURE","DISCRET")
        !~ K=vFas(I)%iSpc
        !~ tNuFas(I,1:nCp)= tNuSp(K,1:nCp)
        !~ tAlfFs(1:nCp,I)= tAlfSp(1:nCp,K)
      !~ !
      !~ CASE("MIXT")
        !~ !
        !~ SF=vMixFas(vFas(I)%iMix) !-> the solution phase
        !~ SM=vMixModel(SF%iModel)  !-> the corresponding model
        !~ NPole=SM%NPole
        !~ !
        !~ WHERE(.NOT.SF%vLPole(1:NPole)) SF%vXPole(1:SM%NPole)=Zero
        !~ !
        !~ DO iCp=1,nCp
          !~ tNuFas(I,iCp)= Zero
          !~ tAlfFs(iCp,I)= Zero
          !~ !
          !~ DO J=1,NPole
            !~ IF(SF%vLPole(J)) THEN
              !~ tNuFas(I,iCp)= tNuFas(I,iCp) &
              !~ &            + tNuSp(SM%vIPole(J),iCp) *SF%vXPole(J)
              !~ tAlfFs(iCp,I)= tAlfFs(iCp,I) &
              !~ &            + tAlfSp(iCp,SM%vIPole(J))*SF%vXPole(J)
            !~ ENDIF
          !~ ENDDO
          !~ !
        !~ ENDDO
      !~ !
    !~ END SELECT
  !~ ENDDO
  !~ !
!~ ENDSUBROUTINE FasTable_Calc

!~ SUBROUTINE Basis_Calc_iW_iH_iOH_iO2( & !
!~ & vSpc, &                  !in
!~ & isW,isH_,isOH,isO2,MWsv) !out
!~ !--
!~ !-- find indexes of species H+, OH-, O2 in vSpc (unchanged for given vSpc)
!~ !--
  !~ USE M_T_Species,ONLY: T_Species,Species_Index
  !~ !
  !~ TYPE(T_Species),INTENT(IN) :: vSpc(:)
  !~ INTEGER,        INTENT(OUT):: isW,isH_,isOH,isO2
  !~ REAL(dp),       INTENT(OUT):: MWsv
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Basis_Calc_iW_iH_iOH_iO2"
  !~ !
  !~ !-------------- ranks of species H+ and OH- in current species list --
  !~ isW=  Species_Index("H2O",vSpc)
  !~ isOH= Species_Index("OH-",  vSpc); IF(isOH==0) isOH=Species_Index("OH[-]",vSpc)
  !~ isH_= Species_Index("H+",   vSpc); IF(isH_==0) isH_=Species_Index("H[+]",vSpc)
  !~ isO2= Species_Index("O2,AQ",vSpc); IF(isO2==0) isO2=Species_Index("O2(AQ)",vSpc)
  !~ isO2= Species_Index("O2_AQ",vSpc)
  !~ !
  !~ IF(isW--0)  CALL Stop_("species H2O not found !!!") !-------------STOP
  !~ IF(isH_--0) CALL Stop_("species H+ not found  !!!") !-------------STOP
  !~ IF(isOH--0) CALL Stop_("species OH- not found !!!") !-------------STOP
  !~ !
  !~ MWsv=vSpc(isW)%WeitKg
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(2(A,1X,I3))') "  isH_",isH_,"  isOH",isOH
  !~ !
  !~ IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Basis_Calc_iW_iH_iOH_iO2"
!~ ENDSUBROUTINE Basis_Calc_iW_iH_iOH_iO2

!~ SUBROUTINE Basis_RedoxAssign(vCpn,vEle,vSpc)
!~ !------------------------------------ assign redox states of elements --
  !~ USE M_T_Element,  ONLY: T_Element
  !~ USE M_T_Component,ONLY: T_Component
  !~ USE M_T_Species,  ONLY: T_Species
  !~ !
  !~ TYPE(T_Component),INTENT(IN) :: vCpn(:)
  !~ TYPE(T_Element),  INTENT(IN) :: vEle(:)
  !~ TYPE(T_Species),  INTENT(IN) :: vSpc(:)
  !~ !
  !~ TYPE(T_Species):: S_
  !~ INTEGER:: nCp,I, ZSp, Z_
  !~ LOGICAL:: fOk
  !~ !
  !~ nCp=SIZE(vCpn)
  !~ !
  !~ DO I=1,nCp
    !~ S_=vSpc(vCpn(I)%iSpc)
    !~ !
    !~ !calculate redox state of element in its Prim'Species -> basis valency
    !~ Z_= DOT_PRODUCT( S_%vStoikio(1:nCp),vEle(1:nCp)%Z ) &
    !~ & - S_%vStoikio(I)*vEle(I)%Z
    !~ !
    !~ !restrict this procedure to ele's with no valency assigned in dtb (Fe,S,...)
    !~ IF(vEle(I)%Redox=="VAR") THEN
      !~ !vEle(I)%Z=(ZSp-Z_)/vStoik(I)
      !~ IF(S_%vStoikio(I)>0) &
      !~ WRITE(fTrc,'(2A,2(A,I3),/)') &
      !~ & "valency of ",vEle(I)%NamEl,&
      !~ & ", Old=",vEle(I)%Z,&
      !~ & ", New=",(ZSp-Z_) /S_%vStoikio(I)
    !~ ENDIF
    !~ !
  !~ ENDDO
  !~ !
!~ ENDSUBROUTINE Basis_RedoxAssign

!~ SUBROUTINE Basis_Calc_ieO_ieH_ieOx(vEle,ieO_,ieH_,ieOx)
!~ != find indexes of elements H and OX in vEle
!~ != -> "static" (stable for a given database)
!~ != used for an aqueous system
  !~ USE M_T_Element,ONLY: T_Element,Element_Index
  !~ !
  !~ TYPE(T_Element),DIMENSION(:),INTENT(IN) :: vEle
  !~ INTEGER,                     INTENT(OUT):: ieO_,ieH_,ieOx
  !~ !
  !~ ieO_= Element_Index("O__",vEle)
  !~ ieH_= Element_Index("H__",vEle)
  !~ ieOx= Element_Index("OX_",vEle)
  !~ !
!~ ENDSUBROUTINE Basis_Calc_ieO_ieH_ieOx

