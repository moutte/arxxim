MODULE M_Basis_Vars

  USE M_Kinds 
  USE M_Trace,ONLY: Stop_
  
  IMPLICIT NONE
  
  PRIVATE
  !
  PUBLIC:: Basis_SpcOrdPrm_Alloc
  PUBLIC:: Basis_SpcTypes_Alloc
  PUBLIC:: Basis_AlfaNu_Alloc
  PUBLIC:: Basis_CleanAlfaNu
  PUBLIC:: Basis_CleanAll
  !~ PUBLIC:: Basis_CpnPrm_Alloc
  !
  !--------------------------------------------- stoichiometry tables --
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: tStoikio
  !-- tStoikio(nCp,nSp)- Stoikio of All Species in Terms of Elements, in Columns
  !-- tStoikio(:,I)- formula of species I, REAL values,
  !-- tStoikio(:,I)- vStoikio(1:nEl)/REAL(div)
  !
  !-- Pr-Primary, As-AqueousSecondary, Ms-MineralNotPrimary
  !-- tAlfXX- Stoikio of All Species in Terms of Components, in Columns
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: & !
  & tAlfSp, & !
  & tAlfPr,tAlfAs,tAlfMs, & !stoichiometric matrices for corresponding species
  & tAlfFs
  !
  !-- tNuSp- stoichio of all species in terms of prim'species
  !-- tNuAs,tNuMs - stoikio of Sec'Species (Aqu.& Min.) in Terms of Prim'Species,
  !-- obtained by extraction from tNuSp
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: & 
  & tNuFas, &
  & tNuSp,  &
  & tNuAs,tNuMs
  !--------------------------------------------/ stoichiometry tables --
  !
  CHARACTER,DIMENSION(:),ALLOCATABLE,PUBLIC:: vSpcCod
  !
  !-- arrays of booleans, tells the species' 'status'
  !-- flags whether component/ species is
  !--   PrimaryInert,AqueousMobile, MineralMobile, etc.
  !-- values will be updated at each basis permutation
  LOGICAL, DIMENSION(:),ALLOCATABLE,PUBLIC:: & 
  & vLCi,vLAx,vLMx,vLAs,vLMs
  !
  !-- permutation arrays for basis changes
  !-- vPrmFw: permutation array from vSpc order to current species order
  INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC:: & !
  !~ & iEleCpn,  & !
  !~ & iCpnEle,  & !
  & vPrmFw,   & ! permutation "forward", from vSpc to vMol,vLnAct,...
  & vPrmBk      ! permutation "bakward"
  !
  INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC:: & !ordering arrays
  & vOrdPr, & ! order of prim'species in vSpc
  & vOrdAq, & ! order of all aqu'species in vSpc
  & vOrdAs, & ! order of aqu'second'species in vSpc
  & vOrdMs    ! order of min'second'species in Spc
  !
  !vOrdPr= order of Prim'species,
  !  -> vSpc(vOrdPr(1:nCp)) gives the Prim'Species array, with Inert'species first
  !vOrdAq= order of Aqu'species,
  !  -> vSpc(vOrdAq(1:nAq)) gives the Aqu'Species aqu, <Prim'Inert><Second'Inert><Prim'Mobile>
  !
  LOGICAL, PUBLIC:: bOH2= .FALSE. !.TRUE. !bOH2=USE OH2 in mole balance on element O
  REAL(dp),PUBLIC:: initMolNu= 1.D-3 !initial mole nrs of species in fluid
  !                                   (-> initial  guess for Newton)
  INTEGER,PUBLIC:: & !
  & nCi,     &  ! nCi= prim'species, "internal"; nPx=prim'species, "mobile" (ReDOx, pH, etc.)
  & nAx,     &  ! nAx= NrOf MOBILE Aqu'Species
  & nMx,     &  ! nMx= Number of Minerals or Gases selected as mobile species
  & nAs,     &  ! nAs= NrOf second'aqu.species, nAs=nAq-nCi-nAx
  & nMs         ! nMs= NrOf "secondary" (i.e. not "mobile") minerals
  !
  !nAq= nCi+nAs+nAx
  !nPx= nAx+nMx (not global)
  !nMn= nMx+nMs 
  !nSp= nAq+nMn=nCi+nAx+nMx+nAs+nMs               
  !
  !species order =
  !blocks= <nCi><nAq-nCi-nAx><nAx><nMx><nMn-nMx>
  !-> all aqu'species on left part, prim'ext'species in middle, then min'int'species
  !-> all intern'species are contiguous, all aqu.species are contiguous,
  !   all mobile species (aqu' + min') are contiguous 
  !-- Species ordering
  !iCi          1:nCi                                = Int'Prim'Aqu'Sp
  !iAs            nCi+1:nCi+nAs                      =      Sec'Aqu'Sp
  !iAx                  nAs+nAx+1:nAq                = eXt'Prim'Aqu'Sp 
  !iMx                            nAq+1:nAq+nMx      = eXt'Prim'Min'Sp
  !iMs                                  nAq+nMx+1:nSp=      Sec'Min'Sp
  !-- Component ordering
  !1:nCi         = Int'Prim'Aqu'Sp
  !nCi+1:nCi+nAx = eXt'Prim'Aqu'Sp 
  !nCi+nAx+1:nCp = eXt'Prim'Min'Sp
  !
  INTEGER,PUBLIC:: & !for aqueous solutions, indexes of components
  & iO_= 0,  &   !index of component assoc'd with element O  
  & iH_= 0,  &   !index of component assoc'd with element H
  & iOx= 0,  &   !index of component assoc'd with OXydation "element"
  & iBal=0       !index of balance component
  !
  INTEGER,PUBLIC:: & !for aqueous solutions, indexes of species
  & isW= 0,  &   !index of solvent species, e.g. H2O
  & isH_=0,  &   !index of H+ species, for pH calculation
  & isOH=0,  &   !index of OH- species, for basis change
  & isO2=0       !index of O2,aq species, used for pE conversion
  !
  REAL(dp),PUBLIC::MWSv !molar weit of solvent
  !
CONTAINS

!~ SUBROUTINE Basis_CpnPrm_Alloc(N)
  !~ INTEGER,INTENT(IN):: N
  !~ IF(ALLOCATED(iEleCpn)) DEALLOCATE(iEleCpn); ALLOCATE(iEleCpn(1:N)) 
  !~ IF(ALLOCATED(iCpnEle)) DEALLOCATE(iCpnEle); ALLOCATE(iCpnEle(1:N)) 
!~ ENDSUBROUTINE Basis_CpnPrm_Alloc

SUBROUTINE Basis_SpcTypes_Alloc(nSp)
  INTEGER,INTENT(IN):: nSp
  !
  IF(nSp==0) CALL Stop_("NO SPECIES ???")
  !
  IF(ALLOCATED(vSpcCod)) DEALLOCATE(vSpcCod); ALLOCATE(vSpcCod(1:nSp)); vSpcCod="Z"
  !
  IF(ALLOCATED(vLCi)) DEALLOCATE(vLCi); ALLOCATE(vLCi(1:nSp)); vLCi=.FALSE.
  IF(ALLOCATED(vLAx)) DEALLOCATE(vLAx); ALLOCATE(vLAx(1:nSp)); vLAx=.FALSE.
  IF(ALLOCATED(vLAs)) DEALLOCATE(vLAs); ALLOCATE(vLAs(1:nSp)); vLAs=.FALSE.
  IF(ALLOCATED(vLMx)) DEALLOCATE(vLMx); ALLOCATE(vLMx(1:nSp)); vLMx=.FALSE.
  IF(ALLOCATED(vLMs)) DEALLOCATE(vLMs); ALLOCATE(vLMs(1:nSp)); vLMs=.FALSE.
  !
ENDSUBROUTINE Basis_SpcTypes_Alloc

SUBROUTINE Basis_SpcOrdPrm_Alloc(nCp_,nSp_,nAq_,nAs_,nMs_)
  INTEGER,INTENT(IN):: nCp_,nSp_,nAq_,nAs_,nMs_
  !
  IF(nCp_==0) CALL Stop_("NO COMPONENTS ???")
  !
  IF(ALLOCATED(vOrdPr)) DEALLOCATE(vOrdPr)  ; ALLOCATE(vOrdPr(1:nCp_))
  !
  IF(ALLOCATED(vOrdAq)) DEALLOCATE(vOrdAq)  ; ALLOCATE(vOrdAq(1:nAq_))
  IF(ALLOCATED(vOrdAs)) DEALLOCATE(vOrdAs)  ; ALLOCATE(vOrdAs(1:nAs_))
  IF(ALLOCATED(vOrdMs)) DEALLOCATE(vOrdMs)  ; ALLOCATE(vOrdMs(1:nMs_))
  !
  IF(ALLOCATED(vPrmFw)) DEALLOCATE(vPrmFw)  ; ALLOCATE(vPrmFw(1:nSp_))
  IF(ALLOCATED(vPrmBk)) DEALLOCATE(vPrmBk)  ; ALLOCATE(vPrmBk(1:nSp_))
  !
ENDSUBROUTINE Basis_SpcOrdPrm_Alloc

SUBROUTINE Basis_AlfaNu_Alloc(nC,nS,nF,nA,nM)
  INTEGER,INTENT(IN):: nC,nS,nF,nA,nM
  !
  CALL Basis_CleanAlfaNu
  !
  ALLOCATE(tAlfSp(1:nC,1:nS))
  ALLOCATE(tAlfPr(1:nC,1:nC))
  ALLOCATE(tAlfAs(1:nC,1:nA))
  ALLOCATE(tAlfMs(1:nC,1:nM))
  !
  ALLOCATE(tNuSp(1:nS,1:nC))
  ALLOCATE(tNuAs(1:nA,1:nC))
  ALLOCATE(tNuMs(1:nM,1:nC))
  !
  ALLOCATE(tNuFas(1:nF,1:nC))
  ALLOCATE(tAlfFs(1:nC,1:nF))
  !
ENDSUBROUTINE Basis_AlfaNu_Alloc

SUBROUTINE Basis_CleanAlfaNu

  IF(ALLOCATED(tAlfSp)) DEALLOCATE(tAlfSp)
  IF(ALLOCATED(tAlfPr)) DEALLOCATE(tAlfPr)
  IF(ALLOCATED(tAlfAs)) DEALLOCATE(tAlfAs)
  IF(ALLOCATED(tAlfMs)) DEALLOCATE(tAlfMs)
  IF(ALLOCATED(tAlfFs)) DEALLOCATE(tAlfFs)
  !
  IF(ALLOCATED(tNuSp))  DEALLOCATE(tNuSp)
  IF(ALLOCATED(tNuAs))  DEALLOCATE(tNuAs) 
  IF(ALLOCATED(tNuMs))  DEALLOCATE(tNuMs)
  IF(ALLOCATED(tNuFas)) DEALLOCATE(tNuFas)
  
ENDSUBROUTINE Basis_CleanAlfaNu

SUBROUTINE Basis_CleanAll

  IF(ALLOCATED(tStoikio)) DEALLOCATE(tStoikio)
  !
  IF(ALLOCATED(vSpcCod)) DEALLOCATE(vSpcCod)
  !
  IF(ALLOCATED(vLCi))   DEALLOCATE(vLCi)
  IF(ALLOCATED(vLAx))   DEALLOCATE(vLAx)
  IF(ALLOCATED(vLAs))   DEALLOCATE(vLAs)
  IF(ALLOCATED(vLMx))   DEALLOCATE(vLMx) !external minerals
  IF(ALLOCATED(vLMs))   DEALLOCATE(vLMs) !internal minerals
  !  
  IF(ALLOCATED(vPrmFw)) DEALLOCATE(vPrmFw)
  IF(ALLOCATED(vPrmBk)) DEALLOCATE(vPrmBk)
  !
  !~ IF(ALLOCATED(iEleCpn)) DEALLOCATE(iEleCpn)
  !~ IF(ALLOCATED(iCpnEle)) DEALLOCATE(iCpnEle)
  !
  IF(ALLOCATED(vOrdPr)) DEALLOCATE(vOrdPr) !
  IF(ALLOCATED(vOrdAq)) DEALLOCATE(vOrdAq) !
  IF(ALLOCATED(vOrdAs)) DEALLOCATE(vOrdAs) !
  IF(ALLOCATED(vOrdMs)) DEALLOCATE(vOrdMs) !
  !
  IF(ALLOCATED(tAlfSp)) DEALLOCATE(tAlfSp)
  IF(ALLOCATED(tAlfPr)) DEALLOCATE(tAlfPr)
  IF(ALLOCATED(tAlfAs)) DEALLOCATE(tAlfAs)
  IF(ALLOCATED(tAlfMs)) DEALLOCATE(tAlfMs)
  !
  IF(ALLOCATED(tNuSp))  DEALLOCATE(tNuSp)
  IF(ALLOCATED(tNuAs))  DEALLOCATE(tNuAs) 
  IF(ALLOCATED(tNuMs))  DEALLOCATE(tNuMs)
  !
  IF(ALLOCATED(tAlfFs)) DEALLOCATE(tAlfFs)
  IF(ALLOCATED(tNuFas))  DEALLOCATE(tNuFas)
  !  
ENDSUBROUTINE Basis_CleanAll

ENDMODULE M_Basis_Vars
