module M_Basis_Vars

  use M_Kinds 
  use M_Trace,only: Stop_
  
  implicit none
  
  private
  !
  public:: Basis_SpcOrdPrm_Alloc
  public:: Basis_SpcTypes_Alloc
  public:: Basis_AlfaNu_Alloc
  public:: Basis_CleanAlfaNu
  public:: Basis_CleanAll
  !
  !--------------------------------------------- stoichiometry tables --
  real(dp),dimension(:,:),allocatable,public:: tStoikio
  !-- tStoikio(nCp,nSp)- Stoikio of All Species in Terms of Elements, in Columns
  !-- tStoikio(:,I)- formula of species I, real values,
  !-- tStoikio(:,I)- vStoikio(1:nEl)/real(div)
  !
  !-- Pr-Primary, As-AqueousSecondary, Ms-MineralNotPrimary
  !-- tAlfXX- Stoikio of All Species in Terms of Components, in Columns
  real(dp),dimension(:,:),allocatable,public:: & !
  & tAlfSp, & !
  & tAlfPr,tAlfAs,tAlfMs, & !stoichiometric matrices for corresponding species
  & tAlfFs
  !
  !-- tNuSp- stoichio of all species in terms of prim'species
  !-- tNuAs,tNuMs - stoikio of Sec'Species (Aqu.& Min.) in Terms of Prim'Species,
  !-- obtained by extraction from tNuSp
  real(dp),dimension(:,:),allocatable,public:: & 
  & tNuFas, &
  & tNuSp,  &
  & tNuAs,tNuMs
  !--------------------------------------------/ stoichiometry tables --
  !
  character,dimension(:),allocatable,public:: vSpcCod
  !
  !-- arrays of booleans, tells the species' 'status'
  !-- flags whether component/ species is
  !-- PrimaryInert,AqueousMobile, MineralMobile, etc.
  !-- values will be updated at each basis permutation
  logical, dimension(:),allocatable,public:: & 
  & vLCi,vLAx,vLMx,vLAs,vLMs
  !
  !-- permutation arrays for basis changes
  !-- vPrmFw: permutation array from vSpc order to current species order
  integer,dimension(:),allocatable,public:: & !
  !  & iEleCpn,  & !
  !  & iCpnEle,  & !
  & vPrmFw,   & ! permutation "forward", from vSpc to vMol,vLnAct,...
  & vPrmBk      ! permutation "bakward"
  !
  integer,dimension(:),allocatable,public:: & !ordering arrays
  & vOrdPr, & ! order of prim'species in vSpc
  & vOrdAq, & ! order of all aqu'species in vSpc
  & vOrdAs, & ! order of aqu'second'species in vSpc
  & vOrdMs    ! order of min'second'species in Spc
  !
  !vOrdPr= order of Prim'species,
  !  -> vSpc(vOrdPr(1:nCp)) gives the Prim'Species array, 
  !     with Inert'species first
  !vOrdAq= order of Aqu'species,
  !  -> vSpc(vOrdAq(1:nAq)) gives the Aqu'Species aqu,
  !     <Prim'Inert><Second'Inert><Prim'Mobile>
  !
  logical, public:: bOH2= .false. !.true. 
  !bOH2=use OH2 in mole balance on element O
  real(dp),public:: initMolNu= 1.D-3 !initial mole nrs of species in fluid
  !                                   (-> initial  guess for Newton)
  integer,public:: & !
  & nCi,     &  ! nCi= prim'species, "internal"; 
  !             ! nPx=prim'species, "mobile" (Redox, pH, etc.)
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
  !-> all aqu'species on left part, prim'ext'species in middle, 
  !   then min'int'species
  !-> all intern'species are contiguous, all aqu.species are contiguous,
  !   all mobile species (aqu' + min') are contiguous 
  !-- Species ordering
  !iCi          1:nCi                                 = Int'Prim'Aqu'Sp
  !iAs            nCi+1:nCi+nAs                       =      Sec'Aqu'Sp
  !iAx                  nAs+nAx+1:nAq                 = eXt'Prim'Aqu'Sp 
  !iMx                            nAq+1:nAq+nMx       = eXt'Prim'Min'Sp
  !iMs                                  nAq+nMx+1:nSp =      Sec'Min'Sp
  !-- Component ordering
  !1:nCi         = Int'Prim'Aqu'Sp
  !nCi+1:nCi+nAx = eXt'Prim'Aqu'Sp 
  !nCi+nAx+1:nCp = eXt'Prim'Min'Sp
  !
  integer,public:: & !for aqueous solutions, indexes of components
  & iO_= 0,  &   !index of component assoc'd with element O  
  & iH_= 0,  &   !index of component assoc'd with element H
  & iOx= 0,  &   !index of component assoc'd with OXydation "element"
  & iBal=0       !index of balance component
  !
  integer,public:: & !for aqueous solutions, indexes of species
  & isW= 0,  &   !index of solvent species, e.g. H2O
  & isH_=0,  &   !index of H+ species, for pH calculation
  & isOH=0,  &   !index of OH- species, for basis change
  & isO2=0       !index of O2,aq species, used for pE conversion
  !
  real(dp),public::MWSv !molar weit of solvent
  !
contains

! subroutine Basis_CpnPrm_Alloc(N)
!   integer,intent(in):: N
!   if(allocated(iEleCpn)) deallocate(iEleCpn); allocate(iEleCpn(1:N)) 
!   if(allocated(iCpnEle)) deallocate(iCpnEle); allocate(iCpnEle(1:N)) 
! end subroutine Basis_CpnPrm_Alloc

subroutine Basis_SpcTypes_Alloc(nSp)
  integer,intent(in):: nSp
  !
  if(nSp==0) call Stop_("NO SPECIES ???")
  !
  if(allocated(vSpcCod)) deallocate(vSpcCod); allocate(vSpcCod(1:nSp)); vSpcCod="Z"
  !
  if(allocated(vLCi)) deallocate(vLCi); allocate(vLCi(1:nSp)); vLCi=.false.
  if(allocated(vLAx)) deallocate(vLAx); allocate(vLAx(1:nSp)); vLAx=.false.
  if(allocated(vLAs)) deallocate(vLAs); allocate(vLAs(1:nSp)); vLAs=.false.
  if(allocated(vLMx)) deallocate(vLMx); allocate(vLMx(1:nSp)); vLMx=.false.
  if(allocated(vLMs)) deallocate(vLMs); allocate(vLMs(1:nSp)); vLMs=.false.
  !
  return
end subroutine Basis_SpcTypes_Alloc

subroutine Basis_SpcOrdPrm_Alloc(nCp_,nSp_,nAq_,nAs_,nMs_)
  integer,intent(in):: nCp_,nSp_,nAq_,nAs_,nMs_
  !
  if(nCp_==0) call Stop_("NO COMPONENTS ???")
  !
  if(allocated(vOrdPr)) deallocate(vOrdPr)  ; allocate(vOrdPr(1:nCp_))
  !
  if(allocated(vOrdAq)) deallocate(vOrdAq)  ; allocate(vOrdAq(1:nAq_))
  if(allocated(vOrdAs)) deallocate(vOrdAs)  ; allocate(vOrdAs(1:nAs_))
  if(allocated(vOrdMs)) deallocate(vOrdMs)  ; allocate(vOrdMs(1:nMs_))
  !
  if(allocated(vPrmFw)) deallocate(vPrmFw)  ; allocate(vPrmFw(1:nSp_))
  if(allocated(vPrmBk)) deallocate(vPrmBk)  ; allocate(vPrmBk(1:nSp_))
  !
  return
end subroutine Basis_SpcOrdPrm_Alloc

subroutine Basis_AlfaNu_Alloc(nC,nS,nF,nA,nM)
  integer,intent(in):: nC,nS,nF,nA,nM
  !
  call Basis_CleanAlfaNu
  !
  allocate(tAlfSp(1:nC,1:nS))
  allocate(tAlfPr(1:nC,1:nC))
  allocate(tAlfAs(1:nC,1:nA))
  allocate(tAlfMs(1:nC,1:nM))
  !
  allocate(tNuSp(1:nS,1:nC))
  allocate(tNuAs(1:nA,1:nC))
  allocate(tNuMs(1:nM,1:nC))
  !
  allocate(tNuFas(1:nF,1:nC))
  allocate(tAlfFs(1:nC,1:nF))
  !
  return
end subroutine Basis_AlfaNu_Alloc

subroutine Basis_CleanAlfaNu

  if(allocated(tAlfSp)) deallocate(tAlfSp)
  if(allocated(tAlfPr)) deallocate(tAlfPr)
  if(allocated(tAlfAs)) deallocate(tAlfAs)
  if(allocated(tAlfMs)) deallocate(tAlfMs)
  if(allocated(tAlfFs)) deallocate(tAlfFs)
  !
  if(allocated(tNuSp))  deallocate(tNuSp)
  if(allocated(tNuAs))  deallocate(tNuAs) 
  if(allocated(tNuMs))  deallocate(tNuMs)
  if(allocated(tNuFas)) deallocate(tNuFas)
  
  return
end subroutine Basis_CleanAlfaNu

subroutine Basis_CleanAll
  !
  if(allocated(tStoikio)) deallocate(tStoikio)
  !
  if(allocated(vSpcCod)) deallocate(vSpcCod)
  !
  if(allocated(vLCi))   deallocate(vLCi)
  if(allocated(vLAx))   deallocate(vLAx)
  if(allocated(vLAs))   deallocate(vLAs)
  if(allocated(vLMx))   deallocate(vLMx) !external minerals
  if(allocated(vLMs))   deallocate(vLMs) !internal minerals
  !  
  if(allocated(vPrmFw)) deallocate(vPrmFw)
  if(allocated(vPrmBk)) deallocate(vPrmBk)
  !
  !~ if(allocated(iEleCpn)) deallocate(iEleCpn)
  !~ if(allocated(iCpnEle)) deallocate(iCpnEle)
  !
  if(allocated(vOrdPr)) deallocate(vOrdPr) !
  if(allocated(vOrdAq)) deallocate(vOrdAq) !
  if(allocated(vOrdAs)) deallocate(vOrdAs) !
  if(allocated(vOrdMs)) deallocate(vOrdMs) !
  !
  if(allocated(tAlfSp)) deallocate(tAlfSp)
  if(allocated(tAlfPr)) deallocate(tAlfPr)
  if(allocated(tAlfAs)) deallocate(tAlfAs)
  if(allocated(tAlfMs)) deallocate(tAlfMs)
  !
  if(allocated(tNuSp))  deallocate(tNuSp)
  if(allocated(tNuAs))  deallocate(tNuAs) 
  if(allocated(tNuMs))  deallocate(tNuMs)
  !
  if(allocated(tAlfFs)) deallocate(tAlfFs)
  if(allocated(tNuFas))  deallocate(tNuFas)
  !  
  return
end subroutine Basis_CleanAll

end module M_Basis_Vars
