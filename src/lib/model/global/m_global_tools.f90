module M_Global_Tools
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,fHtm,T_,Stop_,Pause_
  !
  implicit none
  !
  private
  !
  public:: Global_Zero
  public:: Global_Clean
  public:: Global_Species_Select
  public:: Global_TP_Update
  public:: Global_Show
  public:: Species_TP_Update
  !
contains

subroutine Global_Show
!--
!--
  use M_IoTools
  use M_T_Phase,only: T_Phase
  use M_Global_Vars, only: vSpc,vFas
  use M_Global_Vars, only: vSolModel,vMixModel
  use M_Global_Vars, only: vSolFas,vMixFas
  ! vDiscretModel,vDiscretParam
  !
  integer:: I,J,K
  integer:: F
  type(T_Phase):: fs
  
  call GetUnit(F)
  open(F,file="global_show.txt")

  write(F,'(A)') "NamFs,iSpc,iMix,iSol"
  do I=1,size(vFas)
  
    fs= vFas(I)
    
    write(F,'(A,3I4,3X)',advance="NO") &
    & fs%NamFs,fs%iSpc,fs%iMix,fs%iSol
    
    if(fs%iSpc >0) write(F,'(A)') vSpc(fs%iSpc)%NamSp
    
    if(fs%iMix >0) then
      do J=1,vMixModel(vMixFas(fs%iMix)%iModel)%nPole
        K= vMixModel(vMixFas(fs%iMix)%iModel)%vIPole(J)
        write(F,'(A,1X)',advance="NO") trim(vSpc(K)%NamSp)
      end do
      write(F,*)
    end if
    
    if(fs%iSol >0) then
      K= vSolModel(vSolFas(fs%iSol)%iModel)%iSolvent
      write(F,'(A,1X)',advance="NO") trim(vSpc(K)%NamSp)
      do J=1,vSolModel(vSolFas(fs%iSol)%iModel)%nSolute
        K= vSolModel(vSolFas(fs%iSol)%iModel)%vISolute(J)
        write(F,'(A,1X)',advance="NO") trim(vSpc(K)%NamSp)
      end do
      write(F,*)
    end if
    
  end do
  
  close(F)
  
  return
end subroutine Global_Show

subroutine Global_TP_Update( &
& TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
& vSpc,vMixModel,vMixFas,vFas) !inout
!--
!-- update (T,P) dependent properties of species, sol'models, phases,
!-- to be called everytime T,P changes
!--
  use M_T_Species,  only: T_Species,T_SpeciesDtb
  use M_T_MixModel, only: T_MixModel
  use M_T_MixPhase, only: T_MixPhase,MixPhase_CalcActivs
  use M_T_Phase,    only: T_Phase,Phase_Calc
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  !
  real(dp),            intent(in):: TdgK,Pbar
  type(T_SpeciesDtb),  intent(in):: vSpcDtb(:)
  type(T_DiscretModel),intent(in):: vDiscretModel(:)
  type(T_DiscretParam),intent(in):: vDiscretParam(:)
  !
  type(T_Species), intent(inout):: vSpc(:)
  type(T_MixModel),intent(inout):: vMixModel(:)
  type(T_MixPhase),intent(inout):: vMixFas(:)
  type(T_Phase),   intent(inout):: vFas(:)
  !
  integer :: I
  !
  call Species_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel)
  !
  !-------------------update activ's of end-members in non-aqu'solutions
  do I=1,size(vMixFas) !
    call MixPhase_CalcActivs( &
    & TdgK,Pbar,                    & !IN
    & vMixModel(vMixFas(I)%iModel), & !IN
    & vMixFas(I))                     !INOUT
  end do
  !
  !----------------------------------------------------update all phases
  do I=1,size(vFas)
    call Phase_Calc( & 
    & TdgK,Pbar,vSpc,vMixModel,vMixFas, & !IN
    & vFas(I))                            !OUT
  end do
  !
end subroutine Global_TP_Update

subroutine Species_TP_Update( &
& TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
& vSpc,vMixModel) !inout
!--
!-- update (T,P) dependent properties of species, sol'models, phases,
!-- to be called everytime T,P changes
!--
  use M_T_DtbH2OHkf,only: DtbH2OHkf_Calc,T_DtbH2OHkf
  use M_T_Species,  only: T_Species,T_SpeciesDtb
  use M_T_MixModel, only: T_MixModel
  use M_T_MixModel, only: MixModel_Param_Update
  use M_Dtb_Calc,   only: Species_TP_Update_fromDtb
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  use M_DiscretModel_Tools
  !
  real(dp),            intent(in):: TdgK,Pbar
  type(T_SpeciesDtb),  intent(in):: vSpcDtb(:)
  type(T_DiscretModel),intent(in):: vDiscretModel(:)
  type(T_DiscretParam),intent(in):: vDiscretParam(:)
  !
  type(T_Species), intent(inout):: vSpc(:)
  type(T_MixModel),intent(inout):: vMixModel(:)
  !
  integer :: I,J
  type(T_DtbH2OHkf):: PropsH2O
  !
  ! compute solvent properties, needed for aqu'species
  if(count(vSpc(:)%Typ=='AQU')>0) &
  & call DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
  
  !----------------------- update T,P dependent param's for pure species
  do J=1,size(vSpc)
    if(vSpc(J)%iDtb>0) &
    call Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,vSpc(J))
    !-> update vSpc(J)%G0rt !-> =G/RT
  end do
  
  !-------------------- update T,P dependent parameters in mixing models
  do I=1,size(vMixModel)
    call MixModel_Param_Update( &
    & TdgK,Pbar, &  !in
    & vMixModel(I)) !out
  end do
  
  !------------------------------ update species built by discretization
  if(size(vDiscretModel)>0) then
    call DiscretSpecies_TP_Update( &
    & vMixModel,     & !IN
    & TdgK,Pbar,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
  end if
  
  return
end subroutine Species_TP_Update

subroutine Global_Zero
  use M_Global_Vars
  
  nAq=0; nMn=0; nGs=0
  !
  call Global_Clean
  !
  allocate(vEle(0))
  !
  allocate(vSpc(0))
  allocate(vSpcDtb(0))
  !
  allocate(vFas(0))
  !
  allocate(vDiscretModel(0))
  allocate(vDiscretParam(0))
  !
  allocate(vKinModel(0))
  allocate(vKinFas(0))
  !
  allocate(vMixModel(0))
  allocate(vMixFas(0))
  !
  allocate(vSolModel(0))
  allocate(vSolFas(0))
  !
  allocate(tFormula(0,0))
  !
  !default solvent= water at 25°C/1bar (1 atm ???)
  !! Solvent%Name=        "WATER"
  !! Solvent%ActModel=    "IDEAL"
  !
end subroutine Global_Zero

subroutine Global_Clean
  use M_Global_Vars
  use M_Dtb_Vars
  !
  if(allocated(vEle))    deallocate(vEle)
  !
  if(allocated(vSpc))    deallocate(vSpc)
  if(allocated(vSpcDtb)) deallocate(vSpcDtb)
  !
  if(allocated(vFas))    deallocate(vFas)
  !
  if(allocated(vMixModel)) deallocate(vMixModel)
  if(allocated(vMixFas))   deallocate(vMixFas)
  !
  if(allocated(vSolModel)) deallocate(vSolModel)
  if(allocated(vSolFas))   deallocate(vSolFas)
  !
  if(allocated(vDiscretModel)) deallocate(vDiscretModel)
  if(allocated(vDiscretParam)) deallocate(vDiscretParam)
  !
  if(allocated(vKinModel)) deallocate(vKinModel)
  if(allocated(vKinFas))   deallocate(vKinFas)
  !
  if(allocated(tFormula))  deallocate(tFormula)
  !
  call Dtb_Vars_Clean
  !
end subroutine Global_Clean

subroutine Global_Species_Select
!--
!--read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks,
!--from the vSpc built from database,
!--retrieve species consistent with vEle & redox & vExclude
!--
  use M_T_Species,  only: T_Species
  use M_Global_Vars,only: vSpc
  !
  logical,allocatable:: vExclude(:)
  logical,allocatable:: vInclude(:)
  integer:: iSp,nSp
  !
  type(T_Species),allocatable:: vSpcNew(:)
  !
  !------------------ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  allocate(vExclude(1:size(vSpc)))  ;  vExclude(:)= .false.
  allocate(vInclude(1:size(vSpc)))  ;  vInclude(:)= .true.
  call Species_Read_Excluded(vSpc,vExclude,vInclude)
  !-----------------/ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  !
  !------------------------------------------------ species selection --
  !-- from the vSpc built from database,
  !-- retrieve species consistent
  !-- with vEle & redox & vExclude
  if(idebug>1) write(fTrc,'(A)') &
  & "< Restrict database to species consistent with element list and redox state"
  !
  allocate(vSpcNew(size(vSpc)))
  !
  nSp=0
  do iSp=1,size(vSpc)
    if( vExclude(iSp) .or. (.not. vInclude(iSp)) ) cycle
    nSp=  nSp+1
    vSpcNew(nSp)= vSpc(iSp)
  end do
  !
  deallocate(vExclude)
  deallocate(vInclude)
  !-----------------------------------------------/ species selection --
  !
  !-------------------------------------------- build new sorted vSpc --
  deallocate(vSpc)
  allocate(vSpc(1:nSp))
  vSpc(1:nSp)= vSpcNew(1:nSp)
  !
  return
end subroutine Global_Species_Select

subroutine Species_Read_Excluded(vSpc,vExclude,vInclude)
!--
!-- process the SPECIES.EXCLUDE block
!-- and the SPECIES.INCLUDE block
!--
  use M_Files,    only: NamFInn
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_T_Species,only: T_Species,Species_Index,Species_Rename
  !
  type(T_Species),intent(in) :: vSpc(:)
  logical,intent(out):: vExclude(:),vInclude(:)
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios,I     
  !
  if(idebug>1) write(fTrc,'(/,A)') "< Species_Read_Excluded"
  !
  vExclude= .false.
  vInclude= .true.
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))
  !
  DoFile: do
  
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    
    call AppendToEnd(L,W,EoL)
    select case(W)
    
    case("ENDINPUT"); exit DoFile
    
    case("SPECIES.EXCLUDE")
      !
      vExclude(:)= .false.
      !
      DoLine1: do
      
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoLine1 !skip comment lines
        call AppendToEnd(L,W,EoL)
        
        select case(W)
          case("ENDINPUT"); exit DoFile
          case("END","ENDSPECIES.EXCLUDE"); exit DoLine1
          !-> reading data from several block..end blocks 
          !case("END","ENDSPECIES.EXCLUDE"); exit DoFile
          !!-> reading only the first available block
        end select
        !
        call Species_Rename(W)
        I= Species_Index(trim(W),vSpc)
        !
        !---------------------------------------------------- trace --
        if(I<1) then
          write(fTrc,'(3A)') "Species ",trim(W)," = is not in current vSpc"
        else
          write(fTrc,'(3A)') "Species ",trim(W)," = is excluded"
          vExclude(I)=.true.
        end if
        !---------------------------------------------------/ trace --
        !
      end do DoLine1
    !endcase("SPECIES.EXCLUDE")
    
    case("SPECIES.INCLUDE")
      !
      vInclude(:)= .false.
      !
      DoLine2: do
      
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoLine2 !skip comment lines
        call AppendToEnd(L,W,EoL)
        
        select case(W)
        case("ENDINPUT")                  ;  exit DoFile
        case("END","ENDSPECIES.INCLUDE")  ;  exit DoLine2
        !-> reading data from several block..end blocks 
        !case("END","ENDSPECIES.EXCLUDE"); exit DoFile
        !!-> reading only the first available block
        end select
        
        call Species_Rename(W)
        I= Species_Index(trim(W),vSpc)
        !
        !---------------------------------------------------- trace --
        if(I<1) then
          write(fTrc,'(3A)') "Species ",trim(W)," = is not in current vSpc"
        else
          write(fTrc,'(3A)') "Species ",trim(W)," = is included"
          vInclude(I)=.true.
        end if
        !---------------------------------------------------/ trace --
        !
      end do DoLine2
    !endcase("SPECIES.INCLUDE")
      
    end select
  end do DoFile
  !
  close(F)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ Species_Read_Excluded"
end subroutine Species_Read_Excluded

end module M_Global_Tools
