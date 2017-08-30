module M_DiscretModel_Read
!--
!-- routine for reading block(s) MIXTURE.DISCRETIZE (was SOLUTION.DISCRETIZE)
!-- and build pure species from discretization of mixture phase(s)
!--
  use M_Kinds
  use M_Trace,         only: iDebug,fTrc,T_,Stop_,Pause_,Warning_
  use M_T_DiscretModel,only: T_DiscretModel
  implicit none
  !
  private
  !
  public:: DiscretModel_Read
  !
  !~ type:: T_LnkDiscretModel
    !~ type(T_DiscretModel):: Value
    !~ type(T_LnkDiscretModel),pointer::Next
  !~ end type T_LnkDiscretModel
  !
contains

subroutine DiscretModel_Read(vMixModel)
!--
!-- called by DiscretPhase_Add in module m_global_alloc
!--
  use M_T_MixModel,only: T_MixModel
  !
  use M_Global_Vars,only: vDiscretModel
  !
  type(T_MixModel),intent(in) :: vMixModel(:)  
  !
  type(T_DiscretModel),allocatable:: vTmp(:)
  integer:: N
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_Read"
  !
  allocate(vTmp(100))
  call DiscretModel_ReadFile(vMixModel,vTmp,N)
  !
  if(N>0) then
    deallocate(vDiscretModel)
    allocate(vDiscretModel(N))
    vDiscretModel(1:N)= vTmp(1:N)
  end if
  !
  deallocate(vTmp)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretModel_Read"
  !! if(iDebug>0) print *,"DiscretModel_Read, N=", N
  !
end subroutine DiscretModel_Read

subroutine DiscretModel_ReadFile( & !
!--
!-- reads the parameters for the "discretization" of a mixture
!--
& vMixModel,     & !in,  database of solution models
& vDiscretModel, & !out
& N)               !OUT
  use M_IOTools
  use M_Files,     only: NamFInn
  use M_T_MixModel,only: T_MixModel
  use M_T_DiscretModel
  !
  type(T_MixModel),    intent(in) :: vMixModel(:)
  type(T_DiscretModel),intent(out):: vDiscretModel(:)
  integer,             intent(out):: N
  !
  type(T_DiscretModel):: DsModl
  type(T_MixModel)    :: MxModl
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios
  integer:: K,iMix,DD,i,DimMix
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_ReadFile"
  !
  N= 0
  !
  if(size(vDiscretModel)<1) return
  !
  iMix= 0
  call GetUnit(F)
  open(f,file=trim(NamFInn))
  !
  DoFile: do
    !
    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    !
    call AppendToEnd(L,W,EoL)
    select case(W)
    !
    case("ENDINPUT"); exit DoFile
    !
    case("SOLUTION.DISCRETIZE","MIXTURE.DISCRETIZE") !for reading one "block"
      !! MIXTURE.DISCRETIZE
      !!   NAME   BIOT
      !!   NUMBER 9
      !!   MODEL  BIOTITE_MG_IDEAL
      !! endSOLUTION.DISCRETIZE
      DoBlock: do
      
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines
        call AppendToEnd(L,W,EoL)
        
        select case(W)
          !
          case("ENDINPUT"); exit DoFile
          !
          case default
            call Warning_("Unknown (or Obsolete) Keyword: "//trim(W))
          !
          case("END","ENDSOLUTION.DISCRETIZE","ENDMIXTURE.DISCRETIZE")
            !
            if(iMix>0) then
              DimMix= vMixModel(DsModl%iMix)%NPole
              !
              if(DimMix>3) then
                print *,"only binary or ternary models !!!"
              else
                N= N+1
                !
                if(DimMix==2) DsModl%Dim2= 0
                if(DimMix==3) DsModl%Dim2= DsModl%Dim1
                call DiscretModel_Dim(DsModl%Dim1,DsModl%Dim2,DsModl%DimTot)
                !
                vDiscretModel(N)= DsModl
                !
                if(iDebug>2) then
                print '(3(A,I4),4A)', &
                & "  Dim1=  ", DsModl%Dim1,   &
                & "  Dim2=  ", DsModl%Dim2,   &
                & "  DimTot=", DsModl%DimTot, &
                & "  MxModl=", trim(vMixModel(DsModl%iMix)%Name), &
                & "  DsModl=", trim(DsModl%Name)
                end if
                !
                if(iDebug>0) write(fTrc,'(3(A,I4),4A)') &
                & "  Dim1=  ", DsModl%Dim1,   &
                & "  Dim2=  ", DsModl%Dim2,   &
                & "  DimTot=", DsModl%DimTot, &
                & "  MxModl=", trim(vMixModel(DsModl%iMix)%Name), &
                & "  DsModl=", trim(DsModl%Name)
              end if
              !
            end if
            exit DoBlock
            !
          case("NAME")
            call LinToWrd(L,W,EoL); DsModl%Name=trim(W)
            !
          case("NUMBER")
            call LinToWrd(L,W,EoL); call WrdToInt(W,DD) 
            DD= MIN(DD,MaxDiscret)
            DD= MAX(DD,MinDiscret)
            !
            DsModl%Dim1=   DD -1
            DsModl%DimTot= DsModl%Dim1
            DsModl%Dim2=   0
            !
            if(.not. Eol) then
              call LinToWrd(L,W,EoL); call WrdToInt(W,DD)
              DD= MIN(DD,MaxDiscret)
              DD= MAX(DD,MinDiscret)
              DsModl%Dim2= DD -1
              DsModl%Dim2= DsModl%Dim1
            end if
            !
          case("MODEL")
            call LinToWrd(L,W,EoL)
            ! find the solution name in vMixModel%Name
            iMix=0
            do K=1,size(vMixModel) 
              if(trim(W)==trim(vMixModel(K)%Name)) iMix=K
            end do
            !
            if(iMix>0) then
              if(vMixModel(iMix)%NPole<2) iMix=0
            end if
            if(iMix==0) call Stop_ &
            & ("MIXTURE.DISCRETIZE: "//trim(W)//"-> Mixing model not found ...")
            !
            DsModl%iMix= iMix
            DsModl%P1= 1 !default values
            DsModl%P2= 2 !id
            DsModl%P3= 0 !id
            !
          case("POLE")
            if(iMix>0) then
              MxModl= vMixModel(iMix)
              do
                call LinToWrd(L,W,EoL)
                I=0
                do K=1,MxModl%NPole
                  if(trim(W)==trim(MxModl%vNamPole(K))) I=K
                end do
                if(EoL) exit
              end do
            end if
            !
        end select
        
      end do DoBlock
      
    end select
    
  end do DoFile
  close(F)
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ DiscretModel_ReadFile"
  !
end subroutine DiscretModel_ReadFile

subroutine DiscretModel_Dim(Dim1,Dim2,N)
  integer,intent(in) :: Dim1,Dim2
  integer,intent(out):: N
  !
  integer:: Dims,i0,i,j,k
  !
  Dims= Dim1 +1
  if(Dim2>0) then  ;  i0= 2  ! ternary
  else             ;  i0= 1  ! binary
  end if
  N=0
  i=0
  do
    i= i+1
    k= 0
    do
      if(Dim2>0) k= k +1
      j= Dims -i -k
      N= N+1
      if(Dim2==0) exit
      if(k==Dims-1 -i) exit
    end do
    if(i==Dims -i0) exit
  end do
  !
  return
end subroutine DiscretModel_Dim

!~ subroutine DiscretModel_Read_(vMixModel)
!~ !.called by DiscretPhase_Add in module m_global_alloc
  !~ use M_T_MixModel,only: T_MixModel
  !~ !
  !~ use M_Global_Vars,only: vDiscretModel
  !~ !
  !~ type(T_MixModel),intent(in) :: vMixModel(:)  
  !~ !
  !~ type(T_LnkDiscretModel),pointer:: Lnk
  !~ integer:: N
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_Read_begin"
  !~ !
  !~ call DiscretModel_BuildLnk(vMixModel,Lnk,N)
  !~ !
  !~ if(N>0) then
    !~ deallocate(vDiscretModel)
    !~ allocate(vDiscretModel(N))
    !~ call DiscretModel_LnkToVec(Lnk,vDiscretModel)
  !~ end if
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A,/)') "< DiscretModel_Read_end"
  !~ !! if(iDebug>0) print *,"DiscretModel_Read, N=", N
  !~ !
!~ end subroutine DiscretModel_Read_

!~ subroutine DiscretModel_BuildLnk( &
!~ !.reads the parameters for the "discretization" of a solid solution 
!~ & vMixModel, & !IN,  database of solution models
!~ & Lnk, & !pointer
!~ & N) !OUT
  !~ use M_IOTools
  !~ use M_Files,     only: NamFInn
  !~ use M_T_MixModel,only: T_MixModel
  !~ use M_T_DiscretModel
  !~ !
  !~ type(T_MixModel),intent(in) :: vMixModel(:)
  !~ type(T_LnkDiscretModel),pointer:: Lnk
  !~ integer,         intent(out):: N
  !~ !
  !~ type(T_LnkDiscretModel),pointer:: pCur
  !~ type(T_DiscretModel):: DiscretModel
  !~ !
  !~ character(len=512):: L,W
  !~ logical:: EoL
  !~ integer:: F,ios
  !~ integer:: K,I,DD
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A)') "< DiscretModel_Read_begin"
  !~ !
  !~ N= 0
  !~ I= 0
  !~ nullify(Lnk)
  !~ !
  !~ call GetUnit(F)
  !~ open(f,file=trim(NamFInn))
  !~ DoFile: do 
    !~ read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    !~ call LinToWrd(L,W,EoL)
    !~ if(W(1:1)=='!') cycle DoFile !skip comment lines
    !~ call AppendToEnd(L,W,EoL)
    !~ select case(W)
    !~ !
    !~ case("ENDINPUT"); exit DoFile
    !~ !
    !~ case("SOLUTION.DISCRETIZE") !for reading one "block"
      !~ !! SOLUTION.DISCRETIZE
      !~ !!   NAME   BIOT
      !~ !!   NUMBER 9
      !~ !!   MODEL  BIOTITE_MG_IDEAL
      !~ !! endSOLUTION.DISCRETIZE
      !~ DoBlock: do
        !~ read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        !~ call LinToWrd(L,W,EoL)
        !~ if(W(1:1)=='!') cycle DoBlock !skip comment lines
        !~ call AppendToEnd(L,W,EoL)
        !~ select case(W)
          !~ case("ENDINPUT"); exit DoFile
          !~ case("END","ENDSOLUTION.DISCRETIZE")
            !~ if(I/=0) then
              !~ N= N+1
              !~ !
              !~ if(iDebug>0) write(fTrc,'(A,2(A,I3))') &
              !~ & DiscretModel%Name, &
              !~ & "Dim1=  ", DiscretModel%Dim1, &
              !~ & "Model= ", DiscretModel%iMix
              !~ !
              !~ if(N==1) then !_________"fill" the linked list
                !~ allocate(Lnk); nullify(Lnk%next)
                !~ Lnk%Value=DiscretModel; pCur=>Lnk
              !~ else
                !~ allocate(pCur%next); nullify(pCur%next%next)
                !~ pCur%next%Value=DiscretModel; pCur=>pCur%next
              !~ end if
              !~ !
            !~ end if
            !~ exit DoBlock
          !~ case("NAME")
            !~ call LinToWrd(L,W,EoL); DiscretModel%Name=trim(W)
          !~ case("NUMBER")
            !~ call LinToWrd(L,W,EoL); call WrdToInt(W,DD) 
            !~ DD= MIN(DD,MaxDiscret)
            !~ DD= MAX(DD,MinDiscret)
            !~ !
            !~ DiscretModel%Dim1=   DD -1
            !~ DiscretModel%DimTot= DiscretModel%Dim1
            !~ DiscretModel%Dim2=   0
            !~ !
            !~ if(.not. Eol) then
              !~ call LinToWrd(L,W,EoL); call WrdToInt(W,DD)
              !~ DiscretModel%Dim2= DiscretModel%Dim1
            !~ end if
            !~ !
          !~ case("MODEL")
            !~ call LinToWrd(L,W,EoL)
            !~ I=0
            !~ do K=1,size(vMixModel) ! find the solution name in vMixModel%Name
              !~ if(trim(W)==trim(vMixModel(K)%Name)) I=K
            !~ end do
            !~ !
            !~ if(vMixModel(I)%NPole<2) I=0
            !~ !
            !~ DiscretModel%iMix= I
            !~ ! if(I==0) &
            !~ ! & call Stop_("SOLUTION.DISCRETIZE: "//trim(W)//"-> Mixing model not found ...")
        !~ end select
      !~ end do DoBlock
    !~ end select
  !~ end do DoFile
  !~ close(F)
  !~ !
  !~ if(iDebug>0) write(fTrc,'(A,/)') "< DiscretModel_BuildLnk_end"
  !~ !
!~ end subroutine DiscretModel_BuildLnk

!~ subroutine DiscretModel_LnkToVec(Lnk,vDiscretModel)
  !~ use M_T_DiscretModel,only: T_DiscretModel
  !~ !
  !~ type(T_LnkDiscretModel), pointer    :: Lnk
  !~ type(T_DiscretModel),    intent(out):: vDiscretModel(:)
  !~ !
  !~ type(T_LnkDiscretModel),pointer:: pCur, pPrev
  !~ integer:: K
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A,/)') "< DiscretModel_LnkToVec_begin"
  !~ if(iDebug>1) print '(A)',"DiscretModel_LnkToVec_begin"
  !~ !
  !~ pCur=>Lnk
  !~ K=0
  !~ do while (associateD(pCur))
    !~ K=K+1
    !~ vDiscretModel(K)=pCur%Value
    !~ pPrev=>pCur
    !~ pCur=> pCur%next
    !~ deallocate(pPrev)
    !~ !
    !~ if(iDebug>0) write(fTrc,'(I3,2A,I3)') &
    !~ & k," DiscretModel(i)%Name=",vDiscretModel(K)%Name,vDiscretModel(K)%iMix
  !~ end do
  !~ !
  !~ if(iDebug>0) write(fTrc,'(/,A,/)') "< DiscretModel_LnkToVec_end"
  !~ !
!~ end subroutine DiscretModel_LnkToVec

end module M_DiscretModel_Read
