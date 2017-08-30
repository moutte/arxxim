module M_KinModel_Read
!--
!-- routines for reading kinetic models
!--
  use M_Kinds
  use M_Trace,     only: Stop_,iDebug,fTrc,T_, Fatal_, Warning_,Pause_
  use M_T_Kinmodel,only: T_KinModel
  implicit none
  !
  private
  !
  public:: KinModel_Init
  !
  !public:: KinModel_FileToLnk
  !public:: KinModel_Alloc
  !
  type T_LnkKinModel
    type(T_KinModel)   :: Value
    type(T_LnkKinModel),pointer:: Next
  end type T_LnkKinModel
  type(T_LnkKinModel),pointer:: Lnk
  !
  !type(T_KinModel),dimension(:),allocatable:: vKinModel
  !
contains

subroutine KinModel_Init(vSpc)
  use M_T_Species,  only: T_Species
  use M_Global_Vars,only: vKinModel
  
  type(T_Species),dimension(:),intent(in):: vSpc
  
  type(T_LnkKinModel),pointer:: Lnk
  integer::N, i
  
  call KinModel_FileToLnk(vSpc,N,Lnk) !-> available kinetic data on minerals
  !
  if(N>0) then
    if(allocated(vKinModel)) deallocate(vKinModel)
    allocate(vKinModel(1:N))
    call KinModel_Alloc(Lnk,vKinModel)
  else
    if(iDebug>0) write(fTrc,'(A)') "NO Kinetic Data Found !!!!!!!!!!!!!!"
    if(iDebug>2) print *,"NO Kinetic Database Found !!!"
  end if
  !
  if(iDebug==4) then
    print *,"KinModel_Init >>"
    do i=1,size(vKinModel)
      print *,vKinModel(i)%Name
    enddo
    call pause_
  end if
  
  return
end subroutine KinModel_Init

subroutine BuildLnk(B,E,Lnk,pCur)
  logical            :: B !TRUE=first element
  type(T_KinModel)   :: E
  type(T_LnkKinModel),pointer:: Lnk,pCur
  
  if(B) nullify(Lnk)
  
  if(B) then
    allocate(Lnk)
    nullify(Lnk%next)
    Lnk%Value=E
    pCur => Lnk
  else
    allocate(pCur%next)
    nullify(pCur%next%next)
    pCur%next%Value=E
    pCur => pCur%next
  end if
  
end subroutine BuildLnk

subroutine KinModel_FileToLnk( &
& vSpc, & !IN
& N,Lnk)  !OUT
!--
!-- read kinetic data for minerals (and any other species)
!-- -> build LnkKin
!--
  use M_IOTools !, only:dimV,LinToWrd,ReadRValsV
  use M_Files,   only: NamFKin,DirLog
  use M_Global_Vars,  only: nAq !,vSpc,nSp
  use M_T_Species, only: T_Species,Species_Index,Species_Rename
  use M_T_Kinmodel,only: MaxKinTerm,T_KinModel
  !
  type(T_Species),dimension(:),intent(in):: vSpc
  integer,            intent(out):: N
  type(T_LnkKinModel),pointer    :: Lnk
  !
  type(T_LnkKinModel),pointer::p
  character(len=512):: L,W,W1
  type(T_KinModel)  ::MK
  logical :: EoL,NewM,OldFormat
  logical :: LPrecip,LDissol,ModelIsOk
  integer :: iAq,J
  integer :: f,ios
  real(dp):: X1,X2,X3
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< KinModel_FileToLnk"
  !
  call GetUnit(f)
  open(f,file=trim(NamFKin))
  !
  OldFormat=  .false.
  N= 0
  !
  DoFile: do 
    
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    
    S1: select case(W)
    
    case("ENDINPUT") S1
      exit  DoFile  
    
    case("KINETICS") S1
    !format for kinetic models, a la USGS, new version with possibly two species
      !
      if(.not. EoL) then
        call LinToWrd(L,W1,EoL)
        OldFormat= trim(W1)=="OLD"
      end if
      !
      if(OldFormat) then
      !
      N=0; NewM=.false.
      !
      DoBlockOld: do !build LnkKin of all "kinetic" minerals
      
        read(F,'(A)',iostat=ios) L
        if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlockOld !skip comment lines
        call AppendToEnd(L,W,EoL)
        
        S2_: select case(W)
          case("ENDINPUT") S2_; exit DoFile !to prevent bad file format 
          case("END","ENDKINETICS") S2_; exit DoBlockOld
        end select S2_
        
        if(W(1:1)/='&') then
          !-- if first char is not '&', the line contain a name
          if(NewM) then
            call BuildLnk(N==1,MK,Lnk,p)
            NewM=.false.
            !!if(iDebug>0) write(fTrc,'(A24,I3,A24)') "Min=',M%Name,I," =Spc ",Species_Index(M%Name)
          end if
          !
          !!MK%Special=0
          MK%NTermD=0; MK%pK_d=Zero;  MK%E_d=Zero; MK%N_d=Zero; MK%AlfaD=One; MK%BetaD=One
          MK%NTermP=0; MK%pK_p=Zero;  MK%E_p=Zero; MK%N_p=Zero; MK%AlfaP=One; MK%BetaP=One
          MK%ISpcDiss(:)=0; MK%JSpcDiss(:)=0 
          MK%ISpcPrec(:)=0; MK%JSpcPrec(:)=0 
          !
          MK%Name=trim(W)
          if(.not. EoL) call LinToWrd(L,W,EoL)
          J=Species_Index(trim(W),vSpc)
          !check whether the mineral is also in the thermodynamic dtb
          if(J>0) then !if mineral is Ok, then set NewM true -> will be saved in LnkKin
            N=N+1
            NewM=.true.
          end if
          
        else
          !-- W(1:1)--'&'
          !-- if first char is '&',
          !-- the line is the continuation of preceding, i.e. contains data
          if(NewM) then
            !-- read first word after &, should be either DISSOL or PRECIP
            call LinToWrd(L,W,Eol)
            !
            select case(trim(W))
              case("PRECIP")
                LPrecip=.true.; LDissol=.false.
              case("DISSOL")
                LDissol=.true.; LPrecip=.false.
              case default
                call Stop_(trim(W)//"should have either PRECIP or DISSOL !!!...")
            end select
            !
            call LinToWrd(L,W,Eol) !-> W will contain species name, or "QSK"
            !
            if(trim(W)=="QSK") then
              !--- read Alpha, Beta
              call LinToWrd(L,W,Eol); call WrdToReal(W,X1)
              if(LPrecip) MK%AlfaP=X1
              if(LDissol) MK%AlfaD=X1 
              
              call LinToWrd(L,W,Eol); call WrdToReal(W,X1)
              if(LPrecip) MK%BetaP=X1
              if(LDissol) MK%BetaD=X1 
            
            else
              !--- read activator parameters
              call Species_Rename(W)
              iAq= Species_Index(trim(W),vSpc)
              if (iAq==0) then
                if(iDebug>0) &
                & write(fTrc,'(2A)') &
                & "in KinModel "//trim(MK%Name),", unknown species : "//trim(W)
              end if
              !
              call LinToWrd(L,W,Eol); call WrdToReal(W,X1) !-> read kinetic parameters
              call LinToWrd(L,W,Eol); call WrdToReal(W,X2) !-> read kinetic parameters
              call LinToWrd(L,W,Eol); call WrdToReal(W,X3) !-> read kinetic parameters
              !
              if(iAq>0 .and. iAq<=nAq) then
                !! if(iDebug>0) write(fTrc,'(A,I3,2A)') "iAq=",iAq," <- ",trim(W)
                if(LPrecip .and. MK%NTermP<MaxKinTerm) then
                !activation energies in kiloJoule in input file !!!!!!!!!!!!!
                  MK%NTermP= MK%NTermP+1
                  MK%ISpcPrec(MK%NTermP)= iAq !-> variable, updated in KinModel_UpDate
                  MK%pK_P(MK%NTermP)=     X1
                  MK%E_P(MK%NTermP)=      X2*1.D3
                  MK%N_P(MK%NTermP)=      X3
                end if
                if(LDissol .and. MK%NTermD<MaxKinTerm) then
                  MK%NTermD= MK%NTermD+1
                  MK%ISpcDiss(MK%NTermD)= iAq !-> variable, updated in KinModel_UpDate
                  MK%pK_D(MK%NTermD)=     X1
                  MK%E_D(MK%NTermD)=      X2*1.D3
                  MK%N_D(MK%NTermD)=      X3
                end if
              end if
              
            end if
            !! M%KModel=MK
          end if !if NEWM
        end if
      enddo DoBlockOld
      
      if(NewM) call BuildLnk(N==1,MK,Lnk,p) !save last mineral in list
      !
      else
      !
      N=0; NewM=.false.
      !
      DoBlock: do !build LnkKin of all "kinetic" minerals
        !
        read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoBlock !skip comment lines
        call AppendToEnd(L,W,EoL)
        !
        S2: select case(W)
          case("ENDINPUT") S2; exit DoFile !to prevent bad file format 
          case("END","ENDKINETICS") S2; exit DoBlock
        end select S2
        !
        if(W(1:1)/='&') then
          !-- if first char is not '&', the line contains the model name
          !
          !-- first save model from current buffer to list
          if(NewM .and. ModelIsOk) then
            N=N+1
            call BuildLnk(N==1,MK,Lnk,p)
            NewM=.false.
            !!if(iDebug>0) write(fTrc,'(A24,I3,A24)') "Min=',M%Name,I," =Spc ",Species_Index(M%Name)
          end if
          !
          MK%Name=trim(W)
          !
          MK%NTermD=0; MK%pK_d=Zero;  MK%E_d=Zero; MK%N_d=Zero; MK%AlfaD=One; MK%BetaD=One
          MK%NTermP=0; MK%pK_p=Zero;  MK%E_p=Zero; MK%N_p=Zero; MK%AlfaP=One; MK%BetaP=One
          MK%ISpcDiss(:)=0; MK%JSpcDiss(:)=0 
          MK%ISpcPrec(:)=0; MK%JSpcPrec(:)=0 
          !
          NewM=.true.
          ModelIsOk= .true.
          !
          !QUARTZ USGS2004-1068
          !&  DISSOL  1     1     QsK      
          !&  DISSOL  13.4  90.9  H2O  0    
          !&  PRECIP  1     1     QsK      
          !&  PRECIP  13.4  90.9  H2O  0    
        else
          !-- W(1:1)--'&'
          !-- if first char is '&',
          !-- the line is the continuation of preceding, i.e. contains data
          if(NewM) then
            call LinToWrd(L,W,Eol) !read first word after &, should be either DISSOL or PRECIP
            select case(trim(W))
              case("PRECIP"); LPrecip=.true.; LDissol=.false.
              case("DISSOL"); LDissol=.true.; LPrecip=.false.
              case default
                call Stop_(trim(W)//"should have either PRECIP or DISSOL !!!...")
            end select
            !
            call LinToWrd(L,W,Eol); call WrdToReal(W,X1)
            call LinToWrd(L,W,Eol); call WrdToReal(W,X2)
            !
            call LinToWrd(L,W,Eol) !-> W will contain a species name, or "QSK"
            !
            if(trim(W)=="QSK") then
              !-- read Alpha, Beta
              if(LPrecip) MK%AlfaP=X1
              if(LDissol) MK%AlfaD=X1
              if(LPrecip) MK%BetaP=X2
              if(LDissol) MK%BetaD=X2
              !
            else
              !-- read activator parameters
              call Species_Rename(W)
              iAq=Species_Index(trim(W),vSpc)
              !
              if (iAq==0) then
                if(iDebug>0) &
                & write(fTrc,'(3A)') &
                & "in KinModel "//trim(MK%Name), &
                & ", unknown species : "//trim(W), &
                & ", >> KINETIC MODEL NOT includeD"
                ModelIsOk= .false.
              end if
              !
              call LinToWrd(L,W,Eol)  ;  call WrdToReal(W,X3)
              !
              if(iAq>0 .and. iAq<=nAq) then
                if(iDebug>0) write(fTrc,'(A,I3,2A)') "iAq=",iAq," <- ",vSpc(iAq)%NamSp
                !
                if(LPrecip .and. MK%NTermP<MaxKinTerm) then
                  !
                  MK%NTermP= MK%NTermP +1
                  MK%ISpcPrec(MK%NTermP)=iAq
                  !
                  MK%pK_P(MK%NTermP)= X1
                  MK%E_P(MK%NTermP)=  X2*1.D3 !activation energies in kiloJoule in input file
                  MK%N_P(MK%NTermP)=  X3
                  !
                  MK%JSpcPrec(MK%NTermP)=0
                  !
                  if(.not.EoL) then !case of second species
                    call LinToWrd(L,W,Eol)
                    call Species_Rename(W)
                    iAq=Species_Index(trim(W),vSpc)
                    if(iAq>0 .and. iAq<=nAq) then 
                      MK%JSpcPrec(MK%NTermP)=iAq
                      call LinToWrd(L,W,Eol); call WrdToReal(W,X1)
                      MK%NJ_P(MK%NTermP)=X1
                    else
                      if(iDebug>0) &
                      & write(fTrc,'(3A)') &
                      & "in KinModel "//trim(MK%Name), &
                      & ", unknown species : "//trim(W), &
                      & ", >> KINETIC MODEL NOT includeD"
                      ModelIsOk= .false.
                    end if
                  end if
                  !
                end if
                !
                if(LDissol .and. MK%NTermD<MaxKinTerm) then
                  !
                  MK%NTermD= MK%NTermD +1
                  MK%ISpcDiss(MK%NTermD)=iAq
                  !
                  MK%pK_D(MK%NTermD)= X1
                  MK%E_D(MK%NTermD)=  X2*1.D3
                  !!! activation energies in kiloJoule in input file !!!
                  MK%N_D(MK%NTermD)=  X3
                  !
                  MK%JSpcDiss(MK%NTermD)=0
                  !
                  if(.not.EoL) then !case of second species
                    call LinToWrd(L,W,Eol)
                    iAq=Species_Index(trim(W),vSpc)
                    if(iAq>0 .and. iAq<=nAq) then 
                      MK%JSpcDiss(MK%NTermD)=iAq
                      call LinToWrd(L,W,Eol); call WrdToReal(W,X1)
                      MK%NJ_D(MK%NTermD)=X1
                    end if
                  end if
                end if
                !
              end if
              
            end if
          end if !if NEWM
        end if !if W(1:1)=='&'
      enddo DoBlock
      if(NewM .and. ModelIsOk) then
        N= N+1
        call BuildLnk(N==1,MK,Lnk,p) !save last mineral in list
      end if
      !
      end if
    !endcase("KINETICS")  
    end select S1
     
  enddo DoFile
  !
  close(f)
  !
  if(iDebug>0) write(fTrc,'(A,I3)') "N=", N
  !
  if(iDebug>0) write(fTrc,'(A,/)') "</ KinModel_FileToLnk"
end subroutine KinModel_FileToLnk

subroutine KinModel_Alloc(Lnk,vKinModel)
  !
  type(T_LnkKinModel),pointer:: Lnk
  type(T_KinModel),dimension(:),intent(out):: vKinModel
  !
  type(T_LnkKinModel),pointer::pCur, pPrev
  integer::I
  !
  if(iDebug>0) write(fTrc,'(/,A)') "< KinModel_Alloc"
  I=0
  pCur=>Lnk
  do while (associateD(pCur))
    I= I+1
    vKinModel(I)=pCur%Value 
    if(iDebug>0) write(fTrc,'(I4,A1,A12)') I," ",vKinModel(I)%Name
    pPrev=>pCur; pCur=> pCur%next; deallocate(pPrev)
  enddo
  if(iDebug>0) write(fTrc,'(A,/)') "</ KinModel_Alloc"
end subroutine KinModel_Alloc

end module M_KinModel_Read


