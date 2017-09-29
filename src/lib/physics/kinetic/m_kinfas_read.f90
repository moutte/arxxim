module M_KinFas_Read
!--
!-- routines for reading kinetic models
!--
  use M_Kinds
  use M_Trace,   only: iDebug,fTrc,T_,Stop_,Warning_,Pause_
  use M_T_KinFas,only: T_KinFas
  use M_T_Phase, only: T_Phase
  implicit none
  !
  private
  !
  public:: KinFas_BuildLnk
  public:: KinFas_LnkToVec
  public:: T_LnkKin 
  !
  type T_LnkKin
    type(T_KinFas)         ::Value
    type(T_LnkKin), pointer::Next
  end type T_LnkKin
  !
contains

subroutine KinFas_LnkToVec(LnkKin,vFas,vKinFas)
  !---------------------------------------------------------------------
  type(T_LnkKin),pointer    :: LnkKin
  type(T_Phase), intent(in) :: vFas(:)
  type(T_KinFas),intent(out):: vKinFas(:)
  !---------------------------------------------------------------------
  type(T_LnkKin),pointer::pCur, pPrev
  integer:: I, J
  !
  if(idebug>1) write(fTrc,'(/,A)') "< KinFas_LnkToVec"
  !
  I=0
  J=0
  pCur=> LnkKin
  do while (associated(pCur))
    !
    I= I+1
    vKinFas(I)= pCur%Value
    !
    if(idebug>1) write(fTrc,'(I4,A1,A12)') I," ",vKinFas(I)%NamKF
    !
    pPrev=>pCur
    pCur=> pCur%next
    deallocate(pPrev)
    !
  enddo
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ KinFas_LnkToVec"
  !
end subroutine KinFas_LnkToVec

subroutine BuildLnkKin(B,E,L,P)
  logical               ::B  !.true.=first element
  type(T_KinFas)        ::E
  type(T_LnkKin),pointer::L,P
  if(B) nullify(L)
  if(B) then
    allocate(L);      nullify(L%next);      L%Value=E;       P=>L
  else
    allocate(P%next); nullify(P%next%next); P%next%Value=E;  P=>P%next
  end if
end subroutine BuildLnkKin

subroutine KinFas_BuildLnk(vFas,vKinModel,sModelSurf,N,LnkKin)
!--
!-- read the "ROCK" block from input.file -> build link list --
!--
  use M_IOTools
  use M_Files,      only: NamFInn
  !!use M_T_Species, only: T_Species,Species_Index
  use M_T_Phase,    only: T_Phase,Phase_Index
  use M_T_Kinmodel, only: T_KinModel,KinModel_Index
  use M_T_KinFas,   only: KinFas_Zero,KinFas_Surf_Zero
  !
  type(T_Phase),   intent(in) :: vFas(:)
  type(T_KinModel),intent(in) :: vKinModel(:)
  character(len=7),intent(in) :: sModelSurf
  integer,         intent(out):: N
  type(T_LnkKin),  pointer    :: LnkKin
  !
  real(dp),dimension(dimV)::vX
  type(T_LnkKin),pointer::pKin
  type(T_KinFas)    :: M
  character(len=255):: L
  character(len=80) :: W1,W2,W3
  real(dp):: xGram,xMole
  logical :: EoL,BlockFound,OldFormat
  integer :: I,jKin,mDum,f_,ios
  !
  !dimV: defined in ModIo,
  !dimension of "buffer array" scanned from words of a string
  !
  call GetUnit(f_)
  open(f_,file=trim(NamFInn))
  !
  if(iDebug>2) then
    do i=1,size(vFas)
      print *,"NamFs= ",trim(vFas(i)%NamFs)
    end do
  end if

  ! call pause_
  N= 0
  !
  BlockFound= .false.
  OldFormat=  .false.
  if(idebug>1) write(fTrc,'(/,A)') "< KinFas_BuildLnk"
  DoFile: do 
    !
    read(F_,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W1,EoL)
    if(W1(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W1,EoL)
    !
    select case(W1)
    !
    case("ENDINPUT"); exit DoFile
    !
    !NB: the DYNAMIC.ROCK block is read after KINETICS block
    !
    case("DYNAMIC.ROCK")
      BlockFound=.true.
      N=0
      if(.not. EoL) then
        call LinToWrd(L,W3,EoL)
        OldFormat= trim(W3)=="OLD"
      end if
      if(OldFormat) then
      !------------------------------------------------------ old format
        DoReadOld: do
        
          read(F_,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W1,EoL)
          if(W1(1:1)=='!')   cycle DoReadOld !skip comment lines
          call AppendToEnd(L,W1,EoL)
          !
          select case(W1)
          case("ENDINPUT")                ; exit DoFile
          case("END","ENDDYNAMIC.ROCK")   ; exit DoReadOld
          end select
          !
          call LinToWrd(L,W2,EoL)
          !
          if(trim(W2)/="D" .and. trim(W2)/="R") then
            !-> if next word is not D or R, then it is Species Name
            
            !I= Species_Index(W2,vSpc) !find species named W2 in vSpc
            I= Phase_Index(W2,vFas) !find phase named W2 in vFas
            
            call LinToWrd(L,W2,EoL) !-> next word is Name Of Kinetic Model
            jKin= KinModel_Index(W2,vKinModel) !find species named W2 in vKinModel
            
            call LinToWrd(L,W2,EoL) !read next word
            if(trim(W2)/="D" .and. trim(W2)/="R") call Stop_( &
            & "texture code should be either R (radius) or D (density)" )
            M%cMode=W2(1:1)
          
          else
            I= Phase_Index(W1,vFas)
            !I= Species_Index(W1,vSpc)
            jKin= KinModel_Index(W1,vKinModel)
            M%cMode=W2(1:1)
          
          end if
          !
          if(I==0) write(*,'(A)') "Species not found for "//trim(W1)
          if(jKin==0) write(*,'(A)') "Kinetic model not found for "//trim(W1)
          !
          if(I>0 .and. jKin>0) then
            !
            call KinFas_Zero(M)
            !
            M%NamKF=trim(W1) !-> the Kinetic Mineral Name
            !M%iSpc=I
            M%iFas=I
            M%iKin=jKin
            !
            if(idebug>1) write(fTrc,'(3(A12,1X))') &
            !& M%Name,vSpc(M%iSpc)%Name,vKinModel(M%iKin)%Name
            & M%NamKF,vFas(M%iFas)%NamFs,vKinModel(M%iKin)%Name
            !
            call ReadRValsV(L,mDum,vX) !read radius and relative fraction
            M%Dat%Radius= vX(1)
            M%Dat%PhiM=   vX(2)
            if(vX(3)/=Zero) M%QsKSeuil= vX(3)
            if(vX(4)/=Zero) M%ReacCoef= vX(4)
            !
            !!M%Dat%cSat="1" !"PRIMARY" 
            !!if(M%Dat%PhiM<1.0E-5) M%Dat%cSat="2" !"SECONDA" 
            !!!-> mineral is "secondary",i.e. not really in initial assemblage
            !
            N=N+1; call BuildLnkKin(N==1,M,LnkKin,pKin)
          end if
        enddo DoReadOld
        !
      !------------------------------------------------------/old format
      else
      !------------------------------------------------------ new format
        !
        DoReadRock: do
        
          read(F_,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
          call LinToWrd(L,W1,EoL)
          if(W1(1:1)=='!')   cycle DoReadRock !skip comment lines
          call AppendToEnd(L,W1,EoL)
          select case(W1)
            case("ENDINPUT")              ; exit DoFile
            case("END","ENDDYNAMIC.ROCK") ; exit DoReadRock
          end select
          !
          I=    0
          jKin= 0
          xGram= Zero
          xMole= Zero
          !
          call KinFas_Zero(M)
          !
          M%NamKF= trim(W1) !-> the Kinetic Phase Name
          call LinToWrd(L,W2,EoL)
          !
          do
            
            if(EoL) exit
            !
            select case(trim(W2))
            !
            case default
              call Stop_(trim(W2)//" <<Unknown Keyword in DYNAMIC.ROCK")
            !
            case("SPECIES","SOLUTION","MIXTURE")
              if(trim(W2)=="SOLUTION") &
              & call Warning_(trim(W2)//" soon Obsolete, better use MIXTURE !!!")
              !
              call LinToWrd(L,W2,EoL)
              I= Phase_Index(W2,vFas)
              if(I==0) print '(A)',"SPECIES not found for "//trim(W2)
            !
            case("KINETICS")
              call LinToWrd(L,W2,EoL)
              jKin= KinModel_Index(W2,vKinModel)
              if(jKin==0) &
              & call Stop_("DYNAMIC.ROCK: kinetic model not found for "//trim(W2))
            !
            case("MODE") !=> mode= Radius / Density / Surface/ Equil / Adsorption
              call LinToWrd(L,W2,EoL) ; M%cMode= W2(1:1)
            !
            case("RADIUS") !=> radius (metre)
              call LinToWrd(L,W2,EoL) ; call WrdToReal(W2,M%Dat%Radius)
            !
            case("SURFACE") !=> specific surface m2/kg
              call LinToWrd(L,W2,EoL) ; call WrdToReal(W2,M%Dat%SurfKg) 
            !
            case("VOLUME") !=> volume (relative to other volumes in min'assemblages)
              call LinToWrd(L,W2,EoL) ; call WrdToReal(W2,M%Dat%PhiM)
            !
            case("MOLE") !=> mole number
              call LinToWrd(L,W2,EoL) ; call WrdToReal(W2,xMole)
            !
            case("GRAM") !=> gram number
              call LinToWrd(L,W2,EoL) ; call WrdToReal(W2,xGram)
            !
            end select
            !
            call LinToWrd(L,W2,EoL)
            
          enddo
          !
          if(M%Dat%SurfKg >Zero) M%Dat%Radius= -One
          ! otherwise Radius would take default value set in KinFas_Zero ...
          !
          if(I==0) I= Phase_Index(W1,vFas)
          !-> default species name= name of kinetic phase
          if(I==0) &
          & call Stop_("DYNAMIC.ROCK: Species not found for "//trim(M%NamKF))
          !
          if(M%cMode/="E" .and. M%cMode/="I") then
            if(jKin==0) jKin= KinModel_Index(W1,vKinModel) !default kin'model name= name of kinetic phase
            if(jKin==0) call Stop_("Kinetic model not found for "//trim(M%NamKF))
          end if
          !
          M%iFas= I
          if(xGram>Zero) M%Dat%PhiM= xGram /1.0D3 /vFas(M%iFas)%WeitKg *vFas(M%iFas)%VolM3
          if(xMole>Zero) M%Dat%PhiM= xMole *vFas(M%iFas)%VolM3
          !
          if(M%cMode/="E".and. M%cMode/="I") then
            M%iKin= jKin
          else
            M%iKin= 0
          end if
          !
          call KinFas_Surf_Zero(vFas,M)
          !-> calc' SurfKg from Radius, or Radius from SurfKg
          !
          N=N+1
          call BuildLnkKin(N==1,M,LnkKin,pKin)
          
        enddo DoReadRock
        !
      end if
      !------------------------------------------------------/new format
    !endcase("DYNAMIC.ROCK")
    !
    end select
    !
  enddo DoFile
  close(f_)
  !
  if(.not. BlockFound) print '(A)',"!!!WARNING!!! Block DYNAMIC.ROCK Not Found"
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ KinFas_BuildLnk"
  !
end subroutine KinFas_BuildLnk

subroutine KinFas_Check(vFas_,vKinFas_)
  use M_Numeric_Const,  only: Ln10
  use M_T_Phase,  only: T_Phase,Phase_Index
  !
  type(T_Phase), intent(in):: vFas_(:)
  type(T_KinFas),intent(in):: vKinFas_(:)
  !
  integer ::I,J
  real(dp)::RhoM
  
  do I=1,size(vKinFas_)
    !!J=vKinFas_(I)%iSpc !J= index (in vSpc_) of species named vKinFas(I)%Name
    !!RhoM=vSpc_(J)%WeitKg / vSpc_(J)%V0
    J=vKinFas_(I)%iFas !J= index (in vFas_) of species named vKinFas(I)%Name
    RhoM= vFas_(J)%WeitKg / vFas_(J)%VolM3
    write(fTrc,'(A,I3,A1,A,A15,A1,3(A7,F12.6,A1))') &
    & "iMk=",  I,                                  T_, &
    & "vKinFas=", vFas_(J)%NamFs,                  T_, &
    & "RhoM=", RhoM,                               T_, &
    & "VMol=", vFas_(vKinFas_(I)%iFas)%VolM3,      T_, &
    & "logK=", - vFas_(vKinFas_(I)%iFas)%Grt/Ln10, T_
  enddo
  
end subroutine KinFas_Check

end module M_KinFas_Read


