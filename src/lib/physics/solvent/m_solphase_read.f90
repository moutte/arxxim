module M_SolPhase_Read
!--
!-- tools for reading solution phases
!-- (- phase of variable composition, with assymetric model)
!-- from text files,
!-- and build an array of T_SolPhase (- container for a mixture phase)
!--
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_,Warning_
  !
  implicit none
  !
  private
  !
  public:: SolPhase_Read

contains

subroutine SolPhase_Read(vSpc,vSolModel,Ok,MsgError)
!--
!-- scan the ELECTROLYTE block(s)
!--
  use M_IoTools
  use M_Files,     only: NamFInn
  use M_T_Species, only: T_Species,Species_Index
  use M_T_SolModel,only: T_SolModel
  use M_T_SolPhase,only: T_SolPhase !,SolPhase_ActPole
  !
  use M_Global_Vars,only: vSolFas
  !
  type(T_Species), intent(in) :: vSpc(:)
  type(T_SolModel),intent(in) :: vSolModel(:)
  logical,         intent(out):: Ok
  character(*),    intent(out):: MsgError
  !
  character(len=255):: L !,sList0 !=elements in Database
  character(len=80) :: W !,V1,V2
  logical :: EoL
  integer :: ios,fInn
  !
  type(T_SolModel)  :: SolModel
  type(T_SolPhase)  :: SolFas
  type(T_SolPhase)  :: vSolTmp(10)
  !
  integer :: I,J,K,nS,N
  real(dp):: X1

  if(iDebug>0) write(fTrc,'(/,A)') "< SolPhase_Read"

  Ok= .true.
  MsgError= "Ok"

  call GetUnit(fInn)
  open(fInn,file=trim(NamFInn))

  N=0
  !
  DoFile: do
    read(fInn,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)
    !
    select case(W)

    case("ENDINPUT"); exit DoFile

    case("ELECTROLYTE")

      ! ELECTROLYTE WATER1
      !   MODEL MYAQUEOUS-MODEL
      !   COMPOSITION
      !     H2O   55.55
      !     CA+2  0.001
      !     NA+   0.003
      !     NACL(AQ) 0.002
      !     H+   1.e-7
      !     OH-  1.e-7
      !   end
      ! end

      call LinToWrd(L,W,EoL) !-> phase name

      if(Species_Index(W,vSpc)>0) then
        Ok= .false.
        MsgError= "Solution name should not match with any species name !!!"
        return !--------------------------------------------------return
      end if

      SolFas%Name=trim(W)

      DoSolFas: do
        !
        read(fInn,'(A)',iostat=ios) L  ;  if(ios/=0) exit DoFile
        call LinToWrd(L,W,EoL)
        if(W(1:1)=='!') cycle DoSolFas
        call AppendToEnd(L,W,EoL)
        !
        select case(W)

        case("ENDINPUT")             ; exit DoFile

        case("END","ENDELECTROLYTE") ; exit DoSolFas

        !----------------------------------------------- read model name 
        case("MODEL")
        ! read model name and check it is among the model base vSolModel
        ! -> return its index in array vSolModel
          call LinToWrd(L,W,EoL)
          !
          !find the model name in vSolModel%Name
          K=0
          do I=1,size(vSolModel)
            if(trim(W)==trim(vSolModel(I)%Name)) then
              K= I
              exit
            end if
          end do
          if(K==0) then
            Ok= .false.
            MsgError=             trim(W)//"= ELECTROLYTE.MODEL UNKNOWN"
            return !----------------------------------------------return
          end if

          if(iDebug>0) write(fTrc,'(2A)') "MODEL=", trim(W)
          SolFas%iModel= K

          nS= vSolModel(K)%nSpecies
          allocate(SolFas%vXSpecies(nS))
          SolFas%vXSpecies(:)= Zero !!1.0D-32

        !/case("MODEL")
        !-----------------------------------------------/read model name 

        !-----------------------------------------/read phase compositon
        case("COMPOSITION")
          !
          if(iDebug>0) write(fTrc,'(A)') "!!! !!!Read_COMPOSITION"
          !
          if(SolFas%iModel==0) then
            Ok= .false.
            MsgError=             "must define MODEL before COMPOSITION"
            return !----------------------------------------------return
          end if
          !
          SolModel=vSolModel(SolFas%iModel)
          !-------------------input on several lines, one per end-member
          do

            read(fInn,'(A)') L
            call LinToWrd(L,W,EoL)
            if(W=="!") cycle !-> comment line
            call AppendToEnd(L,W,EoL)

            if(W=="END") exit
            if(W=="ENDCOMPOSITION") exit

            I= Species_Index(trim(W),vSpc)
            !--- search for the species with global index I
            !--- in the species list, vISpecies(:), of solution model SolModel
            !--- and read mole fraction
            J= 0
            !~ print *,"I=",I
            do K=1,SolModel%nSpecies
              !~ print *,"  vISpecies(K)=",SolModel%vISpecies(K)
              if(SolModel%vISpecies(K)==I) then
                J= I
                exit
              end if
            end do

            if(J==0) then
              Ok= .false.
              MsgError=   trim(W)//" NOT FOUND in species list of MODEL"
              return !--------------------------------------------return
            end if

            call LinToWrd(L,W,EoL) !-> rest of line -> composition
            call WrdToReal(W,X1)

            SolFas%vXSpecies(K)= X1

          end do
          !
          N=N+1
          !
          if(iDebug>0) &
          & write(fTrc,'(I3,1X,A24,A24)') N, SolFas%Name, SolModel%Name

          vSolTmp(N)%Name=   SolFas%Name
          vSolTmp(N)%iModel= SolFas%iModel
          nS= vSolModel(SolFas%iModel)%nSpecies
          allocate(vSolTmp(N)%vXSpecies(nS))
          vSolTmp(N)%vXSpecies(:)= SolFas%vXSpecies(:)

          deallocate(SolFas%vXSpecies)

        !/case("COMPOSITION")
        !-----------------------------------------/read phase compositon
        !
        end select !case(W)
        !
      end do DoSolFas
    !/case(ELECTROLYTE)

    end select !case(W)

  end do DoFile
  !
  close(fInn)

  if(N>0) then
    if(allocated(vSolFas)) deallocate(vSolFas)
    allocate(vSolFas(N))
    do I=1,N
      vSolFas(I)%Name= vSolTmp(I)%Name
      vSolFas(I)%iModel= vSolTmp(I)%iModel
      nS= vSolModel(vSolFas(I)%iModel)%nSpecies
      allocate(vSolFas(I)%vXSpecies(nS))
      vSolFas(I)%vXSpecies(:)= vSolTmp(I)%vXSpecies(:)
    end do
  end if

  if(iDebug>0) write(fTrc,'(A,/)') "</ SolPhase_BuildLnk"

  return
end subroutine SolPhase_Read

end module M_SolPhase_Read
