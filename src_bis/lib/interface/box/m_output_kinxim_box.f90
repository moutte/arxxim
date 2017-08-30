module M_output_KINXIM_BOX

  use M_Kinds,only: dp
  use M_iolog
  implicit none
  private

  public :: KINXIM_sorties_COORES
  public :: KINXIM_sorties_tag
  
  character(len=10) :: KINXIM_sorties_tag  = ""

contains
  !--------

  subroutine KINXIM_sorties_COORES

    use M_string
    use M_exchange_KINXIM_BOX
    use M_Box_Lumping, only : Box_Lumping_Output_Simres, Box_Lumping_Output_Coupler, Box_Lumping_Output_LumpCompo
    use M_Box_Lumping_Vars, only : LOption_Lumping
    implicit none

    integer, parameter :: nmax=100
    integer, parameter :: nunit = 1278 ! unite BILAN.dat

    integer, parameter :: nunitR = 1279 ! unite compx-solid.dat
    integer, parameter :: nunitW = 1280 ! unite compx-aqu.dat
    integer, parameter :: nunitP = 1281 ! unite compx-aqu-purewater.dat
    integer, parameter :: nunitG = 1282 ! unite geochemx.dat
    integer, parameter :: nunitL = 1283 ! unite lumping-model.dat
    integer, parameter :: nunitC = 1284 ! unite lumping-model.dat
    
    integer :: nbEAQ, nbMIN, nbC
    integer :: i,j, iligne, imaxligne,ix
    character (len=20), dimension(:), allocatable :: nomEAQ, nomC, nomMIN
    character (len=500) :: lineE, lineM !line, line2, 
    real (dp), dimension(100) :: NE, PHIM, logGAMMA, logQsKMIN
    character, dimension(100) :: contextMIN
    real (dp) :: MASSH2O, PHI, RHOW !NETOTAL, 
    integer :: iH2O, iHplus
    real (dp), dimension(:), allocatable :: NC, MASSMOLC, RHOM, MASSMOLM, MASSMOLE, ZE
    real (dp), dimension(:,:), allocatable :: EAQtoC, MINtoC

    integer , dimension(:), allocatable :: indexEAQ, indexMIN, indexC
    real (dp) :: Ionic, xsalt, ChargeBalance

    character (len=30) :: FileReport
    character (len=10) :: tag
    
    ! ------ instructions
    tag = adjustl(KINXIM_sorties_tag)
    FileReport = 'box-report'//trim(tag)//'.res'
    
!!$    FileReport = 'BILAN.res'
!!$    tag = ''
    
    write(*,*) "KINXIM_OUTPUT_COORES "// tag

    call KINXIM_get_nbEAQ(nbEAQ)
    call KINXIM_get_nbC(nbC)
    call KINXIM_get_nbMIN(nbMIN)
    call KINXIM_get_Ionic(Ionic)

    allocate(nomC(nbC))
    allocate(nomEAQ(nbEAQ))
    allocate(ZE(nbEAQ))
    allocate(nomMIN(nbMIN))

    call KINXIM_get_nomC(nomC)
    call KINXIM_get_nomMIN(nomMIN)
    call KINXIM_get_nomEAQ(nomEAQ)
    call KINXIM_get_ZE(ZE)

    call KINXIM_get_NE(NE(1:nbEAQ))
    call KINXIM_get_logGAMMAEAQ(logGAMMA(1:nbEAQ))
    call KINXIM_get_logQsKMIN(logQsKMIN(1:nbMIN))
    call KINXIM_get_contextMIN(contextMIN(1:nbMIN))
    call KINXIM_get_PHIM(PHIM(1:nbMIN))
    call KINXIM_get_nomMIN(nomMIN(1:nbMIN))

    allocate(EAQtoC(nbEAQ, nbC))
    allocate(MINtoC(nbMIN, nbC))
    allocate(MASSMOLC(nbC))
    allocate(NC(nbC))
    allocate(RHOM(nbMIN))
    allocate(MASSMOLM(nbMIN))
    allocate(MASSMOLE(nbEAQ))

    !write(33,*) '-------OUTPUT---------'
    call KINXIM_get_EAQtoC(EAQtoC, nbEAQ, nbC)

    call KINXIM_get_MINtoC(MINtoC, nbMIN, nbC)
    call KINXIM_get_massmolC(MASSMOLC)

    allocate(indexEAQ(nbEAQ))
    allocate(indexC(nbC))
    allocate(indexMIN(nbMIN))
    
    open(unit=nunit, file=FileReport)
    call iolog_open(nunit, FileReport, ' ')

    open(unit=nunitW, file='compx-aqu'//trim(tag)//'.dat')
    call iolog_open(nunitW, 'compx-aqu'//trim(tag)//'.dat',' ')

    open(unit=nunitP, file='compx-aqu-purewater.dat')
    call iolog_open(nunitP, 'compx-aqu-purewater.dat',' ')
    
    open(unit=nunitG, file='geochemx.dat')
    call iolog_open(nunitW, 'geochemx.dat',' ')

    open(unit=nunitR, file='compx-solid'//trim(tag)//'.dat')
    call iolog_open(nunitR, 'compx-solid'//trim(tag)//'.dat', ' ')

    if (LOption_Lumping) then
      open(unit=nunitL, file='pvt-model.dat')
      call iolog_open(nunitL, 'pvt-model.dat', ' ')
   
      open(unit=nunitC, file='coupler-model.cor')
      call iolog_open(nunitC, 'coupler-model.cor', ' ')
    end if

    ! ----- find iH2O and iH+
    do i=1, nbEAQ
       if (nomEAQ(i)=='h2o') iH2O=i
       if (nomEAQ(i)=='h+') iHplus=i
    end do

    ! ----- index lexicographique
    call string_sort(nbEAQ, nomEAQ, indexEAQ)
    call string_sort(nbC, nomC, indexC)
    call string_sort(nbMIN, nomMIN, indexMIN)

    !----------

    nC(1:nbC)=matmul(transpose(EAQtoC),NE(1:nbEAQ))
!lt    RHOW = sum(NC(1:nbC) * MASSMOLC(1:nbC))/PHI

    call KINXIM_get_RHOM(RHOM(1:nbMIN))
    MASSMOLM=0.d0

    do i=1, nbMIN
       do j=1, nbC
          MASSMOLM(i)=MASSMOLM(i)+MINtoC(i, j)*MASSMOLC(j)
       end do
    end do

!!$    write(*,*) '------------------------------'
!!$    do i=1, nbEAQ
!!$       write(*,*) '---------------'
!!$       write(*,'(A,I6, A10)') "EAQ [",i,"] "//trim(nomEAQ(i))
!!$       do j=1, nbC
!!$          if (abs(EAQtoC(i,j))>1.d-10) then
!!$             write(*,'(A,F8.3,A)') "           ", EAQtoC(i,j), ' '//trim(nomC(j))
!!$          end if
!!$       end do
!!$    end do
!!$    write(*,*) '------------------------------'

    MASSMOLE=0.d0

    do ix=1, nbEAQ
       i=indexEAQ(ix)
       !!write(*,*) '---------------'
       !!write(*,*) "i, EAQ = ", trim(nomEAQ(i))
       do j=1, nbC
          !write(*,*) " --->     ", EAQtoC(i,j), trim(nomC(j)) , MASSMOLC(j)
          MASSMOLE(i)=MASSMOLE(i)+ EAQtoC(i,j)*MASSMOLC(j)
       end do
       !!write(*,*) 'MASSMOL = ', MASSMOLE(i)
    end do
    !write(*,*) '------------------------------'

    write(nunit,'(A)') "<<----------------------------------------------------------"
    write(nunit,'(A)') "<< KINXIM BOX CHEMISTRY REPORT                              "
    write(nunit,'(A)') "<<----------------------------------------------------------"
    write(nunit,'(A)') " "
    write(nunit,*)     " - geochemx.dat         : geochemical system"
    write(nunit,*)     " - compx-aqu.dat        : aqueous phase composition"
    write(nunit,*)     " - compx-solid.dat      : rock mineral composition"
    write(nunit,*)     " "
    write(nunit,'(A)') "<<-----------  geochemx.dat - Geochemical System -----"
    
    iligne=0
    imaxligne=8
    lineM=""
    write(nunit,'(A)') "MINERAL ="
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunit,'(A)') '"'//trim(nomMIN(i))//'"'
       if (iligne>imaxligne) then ! passage a la ligne
          !write(nunit,'(A)') lineM
          lineM=""
          iligne=0
       end if
       iligne = iligne+1
       lineM = trim(lineM)//' "'//trim(nomMIN(i))//'"'
    end do
    !write(nunit,*) trim(lineM)
    write(nunit,*) " "

    iligne=0
    imaxligne=8
    lineE=""
    write(nunit,'(A)') "AQUSPECIES = "

    do ix=1,nbEAQ
       i=indexEAQ(ix)
       if (i==iH2O) then
          write(nunit,'(A)') '< "'//trim(nomEAQ(i))//'"'
       else
          write(nunit,'(A)') '"'//trim(nomEAQ(i))//'"'
       end if
       if (iligne>imaxligne) then ! passage a la ligne
          !write(nunit,*) trim(lineE)
          lineE=""
          iligne=0
       end if
       iligne=iligne+1
       lineE=trim(lineE)//' "'//trim(nomEAQ(i))//'"'
    end do
    ! write(nunit,*) trim(lineE)


    !--------------------------------------------

    write(nunitG,'(A)') "<< KinXim Result "
    write(nunitG,'(A)') "<<-----------  geochemx.dat - Geochemical System -----"
    
    iligne=0
    imaxligne=8
    lineM=""
    write(nunitG,'(A)') "MINERAL ="
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunitG,'(A)') '"'//trim(nomMIN(i))//'"'
       if (iligne>imaxligne) then ! passage a la ligne
          !write(nunit,'(A)') lineM
          lineM=""
          iligne=0
       end if
       iligne = iligne+1
       lineM = trim(lineM)//' "'//trim(nomMIN(i))//'"'
    end do
    !write(nunit,*) trim(lineM)
    write(nunitG,*) " "

    iligne=0
    imaxligne=8
    lineE=""
    write(nunitG,'(A)') "AQUSPECIES = "

    do ix=1,nbEAQ
       i=indexEAQ(ix)
       if (i==iH2O) then
          write(nunitG,'(A)') '< "'//trim(nomEAQ(i))//'"'
       else
          write(nunitG,'(A)') '"'//trim(nomEAQ(i))//'"'
       end if
       if (iligne>imaxligne) then ! passage a la ligne
          !write(nunit,*) trim(lineE)
          lineE=""
          iligne=0
       end if
       iligne=iligne+1
       lineE=trim(lineE)//' "'//trim(nomEAQ(i))//'"'
    end do
    ! write(nunit,*) trim(lineE)

    !-----------------------------------------

    MASSH2O = NE(iH2O)*MASSMOLE(iH2O)
    ChargeBalance = sum(NE(1:nbEAQ)*ZE(1:nbEAQ))/MASSH2O
    call KINXIM_get_SALINITY(xsalt)

    write(nunit,'(A)') "  "
    write(nunit,'(A)') "<<----------- compx-aqu.dat - Water Composition ( mol/kgH2O ) ------------------"
    write(nunit,'(A)') '<BOUND-AQU "BOUND " ='
    write(nunit,'(A)') '<COMP-AQU "ZONE" ='
    
    do ix=1,nbEAQ
       i=indexEAQ(ix)
       if (i==iH2O) then
          write (nunit,'(1A, D20.13, 1A, 1I3, 1A, 1A20, 1F8.0)')  &
               '<', NE(iH2O)/MASSH2O, '   < [',iH2O,' ] ',   nomEAQ(iH2O), ZE(iH2O)
       else
          write(nunit,'(1D20.13, 1A, 1I3,1A,  1A20, 1F8.0)')   &
               NE(i)/MASSH2O,    '   < [', i,' ] ',  nomEAQ(i), ZE(i)
       end if
    end do

    write(nunit,*)
    write(nunit,'(A)') "<===== Salinity ( Equivalent to NaCl g/kgH2O )"
    
    write(nunit,'(1A,1A,1G20.13)') &
         "<",'SALINITY = ', xsalt/MASSH2O   
    
    if (LOption_Lumping) call Box_Lumping_Output_LumpCompo(nunit) 

    !--- copy to compaqu.dat ---
    write(nunitW,'(A)') "<< ArXim Result "
    write(nunitW,'(A)') "<<----------- compx-aqu.dat - Water Composition ( mol/kgH2O ) ------------------"
    write(nunitW,'(A)') '<BOUND-AQU "BOUND " ='
    write(nunitW,'(A)') '<COMP-AQU "ZONE" ='
   
    do ix=1,nbEAQ
       i=indexEAQ(ix)
       if (i==iH2O) then
          write (nunitW,'(1A, D20.13, 1A, 1I3, 1A, 1A20)')  &
               '<', NE(iH2O)/MASSH2O, '   < [',iH2O,' ] ',   nomEAQ(iH2O)
       else
          write(nunitW,'(1D20.13, 1A, 1I3,1A,  1A20)')   &
               NE(i)/MASSH2O,    '   < [', i,' ] ',  nomEAQ(i)
       end if
    end do

    !-------------------------

    write(nunitP,*) " "
    write(nunitP,'(A)') "<<----------- compx-aqu-purewater.dat - Empty Water Composition ( mol/kgH2O ) --------------"
    write(nunitP,'(A)') '<BOUND-AQU "BOUND " ='
    write(nunitP,'(A)') '<COMP-AQU "ZONE" ='

    do ix=1,nbEAQ
       i=indexEAQ(ix)
       if (i==iH2O) then
          write (nunitP,'(1A, D20.13, 1A, 1I3, 1A, 1A20)')  &
               '<', NE(iH2O)/MASSH2O, '   < [',iH2O,' ] ',   nomEAQ(iH2O)
       else
          write(nunitP,'(A20, 1A, 1I3,1A,  1A20)')   &
               '0',    '   < [', i,' ] ',  nomEAQ(i)
       end if
    end do    

    write(nunit,*) " "

    PHI = 1.d0 - sum(PHIM(1:nbMIN))
    write(nunit,'(A)') "<<----------- info Kinxim -------------------------"
    write(nunit,'(1A, 1D20.13)') "<  Porosite :: ", PHI
    write(nunit,*) " "
    write(nunit,'(A)') "<<----------- compx-solid.dat - Rock Mineral Composition ( m3/m3 ) ------"


    write(nunit,'(A)') '<COMP-SOLID "ZONE" ='
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunit,'(1D20.13, 1A, 1I3,1A,  1A20)') &
            PHIM(i)/(1.D0 - PHI),    '   < [', i,' ] ',  nomMIN(i)
    end do

    !--- copy to comprock.dat ---
    write(nunitR,'(A)') "<< ArXim Result "
    write(nunitR,'(A)') "<<----------- compx-solid.dat - Rock Mineral Composition ( m3/m3 ) ------"
    write(nunitR,'(A)') '<COMP-SOLID "ZONE" ='
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunitR,'(1D20.13, 1A, 1I3,1A,  1A20)') &
            PHIM(i)/(1.D0 - PHI),    '   < [', i,' ] ',  nomMIN(i)
    end do
    !----------------------------

    write(nunit,*)
    write(nunit,'(A)') "<<----------- Final Result Informations ------------------------"
    
    write(nunitW,'(A)') 
    write(nunitW,'(A)') "<<----------- Additional Informations ------------------------"
    
    !--- Balance ---
    write(nunit,*)
    write(nunit,'(A)') "<===== Charge Balance ( - )"
    write(nunit,'(1A,1A,1G20.13)') &
         "<",' Balance = ', ChargeBalance

    write(nunitW,'(1A,1A20,1G20.13)') &
         "<",' Charge Balance = ', ChargeBalance

    !--- Ionic Strength ---
    write(nunit,*)
    write(nunit,'(A)') "<===== Ionic Strength ( - )"
    write(nunit,'(1A,1A,1G20.13)') &
         "<",' Ionic Strength = ', Ionic
    
    write(nunitW,'(1A,1A20,1G20.13)') &
         "<",' Ionic Strength = ', Ionic

    !--- Salinity ---
       
    write(nunit,*)
    write(nunit,'(A)') "<===== Salinity ( Equivalent to NaCl g/kgH2O )"
    write(nunit,'(1A,1A,1G20.13)') &
         "<",' SALINITY = ', xsalt/MASSH2O

    write(nunitW,'(1A,1A20,1G20.13, A)') &
         "<",' Salinity = ', xsalt/MASSH2O
    
    !--- pH ---
    write(nunit,*)
    write(nunit,'(A)') "<===== pH ( - )"
    write(nunit,'(1A,1A,1G20.13)') &
         "<",' pH = -log10 (a(H+)) = ', &
         -(log10(NE(iHplus)) - log10(MASSH2O) + logGAMMA(iHplus))

    write(nunitW,'(1A,1A20,1G20.13)') &
         "<",' pH = ', &
         -(log10(NE(iHplus)) - log10(MASSH2O) + logGAMMA(iHplus))
    
    !--- logQSK ---
    write(nunit,*)
    write(nunit,'(A)') "<===== QsK( - ) , logQsK( - ) , mode ( - )"
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunit,'(1A, 1A20,2G20.8, 1A2, 2A)')  "<",nomMIN(i), 10.d0**logQsKMIN(i),  logQsKMIN(i), &
            ' [', contextMIN(i),']'
    end do

    write(nunitW,'(A)') "<<---"
    do ix=1,nbMIN
       i=indexMIN(ix)
       write(nunitW,'(1A, 1A20, 1A, 1G20.8)')  "< logQsK ",nomMIN(i), " = " , logQsKMIN(i) 
    end do

     !--- logGamma ---
    write(nunit,*)
    write(nunit,'(A)') "<===== Gamma (- ) , logGamma ( - )"
    do ix=1,nbEAQ
       i=indexEAQ(ix)
       write(nunit,'(1A, 1A20,2G20.8)')  "< ",nomEAQ(i), 10.d0**logGAMMA(i),  logGAMMA(i)
    end do

    write(nunit,*)
    write(nunit,'(A)') "<<----------- Element Molal Composition ( mol/kgH2O ) --------------"
    do ix=1,nbC
       i=indexC(ix)
       write(nunit,'(1A,1A20,1D20.13)')  "< ",nomC(i) , NC(i)/MASSH2O
    end do

    write(nunitW,*)
    write(nunitW,'(A)') "<<----------- Elemental Composition ----------------------"
    do ix=1,nbC
       i=indexC(ix)
       write(nunitW,'(1A,1A20,1A, 1D20.13)')  "< Molality ",nomC(i) , " = ", NC(i)/MASSH2O
    end do
            
    !------------------------- Lumping 
    
    if (LOption_Lumping) then
      write(nunit,*) " "
      call Box_Lumping_Output_Coupler(nunit) 
      call Box_Lumping_Output_Simres(nunit)
      
      call Box_Lumping_Output_LumpCompo(nunitW)
      call Box_Lumping_Output_Coupler(nunitC) 
      call Box_Lumping_Output_Simres(nunitL)      
    end if

    !------------------------- end FILE

    write(nunit,*) " "
    write(nunit,'(A)') "<------- end "

    call iolog_close(nunit)
    close(nunit)

    call iolog_close(nunitW)
    close(nunitW)
    
    call iolog_close(nunitP)
    close(nunitP)

    call iolog_close(nunitR)
    close(nunitR)

    call iolog_close(nunitG)
    close(nunitG)

    if (LOption_Lumping) then
    call iolog_close(nunitL)
    close(nunitL)
    call iolog_close(nunitC)
    close(nunitC)
    end if
 
    deallocate(MINtoC)
    deallocate(EAQtoC)
    deallocate(MASSMOLC)
    deallocate(MASSMOLM)
    deallocate(NC)
    deallocate(RHOM)
    deallocate(ZE)

  end subroutine KINXIM_sorties_COORES
    
end module M_output_KINXIM_BOX
