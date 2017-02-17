module M_stockvar_KINXIM

  use M_Kinds,only: dp
  
  implicit none
  private

  !
  ! module de stockage des valeurs de flux par espece
  ! sur l'ensemble des pas de temps ( MAXNBTIMESTEP )
  !
  logical, public::  LSTOCK =.false. !.true. !!!JM
  integer :: MAX_NBTIMESTEP
  integer, public:: size_X
  logical, private :: warnout =.false.

  integer, public:: index_stockvar
  real(dp), allocatable, dimension (:,:), public:: stock_X
  real(dp), allocatable, dimension (:), public:: stock_T
  real(dp), allocatable, dimension (:), public:: stock_DT
  real(dp), allocatable, dimension (:), public:: stock_V
  real(dp), allocatable, dimension (:), public:: stock_EPSILON
  character,    allocatable, dimension (:), public:: stock_TAG

  !// public functions
  public:: reset_stockvar
  public:: set_stockvar
  public:: init_stockvar
  public:: del_stockvar
  public:: print_stockvar
  public:: computemean_stockvar

contains

  !----

  subroutine print_stockvar(filename)

    implicit none
    integer :: nunit = 45673

    character(len=*) :: filename
    character(len=30) :: monformat

    integer :: iL, sizetable
    real(dp), dimension(:), allocatable :: table

    real(dp), dimension(:), allocatable :: X
    real(dp) :: T, DT, EPSILON, V !, tol

!!$
!!$    allocate(X(size_X))
!!$
!!$    call computemean_stockvar(X,V, T,DT, EPSILON)
!!$    call set_stockvar(X,V, T,DT, EPSILON,'B')
!!$
!!$    open(unit=nunit, file=filename)
!!$
!!$    sizetable = size_X + 4
!!$    allocate(table(sizetable))
!!$
!!$    do iL=1, index_stockvar-1
!!$       write(monformat, '(1A,1I2,1A)') '(1A, 1I6,',sizetable,'D14.6)'
!!$       !
!!$       table(1) = stock_T (iL)
!!$       table(2) = stock_DT (iL)
!!$       table(3) = stock_EPSILON (iL)
!!$       table(4) = stock_V (iL)
!!$       table(5:size_X + 4)= stock_X (iL,1:size_X)
!!$
!!$       write(nunit,trim(monformat)) stock_TAG(iL), iL, table(1: sizetable)
!!$
!!$    end do
!!$
!!$    deallocate(table)
!!$    deallocate(X)
!!$
!!$    close(nunit)

  end subroutine print_stockvar

  !----

  subroutine reset_stockvar

    implicit none

!!$    stock_X  ( 1:index_stockvar, 1:size_X ) = 0.D0
!!$    stock_T ( 1:index_stockvar )  = 0.D0
!!$    stock_DT ( 1:index_stockvar ) = 0.D0
!!$    stock_EPSILON ( 1:index_stockvar ) = 0.D0
!!$
!!$    index_stockvar = 1
!!$    warnout=.false.

  end subroutine reset_stockvar

  !----

  subroutine set_stockvar(X, V, T, DT, EPSILON, TAG)

    implicit none
    real(dp), dimension (:) :: X
    real(dp) :: T, DT, V
    real(dp) :: EPSILON
    character :: TAG

     !!write(*,*) '-----------------Set stockvar', X, V, T, DT, EPSILON, TAG 
     if (index_stockvar < MAX_NBTIMESTEP) then
       stock_X (index_stockvar,1:size_X) = X( 1:size_X)
       stock_DT (index_stockvar) = DT
       stock_V (index_stockvar) = V
       stock_T (index_stockvar) = T
       stock_EPSILON (index_stockvar) = EPSILON
       stock_TAG (index_stockvar) = TAG

       index_stockvar = index_stockvar + 1
    else
       if (.not.(warnout)) then
          write(*,*) "==============================================="
          write(*,*) "stop stockvar ...  NBTIMESTEP >", MAX_NBTIMESTEP, '!'
          write(*,*) "==============================================="
          warnout=.true.
          stop "GEOCHIMIE : NBTIMESTEP > MAX_NBTIMESTEP"
       end if
    end if

  end subroutine set_stockvar

  !----

  subroutine init_stockvar (maxtstep, nbx)

    implicit none
    integer :: maxtstep
    integer :: nbx

    MAX_NBTIMESTEP = maxtstep
    size_X = nbx
    !!!write(*,*) '-----------------Initialization stockvar ::', size_X
    allocate(stock_X(MAX_NBTIMESTEP,size_X))
    allocate(stock_DT(MAX_NBTIMESTEP))
    allocate(stock_T(MAX_NBTIMESTEP))
    allocate(stock_EPSILON(MAX_NBTIMESTEP))
    allocate(stock_TAG(MAX_NBTIMESTEP))
    allocate(stock_V(MAX_NBTIMESTEP))

    index_stockvar = 1

  end subroutine init_stockvar

  ! ----

  subroutine del_stockvar

    implicit none
    deallocate(stock_X)
    deallocate(stock_DT)
    deallocate(stock_T)
    deallocate(stock_EPSILON)
    deallocate(stock_TAG)
    deallocate(stock_V)
    !!!write(*,*) '-----------------Delete stockvar', size_X
  end subroutine del_stockvar

!---

  subroutine computemean_stockvar(X,V, T,DT, EPSILON)

    implicit none
    integer :: Nstep
    integer :: istep !ivar,
    real(dp), dimension(:) :: X
    real(dp) :: T, tol, DT, EPSILON, V
    !character(len=8) :: strI
    !integer :: iFileStock = 1
    !character(len=40) :: Filename
    !---

    X = 0.D0
    T=0.D0
    tol=0.D0
    DT=0.D0
    EPSILON=0.D0
    Nstep=0
   

    do istep=1, index_stockvar-1
       if (stock_TAG(istep)=='C') then
         
          Nstep = Nstep+1
          DT = stock_DT(istep)
          
          X(1:size_X) = X(1:size_X) + stock_X(istep,1:size_X) * DT
          EPSILON =  EPSILON + stock_EPSILON(istep) * DT
          T = T + DT
          V = V + stock_V(istep) * DT
!!$          write(*,*) 'found data', Nstep
!!$          write(*,*) 'V', V
!!$          write(*,*) 'T', T
!!$          write(*,*) 'EPSILON', EPSILON
!!$          write(*,*) 'X', X
       end if
    end do

    X(1:size_X) = X(1:size_X)/ T
    EPSILON = EPSILON/T
    DT = T/Nstep
    V = V/T
!!$    write(*,*) 'found data', Nstep
!!$    write(*,*) 'V', V
!!$    write(*,*) 'T', T
!!$    write(*,*) 'EPSILON', EPSILON
!!$    write(*,*) 'X', X
!!$    write(*,*) '--------------------'
  end subroutine computemean_stockvar



end module M_stockvar_KINXIM
