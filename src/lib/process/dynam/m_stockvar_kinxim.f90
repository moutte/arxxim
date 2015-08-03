MODULE M_stockvar_KINXIM

  IMPLICIT NONE
  PRIVATE

  !
  ! MODULE de stockage des valeurs de flux par espece
  ! sur l'ensemble des pas de temps ( MAXNBTIMESTEP )
  !
  LOGICAL, PUBLIC::  LSTOCK =.false. !.true. !!!JM
  INTEGER :: MAX_NBTIMESTEP
  INTEGER, PUBLIC:: SIZE_X
  LOGICAL, PRIVATE :: warnout =.false.

  INTEGER, PUBLIC:: index_stockvar
  REAL(kind=8), ALLOCATABLE, DIMENSION (:,:), PUBLIC:: stock_X
  REAL(kind=8), ALLOCATABLE, DIMENSION (:), PUBLIC:: stock_T
  REAL(kind=8), ALLOCATABLE, DIMENSION (:), PUBLIC:: stock_DT
  REAL(kind=8), ALLOCATABLE, DIMENSION (:), PUBLIC:: stock_V
  REAL(kind=8), ALLOCATABLE, DIMENSION (:), PUBLIC:: stock_EPSILON
  CHARACTER,    ALLOCATABLE, DIMENSION (:), PUBLIC:: stock_TAG

  !// PUBLIC FUNCTIONs
  PUBLIC:: reset_stockvar
  PUBLIC:: set_stockvar
  PUBLIC:: init_stockvar
  PUBLIC:: del_stockvar
  PUBLIC:: PRINT_stockvar
  PUBLIC:: computemean_stockvar

CONTAINS

  !----

  SUBROUTINE PRINT_stockvar(fiLEName)

    IMPLICIT NONE
    INTEGER :: nunit = 45673

    CHARACTER(LEN=*) :: fiLEName
    CHARACTER(LEN=30) :: monformat

    INTEGER :: iL, SIZEtable
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: table

    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: X
    REAL(kind=8) :: T, DT, EPSILON, V !, tol

!!$
!!$    ALLOCATE(X(SIZE_X))
!!$
!!$    CALL computemean_stockvar(X,V, T,DT, EPSILON)
!!$    CALL set_stockvar(X,V, T,DT, EPSILON,'B')
!!$
!!$    OPEN(unit=nunit, file=fiLEName)
!!$
!!$    SIZEtable = SIZE_X + 4
!!$    ALLOCATE(table(SIZEtable))
!!$
!!$    DO iL=1, index_stockvar-1
!!$       WRITE(monformat, '(1A,1I2,1A)') '(1A, 1I6,',SIZEtable,'D14.6)'
!!$       !
!!$       table(1) = stock_T (iL)
!!$       table(2) = stock_DT (iL)
!!$       table(3) = stock_EPSILON (iL)
!!$       table(4) = stock_V (iL)
!!$       table(5:SIZE_X + 4)= stock_X (iL,1:SIZE_X)
!!$
!!$       WRITE(nunit,TRIM(monformat)) stock_TAG(iL), iL, table(1: SIZEtable)
!!$
!!$    END DO
!!$
!!$    DEALLOCATE(table)
!!$    DEALLOCATE(X)
!!$
!!$    CLOSE(nunit)

  END SUBROUTINE PRINT_stockvar

  !----

  SUBROUTINE reset_stockvar

    IMPLICIT NONE

!!$    stock_X  ( 1:index_stockvar, 1:SIZE_X ) = 0.D0
!!$    stock_T ( 1:index_stockvar )  = 0.D0
!!$    stock_DT ( 1:index_stockvar ) = 0.D0
!!$    stock_EPSILON ( 1:index_stockvar ) = 0.D0
!!$
!!$    index_stockvar = 1
!!$    warnout=.false.

  END SUBROUTINE reset_stockvar

  !----

  SUBROUTINE set_stockvar(X, V, T, DT, EPSILON, TAG)

    IMPLICIT NONE
    REAL(kind=8), DIMENSION (:) :: X
    REAL(kind=8) :: T, DT, V
    REAL(kind=8) :: EPSILON
    CHARACTER :: TAG

     !!WRITE(*,*) '-----------------Set stockvar', X, V, T, DT, EPSILON, TAG 
     IF (index_stockvar < MAX_NBTIMESTEP) THEN
       stock_X (index_stockvar,1:SIZE_X) = X( 1:SIZE_X)
       stock_DT (index_stockvar) = DT
       stock_V (index_stockvar) = V
       stock_T (index_stockvar) = T
       stock_EPSILON (index_stockvar) = EPSILON
       stock_TAG (index_stockvar) = TAG

       index_stockvar = index_stockvar + 1
    ELSE
       IF (.not.(warnout)) THEN
          WRITE(*,*) "==============================================="
          WRITE(*,*) "STOP stockvar ...  NBTIMESTEP >", MAX_NBTIMESTEP, '!'
          WRITE(*,*) "==============================================="
          warnout=.true.
          stop "GEOCHIMIE : NBTIMESTEP > MAX_NBTIMESTEP"
       END IF
    END IF

  END SUBROUTINE set_stockvar

  !----

  SUBROUTINE init_stockvar (maxtstep, nbx)

    IMPLICIT NONE
    INTEGER :: maxtstep
    INTEGER :: nbx

    MAX_NBTIMESTEP = maxtstep
    SIZE_X = nbx
    !!!WRITE(*,*) '-----------------Initialization stockvar ::', SIZE_X
    ALLOCATE(stock_X(MAX_NBTIMESTEP,SIZE_X))
    ALLOCATE(stock_DT(MAX_NBTIMESTEP))
    ALLOCATE(stock_T(MAX_NBTIMESTEP))
    ALLOCATE(stock_EPSILON(MAX_NBTIMESTEP))
    ALLOCATE(stock_TAG(MAX_NBTIMESTEP))
    ALLOCATE(stock_V(MAX_NBTIMESTEP))

    index_stockvar = 1

  END SUBROUTINE init_stockvar

  ! ----

  SUBROUTINE del_stockvar

    IMPLICIT NONE
    DEALLOCATE(stock_X)
    DEALLOCATE(stock_DT)
    DEALLOCATE(stock_T)
    DEALLOCATE(stock_EPSILON)
    DEALLOCATE(stock_TAG)
    DEALLOCATE(stock_V)
    !!!WRITE(*,*) '-----------------Delete stockvar', SIZE_X
  END SUBROUTINE del_stockvar

!---

  SUBROUTINE computemean_stockvar(X,V, T,DT, EPSILON)

    IMPLICIT NONE
    INTEGER :: Nstep
    INTEGER :: istep !ivar,
    REAL(kind=8), DIMENSION(:) :: X
    REAL(kind=8) :: T, tol, DT, EPSILON, V
    !CHARACTER(LEN=8) :: strI
    !INTEGER :: iFileStock = 1
    !CHARACTER(LEN=40) :: FiLEName
    !---

    X = 0.D0
    T=0.D0
    tol=0.D0
    DT=0.D0
    EPSILON=0.D0
    Nstep=0
   

    DO istep=1, index_stockvar-1
       IF (stock_TAG(istep)=='C') THEN
         
          Nstep = Nstep+1
          DT = stock_DT(istep)
          
          X(1:SIZE_X) = X(1:SIZE_X) + stock_X(istep,1:SIZE_X) * DT
          EPSILON =  EPSILON + stock_EPSILON(istep) * DT
          T = T + DT
          V = V + stock_V(istep) * DT
!!$          WRITE(*,*) 'found DATA', Nstep
!!$          WRITE(*,*) 'V', V
!!$          WRITE(*,*) 'T', T
!!$          WRITE(*,*) 'EPSILON', EPSILON
!!$          WRITE(*,*) 'X', X
       END IF
    END DO

    X(1:SIZE_X) = X(1:SIZE_X)/ T
    EPSILON = EPSILON/T
    DT = T/Nstep
    V = V/T
!!$    WRITE(*,*) 'found DATA', Nstep
!!$    WRITE(*,*) 'V', V
!!$    WRITE(*,*) 'T', T
!!$    WRITE(*,*) 'EPSILON', EPSILON
!!$    WRITE(*,*) 'X', X
!!$    WRITE(*,*) '--------------------'
  END SUBROUTINE computemean_stockvar



ENDMODULE M_stockvar_KINXIM
