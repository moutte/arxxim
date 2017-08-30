!! *******************************************************************
!! * File name:     M_string.f90
!! *
!! * Purpose:       - interface avec une fonction de bas niveau
!! *                  low level string utility from John Burkardt
!! *
!! *                - utilitaires pour trouver un index de mot
!! *                  dans une liste
!! *
!! * Author:        Anthony Michel (anthony.michel@ifp.fr)
!! *
!! * Created:       2005
!! *
!! * Modification:  2005
!! *
!! *
!! ********************************************************************

module M_string

  implicit none
  private

  public :: string_sort  ! tri ordre lexicographique
  public :: find_index
  public :: find_last_index
  public :: string_lowcase ! tout en minuscules
  public :: string_subst   ! substitue char1 by char2

contains

  subroutine ch_low ( c )

    !*******************************************************************************
    !
    !! CH_LOW lowercases a single character.
    !
    !  MODIFIED:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character C, the character to be lowercased.
    !
    implicit none

    character c
    integer itemp

    itemp = ichar ( c )

    if ( 65 <= itemp .and. itemp <= 90 ) then
       c = char ( itemp + 32 )
    end if

    return
  end subroutine ch_low

  subroutine string_subst ( s, c1, c2 )

    !*******************************************************************************
    !
    !! S_LOW replaces all uppercase letters by lowercase ones.
    !
    !  MODIFIED:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be
    !    transformed.  On output, the string is all lowercase.
    !
    implicit none

    integer i
    character ( len = * ) s
    character ( len = 1 ) c1, c2

    do i = 1, len_trim ( s )
       if (s(i:i) == c1)  s(i:i) = c2
    end do

    return
  end subroutine string_subst

  subroutine string_lowcase ( s )

    !*******************************************************************************
    !
    !! S_LOW replaces all uppercase letters by lowercase ones.
    !
    !  MODIFIED:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be
    !    transformed.  On output, the string is all lowercase.
    !
    implicit none

    integer i
    character ( len = * ) s

    do i = 1, len_trim ( s )
       call ch_low ( s(i:i) )
    end do

    return
  end subroutine string_lowcase

  !-----------

  subroutine find_index(listeCONS, nbCONS, nomCONS, iCONS, ierror)

    !----------------
    ! purpose : retrove le premier index iCONS dans une liste
    ! nbfind = 0 si trouve, 1 sinon
    ! iCONS premier index trouve, 0 sinon
    !----------------

    integer :: ierror
    integer :: nbCONS
    character(len=*), dimension(nbCONS)::  listeCONS
    character(len=*)::  nomCONS
    integer :: iCONS

    !--- local
    integer :: i
    !---instructions
    iCONS=0
    ierror=1
    !!write(*,*) "nbCONS=", nbCONS
    do i=1, nbCONS
       if (listecons(i)==nomCONS) then
          iCONS=i
          ierror=0
          return
       end if
    end do

  end subroutine find_index

  !----

  subroutine find_last_index(listeCONS, nbCONS, nomCONS, iCONS, ierror, nbfind)

    !----------------
    ! purpose : retrove index iCONS dans une liste
    ! ierror = 0 si trouve, 1 si echec
    ! nbfind = nb index trouves
    ! iCONS : dernier index trouve, 0 sinon
    !----------------

    integer, intent(out) :: ierror
    integer :: nbCONS
    character(len=*), dimension(nbCONS)::  listeCONS
    character(len=*)::  nomCONS
    integer :: iCONS
    integer, intent(out), optional :: nbfind

    !--- local
    integer :: i
    integer :: nb

    !---instructions
    iCONS=0
    ierror=1
    nb=0
    do i=1, nbCONS
       if (listecons(i)==nomCONS) then
          iCONS=i
          ierror=0
          nb=nb+1
       end if
    end do

    if (present(nbfind)) nbfind=nb

  end subroutine find_last_index

  !----

  subroutine string_sort ( n, sarray, indx)
    !----------------
    ! purpose : tri chaines dans ordre lexicographique
    !----------------

    implicit none

    integer , parameter :: maxlen =64
    integer, intent(in) :: n
    character(len=*), intent(in) :: sarray(:)
    integer, intent(out) :: indx(:)

    !------- local variables
    character(len=maxlen) :: local_sarray(n)
    integer :: local_indx(n)
    integer :: i


    do i=1, n
       local_sarray(i)=sarray(i)
    end do

    do i=1, n
       local_indx(i)=i
    end do

    if (n>1) call svec_sort_heap_a_index ( n, local_sarray, local_indx )

    do i=1, n
       indx(i)=local_indx(i)
    end do


  end subroutine string_sort

  !-----------

  subroutine svec_sort_heap_a_index ( n, sarray, indx )

    !*******************************************************************************
    !
    !! SVEC_SORT_HEAP_A_INDEX does a case-sensitive indexed heap sort of a vector of strings.
    !
    !  Discussion:
    !
    !    The sorting is not actually carried out.
    !    Rather an index array is created which defines the sorting.
    !    This array may be used to sort or index the array, or to sort or
    !    index related arrays keyed on the original array.
    !
    !    The ASCII collating sequence is used, and case is significant.
    !    This means
    !
    !      A < B < C < .... < Y < Z < a < b < .... < z.
    !
    !    Numbers and other symbols may also occur, and will be sorted according to
    !    the ASCII ordering.
    !
    !  MODIFIED:
    !
    !    27 July 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in SARRAY.
    !
    !    Input, character ( len = * ) SARRAY(N), an array to be sorted.
    !
    !    Output, integer INDX(N), contains the sort index.  The
    !    I-th element of the sorted array is SARRAY ( INDX(I) ).
    !
    implicit none

    integer, parameter :: MAX_CHAR = 255
    integer n

    integer i
    integer indx(n)
    integer indxt
    integer ir
    integer j
    integer l
    character ( len = * ) sarray(n)
    character ( len = MAX_CHAR ) string

    do i = 1, n
       indx(i) = i
    end do

    l = n / 2 + 1
    ir = n

    do

       if ( 1 < l ) then

          l = l - 1
          indxt = indx(l)
          string = sarray(indxt)

       else

          indxt = indx(ir)
          string = sarray(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir == 1 ) then
             indx(1) = indxt
             return
          end if

       end if

       i = l
       j = l + l

       do while ( j <= ir )

          if ( j < ir ) then
             if ( llt ( sarray ( indx(j) ), sarray ( indx(j+1) ) ) ) then
                j = j + 1
             end if
          end if

          if ( llt ( string, sarray ( indx(j) ) ) ) then
             indx(i) = indx(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if

       end do

       indx(i) = indxt

    end do

    return

  end subroutine svec_sort_heap_a_index

  !-----------

end module M_string
