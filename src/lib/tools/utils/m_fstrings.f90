!*******************************************************************************
! fstrings.for -- routines for converting between FORTRAN and C character
!                 strings.
!
! Mark Showalter, PDS Rings Node, September 2002
!*******************************************************************************

module M_FStrings

  public :: Fort_CString
  public :: Fort_FString

contains

  !
  !*******************************************************************************
  !$ Component_name:
  !	FORT_CSTRING (fstrings.for)
  !$ Abstract:
  !	Converts a FORTRAN character string to a null-terminated byte array,
  !	for passage to a C function.
  !$ Keywords:
  !	UTILITY, FORTRAN_C
  !	FORTRAN, INTERNAL, subroutine
  !$ Declarations:
  !	subroutine	FORT_CSTRING(string, array, nbytes)
  !	character*(*)	string
  !	integer*1	array(*)
  !	integer*4	nbytes
  !$ Inputs:
  !	string		character string to convert.
  !	nbytes		dimensioned length of byte array.
  !$ Outputs:
  !c	array(1...)	string of bytes with terminal null.
  !$ Returns:
  !	none
  !$ Detailed_description:
  !	This subroutine converts a FORTRAN character string to a null-terminated
  !	byte array, for passage to a C function.  Blank characters at the end of
  !	the character string are not considered significant.  The string is
  !	truncated if necessary to fit into the array.
  !$ External_references:
  !	none
  !$ Examples:
  !	none
  !$ Error_handling:
  !	none
  !$ Limitations:
  !	The dimensioned length of the byte array must be at least one greater
  !	than the effective length of the character string.
  !$ Author_and_institution:
  !	Mark R. Showalter
  !	PDS Rings Node, NASA/Ames Research Center
  !$ Version_and_date:
  !	1.0: January 1994
  !	1.1: September 2002
  !$ Change_history:
  !	1.1: MODIFIED for compatibility with Absoft FORTRAN for Macintosh OS X.
  !*******************************************************************************

  subroutine Fort_CString(string, array, nbytes)
    character*(*), intent(in)    :: string
    integer(kind=4), intent(in)  :: nbytes
    integer(kind=1), intent(out) :: array(*)

    integer :: last, i

    ! Search for the last character actually used.
    do last = len(string), 1, -1
      if (string(last:last) .NE. ' ') exit
    end do

    ! Truncate string if necessary
    if (last .GT. nbytes-1) last = nbytes-1

    ! Copy bytes from character string
    do i = 1, last
      array(i) = ichar( string(i:i) )
    end do

    ! Append null terminator
    array(last+1) = 0

    return
  end subroutine Fort_CString

  !
  !*******************************************************************************
  !$ Component_name:
  !	FORT_FSTRING (fstrings.for)
  !$ Abstract:
  !	Converts a null-terminated byte array returned from a C function to a
  !	FORTRAN character string.
  !$ Keywords:
  !	UTILITY, FORTRAN_C
  !	FORTRAN, INTERNAL, subroutine
  !$ Declarations:
  !	subroutine	FORT_FSTRING(array, string)
  !	integer*1	array(*)
  !	character*(*)	string
  !$ Inputs:
  !	array(1...)	string of bytes with terminal null.
  !$ Outputs:
  !	string		FORTRAN character string.
  !$ Returns:
  !	none
  !$ Detailed_description:
  !	This subroutine converts a null-terminated byte array returned from a C
  !	function to a FORTRAN character string.  The string is truncated if
  !	necessary.
  !$ External_references:
  !	none
  !$ Examples:
  !	none
  !$ Error_handling:
  !	none
  !$ Limitations:
  !	none
  !$ Author_and_institution:
  !	Mark R. Showalter
  !	PDS Rings Node, NASA/Ames Research Center
  !$ Version_and_date:
  !	1.0: January 1994
  !$ Change_history:
  !	none
  !*******************************************************************************

  subroutine Fort_FString(array, string)
    integer(kind=1), intent(in) :: array(*)
    character*(*), intent(out) :: string

    integer :: i

    ! Copy bytes into string, one by one
    do  i = 1, len(string)
      if (array(i) .EQ. 0) then
        ! Pad remainder of string with blanks
        string(i:) = ' '
        return
      end if
      string(i:i) = char(array(i))
    end do
    return

  end subroutine Fort_FString
  !*******************************************************************************
  !
end module M_FStrings
