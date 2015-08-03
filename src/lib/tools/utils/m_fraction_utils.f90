module M_Fraction_Utils
  
   implicit none
   private

   public :: fraction_gcd
   public :: fraction_ppcm
   public :: fraction_reduce
   
 contains 
   
   function fraction_ppcm(a,b) result (p)
     implicit none
     integer, intent(in) :: a, b
     integer :: p
     !--
     p = a*b/gcd(a,b)

   end function fraction_ppcm

   !---
   
   function fraction_gcd(a,b) result (d)
     implicit none
     integer, intent(in) :: a, b
     integer :: d
     !--
     d = gcd(a,b)
     
   end function fraction_gcd

   !---

   subroutine fraction_reduce(a,b)
     implicit none
     integer, intent(inout) :: a, b
     !--
     integer :: d
     d = gcd(a,b)
     a = a/d
     b = b/d
   end subroutine fraction_reduce
   
   recursive function gcd( a, b ) result( divisor )

      ! Note that the function itself is not declared when the result
      ! variable is present.  The type of the function is the type of 
      ! the result variable.  Thus, only the result variable may be 
      ! declared. 

      integer divisor
      integer, intent(in) :: a, b

      integer m, n

      ! Multiple statements may be written on a single source line
      ! provided they are delimited with semicolons.

      m = abs(a); n = abs(b)
      if ( m > n ) call swap( m, n )  ! Insure that m <= n.

      ! When the function invokes itself recursively, the result variable
      ! should be used to store the result of the function.  The function
      ! name is used to invoke the function.  Thus, the function name should
      ! not appear on the left-hand side of an assignment statement.

      if ( m == 0 ) then 
         divisor = n
      else
         divisor = gcd( mod( n, m ), m )
      end if

      ! Unlike internal subprograms, module procedures, such as gcd, may 
      ! have internal subprograms defined within them (to one level only). 
      ! As with module procedures, internal subprograms also have an
      ! explicit interface.  Thus, the swap subroutine is not declared in 
      ! the gcd function above---its interface is explicit.


      contains

      ! The swap subroutine is internal to gcd and is not visible
      ! elsewhere, not even in this module.

      subroutine swap( x, y )
         integer, intent(inout) :: x, y

         integer tmp

         tmp = x; x = y; y = tmp 
      end subroutine swap

   end function gcd

 end module M_Fraction_Utils
