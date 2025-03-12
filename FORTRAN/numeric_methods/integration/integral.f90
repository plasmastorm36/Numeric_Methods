module integral
   implicit none
   contains
      real function midpoint (f, a, b, n) result (anti)
         interface
               real function f(x) result (y)
                  real, intent(in) :: x
            end function f
         end interface
         real, intent(in) :: a, b ! start and endpoint of an integral
         integer, intent(in) :: n ! how many steps to be used
         real :: del, x
         integer :: i
         if (n .lt. 2) then
            STOP "too little steps"
         endif
         del = (b - a) / real(n)
         anti = 0
         x = a + del / 2.0

         do i = 1, n
            anti = anti + f(x) * del
            x = x + del
         end do
      end function midpoint

      real function trapezoid (f, a, b, n) result(anti)
         interface
               real function f(x) result (y)
                  real, intent(in) :: x
            end function f
         end interface
         real, intent(in) :: a, b
         integer, intent(in) :: n
         real :: del, x
         integer :: i

         if (n .lt. 2) then
            stop "not enough steps"
         endif
         del = (b - a) / real(n)
         anti = (f(a) + f(b)) / 2.0
         x = (a + del)
         do i = 1, n - 1
            anti = anti + f(x)
            x = x + del
         end do
         anti = anti * del
      end function trapezoid

      real function simpson (f, a, b, n) result (anti)
         interface
               real function f(x) result (y)
                  real, intent(in) :: x
            end function f
         end interface
         real, intent(in) :: a, b
         integer, intent(in) :: n
         real :: del, x
         integer :: i
         if (n .lt. 2) then ! there must be at least 4 points for simpson's rule to work
            STOP "n is too small"
         else if (mod(n, 2) .ne. 0) then ! there must be even points for simpsons rule to work
            STOP "n is odd, simpsons rule needs an even amount of points to work"
         endif
         del = (b - a) / real(n)
         anti = f(a) + f(b)
         x = a + del
         do i = 1, n - 1
            if (mod(i, 2) .eq. 0) then
               anti = anti + (2.0 * f(x))
            else
               anti = anti + (4.0 * f(x))
            endif
            x = x + del
         end do

         anti = anti * (del / 3)
      end function simpson
end module integral