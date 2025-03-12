program main
   implicit none

    real :: x

   read *, x

   print *, forward_difference(x, 0.01)

   contains
   real function f(x) result (y)
      real, intent(in) :: x
      y = x**2
   end function f

   real function forward_difference(x, h) result (y)
      real, intent(in) :: x, h
      if (h .le. 0.0) then
         STOP
      else
         y = (f(x+h) - f(x)) / h
      endif
   end function forward_difference
end