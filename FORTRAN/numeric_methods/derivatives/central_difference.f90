program main
   real :: x
   read *, x
   print *, central_difference(x, 0.5)
   contains
      real function f(x) result(y)
         real, intent(in) :: x
         ! since I declared the function as real the return is assumed to be real
         y = x**2
      end function f

      real function central_difference(x, h) result (y)
         real, intent(in) :: x, h
         if (h .le. 0.0) then
            STOP
         else
            y = (f(x + h) - f(x - h)) / (2*h)
         endif
      end function central_difference
end
