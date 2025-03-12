program main
   real :: x
   read *, x
   print *, rich(x, 0.5)

   contains
   real function f(x) result (y)
      real, intent(in) :: x
      y = x**2
   end function f

   real function D(x, h) result (dy) 
      real, intent(in) :: x, h
      if (h .le. 0.0) then
         stop "illegal step size"
      else
         dy = (f(x + h) - f(x - h)) / (2*h)
      endif
   end function D

   real function rich (x, h) result (dy)
      real, intent(in) :: x, h
      if (h .le. 0.0) then
         stop "illegal step size"
      else
         dy = (4*D(x, h / 2) - D(x, h)) / 3
      endif
   end function rich
end
