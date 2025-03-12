function equation (x) result (y)
   real :: x
   real :: y
   y = x**2
end function equation

function back_difference(x, h) result (y)
   real, intent(in) :: x, h
   real :: y
   if (h .le. 0.0)  then
      stop
   else
      y = (equation(x) - equation(x - h)) / h
   endif
end function back_difference

program main
   real :: x
   read *, x
   print *, back_difference(x, 0.00005)
end