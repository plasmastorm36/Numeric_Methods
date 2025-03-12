module rk
   abstract interface
      function ODE_Function(t, y) result(dy) ! this is being treated as a system to better solve
         implicit none
         real(8), intent(in) :: t ! t value
         real(8), intent(in) :: y(:) ! initial y value (being ordered in order of derivatives)
         real(8) :: dy(size(y)) ! how y changes
      end function ODE_Function
   end interface

   contains
      subroutine rk1(f, y, t0, h, n)
         implicit none
         procedure(ODE_Function) :: f ! ODE system
         real(8), intent(inout) :: y(:) ! initial values
         real(8), intent(in) :: t0, h ! initial time and step size
         integer, intent(in) :: n ! number of iterations

         integer :: i
         real(8) :: t, dy(size(y))
         open(11, file = "rk1.txt", status="replace") ! i want to overwrite/erase the previous entry
         write(11, "(A, *(A,I0))") "t ", ("y", i, i = 1, size(y)) ! column headers

         t = t0

         do i = 1, n
            dy = f(t, y)
            y = y + h * dy
            t = t + h
            write(11, *) t, y
         end do

         close(11)
      end subroutine rk1

      subroutine rk2(f, y, t0, h, n)
         implicit none
         procedure(ODE_Function) :: f ! ode system
         real(8), intent(inout) :: y(:) ! initial values
         real(8), intent(in) :: t0, h ! initial time and step size
         integer, intent(in) :: n ! number of iterations

         integer :: i
         real(8) :: t, k1(size(y)), k2(size(y))
         open(12, file = "rk2.txt", status="replace")
         write(12, "(A, *(A,I0))") "t ", ("y", i, i = 1, size(y)) ! column headers

         t = t0
         do i = 1, n
            k1 = f(t, y)
            k2 = f(t + 0.5*h, y + 0.5*k1*h)
            y = y +  (0.5 * k1 + 0.5*k2)*h
            t = t + h
            write(12, *) t, y
         end do
      
         close(12)
      end subroutine rk2

      subroutine rk3(f, y, t0, h, n)
         implicit none
         procedure(ODE_Function) :: f ! ode system
         real(8), intent(inout) :: y(:) ! initial values
         real(8), intent(in) :: t0, h ! initial time and step size
         integer, intent(in) :: n ! number of iterations

         integer :: i
         real(8) :: t, k1(size(y)), k2(size(y)), k3(size(y))
         open(13, file = "rk3.txt", status = "replace")
         write(13, "(A, *(A,I0))") "t ", ("y", i, i = 1, size(y)) ! column headers

         t = t0
         do i = 1, n
            k1 = h*f(t, y)
            k2 = h*f(t + 0.5*h, y + 0.5*k1)
            k3 = h*f(t + h, y + h*f(t + h, y + k1))
            y = y + (1.0/6.0)*(k1 + 4*k2 + k3)
            t = t + h
            write(13, *) t, y
         end do
         close(13)
      end subroutine rk3

      subroutine rk4(f, y, t0, h, n)
         implicit none
         procedure(ODE_Function) :: f ! ode system
         real(8), intent(inout) :: y(:) ! initial values
         real(8), intent(in) :: t0, h ! initial time and step size
         integer, intent(in) :: n ! number of iterations

         integer :: i
         real(8) :: t, k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y))
         open(14, file = "rk4.txt", status="replace")
         write(14, "(A, *(A,I0))") "t ", ("y", i, i = 1, size(y))

         t = t0
         do i = 1, n
            k1 = h * f(t, y)
            k2 = h * f(t + 0.5*h, y + 0.5*k1)
            k3 = h * f(t + 0.5*h, y + 0.5*k2)
            k4 = h * f(t + h, y + k3)
            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
            t = t + h
            write(14, *) t, y
         end do
         close(14)
      end subroutine rk4

end module rk
