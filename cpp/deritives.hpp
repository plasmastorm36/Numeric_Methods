#include <iostream>
/**
 * differential method using a backwards function
 * 
 * @param f pointer to any function that returns a double
 * @param x point at which you want to find the derivative
 * @param h step size of the numerical method
 */
double back_difference (double (*f)(double), double x, double h) {
   if (h <= 0.0) {
      perror("step size too small");
      exit(EXIT_FAILURE);
   }

   return (f(x) - f(x - h)) / h;
}

double forward_difference (double (*f)(double), double x, double h) {
   if (h <= 0.0) {
      perror("step size too small");
      exit(EXIT_FAILURE);
   }

   return (f(x + h) - f(x)) / h;
}

double central_difference (double (*f)(double), double x, double h) {
   if (h < 0.0) {
      perror("step size too small");
      exit(EXIT_FAILURE);
   }

   return (f(x + h) - f(x - h)) / (2*h);
}

double richard (double (*f)(double), double x, double h) {
   if (h < 0.0) {
      perror("step size too small");
      exit(EXIT_FAILURE);
   }

   return (4*central_difference(f, x, h/2) - central_difference(f, x, h)) / 3;
}
