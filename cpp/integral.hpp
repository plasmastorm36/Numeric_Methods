#include<iostream>
using namespace std;

double midpoint (double (*f)(double), double a, double b, int n) {
   double del = (b - a) / double(n);
   double x = (a + del) / 2;
   double result = 0.0;
   for (; x < b; x += del) {
      result += f(x);
   }

   return del * result;
}

double trapezoid (double (*f)(double), double a, double b, int n) {
   double del = (b - a) / n;
   double x = a + del;
   double result = (f(a) + f(b)) / 2.0;
   for (; x < b; x += del) {
      result += f(x);
   }

   return del * result;
}

double simpson (double (*f)(double), double a, double b, int n) {
   if (n < 2) {
      cerr << "n must be at least 2"<<endl;
      exit(EXIT_FAILURE);
   } else if (n % 2 != 0) {
      cerr<<"n must be even"<<endl;
      exit(EXIT_FAILURE);
   }

   double del = (b - a) / double(n);
   double x = a + del;
   double result = f(a) + f(b);

   for (int i = 1; i < n; i++) {
      i % 2 == 0 ? result += 2 * f(x) : result += 4 * f(x);
      x += del;
   }

   return del * result / 3.0;
}
