#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

/**
 * rk1 ODE solver, puts all of the values it calculates into a file named rk1.txt
 * 
 * Parameters:
 * f: function respresenting the ODE
 * order: order of the ODE
 * y: an array of ivps starting from lowest order solution to hightest order
 * t0: initial time of the ivps
 * h: step size
 * n: number of iterations
 */
void rk1 (double* (*f)(double, double*), int order, double* y, double t0, double h, int n) {
   ofstream out_file ("rk1.txt");
   double t = t0;
   double* y_calc = new double[order];

   out_file<<"t";
   for (int i = 0; i < order; ++i) {
      y_calc[i] = y[i];
      out_file<<", y"<<i;
   }

   out_file<<endl;

   for (int i = 0; i < n; i++) {
      double* k1 = f(t, y_calc);
      for (int j = 0; j < order; j++) {
         y_calc[j] += h*k1[j];
      }

      t += h;
      out_file<<t;
      for (int j = 0; j < order; j++) {
         out_file<<", "<<y_calc[j];
      }

      out_file<<endl;

      delete[] k1;
   }

   delete[] y_calc;
   out_file.close();
}

void rk2(double* (*f)(double, double*), int order, double* y, double t0, double h, int n) {
   ofstream out_file("rk2.txt");

   double t = t0;
   double* y_calc = new double[order];

   out_file<<"t";
   for (int i = 0; i < order; i++) {
      y_calc[i] = y[i];
      out_file<<", y"<<i;
   }

   out_file<<endl;

   for (int i = 0; i < n; i++) {
      double* k1 = f(t, y_calc);
      double* ypk1 = new double[order];

      for (int j = 0; j < order; j++) {
         ypk1[j] = h*k1[j] + y_calc[j];
      }

      double* k2 = f(t + h, ypk1);

      for (int j = 0; j < order; j++) {
         y_calc[j] += 0.5 * h * (k1[j] + k2[j]);
      }

      t += h;
      out_file<<t;
      for (int j = 0; j < order; j++) {
         out_file<<", "<<y_calc[j];
      }

      out_file<<endl;

      delete[] k1;
      delete[] ypk1;
      delete[] k2;
   }

   delete[] y_calc;
   out_file.close();
}

void rk3(double* (*f)(double, double*), int order, double* y, double t0, double h, int n) {
   ofstream out_file("rk3.txt");
   
   double t = t0;
   double* y_calc = new double[order];

   out_file<<"t";
   for (int i = 0; i < order; i++) {
      out_file<<", y"<<i;
      y_calc[i] = y[i];
   }

   out_file<<endl;

   for (int i = 0; i < n; i++) {
      double* k1 = f(t, y_calc);
      double* ypk1 = new double[order];
      for (int j = 0; j < order; j++) {
         ypk1[j] = y_calc[j] + 0.5*h*k1[j];
      }

      double* k2 = f(t + 0.5*h, ypk1);

      double* ymk1pk2 = new double[order];
      for (int j = 0; j < order; j++) {
         ymk1pk2[j] = y_calc[j] - h*k1[j] + 2*h*k2[j];
      }

      double* k3 = f(t + h, ymk1pk2);

      for (int j = 0; j < order; j++) {
         y_calc[j] += h*(k1[j] + 4*k2[j] + k3[j])/6.0;
      }

      t += h;

      out_file<<t;
      for (int j = 0; j < order; j++) {
         out_file<<", "<<y_calc[j];
      }

      out_file<<endl;

      delete[] k1;
      delete[] ypk1;
      delete[] k2;
      delete[] ymk1pk2;
      delete[] k3;
   }

   delete[] y_calc;
   out_file.close();
}

void rk4(double* (*f)(double, double*), int order, double* y, double t0, double h, int n) {
   ofstream out_file("rk4.txt");

   double t = t0;
   double* y_calc = new double[order];

   out_file<<"t";
   for (int i = 0; i < order; i++) {
      out_file<<", y"<<i;
      y_calc[i] = y[i];
   }

   out_file<<endl;

   for (int i = 0; i < n; i++) {
      double* k1 = f(t, y_calc);
      double* ypk1 = new double[order];
      for (int j = 0; j < order; j++) {
         ypk1[j] = y_calc[j] + 0.5*h*k1[j];
      }

      double* k2 = f(t + 0.5*h, ypk1);
      double* ypk2 = new double[order];
      for (int j = 0; j < order; j++) {
         ypk2[j] = y_calc[j] + 0.5*h*k2[j];
      }

      double* k3 = f(t + 0.5*h, ypk2);
      double* ypk3 = new double[order];
      for (int j = 0; j < order; j++) {
         ypk3[j] = y_calc[j] + h*k3[j];
      }

      double* k4 = f(t + h, ypk3);

      for (int j = 0; j < order; j++) {
         y_calc[j] += h*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6.0;
      }
   
      t += h;
      out_file<<t;
      for (int j = 0; j < order; j++) {
         out_file<<", "<<y_calc[j];
      }

      out_file<<endl;
      delete[] k1;
      delete[] ypk1;
      delete[] k2;
      delete[] ypk2;
      delete[] k3;
      delete[] ypk3;
      delete[] k4;
   }

   delete[] y_calc;
   out_file.close();
}
