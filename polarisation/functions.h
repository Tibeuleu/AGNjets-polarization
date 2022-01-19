#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double angle_princ(double);
double sinc(double);
double Bessel_I(double,double);
double Bessel_K(double,double);
double Bessel_K53(double);
double Bessel_K23(double);
double Bessel_K13(double);
double Fit_Bessel_K53(double);
double Fit_Bessel_K23(double);
double Fit_Bessel_K13(double);
double WL_L(double);
double WL_J(double);
double WL_R(double);
double SimpsonIntegration(double,double,int,double(*f)(double));
double SimpsonIntegration2(double,double,int,double(*f)(double*,double),double*);
double TrapezeIntegration(double,double,int,double(*f)(double));
double TrapezeIntegration2(double,double,int,double(*f)(double*,double),double*);
double GaussIntegration(double,double,int,double(*f)(double));
double GaussIntegration2(double,double,int,double(*f)(double*,double),double*);
double F(double);
double Fit_F(double);
double G(double);

#endif
