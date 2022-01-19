#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H

double L_isol(double,double,double,double,double);
double L_pop(double,double,double,double,double,double*);
double DopplerFactor(double,double);
double phi_rest(double,double);
double I_orth(double,double,double,double,double);
double I_para(double,double,double,double,double);
double Pola_I(double,double,double,double);
double Pola_Q(double,double,double,double);
double Pola_U(double,double,double,double);
double Pola_V(double,double,double,double);

#endif
