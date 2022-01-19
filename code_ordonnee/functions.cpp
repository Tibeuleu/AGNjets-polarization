/* Librairie de fonctions simple utilisées dans plusieurs programmes */

#include "functions.h"
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "boost/math/special_functions/bessel.hpp"
using namespace std;

double angle_princ(double theta){
    while(fabs(theta) > 2.*M_PI){
        if(theta > 0)
            theta -= 2.*M_PI;
        else
            theta += 2.*M_PI;
    }
    if(theta < 0)
        theta += 2.*M_PI;
    return theta;
}

double sinc(double x){
    if(fabs(x) < 1e-5)
        return 1.;
    return sin(x)/x;
}

/*
 * Fonctions de Bessel modifiées
 */


// en utilisant les fonctions bessel de la librairie boost, le fit est prefere ensuite

double Bessel_I(double alpha, double x){        //Simplification de l appel des fonctions pour le code principal
	return boost::math::cyl_bessel_i((double)alpha,(double)x);
}
double Bessel_K(double alpha, double x){
	return boost::math::cyl_bessel_k((double)alpha,(double)x);
}


double Bessel_K53(double x){
	return Bessel_K(5./3.,x);
}
double Bessel_K23(double x){
	return Bessel_K(2./3.,x);
}
double Bessel_K13(double x){
	return Bessel_K(1./3.,x);
}


double AsymBessel_1(double nu, double x){       //Un fit des fonctions de Bessel modifiees mais qui ne me permet pas de calculer Bessel_K13
	return 0.5*tgamma(nu)*pow(0.5*x,-nu);
}
double AsymBessel_2(double nu, double x){
	return sqrt(M_PI/(2.*x))*exp(-x);
}

double Fit_Bessel_K53(double x){            //Voir "Analytical Fits to the Synchrotron Functions", de M. Fouka, S. Ouichaoui sur ResearchGate
	int i;
	double d1,d2, h1=0, h2=0;
	double** a = (double**)malloc(2*sizeof(double*));
	a[0] = (double*)malloc(3*sizeof(double));
	a[0][0] = -1.0194198041210243;
	a[0][1] = 0.28011396300530672;
	a[0][2] = -7.71058491739234908e-2;
	a[1] = (double*)malloc(3*sizeof(double));
	a[1][0] = -15.761577796582387;
	a[1][1] = 0.;
	a[1][2] = 0.;
	for(i=0;i<3;i++){
		h1 += a[0][i]*pow(x,1./(i+1.));
		h2 += a[1][i]*pow(x,1./(i+1.));
	}
	free(a[0]); free(a[1]); free(a);
	d1 = exp(h1);
	d2 = 1-exp(h2);
	return d1*AsymBessel_1(5./3.,x)+d2*AsymBessel_2(5./3.,x);
}
double Fit_Bessel_K23(double x){
	int i;
	double d1,d2, h1=0, h2=0;
	double** a = (double**)malloc(2*sizeof(double*));
	a[0] = (double*)malloc(4*sizeof(double));
	a[0][0] = -1.0010216415582440;
	a[0][1] = 0.88350305221249859;
	a[0][2] = -3.6240174463901829;
	a[0][3] = 0.57393980442916881;
	a[1] = (double*)malloc(4*sizeof(double));
	a[1][0] = -0.2493940736333186;
	a[1][1] = 0.9122693061687756;
	a[1][2] = 1.2051408667145216;
	a[1][3] = -5.5227048291651126;
	for(i=0;i<4;i++){
		h1 += a[0][i]*pow(x,1./(i+1.));
		h2 += a[1][i]*pow(x,1./(i+1.));
	}
	free(a[0]); free(a[1]); free(a);
	d1 = exp(h1);
	d2 = 1-exp(h2);
	return d1*AsymBessel_1(2./3.,x)+d2*AsymBessel_2(2./3.,x);
}
double Fit_Bessel_K13(double x){            //Attention formule valable seulement pour x >> 1
	return Fit_Bessel_K53(x)-4./(3.*x)*Fit_Bessel_K23(x);
}

/*
 * Fonctions Synchrotron
 */

double F(double x){         //Fonction synchrotron introduite par Westfold and Legg (1959)
	double b = 100;
	int n = 8e4;
	return x*SimpsonIntegration(x,b,n,Fit_Bessel_K53);  //Fait ici intervenir une intégration, on preferera le fit trouve dans le meme papier que precedent
}
//un fit de la fonction synchrotron
double AsymF_1(double x){
	double coeff = M_PI*pow(2.,5./3.)/(sqrt(3)*tgamma(1./3.));
	return coeff*pow(x,1./3.);
}
double AsymF_2(double x){
	double coeff = sqrt(0.5*M_PI);
	return coeff*sqrt(x)*exp(-x);
}
double Fit_F(double x){
	int i;
	double d1,d2, h1=0, h2=0;
	double** a = (double**)malloc(2*sizeof(double*));
	a[0] = (double*)malloc(3*sizeof(double));
	a[0][0] = -0.97947838884478688;
	a[0][1] = -0.83333239129525072;
	a[0][2] = 0.15541796026816246;
	a[1] = (double*)malloc(3*sizeof(double));
	a[1][0] = -4.69247165562628882e-2;
	a[1][1] = -0.70055018056462881;
	a[1][2] = 1.03876297841949544e-2;
	for(i=0;i<3;i++){
		h1 += a[0][i]*pow(x,1./(i+1.));
		h2 += a[1][i]*pow(x,1./(i+1.));
	}
	free(a[0]); free(a[1]); free(a);
	d1 = exp(h1);
	d2 = 1-exp(h2);
	return d1*AsymF_1(x)+d2*AsymF_2(x);
}

double G(double x){     //Deuxieme fonction synchrotron
	return x*Bessel_K23(x);
}

/*
 * Fonctions définies sur les fonctions de Bessel pour les calculs des paramètres de Stokes de l'émission synchrotron
 */

double WL_L(double n){
    if(n<=2./3.)
        return 0;
    else
        return pow(2,n-2)*tgamma(0.5*n-1./3.)*tgamma(0.5*n+1./3.);
}

double WL_J(double n){
    if(n<=2./3.)
        return 0;
    else
        return (2./3.+n)/n*WL_L(n);
}

double WL_R(double n){
    if(n<=1./3.)
        return 0;
    else
        return pow(2,n-2)*tgamma(0.5*n-1./6.)*tgamma(0.5*n+1./6.);
}

/*
 * Methodes d'intégration
 */

double SimpsonIntegration(double a, double b, int n, double(*f)(double)){
	double h = (b-a)/n, I = 0, x1, x2;
	int k;
	for(k=0;k<n;k++){
		x1 = a + k*h;
		x2 = a + (k+1)*h;
		I += (x2-x1)/6.0*(f(x1)+4.0*f(0.5*(x1+x2))+f(x2));
	}
	return I;
}

double SimpsonIntegration2(double a, double b, int n, double(*f)(double*,double), double* fix){     //Pour les fonctions de plusieurs variables on stocke celle sur lesquels on integre pas dans le pointeur 'fix'
	double h = (b-a)/n, I = 0, x1, x2;
	int k;
	for(k=0;k<n;k++){
		x1 = a + k*h;
		x2 = a + (k+1)*h;
		I += (x2-x1)/6.0*(f(fix,x1)+4.0*f(fix,0.5*(x1+x2))+f(fix,x2));
	}
	return I;
}

double TrapezeIntegration(double a, double b, int n, double(*f)(double)){
	double h = (b-a)/n, I = 0, x1, x2;
	int k;
	for(k=0;k<n;k++){
		x1 = a+k*h;
		x2 = a+(k+1)*h;
		I += 0.5*(x2-x1)*(f(x1)+f(x2));
	}
	return I;
}

double TrapezeIntegration2(double a, double b, int n, double(*f)(double*,double), double* fix){
	double h = (b-a)/n, I = 0, x1, x2;
	int k;
	for(k=0;k<n;k++){
		x1 = a+k*h;
		x2 = a+(k+1)*h;
		I += 0.5*(x2-x1)*(f(fix,x1)+f(fix,x2));
	}
	return I;
}

double GaussIntegration(double a, double b, int n, double(*f)(double)){
	double h = (b-a)/n, I = 0, x1, xl, xr;
	int k;
	for(k=0;k<n;k++){
		x1 = a+k*h;
		xl = x1+(0.5-sqrt(3.)/6.)*h;
		xr = x1+(0.5+sqrt(3.)/6.)*h;
		I += f(xl)+f(xr);
	}
	return 0.5*h*I;
}

double GaussIntegration2(double a, double b, int n, double(*f)(double*,double), double* fix){
	double h = (b-a)/n, I = 0, x1, xl, xr;
	int k;
	for(k=0;k<n;k++){
		x1 = a+k*h;
		xl = x1+(0.5-sqrt(3.)/6.)*h;
		xr = x1+(0.5+sqrt(3.)/6.)*h;
		I += f(fix,xl)+f(fix,xr);
	}
	return 0.5*h*I;
}
