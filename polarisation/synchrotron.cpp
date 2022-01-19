#include "functions.h"  //librairie perso qui contient quelques fonctions utiles
#include<math.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

const double e = 1.602e-19, me = 9.109e-31, eps0 = 1e-9/(36.*M_PI), mu0 = 4.*M_PI*1e-7, c = 299792458;
double kappa = 1.;
double zeta = e*me/(24.*pow(M_PI,4)*eps0*c);

/*
 * Cas de l electron isole
 */
double beta(double y){      //Calcul du facteur beta (v/c) en fonction du facteur de lorentz (y)
	return sqrt(y*y-1.)/y;
}

double a(double B,double alpha, double y){      //Calcul du rayon de courbure de la trajectoire de l electron de pitch alpha, de facteur de lorentz y dans le champ B
	return sqrt(y*y-1)*me*c/(e*B*sin(alpha));
}

double thetagamma(double theta, double y){      
	return sqrt(1+y*y*theta*theta);
}

double eta(double w, double B, double alpha, double thetgam, double y){
	return (w*a(B,alpha,y)*pow(thetgam,3))/(3.*c*pow(y,3));
}

double L_isol(double w, double B, double phi, double alpha, double y){      //Calcul de l emissivite pour une particule isolee
	double theta = phi-alpha;
	double thetgam = thetagamma(theta,y);
	double eta_w = eta(w,B,alpha,thetgam,y);
	double bet = beta(y);

	double factor = (zeta*bet*bet*w*w*pow(thetgam,4))/(B*pow(y,3)*sin(alpha)*sin(alpha));

	double K23 = Bessel_K23(eta_w), K13 = Bessel_K13(eta_w);	
	return factor*(K23*K23+(thetgam*thetgam-1)/(thetgam*thetgam)*K13*K13);
}

/*
 * Pour une population d electrons
 */
double pas_adapte(int d, double x, double x1, double x2, double y1, double y2){     //permet d adapter le nombre de points d integration (donc le pas) selon une loi puissance negative
	double a = (y1-y2)/(1./pow(x1,d)-1./pow(x2,d));
	return a*(1./pow(x,d)-1./pow(x2,d))+y2;
}


double fonction(double w, double B, double phi, double alpha, double y, double p){  //expression de la fonction a integrer selon alpha et gamma pour obternir l emmissivite
	return 0.5*sin(alpha)*kappa*pow(y,-p)*L_isol(w,B,phi,alpha,y);
}

double fy(double* tab, double y){     //fonction(w,phi,alpha,y) mais dont seul 'y' est variable, pour l integration selon 'y'
	return fonction(tab[0],tab[1],tab[2],tab[4],y,tab[5]);
}
	
double fa(double* tab, double alpha){     //fonction(w,phi,alpha,y) mais dont seul 'alpha' est variable, pour l integration selon 'alpha'
	double w = tab[0], y_p = tab[3];
	tab[4] = alpha;
	int n = (int)pas_adapte(6,log10(w),7.,23.,1e6,50.);
	return GaussIntegration2(1.,y_p,n,fy,tab);      //on stocke les variables de fy qui ne seront pas integrees dans tab
}

double L_pop(double w, double B, double phi, double y_p, double p, double* tab){       //Calcul de l emmissivite pour une population d electrons selon une distribution en energie et en pitch angle
    if(phi>M_PI/2.)
        phi = M_PI-phi;
	double dalpha = 1e-3/sqrt((w/(2*M_PI))*1e-10);
	double alpha1 = phi-dalpha, alpha2 = phi+dalpha;
	int n = 30;//(int)pas_adapte(2,log10(w),7.,23.,1e6,50.);
	tab[0] = w ; tab[1] = B ; tab[2] = phi ; tab[3] = y_p; tab[4] = alpha1; tab[5] = p;
	return GaussIntegration2(alpha1,alpha2,n,fa,tab);   //on stocke les variables de fa qui ne seront pas integrees dans tab
}

/*
 * Boost relativiste
 */

double DopplerFactor(double y, double phi){
	double bet = beta(y);	
	return 1./(y*(1-bet*cos(phi)));
}

double phi_rest(double phi, double y){      //calcul de l angle d observation dans le referentiel du photon relativiste
	double d = DopplerFactor(y,phi);
	return asin(d*sin(phi));
}

/*
 * Polarisation
 */

//On calcule les composantes orthogonales et paralleles de l emissivite

double I_orth(double w, double B, double phi, double alpha, double y){
    double theta = phi-alpha;
    double thetgam = thetagamma(theta,y);
    double eta_w = eta(w,B,alpha,thetgam,y);
    double x = 2*eta_w/pow(thetgam,3);
	double coeff = sqrt(3)*e*e*y*sin(alpha)/(8*M_PI*eps0*c);
	return coeff*(Fit_F(x)+G(x));
}

double I_para(double w, double B, double phi, double alpha, double y){
    double theta = phi-alpha;
    double thetgam = thetagamma(theta,y);
    double eta_w = eta(w,B,alpha,thetgam,y);
    double x = 2*eta_w/pow(thetgam,3);
	double coeff = sqrt(3)*e*e*y*sin(alpha)/(8*M_PI*eps0*c);
	return coeff*(Fit_F(x)-G(x));
}

//Calculs pour les parametres de Stokes d apres Westfold and Legg

double Pola_I(double nu, double B, double alpha, double p){
    double f_B = (e*B)/(2.*M_PI*me*c);
    double factor = (kappa*mu0*e*e*c)/(2.*sqrt(2.))*pow(3./2.,p/2.);
    double J = WL_L((p+1.)/2.);
    return factor*0.5*sin(alpha)*pow(f_B*sin(alpha),(p+1.)/2.)*pow(nu,-(p-1.)/2.)*J;
}

double Pola_Q(double nu, double B, double alpha, double p){
    double f_B = (e*B)/(2.*M_PI*me*c);
    double factor = (kappa*mu0*e*e*c)/(2.*sqrt(2.))*pow(3./2.,p/2.);
    double L = WL_L((p+1.)/2.);
    return factor*0.5*sin(alpha)*pow(f_B*sin(alpha),(p+1.)/2.)*pow(nu,-(p-1.)/2.)*L;
}

double Pola_U(double nu, double B, double alpha, double p){
    return 0;
}

double Pola_V(double nu, double B, double alpha, double p){
    double f_B = (e*B)/(2.*M_PI*me*c);
    double factor = (kappa*mu0*e*e*c)/sqrt(3.)*pow(3./2.,p/2.);
    double L = WL_L(0.5*p);
    double J = WL_J(0.5*p);
    double R = WL_J(0.5*p+1.);
    return factor*0.5*cos(alpha)*pow(f_B*sin(alpha),0.5*p+1.)*pow(nu,-0.5*p)*(R+3*(L-0.5*J));
}
