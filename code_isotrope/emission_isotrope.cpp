/* Calcule et trace les flux en fonction de la longueur d'onde et de la fréquence pour emission isotrope d'un electron isolé et d'une population d'électrons suivant une loi de puissance */

#include "functions.h"
#include<math.h>
#include<fstream>
#include<iostream>
using namespace std;

const double e = 1.602e-19, m = 9.109e-31, eps0 = 1e-9/(36.*M_PI), c = 299792458;
double B = 1e-4, gam = 1e4, pitch = M_PI/4, kappa=1., p=2.0;
double gam2 = gam*gam;
double wg = e*B/m;
double wc = (3./2.)*gam*gam2/sqrt(gam2-1)*wg*sin(pitch);
double a = sqrt(gam2-1)*c/(wg*sin(pitch));
double Tr = 2*M_PI*gam/wg;

double gammafact = (tgamma(p/4.+19./12.)*tgamma(p/4.-1./12.)*tgamma(p/4.+5./4.))/tgamma(p/4.+7./4.);

/*
 * Cas de l electron isole
 */

//On calcule les composantes orthogonales et paralleles de l emissivite

double int_orth(double w, double theta){        //fonction a integrer pour obtenir la composante orthogonale a B de l emissivite
	double thetgam = sqrt(1+gam2*theta*theta);
	double thetgam2 = thetgam*thetgam;
	double eta = w*a*thetgam*thetgam2/(3*c*gam*gam2);
	double K = Fit_Bessel_K23(eta);
	return thetgam2*thetgam2*K*K;
}

double I_orth(double w){
	double b = 500.;
	int n = 5e3;
	double coeff = w*w*e*e*a*a*sin(pitch)/(6*M_PI*M_PI*eps0*c*c*c*gam2*gam2);
	double* m = &w;
	//return coeff*SimpsonIntegral2(-b,b,n,int_orth,m);
	double coeff2 = sqrt(3)*e*e*gam*sin(pitch)/(8*M_PI*eps0*c);
	return coeff2*(Fit_F(w/wc)+G(w/wc));    //on utilise ici la definition avec les fonctions synchrotron de WL
}

double int_para(double w, double theta){        //fonction a integrer pour obtenir la composante parallele a B de l emissivite
	double thetgam = sqrt(1+gam2*theta*theta);
	double thetgam2 = thetgam*thetgam;
	double eta = w*a*thetgam*thetgam2/(3*c*gam*gam2);
	double K = Fit_Bessel_K13(eta);
	return thetgam2*theta*theta*K*K;
}

double I_para(double w){
	double b = 500.;
	int n = 5e3;
	double coeff = w*w*e*e*a*a*sin(pitch)/(6*M_PI*M_PI*eps0*c*c*c*gam2);
	double* m = &w;
	//return coeff*SimpsonIntegral2(-b,b,n,int_para,w);
	double coeff2 = sqrt(3)*e*e*gam*sin(pitch)/(8*M_PI*eps0*c);
	return coeff2*(Fit_F(w/wc)-G(w/wc));    //definition avec lels fonctions synchrotron de WL
}


//On deduit l emissivite totale

double J_isol(double w){        //Emissivite isotrope pour un electron isole, calcule en fonction des composantes orthogonale et parallele
	return (I_orth(w)+I_para(w))/Tr;
}

double L_isol(double w){    //autre expression en utilisant la fonction synchrotron
	double x = w/wc;
	double coeff = sqrt(3)*e*e*e*B*sin(pitch)/(8*M_PI*M_PI*eps0*c*m);
	return coeff*Fit_F(x);
}

/*
 * Pour une population d electrons
 */

double L_pop(double w){     //Emissivite isotrope pour une population d electrons de distribution en energie dE = kappa*E^{-p}
	double coeff = (sqrt(3)*e*e*e*B*kappa*sqrt(M_PI))/(16*M_PI*M_PI*eps0*c*m*(p+1.));
	double depw = pow((w*m*m*m*c*c*c*c)/(3*e*B),-(p-1.)/2.);
	return coeff*gammafact*depw;
	//return pow(B,(p+1.)/2.)*pow(w,-(p-1.)/2.);
}

double pola(double w){      //Calcul du degre de polarisation du rayonnement emis
	return (I_orth(w)-I_para(w))/(I_orth(w)+I_para(w));
}

int main(){
	double w = 1e9;
	fstream data("data_iso.res",ios::out);
	int i, n=5000;
	double h = (17.-9.)/n;
	for(i=0;i<n;i++){
		data << pow(10,9.+h*i) << " " << L_isol(2.*M_PI*pow(10,9.+h*i)) << " " << L_pop(2.*M_PI*pow(10,9.+h*i)) << endl;    //stockage des valeurs obtenues pour une frequence de 10^{9+h*i} (les calculs d emissivite sont en fonction de la pulsation) dans le fichier data_iso.res
	}
	data.close();
	system("python affichage_iso.py &");      //affichage des donnees stockees dans data_iso.res avec python
	return 0;
}
