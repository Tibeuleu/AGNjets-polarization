/* Calcule et trace les flux en fonction de l'angle d'observation pour emission ordonnée et boostée d'une population d'électrons suivant une loi de puissance, pour une valeur fixée de fréquence, d'intensité de champ magnétique et pour plusieurs valeurs de boost : en faisant varier le facteur de lorentz des "blobs" du jet. */

#include "functions.h"
#include<math.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

const double e = 1.602e-19, me = 9.109e-31, eps0 = 1e-9/(36.*M_PI), c = 299792458;
double nu1=1e14, B1 = 1e-4, gam = 1e4, y_p1=pow(10,4.5), phi1 = M_PI/4., pitch = M_PI/4., kappa=1., p1=2.0;       //paramètres physiques par défaut : intensité du champ magnétique, facteur de lorentz des électrons, facteur de lorentz maximal pour la distribution en énergie, angle d'observation, pitch angle de l'électron, paramètres de la distribution en energie dE = kappa*E^{-p}
double gam2 = gam*gam;
double be = sqrt(gam2-1.)/gam, be2 = 1.-1./gam2;
double wg = e*B1/me;
double wc = (3./2.)*gam2/be*wg*sin(pitch);
double ae = be*c/(wg*sin(pitch));
double Tr = 2*M_PI*gam/wg;
double zeta = e*me/(24.*pow(M_PI,4)*eps0*c);        //facteur constant qui intervient dans l'expression de l'émissivité ordonnée pour un électron isolé

double z = 0.158, D_l = 749.*3.085*1e16*1e6;    //paramètres redshift et distance lumière pour 3C273

/*
 * Cas de l'électron isolé
 */
double beta(double y){      //calcul du facteur beta (v/c) en fonction du facteur de lorentz
	return sqrt(y*y-1.)/y;
}

double a(double B,double alpha, double y){      //calcul du rayon de courbure de la trajectoire de l'électron de pitch alpha et de facteur de lorentz y dans le champ magnétique d'intensité B
	return sqrt(y*y-1)*me*c/(e*B*sin(alpha));
}

double thetagamma(double theta, double y){
	return sqrt(1+y*y*theta*theta);
}

double eta(double w, double B, double alpha, double thetgam, double y){
	return (w*a(B,alpha,y)*pow(thetgam,3))/(3.*c*pow(y,3));
}

double L_isol(double w, double B, double phi, double alpha, double y){      //Emissivité ordonnee pour un électron isolé observé sous l'angle phi
	double theta = phi-alpha;
	double thetgam = thetagamma(theta,y);
	double eta_w = eta(w,B,alpha,thetgam,y);
	double bet = beta(y);

	double factor = (zeta*bet*bet*w*w*pow(thetgam,4))/(B*pow(y,3)*sin(alpha)*sin(alpha));

	double K23 = Bessel_K23(eta_w), K13 = Bessel_K13(eta_w);
	return factor*(K23*K23+(thetgam*thetgam-1)/(thetgam*thetgam)*K13*K13);
}


/*
 * Pour une population d'électrons
 */
double pas_adapte(int d, double x, double x1, double x2, double y1, double y2){     //permet d'adapter le nombre de points d'intégration (donc le pas) selon une loi puissance négative
	double a = (y1-y2)/(1./pow(x1,d)-1./pow(x2,d));
	return a*(1./pow(x,d)-1./pow(x2,d))+y2;
}


double fonction(double w, double B, double phi, double alpha, double y, double p){  //expression de la fonction à intégrer selon alpha et gamma pour obternir l'émissivité
	return 0.5*sin(alpha)*kappa*pow(y,-p)*L_isol(w,B,phi,alpha,y);
}

double fy(double* tab, double y){     //fonction(w,phi,alpha,y) mais dont seul 'y' est variable, pour l'intégration selon 'y'
	return fonction(tab[0],tab[1],tab[2],tab[4],y,tab[5]);
}
	
double fa(double* tab, double alpha){     //fonction(w,phi,alpha,y) mais dont seul 'alpha' est variable, pour l'intégration selon 'alpha'
	double w = tab[0], y_p = tab[3];
	tab[4] = alpha;
	int n = (int)pas_adapte(6,log10(w),7.,23.,1e5,50.);
	return GaussIntegration2(1.,y_p,n,fy,tab);
}

double L_pop(double w, double B, double phi, double y_p, double p, double* tab){        //Emissivité ordonnée pour une population d'électrons de distribution en énergie de paramètre p observée sous un angle phi
	double dalpha = 1e-3/sqrt((w/(2*M_PI))*1e-10);
	double alpha1 = phi-dalpha, alpha2 = phi+dalpha;        //seuls les photons émis par les électrons de pitch proche de l'angle d'observation sont observés
	int n = 50;//(int)pas_adapte(5,log10(w),7.,23.,1000.,50.);
	tab[0] = w ; tab[1] = B ; tab[2] = phi ; tab[3] = y_p; tab[4] = alpha1; tab[5] = p;
	return GaussIntegration2(alpha1,alpha2,n,fa,tab);
}

/*
 * Boost relativiste
 */

double DopplerFactor(double y, double phi){
	double bet = beta(y);
	return 1./(y*(1-bet*cos(phi)));
}

double phi_rest(double phi, double y){      //calcul de l'angle d'observation dans le référentiel de l'électron relativiste
	double d = DopplerFactor(1,phi);
    	return acos((cos(phi)-beta(y))/(1-beta(y)*cos(phi)));
	//return asin(d*sin(phi));
}

int main(){
	time_t start = time(nullptr);
	fstream data("data_boost.res",ios::out);
	data << nu1 << " " << B1 << " " << p1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;	//on passe les paramètres physiques à l'affichage python

	double wi = 1e7, wf = 1e15;
	double phi_i = 1e-3, phi_f = M_PI/2.;
	int i,j, n=90;
	double h = (log10(wf)-log10(wi))/n, g = (phi_f-phi_i)/n;
	double* m = (double*)malloc(6*sizeof(double));  //pointeurs dont on va avoir besoin dans L_pop, fy(*m,y) et fa(*m,alpha) pour se passer les arguments fixes
	double* l = (double*)malloc(6*sizeof(double));
    
    	double d = DopplerFactor(1.,M_PI/2.);
	double L_ref = pow(d,2.+(p1-1.)/2.)*L_pop(nu1/d,B1,phi_rest(M_PI/2.,1.),y_p1,p1,m)*pow(1.+z,(p1+1.)/2.)/(4.*M_PI*D_l*D_l);   //Stockage d'une valeur de référence pour une normalisation à l'affichage des données
	data << L_ref;
		
	double *y_bulk = (double*)malloc(9*sizeof(double));     //valeurs de facteur de lorentz pour les "blobs" du jet
	y_bulk[0]=1.;y_bulk[1]=2.;y_bulk[2]=3.;y_bulk[3]=4.;y_bulk[4]=5.;y_bulk[5]=8.;y_bulk[6]=10.;y_bulk[7]=15.;y_bulk[8]=20.;
	for(j=0;j<9;j++)
		data << " " << y_bulk[j];
	data << endl;

	for(i=0;i<n;i++){
		data << phi_i+g*i;
		for(j=0;j<9;j++){
            		d = DopplerFactor(y_bulk[j],phi_i+g*i);
		    	data << " " << pow(d,2.+(p1-1.)/2.)*L_pop(nu1/d,B1,phi_rest(phi_i+g*i,y_bulk[j]),y_p1,p1,m)*pow(1.+z,(p1+1.)/2.)/(4.*M_PI*D_l*D_l);	//On calcule le flux dans le référentiel de l'électron : (nu1/d, phi_rest) sont la fréquence et l'angle d'observation dans le référentiel de l'électron à partir de (nu1, phi) la fréquence et l'angle d'observation dans le référentiel de l'observateur
        	}
		data << endl;
	}
	data.close();
	free(m); free(l); free(y_bulk);
	cout << "Elapsed time : " << difftime(time(nullptr),start) << "s." << endl;
	system("python affichage_boost.py");		//affichage d'une superposition des courbes du flux en fonction de l'angle d'observation pour différentes valeurs de y_bulk
	return 0;
}
