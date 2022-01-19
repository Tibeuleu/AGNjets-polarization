/*Calcule par interpolation linéaire des valeurs de flux stockées dans les fichiers de données "data_ord_p.res" les valeurs obtenues par le modèle pour les mêmes fréquences que celles mesurées dans "photandseds.tbl", le SED du quasar étudié, de sorte à pouvoir comparer les données et choisir le meilleur fit.*/

#include "functions.h"
#include<math.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

const double e = 1.602e-19, me = 9.109e-31, eps0 = 1e-9/(36.*M_PI), c = 299792458;
double B1 = 1e-4, gam = 1e4, y_p1=pow(10,4.5), phi1 = M_PI/4., pitch = M_PI/4., kappa=1., p1=2.0, y_bulk1=4.0;
double gam2 = gam*gam;
double zeta = e*me/(24.*pow(M_PI,4)*eps0*c);

double z = 0.77, D_l = 3407.27e6*3.085e16;    //paramètres redshift et distance lumière pour 3C175

/*
 * Cas de l'électron isolé
 */
double beta(double y){				//calcul du facteur beta (v/c) en fonction du facteur de lorentz y
	return sqrt(y*y-1.)/y;
}

double a(double B,double alpha, double y){	//calcul du rayon de courbure de la trajectoire de l'électron de pitch angle alpha et de facteur de lorentz y dans un champ magnétique d'intensité B
	return sqrt(y*y-1)*me*c/(e*B*sin(alpha));
}

double thetagamma(double theta, double y){	//facteur provenant du développement limité du Flux en theta
	return sqrt(1+y*y*theta*theta);
}

double eta(double w, double B, double alpha, double thetgam, double y){	//facteur provenant d'un changement de variable pour le calcul du flux
	return (w*a(B,alpha,y)*pow(thetgam,3))/(3.*c*pow(y,3));
}

double L_isol(double w, double B, double phi, double alpha, double y){	//Cacul de l'émissivité pour un électron isolé
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
double pas_adapte(int d, double x, double x1, double x2, double y1, double y2){     	//permet d'adapter le nombre de points d'intégration (donc le pas) selon une loi puissance négative
	double a = (y1-y2)/(1./pow(x1,d)-1./pow(x2,d));
	return a*(1./pow(x,d)-1./pow(x2,d))+y2;
}


double fonction(double w, double B, double phi, double alpha, double y, double p){  	//expression de la fonction à intégrer selon alpha et gamma pour obtenir l'émissivité
	return 0.5*sin(alpha)*kappa*pow(y,-p)*L_isol(w,B,phi,alpha,y);
}

double fy(double* tab, double y){     		//fonction(w,phi,alpha,y) mais dont seul 'y' est variable, pour l'intégration selon 'y'
	return fonction(tab[0],tab[1],tab[2],tab[4],y,tab[5]);
}
	
double fa(double* tab, double alpha){     	//fonction(w,phi,alpha,y) mais dont seul 'alpha' est variable, pour l'intégration selon 'alpha'
	double w = tab[0], y_p = tab[3];
	tab[4] = alpha;
	int n = (int)pas_adapte(6,log10(w),7.,23.,1e6,50.);
	return GaussIntegration2(1.,y_p,n,fy,tab);
}

double L_pop(double w, double B, double phi, double y_p, double p, double* tab){	//calcul de l'émissivité pour une population d'électrons
	double dalpha = 1e-3/sqrt((w/(2*M_PI))*1e-10);
	double alpha1 = phi-dalpha, alpha2 = phi+dalpha;
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

double phi_rest(double phi, double y){		//calcul de phi dans le référentiel dans le référentiel au repos d'un électron de facteur de lorentz y dans le référentiel de l'observateur 
	double d = DopplerFactor(y,phi);
	return acos((cos(phi)-beta(y))/(1-beta(y)*cos(phi)));
	//return asin(d*sin(phi));
}

int main(){
	time_t start = time(nullptr);
	fstream data("Flux_modele.res",ios::out);
	double wi = 1e7, wf = 1e15;
	double phi_i = 1e-3, phi_f = M_PI/2.;
	int i,j, n=100;
	double h = (log10(wf)-log10(wi))/n, g = (phi_f-phi_i)/n;
	double* m = (double*)malloc(6*sizeof(double));  	//pointeurs dont on va avoir besoin dans L_pop, fy(*m,y) et fa(*m,alpha) pour se passer les arguments fixes
	double* l = (double*)malloc(6*sizeof(double));
    
	data << phi1;
	
	double *p = (double*)malloc(5*sizeof(double));      	//valeurs de p testées
	for(j=0;j<5;j++){
		p[j] = 2.75+0.125*j;
		data << " " << p[j];
	}	
	data << endl;
	
	double d,w1=0,w2=0,f,u1,u2;
    	double mem[8] = {0.,0.,0.,0.,0.,0.,0.,0.}, fromf1[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.}, fromf2[9] = {1.,0.,0.,0.,0.,0.,0.,0.,0.};	//mem stocke un point interpolé dans l'intervalle [fromf1;fromf2]
    	fstream fich("photandseds.tbl",ios::in);
    	for(j=0;j<4;j++){
        	fich.ignore(1000,'\n');
    	}
    	while(fich >> w1 >> f >> u1 >> u2){		//on récupère les abscisses des points de "photandseds.tbl"
	    	data << w1;
        	if(w1 != w2){
            		fstream forinterpol("data_ord_p_3.res",ios::in);
            		forinterpol.ignore(1000,'\n');
            		while(forinterpol >> fromf2[0] >> fromf2[1] >> fromf2[2] >> fromf2[3] >> fromf2[4] >> fromf2[5]){	//on récupère les valeurs fromf1 et fromf2 du set de données "data_ord_p.res"
                		if(fromf2[0] > w1)
                    			break;
                		for(j=0;j<6;j++)
                    			fromf1[j] = fromf2[j];
            		}
            		forinterpol.close();
            		for(j=0;j<5;j++){		//on effectue l'interpolation linéaire des points de "data_ord_p.res" aux abscisses des points de "photandseds.tbl"
				//d = DopplerFactor(y_bulk1,phi1);
				mem[j] = fromf1[j+1]+(fromf2[j+1]-fromf1[j+1])*(w1-fromf1[0])/(fromf2[0]-fromf1[0]);//pow(d,2.+(p[j]-1.)/2.)*L_pop(2.*M_PI*w1,B1,phi_rest(phi1,y_bulk1),y_p1,p[j],m)*pow(1.+z,(p[j]+1.)/2.)/(4.*M_PI*D_l*D_l);
                		w2 = w1;
            		}
        	}
        	for(j=0;j<5;j++)
            		data << " " << mem[j];
		data << endl;
	}
    	fich.close();
	data.close();
	free(m); free(l); free(p);
	cout << "Elapsed time : " << difftime(time(nullptr),start) << "s." << endl;
	system("python plotfromsed.py");		//tracé du fit avec python
	return 0;
}
