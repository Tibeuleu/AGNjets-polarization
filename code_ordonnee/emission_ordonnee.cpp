/* Calcule et trace les fluxs en fonction de la longueur d'onde et de la fréquence pour émission ordonnée d'une population d'électrons suivant une loi de puissance. On choisit le paramètre que l'on veut faire varier dans le main() en décommentant les lignes et boucles correspondant à l'initialisation et aux calculs des fluxs faisant varier ce paramètre. Il faut penser à ne faire varier qu'un paramètre à la fois, et commenter ceux qui ne serviront pas. */

#include "functions.h"
#include<math.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

const double e = 1.602e-19, me = 9.109e-31, eps0 = 1e-9/(36.*M_PI), c = 299792458;
double B1 = 1e-4, gam = 1e4, y_p1=pow(10,4.5), phi1 = M_PI/4., pitch = M_PI/4., kappa=1., p1=2.0;   //paramètres physiques par défaut : intensité du champ magnétique, facteur de lorentz de électrons, facteur de lorentz maximal, angle d'observation, pitch angle de l'électron, paramètres kappa et p de la distribution en énergie dE = kappa*E^{-p}
double gam2 = gam*gam;
double be = sqrt(gam2-1.)/gam, be2 = 1.-1./gam2;
double wg = e*B1/me;
double wc = (3./2.)*gam2/be*wg*sin(pitch);
double ae = be*c/(wg*sin(pitch));
double Tr = 2*M_PI*gam/wg;
double zeta = e*me/(24.*pow(M_PI,4)*eps0*c);    //facteur constant intervenant dans l'expression de l'émissivité ordonnée pour un électron isolé

double z = 0.158, D_l = 749.*3.085*1e16*1e6;    //paramètres redshift et distance lumière pour 3C273

/*
 * Cas de l electron isole
 */
double beta(double y){      //calcul du facteur beta (v/c) en fonction du facteur de lorentz y
	return sqrt(y*y-1.)/y;
}

double a(double B,double alpha, double y){  //calcul du rayon de courbure de la trajectoire de l'électron de pitch alpha et de facteur de lorentz y dans un champ d'intensité B
	return sqrt(y*y-1)*me*c/(e*B*sin(alpha));
}

double thetagamma(double theta, double y){
	return sqrt(1+y*y*theta*theta);
}

double eta(double w, double B, double alpha, double thetgam, double y){
	return (w*a(B,alpha,y)*pow(thetgam,3))/(3.*c*pow(y,3));
}

double L_isol(double w, double B, double phi, double alpha, double y){  //Calcul de l'émissivité ordonnée pour un électron isolé observé selon en angle phi
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
	int n = (int)pas_adapte(6,log10(w),7.,23.,1e6,50.);
	double I = GaussIntegration2(1.,y_p,n,fy,tab);
    	return I;
}

double L_pop(double w, double B, double phi, double y_p, double p, double* tab){        //Emissivité ordonnée pour une population d'électrons de distribution en energie dE = kappa*E^{-p} observée selon un angle phi (le pointeur tab sert à stocker les variables des fonctions à variables multiples que l'on va intégrer)
	double dalpha = 1e-3/sqrt((w/(2*M_PI))*1e-10);
	double alpha1 = phi-dalpha, alpha2 = phi+dalpha;
	int n = 30;//(int)pas_adapte(5,log10(w),7.,23.,1000.,50.);
	tab[0] = w ; tab[1] = B ; tab[2] = phi ; tab[3] = y_p; tab[4] = alpha1; tab[5] = p;
	double I = GaussIntegration2(alpha1,alpha2,n,fa,tab);       //Intégration de fa sur la variable alpha, les autres variables sont stockées dans tab
    	return I;
}


int main(){
	time_t start = time(nullptr);
	fstream data("data_ord.res",ios::out);
	double wi = 1e7, wf = 1e15;
	int i,j, n=50;
	double h = (log10(wf)-log10(wi))/n;
	double* m = (double*)malloc(6*sizeof(double));  //pointeurs dont on va avoir besoin dans L_pop, fy(*m,y) et fa(*m,alpha) pour passer les variables non integrées
	double* l = (double*)malloc(6*sizeof(double));
	data << y_p1;
	double *phi = (double*)malloc(9*sizeof(double));    //valeurs d angles testées
	phi[0]=M_PI/2.;phi[1]=70.*M_PI/180.;phi[2]=M_PI/4.;phi[3]=M_PI/6.;phi[4]=M_PI/9.;phi[5]=M_PI/12.;phi[6]=M_PI/18.;phi[7]=M_PI/36.;phi[8]=M_PI/180.;
	//for(j=0;j<8;j++)
	//	data << " " << phi[j];
	double *B = (double*)malloc(15*sizeof(double));     //valeurs de B testées
	for(j=0;j<15;j++){
	//	B[j] = pow(10,-(4.+0.5*j));
		data << " " << B[j];
	}
	double *y_p = (double*)malloc(3*sizeof(double));   //valeurs de y_p testées
	for(j=0;j<3;j++){
		y_p[j] = pow(10,1.+2.5*j);
	//	data << " " << y_p[j];
	}
	double *p = (double*)malloc(8*sizeof(double));      //valeurs de p testées
	for(j=0;j<8;j++){
		p[j] = 1.+0.5*j;
		data << " " << p[j];
	}
	data << endl;
	for(i=0;i<n;i++){
		data << pow(10,log10(wi)+h*i);
	//	for(j=0;j<8;j++)
	//		data << " " << L_pop(2.*M_PI*pow(10,log10(wi)+h*i),B1,phi[j],y_p1,p1,m);        //Stockage de la fréquence 10^{initiale+h*i} et des émissivités calculées à cette fréquence en faisant varier un paramètre physique dans le fichier data_ord.res
	//	for(j=0;j<15;j++)
	//		data << " " << L_pop(2.*M_PI*pow(10,log10(wi)+h*i),B[j],phi1,y_p1,p1,m);
    	//  	for(j=0;j<3;j++){
    	//      	m[0]=2.*M_PI*pow(10,log10(wi)+h*i);m[1]=B1;m[2]=phi1;m[3]=y_p[j];m[4]=pitch;m[5]=p1;
	//      	data << " " << L_pop(2.*M_PI*pow(10,log10(wi)+h*i),B1,phi1,y_p[j],p1,m);
    	//  	}
		for(j=0;j<8;j++)
			data << " " << L_pop(2.*M_PI*pow(10,log10(wi)+h*i),B1,phi1,y_p1,p[j],m);
	    	data << endl;
	}
	data.close();
	free(m); free(l); free(phi); free(B); free(y_p); free(p);
	cout << "Elapsed time : " << difftime(time(nullptr),start) << "s." << endl;
	system("python affichage_ord.py &");      //Affichage des SED calculés et stockés dans data_ord.res avec python
	return 0;
}
