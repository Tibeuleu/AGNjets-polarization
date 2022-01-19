/*Calcule la polarisation et le degré de polarisation linéaire de l'émission synchrotron pour différentes valeurs de l'angle d'observation.*/

#include "functions.h"
#include "synchrotron.h"
#include<math.h>
#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<time.h>
using namespace std;

double B1 = 1e-4, gam = 1e4, y_p1=pow(10,4.5), phi1 = M_PI/4., p1=3.0;

int main(){
	time_t start = time(nullptr);
	fstream data("data_pola.res",ios::out);
	double wi = 1e7,wf = 1e15,w=0;
	int i,j,n=100;
	double h=(log10(wf)-log10(wi))/n,alpha=0,para=0,orth=0,d=0;
    data << gam;
    double phi[5] = {0.01,M_PI/6.,M_PI/4.,M_PI/3.,M_PI/2.};
    for(j=0;j<5;j++)
        data << " " << phi[j] << " " << 0 << " " << 0;
    data << endl;
	for(i=0;i<n;i++){
        w = pow(10,log10(wi)+h*i);
        data << w;
        for(j=0;j<5;j++){
            alpha = phi[j]+(2*drand48()-1)*1e-3;
            para = I_para(2.*M_PI*w,B1,phi[j],alpha,gam);
            orth = I_orth(2.*M_PI*w,B1,phi[j],alpha,gam);
            d = (orth-para)/(orth+para);
		    data << " " << para << " " << orth << " " << d;
        }
        data << endl;
    }
	data.close();
	cout << "Elapsed time : " << difftime(time(nullptr),start) << "s." << endl;
	system("python affichage_pola.py &");
	return 0;
}
