/*Génère des électrons dans un espace délimité par un double cône surmonté d'un lobe qui émettent des photons par effet synchrotron (de façon ordonnée boostée ou isotrope).
 * Les photons sont ensuite détectés par des détecteurs placés à différents angles d'observations et qui enregistrent l'intensité reçue, le degré et la direction de polarisation sur chaque pixel.
 * On peut aussi choisir si l'ont prend en compte les effets faradays qui s'appliquent à la polarisation dans le milieu intergalactique.
 */

#include "functions.h"
#include "synchrotron.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<time.h>
using namespace std;

bool const ordonne = true;  	//si "true" l'émission par le jet est ordonnée, si "false" l'émission est isotrope
bool const farad = true;    	//si "true" on prend en compte la rotation et la dépolarisation par effet Faraday
double const B_f = 1.*1e-9;    	//intensité du champ magnétique intergalactique pour les effets Faraday

double const proportion=0.20;   //proportion de photons dans les lobes par rapport au nombre total de photons émis

double const e = 1.602e-19, me = 9.109e-31, eps0 = 1e-9/(36.*M_PI), mu0 = 4.*M_PI*1e-7, n0 = 1.1e2, c = 299792458;  //n0 est la densité d'électrons (en m^-3)
double const h_planck = 6.626e-34;
double const I = 1., Q = 0.1, U = 0.1, V = 0.2; 	//valeurs par défaut des paramètres de Stokes
double const z = 0.77, D_l = 3407.27e6*3.0857e16;   	//redshift et distance lumière de 3C175
double const B_j = 1e-4, y_p_j = pow(10,4.5), p_j = 3.15, y_bulk_j = 4., nu_j = 5e9, theta_j = 10., R_j = 412.01e3*3.0857e16;    //initialisation des paramètres physiques : intensité du champ magnétique, énergie maximale des particules, pente spectrale, facteur de lorentz du jet, fréquence etudiée, demi-angle d'ouverture du jet, taille du jet (pour 3C175)
double const R_l = 0.25*R_j;    //rayon du lobe
double const ep_l = 0.20*R_l, theta_l_max = 3./4.*M_PI;	//épaisseur et ouverture du lobe

int const pas_theta = 10, pas_phi = 10;  		//pas d'observations selon les angles theta et phi
int const n_theta = (int)(180./pas_theta)+1, n_phi = (int)360./pas_phi;

double const theta_det[10] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.}, phi_det = 0.;   			//angles de positionnnement des capteurs
double const x_max = 1.4*R_j, y_max = 1.4*R_j, pas_x = 5.782e3*3.0857e16, pas_y = 5.782e3*3.0857e16;    //Dimensions et résolution des capteurs
int const n_theta_det = 10, n_x = (int)(2*x_max/pas_x), n_y = (int)(2*y_max/pas_y);

int const n_photons = 1e9;

class photon {
    public:
        bool jet;       //"true" si le photon est émis par le jet, "false" s'il est émis par le lobe
        double r_e;     //coordonnées sphériques de l'électron à l'instant de l'émission où le centre de la galaxie est l'origine du repère, le plan galactique est contenu dans le plan xOy et les jets s'étendent le long de l'axe Oz (on suppose le champ magnétique confondu avec Oz)
        double theta_e;
        double phi_e;
        double theta_p; //angle entre la direction d'émission du photon et Oz
        double phi_p;   //angle entre la direction d'émission du photon et Ox
        double E;       //énergie du photon (fréquence)
        double I;       //paramètres de Stokes de l'émission associée au photon
};

/*
 * Pour l'initialisation et la sortie des données des détecteurs
 */

void init(double tab[n_theta][n_phi]){
    for(int i=0;i<n_theta;i++)
        for(int j=0;j<n_phi;j++)
            tab[i][j] = 0;
}

void init_proj(vector<vector<vector<double>>> *tab, int n_det, int n_x, int n_y){
    (*tab).resize(n_det);
    for(int k=0;k<n_det;k++){
        (*tab)[k].resize(n_x);
        for(int i=0;i<n_x;i++){
            (*tab)[k][i].resize(n_y);
            for(int j=0;j<n_y;j++)
                (*tab)[k][i][j] = 0;
        }
    }
}

void sortie(double tab[n_theta][n_phi],char name){
	char file[] = {'d','a','t','a','/','d','a','t','a','_',name,'.','r','e','s','\0'};
	fstream data(file,ios::out);
	for(int i=0;i<n_theta;i++){
		for(int j=0;j<n_phi;j++)
			data << tab[i][j] << " ";
			data << endl;
	}
	data.close();
}

void sortie_proj(vector<vector<double>> tab, int n_x, int n_y, string name, int nbr){
	string filename = "data/detect_"+name+"_"+to_string(nbr)+".res";
	fstream data(filename.c_str(),ios::out);
	for(int j=0;j<n_y;j++){
		for(int i=0;i<n_x;i++)
			data << tab[i][j] << " ";
			data << endl;
	}
	data.close();
}

/*
 * Fonctions de géneration des photons selon la géométrie du problème et des distributions imposées par l'émission synchrotron ordonnée et boostée
 */

void distrib(double f[n_theta], double nu, double B, double y_p, double p, double y_bulk){
    if(ordonne){
        double* m = (double*)malloc(6*sizeof(double));
	    double d = DopplerFactor(1.,M_PI/2.);
	    double L_ref = pow(d,2.+(p-1.)/2.)*L_pop(nu,B,phi_rest(M_PI/2.,1.),y_p,p,m);
	    for(int i=0;i<n_theta;i++){
            try{
                d = DopplerFactor(y_bulk,M_PI/180.*(i*pas_theta));
                f[i] = pow(d,2.+(p-1.)/2.)*L_pop(2.*M_PI*nu,B,phi_rest(M_PI/180.*(i*pas_theta),y_bulk),y_p,p,m)/L_ref;  	//distribution calculée par rapport à la luminosité relative pour une émission boostée et observée sous un certain angle
            } catch(...){       	//pour certaines combinaisons de paramètres (nu,phi) erreur dans la fonction de Bessel calculée en x << 1, donc B_k(x)~0
                f[i] = 0;
            }
        }
        free(m);
    }
    else
        for(int i=0;i<n_theta;i++)  	//pour une émission isotrope la distribution est uniforme sur [0°;180°]
            f[i] = 1.;
}

double max_distrib(double D[n_theta]){
	double M = D[0];
	for(int i=1;i<n_theta;i++)
		if(M<D[i])
			M = D[i];
	return M;
}

void rand_position(photon* p){      	//On place aléatoirement un électron dans l'ensemble jets+lobes (avec une proportion déterminée en debut du programme)
    double prop = proportion/(1.-proportion);
    double r = (2*drand48()-1)*(1.+prop)*R_j;
    if(fabs(r) < R_j){
        (*p).jet = true;
        (*p).theta_e = drand48()*theta_j*M_PI/180.;
        (*p).phi_e = drand48()*2.*M_PI;
        (*p).r_e = R_j/cos((*p).theta_e)*(2*drand48()-1);
    }
    else{
        (*p).jet = false;
        double sgn = r/fabs(r);
        r = R_j+ep_l-R_l;
        double r_s = ((ep_l/R_l)*drand48()+(1.-(ep_l/R_l)))*R_l;
        double theta_s = drand48()*theta_l_max;
        (*p).theta_e = atan(r_s*sin(theta_s)/(r+r_s*cos(theta_s)));
        (*p).phi_e = drand48()*2.*M_PI;
        (*p).r_e = sgn*sqrt(pow(r+r_s*cos(theta_s),2)+pow(r_s*sin(theta_s),2));
    }
}

void rand_direction(double distrib[n_theta],double max, photon* p){	//On tire aléatoirement une direction d'émission du photon
    if((*p).jet){	//Si le photon est émis depuis le jet, émission selon la distribution ordonnée
        double theta = 0, f_theta = 0;
        do{
            theta = drand48()*181.;
            f_theta = drand48()*max;
        }while(f_theta > distrib[(int)(theta/pas_theta)]);
        if((*p).r_e < 0)
            (*p).theta_p = 180.-theta;
        else
            (*p).theta_p = theta;
        (*p).phi_p = drand48()*360.;
    }
    else{               //Si le photon est émis depuis le lobe, émission isotrope
        (*p).theta_p = drand48()*180.;
        (*p).phi_p = drand48()*360.;
    }
}

void proj_plane(double theta, double phi, photon p, double* pos){	//projection du photon sur le plan du capteur où il a été détecté
    double theta_e = p.theta_e, phi_e = p.phi_e, r_e = p.r_e;
    pos[0] = -r_e*sin(theta_e)*sin(phi-phi_e);
    pos[1] = r_e*(sin(theta)*cos(theta_e)-cos(theta)*sin(theta_e)*cos(phi-phi_e));
}

double dist(double D, double theta_d, double phi_d, photon p, double* pos){	//calcul de la distance-luminosité parcourue par le photon jusqu'au capteur
    double r = p.r_e, theta_e = p.theta_e, phi_e = p.phi_e, xp = pos[0], yp = pos[1];
    double ps1 = -D*(sin(theta_e)*sin(theta_d)*cos(phi_e-phi_d)+cos(theta_e)*cos(theta_d));
    double ps2 = -xp*sin(theta_e)*sin(phi_e-phi_d)+yp*(sin(theta_e)*sin(theta_d)*cos(phi_e-phi_d)-cos(theta_e)*sin(theta_d));
    return sqrt(D*D+r*r+(xp*xp+yp*yp)+2.*r*(ps1+ps2));
}

double Faraday_rot(double B, double nu, double d){      //calcul de la rotation faraday de la polarisation par le champ magnétique intergalactique
    double coeff = 0.812;//(e*e*e)/(2.*M_PI*me*me*c*c*c*c);
    double ne = pow(1.+z,3)*1.1e-7;//*n0;
    return coeff*ne*B*1e9*d/3.0857e16*c*c/(nu*nu);
}

double Faraday_dep(double B, double nu, double d, double p0){   //calcul de la dépolarisation faraday par le champ magnétique intergalactique
    double RM = angle_princ(Faraday_rot(B,nu,d));
    return p0*fabs(sinc(RM));
}

/*
 * Appel des fonctions
 */

int main(){
    time_t start = time(nullptr);
    srand48(time(NULL));
    system("mkdir data");
    fstream data("data/photons.res",ios::out);      //création-ouverture des fichiers pour l'affichage 3D des positions des électrons dans l'espace et de la direction des photons émis-reçus par chaque capteur
    fstream data2[n_theta_det];
    for(int k=0;k<n_theta_det;k++){
        string filename = "data/photons_"+to_string(k)+".res";
	    data2[k].open(filename.c_str(),ios::out);
    }

    data << n_theta_det << " " << pas_x << " " << pas_y << " " << nu_j << " " << y_bulk_j << endl;  //transfert des paramètres physiques à python
    data << theta_j << " " << R_j << " " << n_photons << " " << p_j << " " << B_j << endl;
    data << D_l/3.0857e16 << " " << z << " " << n0*pow(1+z,3) << " " << B_f*1e-9 << " " << 0 << endl;

    double* point = (double*)malloc(2*sizeof(double));
    double theta, phi;
    //double D_I[n_theta][n_phi];	                	//Détecteurs pour les paramètres de Stokes pour l'émission totale
    vector<vector<vector<double>>> D_I_proj;	        	//Détecteurs en intensité pour les plans de ciel
    vector<vector<vector<vector<double>>>> D_pola_p;    	//Détecteurs du degré de polarisation pour les plans de ciel
    vector<vector<vector<vector<double>>>> D_pola_theta;
    double distrib_theta[n_theta];	                    	//valeurs pour la distribution du flux en fonction de l'angle d'observation
    distrib(distrib_theta,nu_j,B_j,y_p_j,p_j,y_bulk_j);
    double M = max_distrib(distrib_theta);
    //init(D_I);
    D_pola_p.resize(2); D_pola_theta.resize(2);
    init_proj(&D_I_proj,n_theta_det,n_x,n_y); init_proj(&D_pola_p[0],n_theta_det,n_x,n_y); init_proj(&D_pola_theta[0],n_theta_det,n_x,n_y); init_proj(&D_pola_p[1],n_theta_det,n_x,n_y); init_proj(&D_pola_theta[1],n_theta_det,n_x,n_y); 	//initialisation à 0 des differents capteurs

    photon newphoton;
    int ind[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};      //indice pour le comptage d'électrons à sauter pour ne pas surcharger l'affichage 3D
    int i=0,j=0, x, y;
	double gam=0,para=0,orth=0,p=0,angle=0;

    for(i=0;i<n_photons;i++){
	    newphoton.E = h_planck*nu_j;                    	//émission d'un photon
	    rand_position(&newphoton);
	    rand_direction(distrib_theta,M,&newphoton);
	    theta = newphoton.theta_p;
	    phi = newphoton.phi_p;
	    newphoton.I = I;
	    //D_I[(int)(theta/pas_theta)][(int)(phi/pas_phi)] += newphoton.I;
	    if(n_photons>1e6){                              	//limitation du nombre de photons enregistrés pour la représentation 3D
		    j++;
		    if(j>1000){
			    j=0;
			    data << newphoton.r_e*sin(newphoton.theta_e)*sin(newphoton.phi_e) << " " << newphoton.r_e*sin(newphoton.theta_e)*cos(newphoton.phi_e) << " " << newphoton.r_e*cos(newphoton.theta_e) << " " << newphoton.theta_p << " " << newphoton.phi_p << endl;
		    }
	    }
	    else
		    data << newphoton.r_e*sin(newphoton.theta_e)*sin(newphoton.phi_e) << " " << newphoton.r_e*sin(newphoton.theta_e)*cos(newphoton.phi_e) << " " << newphoton.r_e*cos(newphoton.theta_e) << " " << newphoton.theta_p << " " << newphoton.phi_p << endl;
        for(int k=0;k<n_theta_det;k++){                 	//Détection du photon sur les différents capteurs
	        if((fabs(theta-theta_det[k])<.5)){// & (fabs(phi-phi_det)<.5)){    	// intégration azimutale par des arguments de symmétrie cylindrique
		        //i++;
		        proj_plane(theta_det[k]*M_PI/180.,phi_det*M_PI/180.,newphoton,point);
		        x = n_x-1-(int)((x_max-point[0])/pas_x);    			//projection du photon sur un pixel du capteur
		        y = n_y-1-(int)((y_max-point[1])/pas_y);

                gam = drand48()*(y_p_j-1e0)+1e0;
                para = I_para(2.*M_PI*nu_j,B_j,theta_det[k]*M_PI/180.,theta*M_PI/180.,gam);
                orth = I_orth(2.*M_PI*nu_j,B_j,theta_det[k]*M_PI/180.,theta*M_PI/180.,gam);
                p = (orth-para)/(orth+para);            	//calcul du degré de polarisation

                if(point[1] != 0)
                    angle = angle_princ(atan(point[0]/point[1]));   //calcul de l'angle de polarisation


                else
                    angle = point[0]/fabs(point[0])*M_PI/2.;
                D_I_proj[k][x][y] += newphoton.I;
                D_pola_p[0][k][x][y] += p;              	//ajout du degré de polarisation sur le pixel consideré, qui sera normalisé par l'intensite sur le pixel après l'acquisition
                D_pola_p[1][k][x][y] += Faraday_dep(B_f,nu_j,dist(D_l,theta_det[k],phi_det,newphoton,point),p);     //ajout du degré de polarisation dépolarisé par effet faraday, sera normalisé
                //D_pola_theta[l][k][x][y] += angle;
                D_pola_theta[1][k][x][y] += Faraday_rot(B_f,nu_j,dist(D_l,theta_det[k],phi_det,newphoton,point));   //ajout de l'angle de rotation faraday, sera normalisé
                ind[k]++;
                if(ind[k]>100){
                   		ind[k]=0;
                   		data2[k] << newphoton.r_e*sin(newphoton.theta_e)*sin(newphoton.phi_e) << " " << newphoton.r_e*sin(newphoton.theta_e)*cos(newphoton.phi_e) << " " << newphoton.r_e*cos(newphoton.theta_e) << " " << newphoton.theta_p << " " << phi_det << endl;//newphoton.phi_p << endl;
                }
            }
        }
    }
    double px=0, py=0;
    for(int k=0;k<n_theta_det;k++)
        for(int i=0;i<n_x;i++)
            for(int j=0;j<n_y;j++)
                if(D_I_proj[k][i][j] > 0){
                    D_pola_p[0][k][i][j] = D_pola_p[0][k][i][j]/D_I_proj[k][i][j];      //moyennage de la polarisation sur la ligne de visée correspondant au pixel
                    D_pola_p[1][k][i][j] = D_pola_p[1][k][i][j]/D_I_proj[k][i][j];      //moyennage de la depolarisation sur la ligne de visée correspondant au pixel
                    px = x_max-pas_x*(n_x-1-i);
                    py = y_max-pas_y*(n_y-1-j);
                    //D_pola_theta[0][k][i][j] = D_pola_theta[0][k][i][j]/D_I_proj[k][i][j];      	//moyennage de l'angle de polarisation sur la ligne de visée correspondant au pixel
                    D_pola_theta[1][k][i][j] = D_pola_theta[1][k][i][j]/D_I_proj[k][i][j];      	//moyennage de la rotation faraday sur la ligne de visée correspondant au pixel
                    D_pola_theta[0][k][i][j] = angle_princ(atan(px/py)+D_pola_theta[0][k][i][j]);   	//l'angle de polarisation est calculé en fonction de l'angle projeté entre la direction du champ magnétique et la position de l'électron
                    D_pola_theta[1][k][i][j] = angle_princ(atan(px/py)+D_pola_theta[1][k][i][j]);
                }
    data.close();
    for(int k=0;k<n_theta_det;k++)
        data2[k].close();
    free(point);
    //sortie(D_I,'I');
    for(int k=0;k<n_theta_det;k++)          //sortie des données de capteurs
        sortie_proj(D_I_proj[k],n_x,n_y,"I",k);
    string name;
    for(int l=0;l<2;l++){
        for(int k=0;k<n_theta_det;k++){
            name = to_string(l)+"_pola_p";
            sortie_proj(D_pola_p[l][k],n_x,n_y,name,k);
            name = to_string(l)+"_pola_theta";
            sortie_proj(D_pola_theta[l][k],n_x,n_y,name,k);
        }
        if(not(farad))
            break;
    }
    cout << "Elapsed time : " << difftime(time(nullptr),start) << "s." << endl;
    //system("python affichage_multiple_0.py &");           //affichage 3D+image sans les effets faraday
    if(farad){
        //system("python affichage_multiple_1.py &");       //affichage 3D+image obtenue avec les effets faraday
        system("python affichage.py &");                    //affichage comparé des images obtenues avec et sans effets faraday
    }
    return 0;
}
