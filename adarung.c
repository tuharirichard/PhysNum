#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

const double mf = 5.9736e24;
const double mh = 7.359e22;
const double dA = 405500e3;
const double dP = 363300e3;
const double vA = 964;
const double vP = 1076;
const double G = 6.67384e-11;

void dyp(double* vector, double* deltavec){
	
	deltavec[0] = vector[2];
	deltavec[1] = vector[3];
	deltavec[2] = (G*mf*(-vector[0]))/(pow(vector[0]*vector[0]+vector[1]*vector[1],1.5));
	deltavec[3] = (G*mf*(-vector[1]))/(pow(vector[0]*vector[0]+vector[1]*vector[1],1.5));
}


void runge_kutta(double* vector, double* deltavec, double h, FILE* p){
	double V[4], k[4], k2[4], k3[4], k4[4], avector[4];
	double ah = 2*h, diff;
		
	for(double j=0; j<100000000; j+=h){
		diff = 0;
		dyp(vector, deltavec);
		for(int i=0; i<4; i++){
			k[i] = h * deltavec[i];
			V[i] = vector[i] + 0.5 * k[i]; 
		}

		for(int i=0; i<4; i++){
                        k[i] = h * deltavec[i];
                        V[i] = vector[i] + 0.5 * k[i];
		}	

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k2[i] = h * deltavec[i];
			V[i] = vector[i] + 0.5 * k2[i];
		}

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k3[i] = h * deltavec[i];
			V[i] = vector[i] + k3[i];
		}

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k4[i] = h * deltavec[i];
			vector[i] += ((1./6.)*k[i] + (1./3.)*k2[i] + (1./3.)*k3[i] + (1./6.)*k4[i]);
		}
		//aaaaaaaaaaaaaaaaaaaa
		dyp(vector, deltavec);
		for(int i=0; i<4; i++){
			k[i] = ah * deltavec[i];
			V[i] = vector[i] + 0.5 * k[i]; 
		}

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k2[i] = ah * deltavec[i];
			V[i] = vector[i] + 0.5 * k2[i];
		}

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k3[i] = ah * deltavec[i];
			V[i] = vector[i] + k3[i];
		}

		dyp(V, deltavec);
		for(int i=0; i<4; i++){
			k4[i] = ah * deltavec[i];
			avector[i] += ((1./6.)*k[i] + (1./3.)*k2[i] + (1./3.)*k3[i] + (1./6.)*k4[i]);
		}

		for(int i=0; i<4; i++){
			diff += pow(avector[i]-vector[i],2)/pow(vector[i],2);
		}

		if(diff < 0.05){
			ah *= 2;
			fprintf(p, "%lf %lf %lf %lf \n", avector[0], avector[1], avector[2], avector[3]);
		}
		else{
			ah /=2;
			fprintf(p, "%lf %lf %lf %lf \n", vector[0], vector[1], vector[2], vector[3]);
		}
	}
}

int main(){

double h=2000;

double* vector;
double* deltavec;

vector=(double*)calloc(4,sizeof(double));
deltavec=(double*)calloc(4,sizeof(double));


vector[0]=0;
vector[1]=-dP;
vector[2]=vP;
vector[3]=0;

FILE* p=fopen("adat.dat","w");
runge_kutta(vector,deltavec,h,p);
fclose(p);

return 0;
}
