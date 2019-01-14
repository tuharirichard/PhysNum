#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double M_earth = 5.9736e24;
double M_moon = 7.359e22;
double d_apo = 405500e3;
double d_peri = 363300e3;
double v_apo = 964;
double v_peri = 1076;
double G = 6.67384e-11;

struct vector{
	double x;
	double y;
	double v_x;
	double v_y;
};

void diff(struct vector V, struct vector *dV){
	dV->x = V.v_x;
	dV->y = V.v_y;
	dV->v_x = (G*M_earth * (-V.x))/(pow(pow(V.x,2)+pow(V.y,2),1.5));
	dV->v_y = (G*M_earth * (-V.y))/(pow(pow(V.x,2)+pow(V.y,2),1.5));
}

void euler_step(struct vector* V, struct vector* dV, double dt, FILE* p){	
	for(double i=0; i<10000000; i+=dt){
		diff(*V, dV);
		V->x += dV->x*dt;
		V->y += dV->y*dt;
		V->v_x += dV->v_x*dt;
		V->v_y += dV->v_y*dt;
		fprintf(p, "%lf %lf %lf %lf \n", V->x, V->y, V->v_x, V->v_y);
	}
}


int main(){
	double dt = 1;
	struct vector V={0, -d_peri, v_peri, 0};
	struct vector dV;

	FILE* p=fopen("adat.dat", "w");
	euler_step(&V, &dV, dt, p);
	fclose(p);
	
	return 0;
}
