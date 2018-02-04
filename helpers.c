#include <stdio.h>
#include <math.h>


// SET UP LATTICE

int shift(double betha, int x){
	/*Helpfunctin to calculate the bottom point of each x-column.
	For a lattice */
	double temp = (double)(betha+0.5)*x;
	int shift = - (int)temp;
	if(temp + (double)shift < 0){
		shift = shift+1;
	}
	return shift;
}

int neighbour_index(int nx, int ny, int lx, int ly){
	int k = lx*(2*ny +1) +ly;
	int n = ((2*nx+1)*(2*ny+1)-1); //total number of neighbours, indices go up to n-1!
	int j = -10;

	if(k<((n/2))){
		j = k;
	}
	else if(k>((n/2))){
		j = k-1;
	}
	else if(k==((n/2))){
		printf("Error! self-interaction (lx = nx, ly = ny) is not allowed!");
	}
	return j;
}

void calc_neighbours(int *save_neig, int Ny, int Nx, int ny, int nx, double betha){
	/*Calculate the neighbours of each particle.
	  save_neig has to be filled with -1!!*/
	
	int n = (2*nx+1)*(2*ny+1) -1; //total number of neighbours
	int N = Nx*Ny; //total number of particles

	//fill neighbour array with -1 (hast to be smaller zero) to detect not exixting boundary neighbours correctly!
	int temp_count1, temp_count2,temp_ind;
	for(temp_count1=0;temp_count1<N;temp_count1++){
		for(temp_count2=0;temp_count2<n;temp_count2++){
			temp_ind = n*temp_count1 + temp_count2;
			save_neig[temp_ind] = -1;
		}
	}


	int mx,my,x_temp,y_temp,jIndex,iIndex; /*temporary needed variables*/
	int x,y,lx,ly;
	for(x=0; x<Nx; x++){
		for(y=0; y<Ny; y++){
			iIndex = x*Ny+y;
			for(lx=0; lx<(2*nx+1);lx++){
				mx = -nx+lx;
				x_temp = x+mx;
				if(x_temp < 0){
					continue;
					}
				if(x_temp > (Nx-1)){
					continue;
				}
				for(ly=0; ly<(2*ny+1); ly++){
					my = -ny+ly;
					y_temp = y+my+shift(betha,x)+shift(betha,mx)-shift(betha,x_temp);
					if(y_temp<0){
						continue;
					}
					if(y_temp>(Ny-1)){
						continue;
					}
					if(y_temp == y && x_temp==x){
						continue;
					}
					jIndex = neighbour_index(nx,ny,lx,ly);
					save_neig[(iIndex*n + jIndex)] = x_temp*Ny + y_temp;
				}
			}
		}
	}
}

void calc_lattice(double *save_x, double *save_y, int Nx, int Ny, double alpha, double betha, double rescaled_lat_const){
	int index;
	int x,y;
	for(x=0;x<Nx;x++){
		for(y=0;y<Ny;y++){
			index = x*Ny+y;
			save_x[index] = alpha*(double)x*rescaled_lat_const;
			save_y[index] = ((betha+0.5)*(double)x + (double)(y+shift(betha,x)))*rescaled_lat_const;
		}
	}
}

void calc_ref_lattice(double *save_x, double *save_y, int nx, int ny, double alpha, double betha, double rescaled_lat_const){
	/*Calculate the reference lattice for all neighbours. Place in an array on position j the relative position of the
	  neighbour j.*/
	int index,my,mx;
	int lx,ly;
	for(lx=0; lx<(2*nx+1);lx++){
		mx = -nx+lx;
		for(ly=0; ly<(2*ny+1); ly++){
			my = -ny+ly;
			if( my==0 && mx==0){
				continue;
			}
			index = neighbour_index(nx,ny,lx,ly);
			save_x[index] = alpha*(double)mx*rescaled_lat_const;
			save_y[index] = ((betha+0.5)*(double)mx + (double)(my+shift(betha,mx)))*rescaled_lat_const;
		}
	}
}

// INTERACTION 

/*CALCULATE JACOBIAN, AND HESSIANS*/

double h1(double r, double Q, double l){
	/*Helpfunction 1 for the calculation of the Jacobian.*/
	// SI factor 1/(4*pi*epsilon_0)
	const double SI = 8.987551787997912E+9;
	return ((SI*Q)/(r*r*r)) * (1.0+(r/l)) * exp((-r)/l);
}

double h2(double r, double Q, double l){
	/*Helpfunction 2 for the calculation of the jacobian.
	h2 is the first derivative of h1 with respect to r*/
	// SI factor 1/(4*pi*epsilon_0)
	const double SI = 8.987551787997912E+9;
	return ((-SI*Q)/(r*r*r*r)) * exp((-r)/l) * (3.0*(1.0+(r/l)) + (r*r)/(l*l));	
}

double h3(double r, double Q, double l){
	/*Helpfunction 3 for the calculation of the jacobian.
	h3 is the second derivative of h1 with respect to r*/
	// SI factor 1/(4*pi*epsilon_0)
	const double SI = 8.987551787997912E+9;
	return (SI*Q / (r*r*r*r*r)) * exp((-r)/l) * ((r*r*r)/(l*l*l) + 5.0*(r*r)/(l*l) + 12.0*(1+(r/l)));	
}

double A(double r, double Q, double l, double q, double d){
	/*A = Q*h1(r) + q*h1(rd), factor which comes up several times if 
	  caclulating the jacobian.*/
	double rd = sqrt(r*r + d*d);
	return Q*h1(r,Q,l) + q*h1(rd,Q,l);
}

double B(double r, double Q, double l, double q, double d){
	/*B = Q*h2(r) + q*h2(rd), factor which comes up several times if 
	  caclulating the jacobian.*/
	double rd = sqrt(r*r + d*d);
	return Q*h2(r,Q,l)/r + q*h2(rd,Q,l)/rd;
}

double C(double r, double Q, double l, double q, double d){
	/*C = q*h2(rd), factor which comes up several times if 
	  caclulating the jacobian.*/
	double rd = sqrt(r*r + d*d);
	return q*h2(rd,Q,l)/rd;
}

double D(double r, double Q, double l, double q, double d){
	double rd = sqrt(r*r + d*d);
	return Q*(h3(r,Q,l)/(r*r) - h2(r,Q,l)/(r*r*r)) + q*(h3(rd,Q,l)/(rd*rd) - h2(rd,Q,l)/(rd*rd*rd));
}

double E(double r, double Q, double l, double q, double d){
	double rd = sqrt(r*r + d*d);
	return q*(h3(rd,Q,l)/(rd*rd) - h2(rd,Q,l)/(rd*rd*rd));
}

void calc_jacobian_m(double save_jac[3][3], int j, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double d, double mass){
	/*Calculate the jacobian devided by mass, for the interaction of a central particle with all neighbour particles j, such that the force on  the particle 
	 from neighbour j is given by jac*(di-dj)*/
	double rx = - x_ref_pos[j];
	double ry = - y_ref_pos[j];
	double r = sqrt(rx*rx + ry*ry);

	double a = A(r,Q,l,q,d);
	double b = B(r,Q,l,q,d);
	double c = C(r,Q,l,q,d);

	// save the components
	save_jac[0][0] =  (a + rx*rx*b)/mass;
	save_jac[0][1] =  (ry*rx*b)/mass;
	save_jac[0][2] =  (rx*d*c)/mass;
	save_jac[1][0] =  (ry*rx*b)/mass;
	save_jac[1][1] =  (a + ry*ry*b)/mass;
	save_jac[1][2] =  (ry*d*c)/mass;
	save_jac[2][0] =  (rx*d*c)/mass;
	save_jac[2][1] =  (ry*d*c)/mass;
	save_jac[2][2] =  (a + d*d*c)/mass;
}

void add_self_interaction(double self_jac_m[3][3],double jac_one_m[3][3]){
	/*updates the selfinteraction part by adding the influence of one neighbour*/
	self_jac_m[0][0] += jac_one_m[0][0];
	self_jac_m[0][1] += jac_one_m[0][1];
	self_jac_m[0][2] += jac_one_m[0][2];
	self_jac_m[1][0] += jac_one_m[1][0];
	self_jac_m[1][1] += jac_one_m[1][1];
	self_jac_m[1][2] += jac_one_m[1][2];
	self_jac_m[2][0] += jac_one_m[2][0];
	self_jac_m[2][1] += jac_one_m[2][1];
	self_jac_m[2][2] += jac_one_m[2][2];
}

void calc_hesse_x_m(double save[3][3], int j, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double dist, double mass){
	/*Calculate the hessian for the x component of the force already devided by mass,
	 for the interaction of a central particle with a neighbour particles j.*/
	double rx = - x_ref_pos[j];
	double ry = - y_ref_pos[j];
	double r = sqrt(rx*rx + ry*ry);

	double b = B(r,Q,l,q,dist);
	double d = D(r,Q,l,q,dist);
	double c = C(r,Q,l,q,dist);
	double e = E(r,Q,l,q,dist);

	// save the components
	save[0][0] =  rx * (3*b + rx*rx*d) / mass;
	save[0][1] =  ry * (b + rx*rx*d) / mass;
	save[0][2] =  dist * (c + rx*rx*e) / mass;
	save[1][0] =  ry * (b + rx*rx*d) / mass;
	save[1][1] =  rx * (b + ry*ry*d) / mass;
	save[1][2] =  dist*rx*ry*e / mass;
	save[2][0] =  dist * (c + rx*rx*e) / mass;
	save[2][1] =  dist*rx*ry*e / mass;
	save[2][2] =  rx * (b + dist*dist*e) / mass;
}

void calc_hesse_y_m(double save[3][3], int j, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double dist, double mass){
	/*Calculate the hessian for the x component of the force already devided by mass,
	 for the interaction of a central particle with a neighbour particles j.*/
	double rx = - x_ref_pos[j];
	double ry = - y_ref_pos[j];
	double r = sqrt(rx*rx + ry*ry);

	double b = B(r,Q,l,q,dist);
	double d = D(r,Q,l,q,dist);
	double c = C(r,Q,l,q,dist);
	double e = E(r,Q,l,q,dist);

	// save the components
	save[0][0] =  ry * (b + rx*rx*d) / mass;
	save[0][1] =  rx * (b + ry*ry*d) / mass;
	save[0][2] =  dist*rx*ry*e / mass;
	save[1][0] =  rx * (b + ry*ry*d) / mass;
	save[1][1] =  ry * (3*b + ry*ry*d) / mass;
	save[1][2] =  dist * (c + ry*ry*e) / mass;
	save[2][0] =  dist*rx*ry*e / mass;
	save[2][1] =  dist * (c + ry*ry*e) / mass;
	save[2][2] =  ry * (b + dist*dist*e) / mass;
}

void calc_hesse_z_m(double save[3][3], int j, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double dist, double mass){
	/*Calculate the hessian for the x component of the force already devided by mass,
	 for the interaction of a central particle with a neighbour particles j.*/
	double rx = - x_ref_pos[j];
	double ry = - y_ref_pos[j];
	double r = sqrt(rx*rx + ry*ry);

	double b = B(r,Q,l,q,dist);
	double d = D(r,Q,l,q,dist);
	double c = C(r,Q,l,q,dist);
	double e = E(r,Q,l,q,dist);

	// save the components
	save[0][0] =  dist * (c + rx*rx*e) / mass;
	save[0][1] =  dist*rx*ry*e / mass;
	save[0][2] =  rx * (b + dist*dist*e) / mass;
	save[1][0] =  dist*rx*ry*e / mass;
	save[1][1] =  dist * (c + ry*ry*d) / mass;
	save[1][2] =  ry * (b + dist*dist*e) / mass;
	save[2][0] =  rx * (b + dist*dist*e) / mass;
	save[2][1] =  ry * (b + dist*dist*e) / mass;
	save[2][2] =  dist * (3*c + dist*dist*e) / mass;
}


void set_up_particle_interaction(double save_jacs_m[][3][3], double save_jac_self_m[3][3], int n, double *x_ref_pos, double *y_ref_pos, 
								 double Q, double l, double q, double d, double mass){
	/*Calculate the jacobians for all neighbour positions and already add up the self interaction part*/
	for(int i=0;i<n;i++){
		calc_jacobian_m(save_jacs_m[i], i, x_ref_pos, y_ref_pos, Q, l, q, d,mass);
		add_self_interaction(save_jac_self_m,save_jacs_m[i]);
	}
}

void set_up_particle_interaction_quadratic(double h_x_m[][3][3], double h_x_self_m[3][3], double h_y_m[][3][3], double h_y_self_m[3][3], 
										   double h_z_m[][3][3], double h_z_self_m[3][3],int n, double *x_ref_pos, double *y_ref_pos, 
								 		   double Q, double l, double q, double d, double mass){
	/*Calculate the hessians for all neighbour positions and already add up the self interaction part*/
	for(int i=0;i<n;i++){
		calc_hesse_x_m(h_x_m[i], i, x_ref_pos, y_ref_pos, Q, l, q, d, mass);
		add_self_interaction(h_x_self_m, h_x_m[i]);

		calc_hesse_y_m(h_y_m[i], i, x_ref_pos, y_ref_pos, Q, l, q, d, mass);
		add_self_interaction(h_y_self_m, h_y_m[i]);

		calc_hesse_z_m(h_z_m[i], i, x_ref_pos, y_ref_pos, Q, l, q, d, mass);
		add_self_interaction(h_z_self_m, h_z_m[i]);
	}
}


void calc_accelaration(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], double *vx, double *vy, double *vz, 
					   double fric, double OMEGA_Z_CONF, double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours){
	int ng;
	for(int i=0;i<N;i++){
		ax[i] = self_jac_m[0][0] * x[i] + self_jac_m[0][1] * y[i] + self_jac_m[0][2] * z[i] - fric * vx[i];
		ay[i] = self_jac_m[1][0] * x[i] + self_jac_m[1][1] * y[i] + self_jac_m[1][2] * z[i] - fric * vy[i];
		az[i] = self_jac_m[2][0] * x[i] + self_jac_m[2][1] * y[i] + (self_jac_m[2][2] - OMEGA_Z_CONF*OMEGA_Z_CONF) * z[i] - fric * vz[i];
		for(int j=0;j<n;j++){
			ng = neighbours[i*n+j];
			if(ng < 0){
				continue;
			}
			ax[i] = ax[i] - (jacs_m[j][0][0]*x[ng] + jacs_m[j][0][1]*y[ng] + jacs_m[j][0][2]*z[ng]);
			ay[i] = ay[i] - (jacs_m[j][1][0]*x[ng] + jacs_m[j][1][1]*y[ng] + jacs_m[j][1][2]*z[ng]);
			az[i] = az[i] - (jacs_m[j][2][0]*x[ng] + jacs_m[j][2][1]*y[ng] + jacs_m[j][2][2]*z[ng]);
		}
	}
}

void calc_acc_verlet(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours){
	int ng;
	for(int i=0;i<N;i++){
		ax[i] = self_jac_m[0][0] * x[i] + self_jac_m[0][1] * y[i] + self_jac_m[0][2] * z[i];
		ay[i] = self_jac_m[1][0] * x[i] + self_jac_m[1][1] * y[i] + self_jac_m[1][2] * z[i];
		az[i] = self_jac_m[2][0] * x[i] + self_jac_m[2][1] * y[i] + (self_jac_m[2][2] - OMEGA_Z_CONF*OMEGA_Z_CONF) * z[i];
		for(int j=0;j<n;j++){
			ng = neighbours[i*n+j];
			if(ng < 0){
				continue;
			}
			ax[i] = ax[i] - (jacs_m[j][0][0]*x[ng] + jacs_m[j][0][1]*y[ng] + jacs_m[j][0][2]*z[ng]);
			ay[i] = ay[i] - (jacs_m[j][1][0]*x[ng] + jacs_m[j][1][1]*y[ng] + jacs_m[j][1][2]*z[ng]);
			az[i] = az[i] - (jacs_m[j][2][0]*x[ng] + jacs_m[j][2][1]*y[ng] + jacs_m[j][2][2]*z[ng]);
		}
	}
}


void calc_acc_verlet_complete(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double d, double mass, int *neighbours){
	int ng;
	double xng = 0.;
	double yng = 0.;
	double zng = 0.;
	for(int i=0;i<N;i++){
		ax[i] = 0;
		ay[i] = 0;
		az[i] = (0.118705386) - OMEGA_Z_CONF*OMEGA_Z_CONF*z[i];
		for(int j=0;j<n;j++){
			ng = neighbours[i*n+j];
			xng = 0.0;
			yng = 0.0;
			zng = 0.0;
			if(ng > 0){
				xng = x[ng];
				yng = y[ng];
				zng = z[ng];
			}
			
			double rx = - x_ref_pos[j] - xng + x[i];
			double ry = - y_ref_pos[j] - yng + y[i];
			double rz = -zng + z[i];
			double r = sqrt(rx*rx + ry*ry + rz*rz);
			double rd = sqrt(rx*rx + ry*ry + (rz+d)*(rz+d));
			
			double fr = h1(r,Q,l)/mass;
			double frd= h1(rd,Q,l)/mass;
 
			ax[i] = ax[i] + (Q*fr + q*frd)*rx;
			ay[i] = ay[i] + (Q*fr + q*frd)*ry;
			az[i] = az[i] + Q*fr*rz + q*frd*(rz + d);
		}
	}
}

void calc_acc_verlet_nonlinear(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
							   double vz[], double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours){
	// calculate avarage quadratic z velocity
	double OMEGA_Z_CONF = 0.0;
	double vz_2 = 0.0;
	for(int i=0;i<N;i++){
		vz_2 += vz[i]*vz[i];
	}
	vz_2 = vz_2/((double) N);
	// adjust confinment
	OMEGA_Z_CONF = 125.53*(1 + vz_2*145603);
	
	int ng;
	for(int i=0;i<N;i++){
		ax[i] = self_jac_m[0][0] * x[i] + self_jac_m[0][1] * y[i] + self_jac_m[0][2] * z[i];
		ay[i] = self_jac_m[1][0] * x[i] + self_jac_m[1][1] * y[i] + self_jac_m[1][2] * z[i];
		az[i] = self_jac_m[2][0] * x[i] + self_jac_m[2][1] * y[i] + (self_jac_m[2][2] - OMEGA_Z_CONF*OMEGA_Z_CONF) * z[i];
		for(int j=0;j<n;j++){
			ng = neighbours[i*n+j];
			if(ng < 0){
				continue;
			}
			ax[i] = ax[i] - (jacs_m[j][0][0]*x[ng] + jacs_m[j][0][1]*y[ng] + jacs_m[j][0][2]*z[ng]);
			ay[i] = ay[i] - (jacs_m[j][1][0]*x[ng] + jacs_m[j][1][1]*y[ng] + jacs_m[j][1][2]*z[ng]);
			az[i] = az[i] - (jacs_m[j][2][0]*x[ng] + jacs_m[j][2][1]*y[ng] + jacs_m[j][2][2]*z[ng]);
		}
	}
}

void calc_acc_verlet_quadratic(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double jacs_m[][3][3], double h_x_m[][3][3], double h_y_m[][3][3], 
					   double h_z_m[][3][3], int *neighbours){
	int ng;
	
	// save current vectors
	double di[3];
	double dj[3];
	double dji[3];
	
	for(int i=0;i<N;i++){
		di[0] = x[i];
		di[1] = y[i];
		di[2] = z[i];
		ax[i] = 0.0;
		ay[i] = 0.0;
		az[i] = - OMEGA_Z_CONF*OMEGA_Z_CONF * z[i];
		for(int j=0;j<n;j++){
			ng = neighbours[i*n+j];
			if(ng < 0){
				dj[0] = 0.0;
				dj[1] = 0.0;
				dj[2] = 0.0;
			}
			else{
				dj[0] = x[ng];
				dj[1] = y[ng];
				dj[2] = z[ng];
			}
			dji[0] = di[0]-dj[0];
			dji[1] = di[1]-dj[1];
			dji[2] = di[2]-dj[2];
			
			ax[i] += dji[0]*(jacs_m[j][0][0] + h_x_m[j][0][0]*dji[0] + h_x_m[j][0][1]*dji[1] + h_x_m[j][0][2]*dji[2]);
			ax[i] += dji[1]*(jacs_m[j][0][1] + h_x_m[j][1][0]*dji[0] + h_x_m[j][1][1]*dji[1] + h_x_m[j][1][2]*dji[2]);
			ax[i] += dji[2]*(jacs_m[j][0][2] + h_x_m[j][2][0]*dji[0] + h_x_m[j][2][1]*dji[1] + h_x_m[j][2][2]*dji[2]);
			
			ay[i] += dji[0]*(jacs_m[j][1][0] + h_y_m[j][0][0]*dji[0] + h_y_m[j][0][1]*dji[1] + h_y_m[j][0][2]*dji[2]);
			ay[i] += dji[1]*(jacs_m[j][1][1] + h_y_m[j][1][0]*dji[0] + h_y_m[j][1][1]*dji[1] + h_y_m[j][1][2]*dji[2]);
			ay[i] += dji[2]*(jacs_m[j][1][2] + h_y_m[j][2][0]*dji[0] + h_y_m[j][2][1]*dji[1] + h_y_m[j][2][2]*dji[2]);
			
			az[i] += dji[0]*(jacs_m[j][2][0] + h_z_m[j][0][0]*dji[0] + h_z_m[j][0][1]*dji[1] + h_z_m[j][0][2]*dji[2]);
			az[i] += dji[1]*(jacs_m[j][2][1] + h_z_m[j][1][0]*dji[0] + h_z_m[j][1][1]*dji[1] + h_z_m[j][1][2]*dji[2]);
			az[i] += dji[2]*(jacs_m[j][2][2] + h_z_m[j][2][0]*dji[0] + h_z_m[j][2][1]*dji[1] + h_z_m[j][2][2]*dji[2]);
		}
	}
}

void update_step(int N, double *start_x, double *start_y, double *start_z, double *add_x, double *add_y, double *add_z, double factor,
	double *save_x, double *save_y, double *save_z){
	for(int i=0;i<N;i++){
		save_x[i] = start_x[i] + add_x[i]*factor;
		save_y[i] = start_y[i] + add_y[i]*factor;
		save_z[i] = start_z[i] + add_z[i]*factor;
	}
}

void update_final(int N, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *bx, double *cx, double *dx,
	double *ay, double *by, double *cy, double *dy, double *az, double *bz, double *cz, double *dz, double time_step){
	for(int i=0;i<N;i++){
		x[i] = x[i] + vx[i]*time_step + time_step*time_step*(ax[i]+bx[i]+cx[i])/6;
		y[i] = y[i] + vy[i]*time_step + time_step*time_step*(ay[i]+bx[i]+cx[i])/6;
		z[i] = z[i] + vz[i]*time_step + time_step*time_step*(az[i]+bz[i]+cz[i])/6;

		vx[i] = vx[i] + time_step*(ax[i] + bx[i]*2 + cx[i]*2 + dx[i])/6;
		vy[i] = vy[i] + time_step*(ay[i] + by[i]*2 + cy[i]*2 + dy[i])/6;
		vz[i] = vz[i] + time_step*(az[i] + bz[i]*2 + cz[i]*2 + dz[i])/6;
	}

}
