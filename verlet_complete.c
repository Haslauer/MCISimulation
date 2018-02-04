/* numerical solver for the linear equations of motion for 2D hexagonal plasma crystal*/

#include <stdio.h>
#include <math.h>
#include "helpers.h"

void set_initial_values_to_zeros(int N,double x[],double y[],double z[]){
	/*for the moment set everything to zero*/
	for(int i=0;i<N;i++){
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
	}
}


void solve(const double T_STEP, const double T_MAX, const int SAVE_PERIOD, double *para, double **save_xpos, double **save_ypos, double **save_zpos,
double **save_xvel, double **save_yvel, double **save_zvel, double *save_t, double *vec_initial){

	/*System is fully characterized by the following System parameters*/
	const int NY =  (int)para[0];
	const double INTER_RANGE =  para[1];
	const double FRICTION_PARA = para[2];
	const double MASS = para[3];
	const double TEMP = para[4];
	const double COMPRESSION_FACTOR = para[5];
	const double SHARE_SLOPE = para[6];
	const double LAT_CONST = para[7];
	const double CHARGE_Q = para[8];
	const double CHARGE_q = para[9];
	const double SCREEN_LENGTH =  para[10];
	const double WAKE_DIST = para[11];
	const double OMEGA_Z_CONF = para[12];
	


	// CALCULATE STDDEV FOR NOISE
	const double BOLTZMANN_CONST = 1.38064852E-23;
	const double NOISE_FACTOR = sqrt((2.0*BOLTZMANN_CONST*TEMP*FRICTION_PARA*T_STEP)/MASS);
	// CALCULATE VERLET FACTOR
	const double VERLET_FACTOR = 2.0/(2.0 + FRICTION_PARA*T_STEP);
	// CALCULATE HELP FACTORS
	const double ALPHA = (sqrt(3.0)/2.0)/(COMPRESSION_FACTOR*COMPRESSION_FACTOR);
	const double BETA = (sqrt(3.0)/2.0)*(SHARE_SLOPE/COMPRESSION_FACTOR);

	// NUMBER OF PARTICLES
	const int NX = (int)((NY-1) / ALPHA) + 1;
	const int N = NX * NY;

	// NUMBER OF NEIGHBOURS
	const int ny = (int)(INTER_RANGE/(LAT_CONST*COMPRESSION_FACTOR)) + 1;
	const int nx = (int)(ny/ALPHA);
	const int n = (2*ny+1)*(2*nx+1)-1;

	// CALCULATE NEIGHBOURS
	/*
	each particle has n neighbours where the index j = 0,...,n-1 reffers always
	to the neighbour with identical relative position. If now i ist the particle 
	index, and j is the neighbour index we have the neighbour j of particle i on
	positon i*n+j of the neighbours array.
	*/

	int neighbours[n*N]; // has to be implemented as dynamic  memory allocation to make big systems possible!
	calc_neighbours(neighbours,NY,NX,ny,nx,BETA);

	// CALCULATE REFERENCE LATTICE FOR THE CALCULATION OF THE JACOBIANS
	double x_ref_pos[n];
	double y_ref_pos[n];
	calc_ref_lattice(x_ref_pos,y_ref_pos,nx,ny,ALPHA,BETA,LAT_CONST*COMPRESSION_FACTOR);

	// CALCULATE JACOBIANS AND DIVIDE ALREADY BY MASS
	double jacobians_m[n][3][3];
	double self_jacobian_m[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // save self interaction part
	set_up_particle_interaction(jacobians_m, self_jacobian_m, n, x_ref_pos, y_ref_pos, CHARGE_Q, SCREEN_LENGTH, CHARGE_q, WAKE_DIST, MASS);


	// SET UP POSITION AND VELOCITY VECTORS
	// remark positions are just the deviations from the equlibrium lattice positions
	double x[N];
	double y[N];
	double z[N];
	set_initial_values_to_zeros(N,x,y,z);
	double x_new[N];
	double y_new[N];
	double z_new[N];
	set_initial_values_to_zeros(N,x_new,y_new,z_new);
	double vx[N];
	double vy[N];
	double vz[N];
	set_initial_values_to_zeros(N,vx,vy,vz);

	// save accelerations
	double ax[N];
	double ay[N];
	double az[N];
	set_initial_values_to_zeros(N, ax, ay, az);
	double ax_new[N];
	double ay_new[N];
	double az_new[N];
	set_initial_values_to_zeros(N, ax_new, ay_new, az_new);

	// save noise
	double noise_x[N];
	double noise_y[N];
	double noise_z[N];
	
	// warm up random generator
	double just_safe_without_sense = 0;
	for(int s=0;  s<5000; s++){
		just_safe_without_sense = gaussrnd();
	}


	/*integration*/
	
	int max_iter = (int)(T_MAX / T_STEP);
	double t = 0.0;
	int count = 0;
	int save_count = 0;
	int save_index = 0;
	int progress = 0;

	while(count<max_iter){

		
		// update positions 
		for(int i=0; i<N; i++){
			noise_x[i] = gaussrnd()*NOISE_FACTOR;
			noise_y[i] = gaussrnd()*NOISE_FACTOR;
			noise_z[i] = gaussrnd()*NOISE_FACTOR;
			x_new[i] = x[i] + (2.0*vx[i] + ax[i]*T_STEP + noise_x[i])*VERLET_FACTOR*T_STEP/2.0;
			y_new[i] = y[i] + (2.0*vy[i] + ay[i]*T_STEP + noise_y[i])*VERLET_FACTOR*T_STEP/2.0;
			z_new[i] = z[i] + (2.0*vz[i] + az[i]*T_STEP + noise_z[i])*VERLET_FACTOR*T_STEP/2.0;
		}

		// update accelaration
		calc_acc_verlet_complete(N, n, ax_new, ay_new, az_new, x_new, y_new, z_new, OMEGA_Z_CONF, x_ref_pos, y_ref_pos, CHARGE_Q, SCREEN_LENGTH, CHARGE_q, WAKE_DIST, MASS, neighbours);

		//update velocity
		for(int i=0; i<N; i++){
			vx[i] = vx[i] + (ax[i] + ax_new[i])*T_STEP/2.0 - FRICTION_PARA*(x_new[i] - x[i]) + noise_x[i];
			vy[i] = vy[i] + (ay[i] + ay_new[i])*T_STEP/2.0 - FRICTION_PARA*(y_new[i] - y[i]) + noise_y[i];
			vz[i] = vz[i] + (az[i] + az_new[i])*T_STEP/2.0 - FRICTION_PARA*(z_new[i] - z[i]) + noise_z[i];

			// n+1 -> n, (save new in old)
			x[i] = x_new[i];
			y[i] = y_new[i];
			z[i] = z_new[i];
			ax[i] = ax_new[i];
			ay[i] = ay_new[i];
			az[i] = az_new[i];
		}
		
		// show progress
		if(count % (max_iter/10) == 0){
			progress++;
			printf("Progress: %i %% \n", (int)(progress*10));
		}
		// update counter and etc.
		t = t + T_STEP;
		count++;
		save_count++;

		
		if(save_count == SAVE_PERIOD){
			save_t[save_index] = t;
			for(int part_nr=0; part_nr<N; part_nr++){
			save_xpos[part_nr][save_index] = x[part_nr];
			save_ypos[part_nr][save_index] = y[part_nr];
			save_zpos[part_nr][save_index] = z[part_nr];

			save_xvel[part_nr][save_index] = vx[part_nr];
			save_yvel[part_nr][save_index] = vy[part_nr];
			save_zvel[part_nr][save_index] = vz[part_nr];
			}
			save_index++;
			save_count = 0;
		}
	}
}
