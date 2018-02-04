
#ifndef HEADER_FILE
#define HEADER_FILE

// set_up_lattice.c
int shift(double betha, int x);
int neighbour_index(int nx, int ny, int lx, int ly);
void calc_neighbours(int *save_neig, int Ny, int Nx, int ny, int nx, double betha);
void calc_lattice(double *save_x, double *save_y, int Nx, int Ny, double alpha, double betha, double rescaled_lat_const);
void calc_ref_lattice(double *save_x, double *save_y, int nx, int ny, double alpha, double betha, double rescaled_lat_const);

// interactions.c
void set_up_particle_interaction(double save_jacs_m[][3][3], double save_jac_self_m[3][3], int n, double *x_ref_pos,
								 double *y_ref_pos, double Q, double l, double q, double d, double mass);

void set_up_particle_interaction_quadratic(double h_x_m[][3][3], double h_x_self_m[3][3], double h_y_m[][3][3], double h_y_self_m[3][3], 
										   double h_z_m[][3][3], double h_z_self_m[3][3],int n, double *x_ref_pos, double *y_ref_pos, 
								 		   double Q, double l, double q, double d, double mass);

void calc_accelaration(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], double *vx, double *vy, double *vz, 
					   double fric, double OMEGA_Z_CONF, double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours);

void calc_acc_verlet(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours);
					   
void calc_acc_verlet_nonlinear(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
							   double vz[], double self_jac_m[3][3], double jacs_m[][3][3], int *neighbours);

void calc_acc_verlet_quadratic(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double jacs_m[][3][3], double h_x_m[][3][3], double h_y_m[][3][3], 
					   double h_z_m[][3][3], int *neighbours);

void calc_acc_verlet_complete(int N, int n, double ax[], double ay[], double az[], double x[], double y[], double z[], 
					   double OMEGA_Z_CONF, double *x_ref_pos, double *y_ref_pos, double Q, double l, double q, double d, double mass, int *neighbours);

void update_step(int N, double *start_x, double *start_y, double *start_z, double *add_x, double *add_y, double *add_z, double factor,
				 double *save_x, double *save_y, double *save_z);

void update_final(int N, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *bx, double *cx, double *dx,
				  double *ay, double *by, double *cy, double *dy, double *az, double *bz, double *cz, double *dz, double time_step);



//RNG.c
double gaussrnd();


#endif
