import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
from lattice import*
from wraper_for_c_simulation import*

## path to compiled c-code for the used solver
FOLDER_PATH = '/home/mh/Desktop/RunMCISimulation'

SOLVER_FILE_NAME = 'linear_solver.so'
#SOLVER_FILE_NAME = 'complete_solver.so'

SOLVER_PATH = FOLDER_PATH + '/' + SOLVER_FILE_NAME

load_solver(SOLVER_PATH)

## PHYCISAL CONSTANTS
EL_CHARGE = 1.602176487e-19        # electron charge
BOLTZMANN_CONST = 1.38064852E-23

## SYSTEM PARAMETER
LAT_CONST = 4.804e-04              # lattice constant
NY = 5                             # number of particles in y-direction
# the number of particles in x-direction (NX) will be choosen automatically
# in such  a way that a quadratic lattice will be simulated. 
# (remark: in case of a deformed lattice NY = NX not true) 
INTER_RANGE = 10.0 * LAT_CONST     # finite interaction range
CHARGE_Q = -19000 * EL_CHARGE      # particle charge
CHARGE_q = 0.2 * abs(CHARGE_Q)     # wake charge
SCREEN_LENGHT = 3.800e-04          # screening length of potential
WAKE_DIST = 0.3 * SCREEN_LENGHT    # wake distance from particle
MASS = 6.1e-13                     # particle mass
OMEGA_Z_CONF = 19.8 * 2 * np.pi    # confinement strength
FRICTION_PARA = 1.26               # friction parameter
TEMP = 300.0                       # temperatur in Kelvin

## SET TIME ASPECTS OF SIMULATION
T_MAX = 0.1                       # simulation starts at t=0 and ends at t=T_MAX
T_STEP = 1.0e-04                  # integration step size h
SAVE_INTERVAL = 1.0/1000          # resulution: time difference between to saved points
#---------- help stuff -----------#
SAVE_PERIOD = int(SAVE_INTERVAL / T_STEP)   # save every SAVE_PERIOD step...
SAVE_POINTS = int((T_MAX / T_STEP) / SAVE_PERIOD)   # how many points to save?
#---------------------------------#

## CHOOSE LATTICE CONFIGURATION
# set not deformed values as default
mu = 1.0                           # compresion factor
m = 0.0                            # share slope

# values representing Ingos simulation for compression under 0 degree
H0 = 0                             
MU0 = 1.04798011295                
# vlaues representing Ingox simulation fo compression under 30 degree
H30 = -0.268885015108e-04         
MU30 = 1.01429243175               

angle = None # choose 0 or 30 to simulate one of Ingos configuations

if angle == 0:
    # compression under 0 degree lattice situation
    mu = MU0
    m = (2.0 * H0) / (LAT_CONST * np.sqrt(3))
elif angle == 30:
    #compression under 30 degree lattice situation
    mu = MU30
    m = (2.0 * H30) / (LAT_CONST * np.sqrt(3))  

#---------- help stuff -----------#
# CALCULATE HELP FACTORS
alpha = (np.sqrt(3.0) / (2.0 * mu**2))
betha = np.sqrt(3.0) * m / mu

# CALCULATE SYSTEM SIZE
NX = int((NY - 1) / alpha) + 1  # number of particles in x direction
N = NY * NX                     # total number of particles

# CALCULATE INTERACTION SIZE
ny = int(INTER_RANGE / (LAT_CONST * mu)) + 1
nx = int(ny / alpha)
n = (2 * ny + 1) * (2 * nx + 1) - 1     # total number of neighbours
#----------------------------------#

## LATTICE UNDER INVESTIGATION
L = Lattice(NY, LAT_CONST, mu, m)

# neighbours
neighbour_indices = L.get_neighbour_indices(INTER_RANGE)

# lattice
equilibrium_lattice = L.get_lattice_vectors()

# next neighbours
next_neighbour_indices = L.get_nn_indices()


## SET INITIAL POSITIONS
# let find equlibrium state from zero movement
pos_initial = np.zeros(3*N)
vel_initial = np.zeros(3*N)
vec_initial = np.concatenate((pos_initial, vel_initial))


## INTEGRATE SYSTEM IN C

# set up arrays to save results
save_t = np.zeros(SAVE_POINTS, dtype=np.float64)
x_pos = np.zeros((N, SAVE_POINTS), dtype=np.float64)
y_pos = np.zeros((N, SAVE_POINTS), dtype=np.float64)
z_pos = np.zeros((N, SAVE_POINTS), dtype=np.float64)

x_vel = np.zeros((N, SAVE_POINTS), dtype=np.float64)
y_vel = np.zeros((N, SAVE_POINTS), dtype=np.float64)
z_vel = np.zeros((N, SAVE_POINTS), dtype=np.float64)

PARAS = np.array([NY, INTER_RANGE, FRICTION_PARA, MASS, TEMP, mu, m,
                  LAT_CONST, CHARGE_Q, CHARGE_q, SCREEN_LENGHT, WAKE_DIST,
                  OMEGA_Z_CONF], dtype=np.float64)

# INTEGRATE SYSTEM
solve(T_STEP, T_MAX, SAVE_PERIOD, PARAS, x_pos, y_pos,
      z_pos, x_vel, y_vel, z_vel, save_t, vec_initial)

# SAVE THE RESULT
time_stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

deviations = [x_pos, y_pos, z_pos]
velocities = [x_vel, y_vel, z_vel]

save_path = FOLDER_PATH

os.mkdir(save_path + '/' + time_stamp + '-simulation_result')
np.save(save_path + '/' + time_stamp + '-simulation_result/time', save_t)
np.save(save_path + '/' + time_stamp + '-simulation_result/deviations', deviations)
np.save(save_path + '/' + time_stamp + '-simulation_result/velocities', velocities)

np.save(save_path + '/' + time_stamp + '-simulation_result/equilibrium_lattice', equilibrium_lattice)
np.save(save_path + '/' + time_stamp + '-simulation_result/neighbour_indices', neighbour_indices)
np.save(save_path + '/' + time_stamp + '-simulation_result/next_neighbour_indices', next_neighbour_indices)

plt.figure()
plt.plot(save_t, x_pos[4])
plt.plot(save_t, y_pos[4])
plt.plot(save_t, z_pos[4])
plt.show()

EkC_x = np.mean(x_vel**2, axis=0) * 0.5 * MASS
EkC_y = np.mean(y_vel**2, axis=0) * 0.5 * MASS
EkC_z = np.mean(z_vel**2, axis=0) * 0.5 * MASS
Ekin_C = (EkC_x + EkC_y + EkC_z) / 3.0

E_equpart = 0.5 * BOLTZMANN_CONST * TEMP
# Energy per degree of freedom according to equpartition theorem

plt.figure()
plt.title('Ekin for different dimensions C')
plt.plot(save_t, EkC_x, label='x')
plt.plot(save_t, EkC_y, label='y')
plt.plot(save_t, EkC_z, label='z')
plt.axhline(E_equpart)
plt.legend(loc='best')


plt.figure()
plt.title('Ekin total')
plt.plot(save_t, Ekin_C, label='Ekin c')
plt.axhline(E_equpart)
plt.legend(loc='best')
plt.ylim([0, 2e-20])

plt.show()

