import ctypes as ct
import numpy as np

################################################################

# TO USE 2D ARRAYS WE HAVE TO CONVERT IT INTO AN ARRAY OF POINTERS
# WHICH CAN BE HANDLED BY C


def convert_2D_array_to_array_of_pointers(arry):

    arrypp = (arry.__array_interface__['data'][0] +
              np.arange(arry.shape[0]) * arry.strides[0]).astype(np.uintp)
    return arrypp


INT = ct.c_int
DOUBLE = ct.c_double
ARRAY = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
ARRAYOFPOINTER = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C')


def load_solver(solver_path):
	global _solver_minimal
	_solver_minimal = ct.CDLL(solver_path)
	_solver_minimal.solve.argtypes = (DOUBLE, DOUBLE, INT, ARRAY, ARRAYOFPOINTER,
		                          ARRAYOFPOINTER, ARRAYOFPOINTER,
		                          ARRAYOFPOINTER, ARRAYOFPOINTER,
		                          ARRAYOFPOINTER, ARRAY, ARRAY)


def solve(t_step, t_max, save_interval, paras, xpos, ypos, zpos,
          xvel, yvel, zvel, save_t, vec_initial):

    # convert position arrays to arrays of pointers
    # to save the trajectories
    xpospp = convert_2D_array_to_array_of_pointers(xpos)
    ypospp = convert_2D_array_to_array_of_pointers(ypos)
    zpospp = convert_2D_array_to_array_of_pointers(zpos)

    xvelpp = convert_2D_array_to_array_of_pointers(xvel)
    yvelpp = convert_2D_array_to_array_of_pointers(yvel)
    zvelpp = convert_2D_array_to_array_of_pointers(zvel)

    global _solver_minimal

    # integrate the system
    _solver_minimal.solve(t_step, t_max, save_interval, paras, xpospp, ypospp,
                          zpospp, xvelpp, yvelpp, zvelpp, save_t, vec_initial)
