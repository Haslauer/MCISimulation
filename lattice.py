# python class which defines a deformed hexagonal lattice
# it  will allways be calculated such that we get a cubic
# lattice wich is equally extended in x and y direction

import numpy as np


def calc_y_shift(betha, x):
    "Helpfunction to calculate the bottom point of each x-column."
    temp = (betha + 0.5) * x
    shift = -int(temp)
    if temp + shift < 0:
        shift += 1
    return shift


def calc_neigh_index(nx, ny, lx, ly):
    k = lx * (2 * ny + 1) + ly
    n = (2 * nx + 1) * (2 * ny + 1) - 1

    j = -10
    if(k < n / 2):
        j = k
    elif(k > n / 2):
        j = k - 1
    return j


class Lattice:
    "Describes a quadratic 2D plasma crystal lattice."

    def __init__(self, N, a, mu, m):

        self.max_index_y = N - 1
        # automatically adjuste x-direction to get a quadratic lattice
        # dependend on the compression strength
        self.max_index_x = int((2 / np.sqrt(3)) * (N - 1) * mu**2)
        self.lattice_const = a
        self.compression_factor = mu
        self.share_slope = m

    def get_deformation_matrix(self):
        "returns the 2d transformation matrix which transforms the lattice \
        vectors of a perfect hexagonal lattice into the deformed ones by \
        matrix multiplication. lat_vec_deformed = defrom_matrix*lat_vec."

        deform_matrix = np.zeros((2, 2))
        deform_matrix[0, 0] = 1.0 / self.compression_factor
        deform_matrix[1, 1] = self.compression_factor
        deform_matrix[1, 0] = self.share_slope

        return deform_matrix

    def get_lattice_vectors(self):
        "returns the 2d lattice voctors of the lattice.\
        return: [particle index][dimension]"

        # total number of particles
        number_particles_total = (
            self.max_index_x + 1) * (self.max_index_y + 1)
        # save lattice vectors
        lattice_vecs = np.zeros((number_particles_total, 2))

        # set up basis vectors for not deformed lattice
        a1 = np.array([0, 1]) * self.lattice_const
        a2 = np.array([np.sqrt(3), 1]) * self.lattice_const / 2
        # transformation matrix
        T = self.get_deformation_matrix()
        # transform basis vectors
        a1 = np.dot(T, a1)
        a2 = np.dot(T, a2)

        # just a helping factor for the calculatin of the quadratic lattice
        beta = (np.sqrt(3) / 2) * (self.share_slope / self.compression_factor)

        # set up lattice with 1d index convention, index = x*(max_index_y+1)+y
        # where x,y are the x and y position 2d index.

        Nx = self.max_index_x + 1
        Ny = self.max_index_y + 1

        for y in range(Ny):
            for x in range(Nx):
                lattice_vecs[x * Ny + y] = x * a2 + \
                    (y + calc_y_shift(beta, x)) * a1

        return lattice_vecs

    def get_neighbour_indices(self, interaction_range):
        "returns the indices of the neighbours for each particle.\
        return: [particle index][neighbour index] = particle index \
        of neighbour index."

        ny = int(interaction_range / (
                 self.lattice_const * self.compression_factor)) + 1
        nx = int((2 * ny * self.compression_factor**2) / np.sqrt(3))

        number_neighbours_total = (2 * nx + 1) * (2 * ny + 1)

        number_particles_total = (self.max_index_x + 1) * (
            self.max_index_y + 1)

        # save the indices
        neighbours = np.ones((number_particles_total,
                              number_neighbours_total), dtype=int) * (-1)

        Nx = self.max_index_x + 1
        Ny = self.max_index_y + 1

        beta = (np.sqrt(3) / 2) * (self.share_slope / self.compression_factor)

        for x in range(Nx):
            for y in range(Ny):
                temp = x * Ny + y
                for lx in range(2 * nx + 1):
                    mx = -nx + lx
                    x_temp = x + mx
                    if x_temp < 0:
                        continue
                    elif x_temp > (Nx - 1):
                        continue
                    for ly in range(2 * ny + 1):
                        my = -ny + ly
                        y_temp = y + my + \
                            calc_y_shift(beta, x) + calc_y_shift(beta, mx) -\
                            calc_y_shift(beta, x_temp)
                        if y_temp < 0:
                            continue
                        elif y_temp > (Ny - 1):
                            continue

                        temp2 = lx * (2 * ny + 1) + ly
                        neighbours[temp][temp2] = (x_temp) * Ny + (y_temp)

        # delte the identity index and return (i is not a neighbour of itself)
        return np.delete(neighbours, (2 * ny + 1) * nx + ny, 1)

    def get_nn_indices(self):
        """
        Calculates the particle index for the next neighbours of each
        particle. There are 6 next neighbours for each particle.
        return: [particle index][neighbour index]
        """

        Nx = self.max_index_x + 1
        Ny = self.max_index_y + 1
        N = Nx * Ny
        next_neighbours = np.ones((N, 6))*(-1)

        beta = (np.sqrt(3) / 2) * (self.share_slope / self.compression_factor)

        for x in range(Nx):
            for y in range(Ny):
                index = x * Ny + y

                nx = +0
                ny = +1
                xi = x + nx
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][0] = xi*Ny + yi

                nx = +1
                ny = +0
                xi = x + nx
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][1] = xi*Ny + yi

                nx = +1
                ny = -1
                xi = x + nx
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][2] = xi*Ny + yi

                nx = +0
                ny = -1
                xi = x + nx
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][3] = xi*Ny + yi

                nx = -1  # for nx < 0 we have to subtract -1
                ny = +0  # additionally to compensate the calc_y_shift
                xi = x + nx  # function.
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - 1 - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][4] = xi*Ny + yi

                nx = -1  # for nx < 0 we have to subtract -1
                ny = +1  # additionaly to compensate the calc_y_shift
                xi = x + nx  # function.
                yi = y + ny + calc_y_shift(beta, x) + (
                              calc_y_shift(beta, nx)) - 1 - (
                              calc_y_shift(beta, xi))

                if(xi >= 0 and xi < Nx):
                    if(yi >= 0 and yi < Ny):
                        next_neighbours[index][5] = xi*Ny + yi

        return next_neighbours

    def get_reference_lattice(self, interaction_range):
        "returns the distance for each neighbour with index j."

        # set up basis vectors for not deformed lattice
        a1 = np.array([0, 1]) * self.lattice_const
        a2 = np.array([np.sqrt(3), 1]) * self.lattice_const / 2
        # transformation matrix
        T = self.get_deformation_matrix()
        # transform basis vectors
        a1 = np.dot(T, a1)
        a2 = np.dot(T, a2)

        ny = int(interaction_range / (
                 self.lattice_const * self.compression_factor)) + 1
        nx = int((2 * ny * self.compression_factor**2) / np.sqrt(3))
        n = (2 * nx + 1) * (2 * ny + 1) - 1
        rel_pos = np.zeros((n, 3))

        betha = (np.sqrt(3) / 2) * (self.share_slope / self.compression_factor)

        for lx in range(2 * nx + 1):
            mx = -nx + lx
            for ly in range(2 * ny + 1):
                my = -ny + ly
                if(my == 0 and mx == 0):
                    continue
                temp = mx * a2 + (my + calc_y_shift(betha, mx)) * a1
                index = calc_neigh_index(nx, ny, lx, ly)
                rel_pos[index][0] = temp[0]
                rel_pos[index][1] = temp[1]

        return rel_pos
